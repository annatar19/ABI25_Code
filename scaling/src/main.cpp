#include <iostream>
#include <ROOT/RNTupleReader.hxx>

#include <ROOT/RNTupleUtil.hxx>
#include <ROOT/RNTupleView.hxx>
#include <TCanvas.h>
#include <TH1I.h>
#include <TROOT.h>
#include <TString.h>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <ostream>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>
using json = nlohmann::json;

using FVec = std::variant<std::vector<std::string>, std::vector<std::int32_t>,
                          std::vector<std::uint32_t>, std::vector<double>>;
using VVec =
    std::variant<ROOT::RNTupleView<std::string>,
                 ROOT::RNTupleView<std::int32_t>,
                 ROOT::RNTupleView<std::uint32_t>, ROOT::RNTupleView<double>>;

std::vector<std::pair<std::string, std::string>>
GetFieldNamesAndTypes(const ROOT::REntry &entry) {
  std::vector<std::pair<std::string, std::string>> fields;

  for (const auto &val : entry) {
    fields.emplace_back(val.GetField().GetFieldName(),
                        val.GetField().GetTypeName());
  }

  std::sort(fields.begin(), fields.end(),
            [](auto const &a, auto const &b) { return a.first < b.first; });

  return fields;
}

template <typename T>
inline void ProcessNumberField(std::ostream &out, const std::string &fieldName,
                               ROOT::RNTupleView<T> &view,
                               std::uint64_t entryId) {
  out << '"' << fieldName << "\":" << '"' << view(entryId) << '"';
}

inline void ProcessStringField(std::ostream &out, std::string_view fieldName,
                               ROOT::RNTupleView<std::string> &view,
                               std::uint64_t entryId) {
  out << '"' << fieldName << "\":\"" << view(entryId) << '"';
}

enum FieldTypes { String, Int32, Uint32, Double };

int main(int argc, char **argv) {
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " <input.ntuple.root> <RNTuple name>"
              << std::endl;
    return 1;
  }
  const char *kNTupleFileName = argv[1];
  const char *kNTupleName = argv[2];

  auto reader = ROOT::RNTupleReader::Open(kNTupleName, kNTupleFileName);
  auto fields = GetFieldNamesAndTypes(reader->GetModel().GetDefaultEntry());

  std::vector<std::pair<std::string, int>> fieldMap;
  std::vector<FVec> fieldsVec;
  std::vector<VVec> viewVec;
  std::cout << "Initializing arrays…" << std::endl;

  for (const auto &[fieldName, fieldType] : fields) {
    if (fieldType == "std::string") {
      fieldsVec.emplace_back(std::in_place_type<std::vector<std::string>>,
                             reader->GetNEntries());
      viewVec.emplace_back(reader->GetView<std::string>(fieldName));
      fieldMap.emplace_back(fieldName, FieldTypes::String);
    } else if (fieldType == "std::int32_t") {
      fieldsVec.emplace_back(std::in_place_type<std::vector<std::int32_t>>,
                             reader->GetNEntries());
      viewVec.emplace_back(reader->GetView<std::int32_t>(fieldName));
      fieldMap.emplace_back(fieldName, FieldTypes::Int32);
    } else if (fieldType == "std::uint32_t") {
      std::cout << "uint!\n";
      fieldsVec.emplace_back(std::in_place_type<std::vector<std::uint32_t>>,
                             reader->GetNEntries());
      viewVec.emplace_back(reader->GetView<std::uint32_t>(fieldName));
      fieldMap.emplace_back(fieldName, FieldTypes::Uint32);
    } else if (fieldType == "double") {
      fieldsVec.emplace_back(std::in_place_type<std::vector<double>>,
                             reader->GetNEntries());
      viewVec.emplace_back(reader->GetView<double>(fieldName));
      fieldMap.emplace_back(fieldName, FieldTypes::Double);
    } else {
      std::cerr << "Found an unsupported fieldtype: " << fieldType << std::endl;
      return 1;
    }
  }

  std::cout << "Reading RNTuple into memory…" << std::endl;

  for (size_t field = 0; field < fieldsVec.size(); ++field) {
    const auto &[fieldName, fieldType] = fields[field];
    std::visit(
        [&](auto &dstVec, auto &srcView) {
          using D = std::decay_t<decltype(dstVec)>;
          using Elem =
              std::remove_const_t<std::remove_reference_t<decltype(srcView(
                  (ROOT::NTupleSize_t)0))>>;
          if constexpr (std::is_same_v<typename D::value_type, Elem>) {
            for (size_t i = 0; i < reader->GetNEntries(); ++i) {
              dstVec[i] = srcView(static_cast<ROOT::NTupleSize_t>(i));
            }
          }
        },
        fieldsVec[field], viewVec[field]);
  }

  std::cout << "Converting to JSON…" << std::endl;

  std::vector<std::string> fieldNames;
  for (const auto &[fieldName, fieldType] : fields) {
    fieldNames.emplace_back(fieldName);
  }
  std::ofstream ofs(std::string(kNTupleName) + ".ndjson");
  for (size_t i = 0; i < reader->GetNEntries(); ++i) {
    json obj = json::object();
    for (size_t f = 0; f < fieldNames.size(); ++f) {
      const auto &key = fieldNames[f];
      std::visit([&](auto &vec) { obj[key] = vec[i]; }, fieldsVec[f]);
    }
    ofs << obj.dump() << "\n";
  }

  std::cout << "Done converting RNTuple to JSON!" << std::endl;

  return 0;
}

