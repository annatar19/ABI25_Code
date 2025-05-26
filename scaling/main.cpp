#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleWriteOptions.hxx>
#include <ROOT/RNTupleWriter.hxx>
#include <ROOT/RNTuple.hxx>
// #include <ROOT/RNTupleInspector.hxx>
#include <iostream>

#include <ROOT/RNTupleUtil.hxx>
#include <ROOT/RNTupleView.hxx>
#include <TCanvas.h>
#include <TH1I.h>
#include <TROOT.h>
#include <TString.h>

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
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

using FVec = std::variant<std::vector<std::string>, std::vector<std::int32_t>,
                          std::vector<std::uint32_t>, std::vector<double>>;
using VVec =
    std::variant<ROOT::RNTupleView<std::string>,
                 ROOT::RNTupleView<std::int32_t>,
                 ROOT::RNTupleView<std::uint32_t>, ROOT::RNTupleView<double>>;

enum FieldTypes { String, Int32, Double };

std::vector<std::pair<std::string, enum FieldTypes>>
GetFieldNamesAndTypes(const ROOT::REntry &entry) {
  std::vector<std::pair<std::string, enum FieldTypes>> fields;

  for (const auto &val : entry) {
    std::string typeName = val.GetField().GetTypeName();
    if (typeName == "std::string") {
      fields.emplace_back(val.GetField().GetFieldName(), FieldTypes::String);
    } else if (typeName == "std::int32_t") {
      fields.emplace_back(val.GetField().GetFieldName(), FieldTypes::Int32);
    } else if (typeName == "double") {
      fields.emplace_back(val.GetField().GetFieldName(), FieldTypes::Double);
    } else {
      throw std::runtime_error("Unsupported fieldtype: " + typeName);
    }
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

int main(int argc, char **argv) {
  if (argc < 4) {
    std::cerr << "Usage: " << argv[0]
              << " <input.ntuple.root> <RNTuple name> <scale>" << std::endl;
    return 1;
  }
  const char *kNTupleFileName = argv[1];
  const char *kNTupleName = argv[2];

  auto reader = ROOT::RNTupleReader::Open(kNTupleName, kNTupleFileName);
  auto fields = GetFieldNamesAndTypes(reader->GetModel().GetDefaultEntry());
  reader->GetModel();

  std::vector<std::pair<std::string, int>> fieldMap;
  std::vector<FVec> fieldsVec;
  std::vector<VVec> viewVec;
  std::cout << "Initializing the view arrays…" << std::endl;

  for (const auto &[fieldName, fieldType] : fields) {
    switch (fieldType) {
    case String:
      fieldsVec.emplace_back(std::in_place_type<std::vector<std::string>>,
                             reader->GetNEntries());
      viewVec.emplace_back(reader->GetView<std::string>(fieldName));
      fieldMap.emplace_back(fieldName, FieldTypes::String);
      break;
    case Int32:
      fieldsVec.emplace_back(std::in_place_type<std::vector<std::int32_t>>,
                             reader->GetNEntries());
      viewVec.emplace_back(reader->GetView<std::int32_t>(fieldName));
      fieldMap.emplace_back(fieldName, FieldTypes::Int32);
    case Double:
      fieldsVec.emplace_back(std::in_place_type<std::vector<double>>,
                             reader->GetNEntries());
      viewVec.emplace_back(reader->GetView<double>(fieldName));
      fieldMap.emplace_back(fieldName, FieldTypes::Double);
      break;
    default:
      throw std::runtime_error("Found an unsupported fieldtype.");
    }
  }

  std::cout << "Initializing the writer…" << std::endl;
  auto model = ROOT::RNTupleModel::Create();

  for (const auto &[fieldName, fieldType] : fields) {
    switch (fieldType) {
    case String:
      model->MakeField<std::string>(fieldName);
      break;
    case Int32:
      model->MakeField<std::int32_t>(fieldName);
      break;
    case Double:
      model->MakeField<double>(fieldName);
      break;
    default:
      throw std::runtime_error("Found an unsupported fieldtype.");
    }
  }
  auto options = ROOT::RNTupleWriteOptions();
  auto writer = ROOT::RNTupleWriter::Recreate(std::move(model), kNTupleName,
                                              kNTupleFileName, options);

  // auto writer = ROOT::RNTupleWriter::Recreate(std::move(reader->GetModel()),
  // kNTupleName, "test.root.rntuple");

  // for (size_t field = 0; field < fieldsVec.size(); ++field) {
  //   const auto &[fieldName, fieldType] = fields[field];
  //   std::visit(
  //       [&](auto &dstVec, auto &srcView) {
  //         using D = std::decay_t<decltype(dstVec)>;
  //         using Elem =
  //             std::remove_const_t<std::remove_reference_t<decltype(srcView(
  //                 (ROOT::NTupleSize_t)0))>>;
  //         if constexpr (std::is_same_v<typename D::value_type, Elem>) {
  //           for (size_t i = 0; i < reader->GetNEntries(); ++i) {
  //             dstVec[i] = srcView(static_cast<ROOT::NTupleSize_t>(i));
  //           }
  //         }
  //       },
  //       fieldsVec[field], viewVec[field]);
  // }
}
