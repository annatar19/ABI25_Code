#include <ROOT/RNTuple.hxx>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleWriteOptions.hxx>
#include <ROOT/RNTupleWriter.hxx>
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
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

using FVec = std::variant<std::vector<std::string>, std::vector<std::int32_t>,
                          std::vector<double>>;
using VVec =
    std::variant<ROOT::RNTupleView<std::string>,
                 ROOT::RNTupleView<std::int32_t>, ROOT::RNTupleView<double>>;

// The results of model->MakeField.
using OVec = std::variant<std::shared_ptr<double>, std::shared_ptr<int>,
                          std::shared_ptr<std::string>>;

enum FieldTypes { String, Int32, Double };

template <class... Ts> struct overloaded : Ts... {
  using Ts::operator()...;
};
template <class... Ts> overloaded(Ts...) -> overloaded<Ts...>;

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

  // std::sort(fields.begin(), fields.end(),
  //           [](auto const &a, auto const &b) { return a.first < b.first; });

  return fields;
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
  // reader->GetModel();

  std::vector<std::pair<std::string, int>> fieldMap;
  std::vector<FVec> inputFieldsVec;
  std::vector<OVec> outputFieldsVec;
  std::vector<VVec> inputViewVec;
  std::cout << "Initializing the view arrays…" << std::endl;

  for (const auto &[fieldName, fieldType] : fields) {
    switch (fieldType) {
    case String:
      inputFieldsVec.emplace_back(std::in_place_type<std::vector<std::string>>,
                                  reader->GetNEntries());
      inputViewVec.emplace_back(reader->GetView<std::string>(fieldName));
      fieldMap.emplace_back(fieldName, FieldTypes::String);
      break;
    case Int32:
      inputFieldsVec.emplace_back(std::in_place_type<std::vector<std::int32_t>>,
                                  reader->GetNEntries());
      inputViewVec.emplace_back(reader->GetView<std::int32_t>(fieldName));
      fieldMap.emplace_back(fieldName, FieldTypes::Int32);
      break;
    case Double:
      inputFieldsVec.emplace_back(std::in_place_type<std::vector<double>>,
                                  reader->GetNEntries());
      inputViewVec.emplace_back(reader->GetView<double>(fieldName));
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
      outputFieldsVec.emplace_back(model->MakeField<std::string>(fieldName));
      break;
    case Int32:
      outputFieldsVec.emplace_back(model->MakeField<std::int32_t>(fieldName));
      break;
    case Double:
      outputFieldsVec.emplace_back(model->MakeField<double>(fieldName));
      break;
    default:
      throw std::runtime_error("Found an unsupported fieldtype.");
    }
  }
  auto options = ROOT::RNTupleWriteOptions();
  options.SetCompression(ROOT::RCompressionSetting::EAlgorithm::EValues::kZSTD,
                         5);

  std::string new_name = kNTupleFileName + std::string(".out");
  auto writer = ROOT::RNTupleWriter::Recreate(std::move(model), kNTupleName,
                                              "B2HHH.ntuple.root", options);

  for (uint64_t i = 0; i < reader->GetNEntries(); ++i) {
    for (size_t field = 0; field < inputFieldsVec.size(); ++field) {
      std::visit(
          overloaded{
                     [&i](std::shared_ptr<std::string> const &out,
                         ROOT::RNTupleView<std::string> &in) { *out = in(i); },
                     [&i](std::shared_ptr<std::int32_t> const &out,
                         ROOT::RNTupleView<std::int32_t> &in) { *out = in(i); },
                     [&i](std::shared_ptr<double> const &out,
                         ROOT::RNTupleView<double> &in) { *out = in(i); },
                     [&](auto const &a, auto &b) {
                       throw std::runtime_error("The data got corrupted!");
                     }},
          outputFieldsVec[field], inputViewVec[field]);
    }
    writer->Fill();
  }
}
