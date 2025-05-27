#include <ROOT/RNTuple.hxx>
#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleWriteOptions.hxx>
#include <ROOT/RNTupleWriter.hxx>
// #include <ROOT/RNTupleInspector.hxx>
#include <cmath>
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

enum FieldTypes { String, Int32, Double };

using VVec =
    std::variant<ROOT::RNTupleView<std::string>,
                 ROOT::RNTupleView<std::int32_t>, ROOT::RNTupleView<double>>;

// The results of model->MakeField.
using OVec = std::variant<std::shared_ptr<double>, std::shared_ptr<int>,
                          std::shared_ptr<std::string>>;

// Needed for the visit.
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

  return fields;
}

void initializeInputViews(
    ROOT::RNTupleReader &reader,
    const std::vector<std::pair<std::string, enum FieldTypes>> &fields,
    std::vector<VVec> &inputViewVec) {
  for (const auto &[fieldName, fieldType] : fields) {
    switch (fieldType) {
    case String:
      inputViewVec.emplace_back(reader.GetView<std::string>(fieldName));
      break;
    case Int32:
      inputViewVec.emplace_back(reader.GetView<std::int32_t>(fieldName));
      break;
    case Double:
      inputViewVec.emplace_back(reader.GetView<double>(fieldName));
      break;
    default:
      throw std::runtime_error("Found an unsupported fieldtype.");
    }
  }
}

void initializeOutputFields(
    ROOT::RNTupleModel &model,
    const std::vector<std::pair<std::string, enum FieldTypes>> &fields,
    std::vector<OVec> &outputFieldsVec) {
  for (const auto &[fieldName, fieldType] : fields) {
    switch (fieldType) {
    case String:
      outputFieldsVec.emplace_back(model.MakeField<std::string>(fieldName));
      break;
    case Int32:
      outputFieldsVec.emplace_back(model.MakeField<std::int32_t>(fieldName));
      break;
    case Double:
      outputFieldsVec.emplace_back(model.MakeField<double>(fieldName));
      break;
    default:
      throw std::runtime_error("Found an unsupported fieldtype.");
    }
  }
}

void scaleEntries(const std::unique_ptr<ROOT::RNTupleWriter> writer,
                  const std::vector<OVec> &outputFieldsVec,
                  std::vector<VVec> &inputViewVec,
                  const uint64_t &unscaled_entries,
                  const uint64_t &scaled_entries) {
  uint64_t i = 0;
  uint64_t curUpperLimit = 0;
  while (i < scaled_entries) {
    // Say unscaled_entries is 10000 and scaled_entries is 23450. This tenary
    // will yield 10000, 10000, and then finally 3450. Perhaps it can be done
    // brachless, but probably not worth the effort and the compiler might do
    // that in that case.
    curUpperLimit = (i + unscaled_entries <= scaled_entries)
                        ? unscaled_entries
                        : scaled_entries % unscaled_entries;
    for (uint64_t j = 0; j < curUpperLimit; ++j) {
      for (size_t field = 0; field < inputViewVec.size(); ++field) {
        std::visit(
            overloaded{
                [&j](const std::shared_ptr<std::string> &out,
                     ROOT::RNTupleView<std::string> &in) { *out = in(j); },
                [&j](const std::shared_ptr<std::int32_t> &out,
                     ROOT::RNTupleView<std::int32_t> &in) { *out = in(j); },
                [&j](const std::shared_ptr<double> &out,
                     ROOT::RNTupleView<double> &in) { *out = in(j); },
                [&](auto const &a, auto &b) {
                  throw std::runtime_error("The data got corrupted!");
                }},
            outputFieldsVec[field], inputViewVec[field]);
      }
      writer->Fill();
    }
    i += curUpperLimit;
  }
}
int main(int argc, char **argv) {
  if (argc < 4) {
    std::cerr << "Usage: " << argv[0]
              << " <input.ntuple.root> <RNTuple name> <scale>" << std::endl;
    return 1;
  }

  std::cout << "Initializing the reader…" << std::endl;

  const char *kNTupleFileName = argv[1];
  const char *kNTupleName = argv[2];
  std::unique_ptr<ROOT::RNTupleReader> reader =
      ROOT::RNTupleReader::Open(kNTupleName, kNTupleFileName);

  std::vector<std::pair<std::string, enum FieldTypes>> fields =
      GetFieldNamesAndTypes(reader->GetModel().GetDefaultEntry());
  std::vector<VVec> inputViewVec;

  initializeInputViews(*reader, fields, inputViewVec);

  std::cout << "Initializing the writer…" << std::endl;
  std::unique_ptr<ROOT::RNTupleModel> model = ROOT::RNTupleModel::Create();
  std::vector<OVec> outputFieldsVec;

  initializeOutputFields(*model, fields, outputFieldsVec);

  ROOT::RNTupleWriteOptions options = ROOT::RNTupleWriteOptions();
  options.SetCompression(ROOT::RCompressionSetting::EAlgorithm::EValues::kZSTD,
                         5);

  // The scaled file is outputted to the output dir with the scaling factor
  // prepended.
  std::string pathName =
      std::string("../output/") + argv[3] + std::string("x_") + kNTupleFileName;
  std::unique_ptr<ROOT::RNTupleWriter> writer = ROOT::RNTupleWriter::Recreate(
      std::move(model), kNTupleName, pathName, options);

  std::cout << "Starting the scaling…" << std::endl;
  const long double scale = std::stold(argv[3]);
  const uint64_t unscaled_entries = reader->GetNEntries();
  const uint64_t scaled_entries =
      static_cast<uint64_t>(std::round(scale * unscaled_entries));

  scaleEntries(std::move(writer), outputFieldsVec, inputViewVec,
               unscaled_entries, scaled_entries);

  std::cout << "Finished the scaling!" << std::endl;
}
