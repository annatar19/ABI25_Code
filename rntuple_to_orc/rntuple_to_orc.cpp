#include <ROOT/RNTupleReader.hxx>

#include <ROOT/RNTupleUtil.hxx>
#include <ROOT/RNTupleView.hxx>

#include <TString.h>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <orc/OrcFile.hh>
#include <orc/Type.hh>
#include <orc/Vector.hh>
#include <orc/Writer.hh>

#include <string>
#include <variant>

using FVec = std::variant<std::vector<std::string>, std::vector<std::int32_t>,
                          std::vector<std::uint32_t>, std::vector<double>>;
using VVec =
    std::variant<ROOT::RNTupleView<std::string>,
                 ROOT::RNTupleView<std::int32_t>,
                 ROOT::RNTupleView<std::uint32_t>, ROOT::RNTupleView<double>>;
using BVec = std::variant<orc::LongVectorBatch *, orc::DoubleVectorBatch *,
                          orc::StringVectorBatch *>;

std::vector<std::pair<std::string, std::string>>
GetFieldNamesAndTypes(const ROOT::REntry &entry) {
  std::vector<std::pair<std::string, std::string>> fields;

  for (const auto &val : entry) {
    fields.emplace_back(val.GetField().GetFieldName(),
                        val.GetField().GetTypeName());
  }

  return fields;
}

// See https://medium.com/@nerudaj/std-visit-is-awesome-heres-why-f183f6437932
// for why this is needed.
template <class... Ts> struct overloaded : Ts... {
  using Ts::operator()...;
};
template <class... Ts> overloaded(Ts...) -> overloaded<Ts...>;

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

  std::cout << "Starting conversion to Apache ORC…" << std::endl;

  std::unique_ptr<orc::Type> schema;
  orc::WriterOptions options;
  std::unique_ptr<orc::Writer> writer;
  std::unique_ptr<orc::OutputStream> outStream =
      orc::writeLocalFile(std::string(kNTupleName) + std::string(".orc"));
  // If the batchSize is reached it is written to the file. The exact value was
  // copied from the example code, it is likely not optimal.
  uint64_t batchSize = 1024;

  // Used to keep track of where in the current batch the loop is;
  uint64_t rows = 0;
  // Used later for itterating over fields.
  uint64_t field_count;

  // These 2 will always exist, regardless of the fields used, so they don't
  // need to be but inside a std::variant.
  std::unique_ptr<orc::ColumnVectorBatch> batch;
  orc::StructVectorBatch *root;

  std::vector<BVec> batches;

  std::cout << "Initializing ORC structures…\n";

  // Orc can use a string for its 'model'.
  std::string schema_str = "struct<";
  bool first_field = true;
  for (const auto &[fieldName, fieldType] : fieldMap) {
    if (!first_field) {
      schema_str += ",";
    }
    first_field = false;
    switch (fieldType) {
    case FieldTypes::String: {
      schema_str += fieldName + ":string";
      break;
    }
    case FieldTypes::Int32: {
      schema_str += fieldName + ":int";
      break;
    }
    case FieldTypes::Uint32: {
      // Orc does not feature unsigned integers. TODO to check for overflow.
      schema_str += fieldName + ":int";
      break;
    }
    case FieldTypes::Double: {
      schema_str += fieldName + ":double";
      break;
    }
    }
  }
  schema_str += ">";
  schema = orc::Type::buildTypeFromString(schema_str);
  writer = orc::createWriter(*schema, outStream.get(), options);
  batch = writer->createRowBatch(batchSize);
  root = dynamic_cast<orc::StructVectorBatch *>(batch.get());
  // Need to loop a second time to initialize the batches.
  int i = 0;
  for (const auto &[fieldName, fieldType] : fieldMap) {
    switch (fieldType) {
    case FieldTypes::String: {
      batches.emplace_back(
          dynamic_cast<orc::StringVectorBatch *>(root->fields[i]));
      break;
    }
    case FieldTypes::Int32: {
      batches.emplace_back(
          dynamic_cast<orc::LongVectorBatch *>(root->fields[i]));

      break;
    }
    case FieldTypes::Uint32: {
      batches.emplace_back(
          dynamic_cast<orc::LongVectorBatch *>(root->fields[i]));

      break;
    }
    case FieldTypes::Double: {
      batches.emplace_back(
          dynamic_cast<orc::DoubleVectorBatch *>(root->fields[i]));

      break;
    }
    }
    ++i;
  }
  std::cout << "Writing ORC…\n";

  for (uint64_t cur_entry = 0; cur_entry < reader->GetNEntries();
       cur_entry += batchSize) {
    uint64_t remaining_rows = reader->GetNEntries() - cur_entry;
    uint64_t rows_to_add =
        remaining_rows >= batchSize ? batchSize : remaining_rows;
    for (uint64_t field = 0; field < fieldsVec.size(); ++field) {
      std::visit(
          overloaded{
              [&](orc::LongVectorBatch *b, std::vector<std::int32_t> &f) {
                for (uint64_t row = 0; row < rows_to_add; row++) {
                  b->data[row] = f[row + cur_entry];
                }
                b->numElements = rows_to_add;
              },
              [&](orc::LongVectorBatch *b, std::vector<std::uint32_t> &f) {
                for (uint64_t row = 0; row < rows_to_add; row++) {
                  b->data[row] = f[row + cur_entry];
                }
                b->numElements = rows_to_add;
              },
              [&](orc::DoubleVectorBatch *b, std::vector<double> &f) {
                for (uint64_t row = 0; row < rows_to_add; row++) {
                  b->data[row] = f[row + cur_entry];
                }
                b->numElements = rows_to_add;
              },
              [&](orc::StringVectorBatch *b, std::vector<std::string> &f) {
                for (uint64_t row = 0; row < rows_to_add; row++) {
                  const auto &s = f[row + cur_entry];
                  char *copy = strdup(s.c_str());
                  b->data[row] = copy;
                  b->length[row] = static_cast<uint64_t>(s.size());
                }
                b->numElements = rows_to_add;
              },
              // This should never happen, this is mostly here because the
              // compiler expects functions for all possible variant
              // combinations.
              [](auto a, auto b) {
                throw std::runtime_error("The data got corrupted!");
              },
          },
          batches[field], fieldsVec[field]);
    }
    root->numElements = rows_to_add;
    writer->add(*batch);
  }

  // Needed in the case of ORC, no standard flushing when going out of scope.
  writer->close();

  std::cout << "Conversion done!\n";

  return 0;
}
