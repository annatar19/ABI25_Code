#include <ROOT/RNTupleReader.hxx>

#include <ROOT/RNTupleUtil.hxx>
#include <ROOT/RNTupleView.hxx>

#include <TString.h>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <orc/Common.hh>
#include <orc/OrcFile.hh>
#include <orc/Type.hh>
#include <orc/Vector.hh>
#include <orc/Writer.hh>

#include <string>
#include <variant>

enum FieldTypes { String, Int32, Double };

using VVec =
    std::variant<ROOT::RNTupleView<std::string>,
                 ROOT::RNTupleView<std::int32_t>,
                 ROOT::RNTupleView<std::uint32_t>, ROOT::RNTupleView<double>>;

using BVec = std::variant<orc::LongVectorBatch *, orc::DoubleVectorBatch *,
                          orc::StringVectorBatch *>;
// using BVec =
//     std::variant<orc::IntegerVectorBatch<long int> *,
//                  orc::FloatingVectorBatch<double> *, orc::StringVectorBatch
//                  *>;

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

// See https://medium.com/@nerudaj/std-visit-is-awesome-heres-why-f183f6437932
// for why this is needed.
template <class... Ts> struct overloaded : Ts... {
  using Ts::operator()...;
};
template <class... Ts> overloaded(Ts...) -> overloaded<Ts...>;

int main(int argc, char **argv) {
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " <input.ntuple.root> <RNTuple name>"
              << std::endl;
    return 1;
  }

  const char *kNTupleFileName = argv[1];
  const char *kNTupleName = argv[2];

  std::cout << "Initializing RNTuple reader…" << std::endl;

  auto reader = ROOT::RNTupleReader::Open(kNTupleName, kNTupleFileName);
  auto fields = GetFieldNamesAndTypes(reader->GetModel().GetDefaultEntry());

  std::vector<VVec> viewVec;

  for (const auto &[fieldName, fieldType] : fields) {
    switch (fieldType) {
    case String:
      viewVec.emplace_back(reader->GetView<std::string>(fieldName));
      break;
    case Int32:
      viewVec.emplace_back(reader->GetView<std::int32_t>(fieldName));
      break;
    case Double:
      viewVec.emplace_back(reader->GetView<double>(fieldName));
      break;
    }
  }

  std::cout << "Initializing ORC Writer…" << std::endl;

  std::unique_ptr<orc::Type> schema;
  orc::WriterOptions options;
  // To set the level equal to RNTuple either that writer has to be changed,
  // the ORC source code has to be changed or the java convert tool has to be
  // used. Convert tool slow, but probably the best option.
  options.setCompression(orc::CompressionKind_ZSTD);
  std::unique_ptr<orc::Writer> writer;
  std::unique_ptr<orc::OutputStream> outStream =
      orc::writeLocalFile( // std::string("../output/") +
          std::string(kNTupleName) + std::string(".orc"));
  // If the batchSize is reached it is written to the file. The exact value
  // was copied from the example code, it is likely not optimal.
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

  // Orc can use a string for its 'model'.
  std::string schema_str = "struct<";
  bool first_field = true;
  for (const auto &[fieldName, fieldType] : fields) {
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
  for (const auto &[fieldName, fieldType] : fields) {
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
    for (uint64_t field = 0; field < fields.size(); ++field) {
      std::visit(
          overloaded{
              [&cur_entry, &rows_to_add](orc::LongVectorBatch *out,
                                         ROOT::RNTupleView<std::int32_t> &in) {
                std::cout << "long int!" << std::endl;
                for (uint64_t row = 0; row < rows_to_add; ++row) {
                  out->data[row] = in(row + cur_entry);
                }
                out->numElements = rows_to_add;
              },
              // orc::FloatingVectorBatch<double>
              [&cur_entry, &rows_to_add](orc::DoubleVectorBatch *out,
                                         ROOT::RNTupleView<double> &in) {
                std::cout << "double!" << std::endl;
                for (uint64_t row = 0; row < rows_to_add; row++) {
                  out->data[row] = in(row + cur_entry);
                }
                out->numElements = rows_to_add;
              },
              [&cur_entry, &rows_to_add](orc::StringVectorBatch *out,
                                         ROOT::RNTupleView<std::string> &in) {
                std::cout << "String!" << std::endl;
                for (uint64_t row = 0; row < rows_to_add; row++) {
                  const auto &s = in(row + cur_entry);
                  char *copy = strdup(s.c_str());
                  out->data[row] = copy;
                  out->length[row] = static_cast<uint64_t>(s.size());
                }
                out->numElements = rows_to_add;
              },
              // This should never happen, this is mostly here because
              // the compiler expects functions for all possible
              // variant combinations.
              [](auto *a, auto &b) {
                std::cout << "a is: " << typeid(a).name() << std::endl;
                std::cout << "b is: " << typeid(b).name() << std::endl;
                throw std::runtime_error("The data got corrupted!");
              }},
          batches[field], viewVec[field]);
    }
    root->numElements = rows_to_add;
    writer->add(*batch);
  }

  // Needed in the case of ORC, no standard flushing when going out of scope.
  writer->close();

  std::cout << "Conversion done!\n";

  return 0;
}
