#include <ROOT/RNTupleReader.hxx>

#include <ROOT/RNTupleUtil.hxx>
#include <ROOT/RNTupleView.hxx>

#include <TH1D.h>
#include <TString.h>

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <list>
#include <orc/OrcFile.hh>
#include <orc/Reader.hh>
#include <orc/Type.hh>
#include <orc/Vector.hh>
#include <orc/Writer.hh>

#include <ostream>
#include <string>

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <tuple>
#include <unistd.h>

#include <cassert>
#include <chrono>
#include <cstdio>
#include <iostream>
#include <memory>
#include <string>

#include <Compression.h>
#include <ROOT/RNTupleReadOptions.hxx>
#include <ROOT/RNTupleReader.hxx>
#include <TApplication.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TClassTable.h>
#include <TFile.h>
#include <TH1D.h>
#include <TROOT.h>
#include <TRootCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreePerfStats.h>
#include <TTreeReader.h>
#include <vector>

enum selected_indices {
  H1_PX,
  H1_PY,
  H1_PZ,
  H1_ProbK,
  H1_ProbPi,
  H1_isMuon,
  H2_PX,
  H2_PY,
  H2_PZ,
  H2_ProbK,
  H2_ProbPi,
  H2_isMuon,
  H3_PX,
  H3_PY,
  H3_PZ,
  H3_ProbK,
  H3_ProbPi,
  H3_isMuon,
};

constexpr double kKaonMassMeV = 493.677;

void clear_cache() {
  sync();
  char signal = '1';
  int fd = open("/proc/sys/vm/drop_caches", O_WRONLY);
  if (fd < 0) {
    perror("open(/proc/sys/vm/drop_caches)");
    return;
  }
  if (write(fd, &signal, 1) != 1) {
    perror("write(/proc/sys/vm/drop_caches)");
  }
  close(fd);
}

static double GetP2(double px, double py, double pz) {
  return px * px + py * py + pz * pz;
}

static double GetKE(double px, double py, double pz) {
  double p2 = GetP2(px, py, pz);
  return sqrt(p2 + kKaonMassMeV * kKaonMassMeV);
}

static void Show(TH1D *h) {
  auto app = TApplication("", nullptr, nullptr);

  gStyle->SetTextFont(42);
  auto c = TCanvas("c", "", 800, 700);
  h->GetXaxis()->SetTitle("m_{KKK} [MeV/c^{2}]");
  h->DrawCopy();
  c.Modified();
  c.Update();
  app.Run();
}

// static void NTupleDirect(const std::string &path) {
// static std::pair<uint64_t, uint64_t> NTupleDirect(const std::string &path) {
static std::tuple<uint64_t, uint64_t, double>
NTupleDirect(const std::string &path) {
  using RNTupleReader = ROOT::RNTupleReader;
  using RNTupleModel = ROOT::RNTupleModel;

  auto ts_init = std::chrono::steady_clock::now();

  auto model = RNTupleModel::Create();
  // auto options = GetRNTupleOptions();
  ROOT::RNTupleReadOptions options;
  auto ntuple =
      RNTupleReader::Open(std::move(model), "DecayTree", path, options);
  // ntuple->EnableMetrics();

  auto viewH1IsMuon = ntuple->GetView<int>("H1_isMuon");
  auto viewH2IsMuon = ntuple->GetView<int>("H2_isMuon");
  auto viewH3IsMuon = ntuple->GetView<int>("H3_isMuon");

  auto viewH1PX = ntuple->GetView<double>("H1_PX");
  auto viewH1PY = ntuple->GetView<double>("H1_PY");
  auto viewH1PZ = ntuple->GetView<double>("H1_PZ");
  auto viewH1ProbK = ntuple->GetView<double>("H1_ProbK");
  auto viewH1ProbPi = ntuple->GetView<double>("H1_ProbPi");

  auto viewH2PX = ntuple->GetView<double>("H2_PX");
  auto viewH2PY = ntuple->GetView<double>("H2_PY");
  auto viewH2PZ = ntuple->GetView<double>("H2_PZ");
  auto viewH2ProbK = ntuple->GetView<double>("H2_ProbK");
  auto viewH2ProbPi = ntuple->GetView<double>("H2_ProbPi");

  auto viewH3PX = ntuple->GetView<double>("H3_PX");
  auto viewH3PY = ntuple->GetView<double>("H3_PY");
  auto viewH3PZ = ntuple->GetView<double>("H3_PZ");
  auto viewH3ProbK = ntuple->GetView<double>("H3_ProbK");
  auto viewH3ProbPi = ntuple->GetView<double>("H3_ProbPi");

  auto hMass = new TH1D("B_mass", "", 500, 5050, 5500);

  unsigned nevents = 0;
  std::chrono::steady_clock::time_point ts_first =
      std::chrono::steady_clock::now();
  for (auto i : ntuple->GetEntryRange()) {

    if (viewH1IsMuon(i) || viewH2IsMuon(i) || viewH3IsMuon(i)) {
      continue;
    }

    constexpr double prob_k_cut = 0.5;
    if (viewH1ProbK(i) < prob_k_cut)
      continue;
    if (viewH2ProbK(i) < prob_k_cut)
      continue;
    if (viewH3ProbK(i) < prob_k_cut)
      continue;

    constexpr double prob_pi_cut = 0.5;
    if (viewH1ProbPi(i) > prob_pi_cut)
      continue;
    if (viewH2ProbPi(i) > prob_pi_cut)
      continue;
    if (viewH3ProbPi(i) > prob_pi_cut)
      continue;

    double b_px = viewH1PX(i) + viewH2PX(i) + viewH3PX(i);
    double b_py = viewH1PY(i) + viewH2PY(i) + viewH3PY(i);
    double b_pz = viewH1PZ(i) + viewH2PZ(i) + viewH3PZ(i);
    double b_p2 = GetP2(b_px, b_py, b_pz);
    double k1_E = GetKE(viewH1PX(i), viewH1PY(i), viewH1PZ(i));
    double k2_E = GetKE(viewH2PX(i), viewH2PY(i), viewH2PZ(i));
    double k3_E = GetKE(viewH3PX(i), viewH3PY(i), viewH3PZ(i));
    double b_E = k1_E + k2_E + k3_E;
    double b_mass = sqrt(b_E * b_E - b_p2);
    hMass->Fill(b_mass);
  }
  auto ts_end = std::chrono::steady_clock::now();
  auto runtime_init =
      std::chrono::duration_cast<std::chrono::microseconds>(ts_first - ts_init)
          .count();
  auto runtime_analyze =
      std::chrono::duration_cast<std::chrono::microseconds>(ts_end - ts_first)
          .count();

  // ntuple->PrintInfo(ROOT::ENTupleInfo::kMetrics);
  Show(hMass);

  // Total events divided by total analysis time. Time is converted from
  // microseconds to seconds.
  double eventsPerSec = static_cast<double>(ntuple->GetNEntries()) /
                        (runtime_analyze / 1'000'000.0);

  delete hMass;
  std::tuple<uint64_t, uint64_t, double> times = {
      static_cast<uint64_t>(runtime_init),
      static_cast<uint64_t>(runtime_analyze), eventsPerSec};

  return times;
}

// static void ORCDirect(const std::string &path) {
static std::tuple<uint64_t, uint64_t, double>
ORCDirect(const std::string &path) {
  auto ts_init = std::chrono::steady_clock::now();

  std::unique_ptr<orc::InputStream> inStream = orc::readLocalFile(path);
  orc::ReaderOptions options;
  orc::ReaderMetrics metrics;
  // Should be the ORC equivalent of ntuple->EnableMetrics();
  // options.setReaderMetrics(&metrics);
  std::unique_ptr<orc::Reader> reader =
      orc::createReader(std::move(inStream), options);

  orc::RowReaderOptions rowReaderOptions;
  std::list<std::string> includedFieldsList = {
      "H1_isMuon", "H2_isMuon", "H3_isMuon", "H1_PX", "H1_PY",    "H1_PZ",
      "H1_ProbK",  "H1_ProbPi", "H2_PX",     "H2_PY", "H2_PZ",    "H2_ProbK",
      "H2_ProbPi", "H3_PX",     "H3_PY",     "H3_PZ", "H3_ProbK", "H3_ProbPi"};
  rowReaderOptions.include(includedFieldsList);

  std::unique_ptr<orc::RowReader> rowReader =
      reader->createRowReader(rowReaderOptions);
  // RNTuple reads a cluster, which is generally bigger than 1024. Trying to
  // read ORCs equivalent, a stripe, yielded slower times however.
  std::unique_ptr<orc::ColumnVectorBatch> batch =
      rowReader->createRowBatch(1024);
  // std::unique_ptr<orc::ColumnVectorBatch> batch =
  //     rowReader->createRowBatch(reader->getRowIndexStride());

  auto hMass = new TH1D("B_mass", "", 500, 5050, 5500);

  auto *root = dynamic_cast<orc::StructVectorBatch *>(batch.get());
  unsigned nevents = 0;
  auto *H1_isMuon = dynamic_cast<orc::LongVectorBatch *>(
      root->fields[selected_indices::H1_isMuon]);
  auto *H2_isMuon = dynamic_cast<orc::LongVectorBatch *>(
      root->fields[selected_indices::H2_isMuon]);
  auto *H3_isMuon = dynamic_cast<orc::LongVectorBatch *>(
      root->fields[selected_indices::H3_isMuon]);

  auto *H1_ProbK = dynamic_cast<orc::DoubleVectorBatch *>(
      root->fields[selected_indices::H1_ProbK]);
  auto *H2_ProbK = dynamic_cast<orc::DoubleVectorBatch *>(
      root->fields[selected_indices::H2_ProbK]);
  auto *H3_ProbK = dynamic_cast<orc::DoubleVectorBatch *>(
      root->fields[selected_indices::H3_ProbK]);
  auto *H1_ProbPi = dynamic_cast<orc::DoubleVectorBatch *>(
      root->fields[selected_indices::H1_ProbPi]);
  auto *H2_ProbPi = dynamic_cast<orc::DoubleVectorBatch *>(
      root->fields[selected_indices::H2_ProbPi]);
  auto *H3_ProbPi = dynamic_cast<orc::DoubleVectorBatch *>(
      root->fields[selected_indices::H3_ProbPi]);

  auto *H1_PX = dynamic_cast<orc::DoubleVectorBatch *>(
      root->fields[selected_indices::H1_PX]);
  auto *H1_PY = dynamic_cast<orc::DoubleVectorBatch *>(
      root->fields[selected_indices::H1_PY]);
  auto *H1_PZ = dynamic_cast<orc::DoubleVectorBatch *>(
      root->fields[selected_indices::H1_PZ]);

  auto *H2_PX = dynamic_cast<orc::DoubleVectorBatch *>(
      root->fields[selected_indices::H2_PX]);
  auto *H2_PY = dynamic_cast<orc::DoubleVectorBatch *>(
      root->fields[selected_indices::H2_PY]);
  auto *H2_PZ = dynamic_cast<orc::DoubleVectorBatch *>(
      root->fields[selected_indices::H2_PZ]);

  auto *H3_PX = dynamic_cast<orc::DoubleVectorBatch *>(
      root->fields[selected_indices::H3_PX]);
  auto *H3_PY = dynamic_cast<orc::DoubleVectorBatch *>(
      root->fields[selected_indices::H3_PY]);
  auto *H3_PZ = dynamic_cast<orc::DoubleVectorBatch *>(
      root->fields[selected_indices::H3_PZ]);

  std::chrono::steady_clock::time_point ts_first =
      std::chrono::steady_clock::now();

  while (rowReader->next(*batch)) {
    for (uint64_t r = 0; r < batch->numElements; ++r) {
      if (H1_isMuon->data[r] || H2_isMuon->data[r] || H3_isMuon->data[r]) {
        continue;
      }
      constexpr double prob_k_cut = 0.5;
      if (H1_ProbK->data[r] < prob_k_cut)
        continue;
      if (H2_ProbK->data[r] < prob_k_cut)
        continue;
      if (H3_ProbK->data[r] < prob_k_cut)
        continue;

      constexpr double prob_pi_cut = 0.5;
      if (H1_ProbPi->data[r] > prob_pi_cut)
        continue;
      if (H2_ProbPi->data[r] > prob_pi_cut)
        continue;
      if (H3_ProbPi->data[r] > prob_pi_cut)
        continue;

      double b_px = H1_PX->data[r] + H2_PX->data[r] + H3_PX->data[r];
      double b_py = H1_PY->data[r] + H2_PY->data[r] + H3_PY->data[r];
      double b_pz = H1_PZ->data[r] + H2_PZ->data[r] + H3_PZ->data[r];
      double b_p2 = GetP2(b_px, b_py, b_pz);
      double k1_E = GetKE(H1_PX->data[r], H1_PY->data[r], H1_PZ->data[r]);
      double k2_E = GetKE(H2_PX->data[r], H2_PY->data[r], H2_PZ->data[r]);
      double k3_E = GetKE(H3_PX->data[r], H3_PY->data[r], H3_PZ->data[r]);
      double b_E = k1_E + k2_E + k3_E;
      double b_mass = sqrt(b_E * b_E - b_p2);
      hMass->Fill(b_mass);
    }
  }
  auto ts_end = std::chrono::steady_clock::now();
  auto runtime_init =
      std::chrono::duration_cast<std::chrono::microseconds>(ts_first - ts_init)
          .count();
  auto runtime_analyze =
      std::chrono::duration_cast<std::chrono::microseconds>(ts_end - ts_first)
          .count();

  // ntuple->PrintInfo(ROOT::ENTupleInfo::kMetrics);
  // Show(hMass);

  // Total events divided by total analysis time. Time is converted from
  // microseconds to seconds.
  double eventsPerSec = static_cast<double>(reader->getNumberOfRows()) /
                        (runtime_analyze / 1'000'000.0);

  delete hMass;
  std::tuple<uint64_t, uint64_t, double> times = {
      static_cast<uint64_t>(runtime_init),
      static_cast<uint64_t>(runtime_analyze), eventsPerSec};
  return times;
}

std::tuple<double, double, double>
mean(const std::vector<std::tuple<uint64_t, uint64_t, double>> &timeVec,
     uint64_t runs) {
  double initMean = 0.0, analysisMean = 0.0, throughputMean = 0.0;
  for (uint64_t i = 0; i < runs; ++i) {
    const auto &[init, analysis, throughput] = timeVec[i];
    initMean += static_cast<double>(init);
    analysisMean += static_cast<double>(analysis);
    throughputMean += throughput;
  }
  initMean /= runs;
  analysisMean /= runs;
  throughputMean /= runs;
  return std::tuple<double, double, double>(initMean, analysisMean,
                                            throughputMean);
}

std::tuple<double, double, double>
STD(const std::vector<std::tuple<uint64_t, uint64_t, double>> &timeVec,
    uint64_t runs, double initMean, double analysisMean,
    double throughputMean) {
  double initSTD = 0.0, analysisSTD = 0.0, throughputSTD = 0.0;
  for (uint64_t i = 0; i < runs; ++i) {
    const auto &[init, analysis, throughput] = timeVec[i];
    initSTD += pow(static_cast<double>(init) - initMean, 2.0);
    analysisSTD += pow(static_cast<double>(analysis) - analysisMean, 2.0);
    throughputSTD += pow(throughput - throughputMean, 2.0);
  }
  initSTD /= runs;
  initSTD = sqrt(initSTD);
  analysisSTD /= runs;
  analysisSTD = sqrt(analysisSTD);
  throughputSTD /= runs;
  throughputSTD = sqrt(throughputSTD);
  return std::tuple<double, double, double>(initSTD, analysisSTD,
                                            throughputSTD);
}

std::vector<std::tuple<uint64_t, uint64_t, double>>
runRNTupleBenchmarks(const std::string &fn, uint64_t runs, bool hotCacheRun) {
  std::vector<std::tuple<uint64_t, uint64_t, double>> timeVec;
  // Even with cold cache it needs a warmup.
  NTupleDirect(fn);
  if (!hotCacheRun) {
    clear_cache();
  }

  for (uint64_t i = 0; i < runs; ++i) {
    if (!hotCacheRun) {
      clear_cache();
    }
    auto result = NTupleDirect(fn);
    timeVec.emplace_back(result);
  }

  return timeVec;
}

std::vector<std::tuple<uint64_t, uint64_t, double>>
runOrcBenchmarks(const std::string &fn, uint64_t runs, bool hotCacheRun) {
  std::vector<std::tuple<uint64_t, uint64_t, double>> timeVec;
  // Even with cold cache it needs a warmup.
  ORCDirect(fn);
  if (!hotCacheRun) {
    clear_cache();
  }

  for (uint64_t i = 0; i < runs; ++i) {
    if (!hotCacheRun) {
      clear_cache();
    }
    auto result = ORCDirect(fn);
    timeVec.emplace_back(result);
  }

  return timeVec;
}

std::string joinJson(const std::vector<uint64_t> &v) {
  std::ostringstream oss;
  oss << '[';
  for (size_t i = 0; i < v.size(); ++i) {
    if (i)
      oss << ',';
    oss << v[i];
  }
  oss << ']';
  return oss.str();
}

std::string joinJson(const std::vector<double> &v) {
  std::ostringstream oss;
  oss << '[';
  for (size_t i = 0; i < v.size(); ++i) {
    if (i)
      oss << ',';
    oss << v[i];
  }
  oss << ']';
  return oss.str();
}

// For the scatterplot.
void writeTimeVecCSV(
    std::ofstream &out, const std::string &fn,
    const std::vector<std::tuple<uint64_t, uint64_t, double>> &timeVec) {
  std::vector<uint64_t> initTimes, analysisTimes;
  std::vector<double> throughputs;
  initTimes.reserve(timeVec.size());
  analysisTimes.reserve(timeVec.size());
  throughputs.reserve(timeVec.size());

  for (const auto &t : timeVec) {
    initTimes.push_back(std::get<0>(t));
    analysisTimes.push_back(std::get<1>(t));
    throughputs.push_back(std::get<2>(t));
  }

  out << '"' << joinJson(initTimes) << '"' << ',' << '"'
      << joinJson(analysisTimes) << '"' << ',' << '"' << joinJson(throughputs)
      << '"' << "\n";
}

void writeTimeStatsCSV(
    std::ofstream &out, const uint64_t runs,
    const std::vector<std::tuple<uint64_t, uint64_t, double>> &timeVec) {
  const auto &[initMean, analysisMean, throughputMean] = mean(timeVec, runs);
  const auto &[initSTD, analysisSTD, throughputSTD] =
      STD(timeVec, runs, initMean, analysisMean, throughputMean);

  std::string initMeanSTR = std::to_string(initMean);
  std::string analysisMeanSTR = std::to_string(analysisMean);
  std::string throughputMeanSTR = std::to_string(throughputMean);
  std::string initSTDSTR = std::to_string(initSTD);
  std::string analysisSTDSTR = std::to_string(analysisSTD);
  std::string throughputSTDSTR = std::to_string(throughputSTD);

  out << initMeanSTR << "," << initSTDSTR << "," << analysisMeanSTR << ","
      << analysisSTDSTR << "," << throughputMeanSTR << "," << throughputSTDSTR
      << "\n";
}

void writeHeadersCSV(std::ofstream &out, const bool meanRun) {
  out << "Filename,Number of Runs,Data Format,Cache,";
  if (meanRun) {
    out << "Initialization Mean,Initialization STD,Analysis Mean,Analysis "
           "STD,Events/sec Mean,Events/sec STD\n";
  } else {
    out << "Initialization times,Analysis Times, Events/sec\n";
  }
}

void writeArgumentsCSV(std::ofstream &out, std::string fn, const uint64_t runs,
                       const bool orcRun, const bool hotCacheRun,
                       const bool meanRun) {
  namespace fs = std::filesystem;
  fs::path p(fn);
  // This is so the relative path to the file is not included in the .csv.
  fn = p.filename().string();

  std::string runsSTR = std::to_string(runs);
  std::string dataFormat = (orcRun ? "ORC" : "RNTuple");
  std::string cacheSetting = (hotCacheRun ? "Hot" : "Cold");
  out << fn << "," << runsSTR << "," << dataFormat << "," << cacheSetting
      << ",";
}

void writeCSV(
    const std::string &fn, const uint64_t runs, const bool orcRun,
    const bool hotCacheRun, const bool meanRun,
    const std::vector<std::tuple<uint64_t, uint64_t, double>> &timeVec) {
  std::string outputFN = "output.csv";
  if (hotCacheRun) {
    outputFN = "hot_" + outputFN;
  } else {
    outputFN = "cold_" + outputFN;
  }
  if (meanRun) {
    outputFN = "stats_" + outputFN;
  } else {
    outputFN = "scatter_" + outputFN;
  }
  std::ofstream out(outputFN, std::ios::app);

  namespace fs = std::filesystem;
  if (!fs::exists(outputFN) || fs::file_size(outputFN) == 0) {
    writeHeadersCSV(out, meanRun);
  }
  writeArgumentsCSV(out, fn, runs, orcRun, hotCacheRun, meanRun);
  if (meanRun) {
    writeTimeStatsCSV(out, runs, timeVec);
  } else {
    writeTimeVecCSV(out, fn, timeVec);
  }
}

int main(int argc, char **argv) {
  if (argc != 6 || (std::atoi(argv[3]) != 0 && std::atoi(argv[3]) != 1) ||
      (std::atoi(argv[4]) != 0 && std::atoi(argv[4]) != 1) ||
      (std::atoi(argv[5]) != 0 && std::atoi(argv[5]) != 1)) {
    std::cerr << "Usage: " << argv[0]
              << " <filename> <# runs> <0: RNTuple, 1: ORC> <0: cold cache, 1: "
                 "hot cache> "
                 "<0: Scatter, 1: Mean and std>"
              << std::endl;
    return 1;
  }
  // The data should always be DecayTree, or something with the same column
  // names and types, but the file name can change a little if scaled versions
  // are used.
  const std::string fn = std::string(argv[1]);
  const uint64_t runs = std::atoi(argv[2]);
  const bool orcRun = static_cast<bool>(std::atoi(argv[3]));
  const bool hotCacheRun = static_cast<bool>(std::atoi(argv[4]));
  const bool meanRun = static_cast<bool>(std::atoi(argv[5]));

  auto timeVec = orcRun ? runOrcBenchmarks(fn, runs, hotCacheRun)
                        : runRNTupleBenchmarks(fn, runs, hotCacheRun);
  writeCSV(fn, runs, orcRun, hotCacheRun, meanRun, timeVec);

  return 0;
}
