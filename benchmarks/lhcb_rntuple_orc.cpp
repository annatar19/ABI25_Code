#include <ROOT/RNTupleReader.hxx>

#include <ROOT/RNTupleUtil.hxx>
#include <ROOT/RNTupleView.hxx>

#include <TH1D.h>
#include <TString.h>

#include <cstddef>
#include <cstdlib>
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
#include <unordered_map>
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
static std::pair<uint64_t, uint64_t> NTupleDirect(const std::string &path) {
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
    nevents++;
    if ((nevents % 100000) == 0) {
      printf("processed %u k events\n", nevents / 1000);
      // printf("dummy is %lf\n", dummy); abort();
    }

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

  std::cout << "Runtime-Initialization: " << runtime_init << "us" << std::endl;
  std::cout << "Runtime-Analysis: " << runtime_analyze << "us" << std::endl;

  // ntuple->PrintInfo(ROOT::ENTupleInfo::kMetrics);
  // Show(hMass);

  delete hMass;
  std::pair<uint64_t, uint64_t> times = {
      static_cast<uint64_t>(runtime_init),
      static_cast<uint64_t>(runtime_analyze)};
  return times;
}

// static void ORCDirect(const std::string &path) {
static std::pair<uint64_t, uint64_t> ORCDirect(const std::string &path) {

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

  std::cout << "Runtime-Initialization: " << runtime_init << "us" << std::endl;
  std::cout << "Runtime-Analysis: " << runtime_analyze << "us" << std::endl;

  // ntuple->PrintInfo(ROOT::ENTupleInfo::kMetrics);
  // Show(hMass);

  delete hMass;
  std::pair<uint64_t, uint64_t> times = {
      static_cast<uint64_t>(runtime_init),
      static_cast<uint64_t>(runtime_analyze)};
  return times;
}

std::pair<double, double>
mean(const std::vector<std::pair<uint64_t, uint64_t>> &timeVec, uint64_t runs) {
  double initMean = 0.0, analysisMean = 0.0;
  for (uint64_t i = 0; i < runs; ++i) {
    const auto &[init, analysis] = timeVec[i];
    initMean += static_cast<double>(init);
    analysisMean += static_cast<double>(analysis);
  }
  initMean /= runs;
  analysisMean /= runs;
  return std::pair<double, double>(initMean, analysisMean);
}

std::pair<double, double>
STD(const std::vector<std::pair<uint64_t, uint64_t>> &timeVec, uint64_t runs,
    double initMean, double analysisMean) {
  double initSTD = 0.0, analysisSTD = 0.0;
  for (uint64_t i = 0; i < runs; ++i) {
    const auto &[init, analysis] = timeVec[i];
    initSTD += pow(static_cast<double>(init) - initMean, 2.0);
    analysisSTD += pow(static_cast<double>(analysis) - analysisMean, 2.0);
  }
  initSTD /= runs;
  initSTD = sqrt(initSTD);
  analysisSTD /= runs;
  analysisSTD = sqrt(analysisSTD);
  return std::pair<double, double>(initSTD, analysisSTD);
}

int main(int argc, char **argv) {
  if (argc < 3 || (std::atoi(argv[2]) != 0 && std::atoi(argv[2]) != 1)) {
    std::cerr << "Usage: " << argv[0]
              << " <# runs> <0: cold cache, 1: hot cache>" << std::endl;
    return 1;
  }
  // -------------------------------------------------------------------------
  // RNTUPLE
  // -------------------------------------------------------------------------
  uint64_t runs = std::atoi(argv[1]);
  bool hot_cache = static_cast<bool>(std::atoi(argv[2]));
  // Even with cold cache it needs a warmup.
  NTupleDirect("B2HHH.ntuple.root");
  if (!hot_cache) {
    clear_cache();
  }
  std::vector<std::pair<uint64_t, uint64_t>> timeVec;
  std::cout << "Running RNTuple benchmarks…\n";
  for (uint64_t i = 0; i < runs; ++i) {
    std::cout << "Running itteration " << i + 1 << "…\n";
    if (!hot_cache) {
      clear_cache();
    }
    auto result = NTupleDirect("B2HHH.ntuple.root");
    timeVec.emplace_back(result);
    // i-1 to account for the skipped initial cold run.
    const auto &[init, analyze] = timeVec[i];
    std::cout << "\tItteration " << i + 1 << " initialization took: " << init
              << " us and analysis took: " << analyze << " us\n";
  }
  const auto &[initMean, analysisMean] = mean(timeVec, runs);
  const auto &[initSTD, analysisSTD] =
      STD(timeVec, runs, initMean, analysisMean);

  std::cout << "On average initialization took: " << initMean << " ± "
            << initSTD << " us\n";
  std::cout << "On average analysis took: " << analysisMean << " ± "
            << analysisSTD << " us\n";

  // -------------------------------------------------------------------------
  // ORC
  // -------------------------------------------------------------------------

  uint64_t runs = std::atoi(argv[1]);
  bool hot_cache = static_cast<bool>(std::atoi(argv[2]));
  // Even with cold cache it needs a warmup.
  ORCDirect("DecayTree.orc");
  if (!hot_cache) {
    clear_cache();
  }

  std::cout << "Running ORC benchmarks…\n";
  std::vector<std::pair<uint64_t, uint64_t>> timeVec;
  std::cout << "Running RNTuple benchmarks…\n";
  for (uint64_t i = 0; i < runs; ++i) {
    std::cout << "Running itteration " << i + 1 << "…\n";
    if (!hot_cache) {
      clear_cache();
    }
    auto result = ORCDirect("DecayTree.orc");
    timeVec.emplace_back(result);
    // i-1 to account for the skipped initial cold run.
    const auto &[init, analyze] = timeVec[i];
    std::cout << "\tItteration " << i + 1 << " initialization took: " << init
              << " us and analysis took: " << analyze << " us\n";
  }

  const auto &[initMean, analysisMean] = mean(timeVec, runs);
  const auto &[initSTD, analysisSTD] =
      STD(timeVec, runs, initMean, analysisMean);

  std::cout << "On average initialization took: " << initMean << " ± "
            << initSTD << " us\n";
  std::cout << "On average analysis took: " << analysisMean << " ± "
            << analysisSTD << " us\n";

  return 0;
}
