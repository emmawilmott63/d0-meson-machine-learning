// train_bdt_jetbins_feature_importance.C
// Six BDTs: inclusive and OS-only, each in three shared-jet-pT bins.
// Split by event number: event % 5 == 0 is test; all other events train.
// mass/isSignal/event/weight/xsWeight/sameJet are not input variables.
//
// In addition to the normal TMVA output, each model gets:
//   tmva/plots/<model>/variable_importance.png
//   A top-level TH1D named "VariableImportance" in its TMVA ROOT file.

#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TH1D.h"

#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <regex>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace {

struct JetBin { const char* tag; float low, high; };

const JetBin bins[] = {
    {"jetPt1_20", 1.f, 20.f},
    {"jetPt20_30", 20.f, 30.f},
    {"jetPt30plus", 30.f, 1.e30f}
};

const char* vars[] = {
    "pxK","pyK","pzK","pxPi","pyPi","pzPi",
    "ptK","ptPi","ptHard","ptSoft","ptSum","ptBalance",
    "etaK","etaPi","dEta","dPhi","dR","openAngle",
    "pairPt","pairEta","cosThetaStar",
    "jetPt","jetEta","jetMass","zJet","dRJet",
    "isoRel02","isoRel03","isoRel","isoRel06","isoNeutralRel",
    "nCharged","nJets","sphericityT","thrustT"
};

// Build a split-frequency ranking from the TMVA XML weights file.
// This deliberately avoids DecisionTree::GetVariableImportance(), which
// crashes in the user's ROOT 6.40.04 run. Each internal split contributes
// one count to the variable selected at that node; scores are normalized
// to sum to one. This is split-frequency importance, not split-gain.
void saveVariableImportance(
    const std::vector<std::string>& variableNames,
    const std::string& model,
    const std::string& outputDir,
    TFile* outputFile
) {
    // TMVA prefixes the weights XML with the full method name. Since the
    // booked method is named "BDT", the actual filename is <model>_BDT.weights.xml.
    const std::string weightsFile =
        "data/dataset_" + model + "/weights/" + model + "_BDT.weights.xml";
    std::ifstream xml(weightsFile);
    if (!xml) {
        std::cerr << "WARNING: cannot open TMVA weights file " << weightsFile
                  << "; skipping variable-importance plot.\n";
        return;
    }

    std::vector<Double_t> scores(variableNames.size(), 0.0);
    const std::regex selectorPattern{R"REGEX(IVar="(-?[0-9]+)")REGEX"};
    std::string line;
    unsigned int splitCount = 0;
    while (std::getline(xml, line)) {
        for (std::sregex_iterator it(line.begin(), line.end(), selectorPattern),
                                  last; it != last; ++it) {
            const int index = std::stoi((*it)[1].str());
            if (index >= 0 && static_cast<std::size_t>(index) < scores.size()) {
                scores[static_cast<std::size_t>(index)] += 1.0;
                ++splitCount;
            }
        }
    }

    if (splitCount == 0) {
        std::cerr << "WARNING: no internal-node variable selectors found in "
                  << weightsFile << "; skipping importance plot.\n";
        return;
    }

    for (Double_t& score : scores)
        score /= static_cast<Double_t>(splitCount);

    std::cout << "Variable importance: counted " << splitCount
              << " internal-node splits from TMVA XML (split frequency).\n";

    std::vector<std::pair<std::string, Double_t>> ranked;
    ranked.reserve(scores.size());
    for (std::size_t i = 0; i < scores.size(); ++i)
        ranked.emplace_back(variableNames[i], scores[i]);
    std::sort(ranked.begin(), ranked.end(),
        [](const auto& a, const auto& b) { return a.second < b.second; });

    const std::string plotDir = outputDir + "/plots/" + model;
    gSystem->mkdir(plotDir.c_str(), true);

    const int n = static_cast<int>(ranked.size());
    TH1D importanceHist(
        "VariableImportance",
        (model + ": BDT split frequency;Input variable;Fraction of splits").c_str(),
        n, 0.0, static_cast<Double_t>(n)
    );
    importanceHist.SetDirectory(nullptr);
    importanceHist.SetStats(false);
    importanceHist.SetFillColor(kAzure - 4);
    importanceHist.SetLineColor(kAzure + 2);
    importanceHist.SetLineWidth(1);
    importanceHist.SetBarWidth(0.8);
    importanceHist.SetBarOffset(0.1);

    Double_t maxScore = 0.0;
    for (int i = 0; i < n; ++i) {
        importanceHist.GetXaxis()->SetBinLabel(i + 1, ranked[i].first.c_str());
        importanceHist.SetBinContent(i + 1, ranked[i].second);
        maxScore = std::max(maxScore, ranked[i].second);
    }
    importanceHist.SetMinimum(0.0);
    importanceHist.SetMaximum(maxScore > 0.0 ? maxScore * 1.15 : 1.0);

    outputFile->cd();
    importanceHist.Write("VariableImportance", TObject::kOverwrite);

    const int canvasHeight = std::max(750, 24 * n + 120);
    TCanvas canvas(("c_variable_importance_" + model).c_str(),
                   "BDT variable split frequency", 1100, canvasHeight);
    canvas.SetLeftMargin(0.28);
    canvas.SetRightMargin(0.05);
    canvas.SetTopMargin(0.08);
    canvas.SetBottomMargin(0.10);
    importanceHist.GetXaxis()->SetLabelSize(0.028);
    importanceHist.GetYaxis()->SetLabelSize(0.035);
    importanceHist.Draw("hbar");
    canvas.SaveAs((plotDir + "/variable_importance.png").c_str());

    std::cout << "\nVariable split frequency for " << model << "\n"
              << "---------------------------------------------\n";
    for (int i = n - 1; i >= 0; --i)
        std::cout << ranked[i].first << " : " << ranked[i].second << "\n";
    std::cout << "Saved: " << plotDir << "/variable_importance.png\n";
}

template<class T>
bool bind(TTree* t, const char* n, T* p) {
    if (!t->GetBranch(n) || t->SetBranchAddress(n, p) < 0) {
        std::cerr << "ERROR: missing/incompatible branch " << n << "\n";
        return false;
    }
    return true;
}

bool trainOne(TTree* src, const std::string& outDir,
              const char* mode, const JetBin& bin, bool osOnly) {
    const std::string model = std::string("BDT_") + mode + "_" + bin.tag;
    const std::string outName = outDir + "/" + model + ".root";
    std::unique_ptr<TFile> out(TFile::Open(outName.c_str(), "RECREATE"));
    if (!out || out->IsZombie()) {
        std::cerr << "ERROR: cannot create " << outName << "\n";
        return false;
    }

    TTree* st = src->CloneTree(0);
    TTree* sv = src->CloneTree(0);
    TTree* bt = src->CloneTree(0);
    TTree* bv = src->CloneTree(0);
    st->SetName("SignalTrain"); sv->SetName("SignalTest");
    bt->SetName("BackgroundTrain"); bv->SetName("BackgroundTest");

    int event = 0, label = 0;
    float weight = 1.f, jetPt = 0.f;
    bool sameSign = false;
    if (!bind(src, "event", &event) || !bind(src, "isSignal", &label) ||
        !bind(src, "weight", &weight) || !bind(src, "jetPt", &jetPt) ||
        !bind(src, "sameSign", &sameSign)) {
        delete st; delete sv; delete bt; delete bv;
        return false;
    }

    Long64_t nst = 0, nsv = 0, nbt = 0, nbv = 0;
    for (Long64_t i = 0; i < src->GetEntries(); ++i) {
        src->GetEntry(i);
        if (jetPt < bin.low || jetPt >= bin.high || (osOnly && sameSign))
            continue;

        const bool test = (event % 5 == 0);
        if (label) {
            (test ? sv : st)->Fill();
            if (test) ++nsv; else ++nst;
        } else {
            (test ? bv : bt)->Fill();
            if (test) ++nbv; else ++nbt;
        }
    }

    if (!nst || !nsv || !nbt || !nbv) {
        std::cerr << "ERROR " << model << ": insufficient class entries. "
                  << "Train S/B=" << nst << "/" << nbt
                  << ", test S/B=" << nsv << "/" << nbv << "\n";
        delete st; delete sv; delete bt; delete bv;
        return false;
    }

    std::cout << "\n" << model << "\nTrain S/B: " << nst << "/" << nbt
              << "  Test S/B: " << nsv << "/" << nbv << "\n";

    const std::string dataName = "data/dataset_" + model;
    TMVA::DataLoader loader(dataName.c_str());

    std::vector<std::string> variableNames;
    for (const char* v : vars) {
        loader.AddVariable(v, 'F');
        variableNames.emplace_back(v);
    }
    if (!osOnly) {
        loader.AddVariable("sameSign", 'F');
        variableNames.emplace_back("sameSign");
    }

    loader.SetWeightExpression("weight");
    loader.AddSignalTree(st, 1.0, TMVA::Types::kTraining);
    loader.AddSignalTree(sv, 1.0, TMVA::Types::kTesting);
    loader.AddBackgroundTree(bt, 1.0, TMVA::Types::kTraining);
    loader.AddBackgroundTree(bv, 1.0, TMVA::Types::kTesting);
    loader.PrepareTrainingAndTestTree(
        TCut(""), TCut(""), "NormMode=NumEvents:!V"
    );

    TMVA::Factory factory(
        model.c_str(), out.get(),
        "!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification"
    );

    factory.BookMethod(
        &loader, TMVA::Types::kBDT, "BDT",
        "!H:!V:NTrees=600:MinNodeSize=2.5%:MaxDepth=3:"
        "BoostType=Grad:Shrinkage=0.10:nCuts=30"
    );

    factory.TrainAllMethods();
    factory.TestAllMethods();
    factory.EvaluateAllMethods();

    // Read split selectors from the exported XML instead of querying the
    // in-memory DecisionTree importance API (which crashes in ROOT 6.40.04).
    saveVariableImportance(variableNames, model, outDir, out.get());

    out->cd();
    st->Write(); sv->Write(); bt->Write(); bv->Write();
    delete st; delete sv; delete bt; delete bv;
    out->Close();

    std::cout << "Wrote " << outName << "\n";
    return true;
}

} // namespace

void train_bdt_jetbins(const char* inputName = "data/candidates_inc_5_100.root",
                       const char* outputDir = "tmva") {
    TMVA::Tools::Instance();
    TFile input(inputName, "READ");
    if (input.IsZombie()) {
        std::cerr << "ERROR: cannot open " << inputName << "\n";
        return;
    }

    TTree* t = nullptr;
    input.GetObject("Candidates", t);
    if (!t) {
        std::cerr << "ERROR: no Candidates tree\n";
        return;
    }

    for (const char* v : vars) {
        if (!t->GetBranch(v)) {
            std::cerr << "ERROR: missing feature branch " << v << "\n";
            return;
        }
    }

    gSystem->mkdir(outputDir, true);
    bool ok = true;
    for (const auto& b : bins) {
        ok = trainOne(t, outputDir, "inclusive", b, false) && ok;
        ok = trainOne(t, outputDir, "OS", b, true) && ok;
    }

    std::cout << "\nTMVA training " << (ok ? "completed." : "had errors.")
              << " Outputs in " << outputDir << "/\n"
              << "Variable-importance plots in " << outputDir << "/plots/<model>/\n";
}
