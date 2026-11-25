#include "analysisUtils.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// Compare signal and background for opposite-sign, same-jet K pi candidates,
// separately in shared-jet pT bins: 1-10, 10-20, and 20+ GeV/c.
// This OS-only companion macro intentionally excludes sameSign as a feature.

namespace {

struct Feature {
    const char* name;
    const char* label;
    int nBins;
    double lo, hi;
};

struct JetPtBin {
    const char* label;
    double lo, hi;
    const char* fileTag;
};

}  // namespace

int main(int argc, char* argv[]) {

    const std::string inputName =
        getInputFile(argc, argv, "data/candidates_os_inclusive.root");
    const std::string tag = getTag(argc, argv, inputName);

    gSystem->mkdir("plots", true);

    const std::vector<Feature> features = {
        {"pxK",          "K p_{x} [GeV/c]", 50, -10.0, 10.0},
        {"pyK",          "K p_{y} [GeV/c]", 50, -10.0, 10.0},
        {"pzK",          "K p_{z} [GeV/c]", 50, -10.0, 10.0},
        {"pxPi",         "#pi p_{x} [GeV/c]", 50, -10.0, 10.0},
        {"pyPi",         "#pi p_{y} [GeV/c]", 50, -10.0, 10.0},
        {"pzPi",         "#pi p_{z} [GeV/c]", 50, -10.0, 10.0},
        {"ptK",          "K p_{T} [GeV/c]", 50, 0.0, 10.0},
        {"ptPi",         "#pi p_{T} [GeV/c]", 50, 0.0, 10.0},
        {"etaK",         "K #eta", 40, -1.1, 1.1},
        {"etaPi",        "#pi #eta", 40, -1.1, 1.1},
        {"dEta",         "#Delta#eta (K - #pi)", 50, -2.2, 2.2},
        {"dPhi",         "#Delta#phi (K - #pi)", 50, -3.2, 3.2},
        {"dR",           "#DeltaR (K, #pi)", 50, 0.0, 4.0},
        {"openAngle",    "Opening angle [rad]", 50, 0.0, 3.2},
        {"ptBalance",    "|pT_{K} - pT_{#pi}| / (pT_{K} + pT_{#pi})", 50, 0.0, 1.0},
        {"ptHard",       "Harder daughter p_{T} [GeV/c]", 50, 0.0, 10.0},
        {"ptSoft",       "Softer daughter p_{T} [GeV/c]", 50, 0.0, 6.0},
        {"ptSum",        "Scalar p_{T} sum (K + #pi) [GeV/c]", 50, 0.0, 15.0},
        {"pairPt",       "Pair p_{T} [GeV/c]", 50, 0.0, 15.0},
        {"pairEta",      "Pair #eta", 40, -1.1, 1.1},
        {"cosThetaStar", "cos#theta* (K, pair rest frame)", 40, -1.0, 1.0},

        // sameJet is deliberately omitted: after the new selection it is
        // always 1 and therefore carries no separation information.
        {"jetPt",        "Shared-jet p_{T} [GeV/c]", 50, 0.0, 40.0},
        {"zJet",         "Pair p_{T} / jet p_{T}", 50, 0.0, 1.5},
        {"dRJet",        "#DeltaR (pair, shared jet axis)", 50, 0.0, 5.0},
        {"jetMass",      "Shared-jet mass [GeV/c^{2}]", 50, 0.0, 8.0},
        {"jetEta",       "Shared-jet #eta", 40, -2.0, 2.0},

        {"isoRel02",     "Charged isolation, R = 0.2", 50, 0.0, 5.0},
        {"isoRel03",     "Charged isolation, R = 0.3", 50, 0.0, 5.0},
        {"isoRel",       "Charged isolation, R = 0.4", 50, 0.0, 5.0},
        {"isoRel06",     "Charged isolation, R = 0.6", 50, 0.0, 5.0},
        {"isoNeutralRel","Neutral isolation, R = 0.4", 50, 0.0, 5.0},

        {"nCharged",     "Charged multiplicity (|#eta|<1)", 50, 0.0, 50.0},
        {"nJets",        "Number of jets (|#eta|<1)", 8, 0.0, 8.0},
        {"sphericityT",  "Transverse sphericity", 40, 0.0, 1.0},
        {"thrustT",      "Transverse thrust", 40, 0.6, 1.0},
        {"mass",         "M_{K#pi} [GeV/c^{2}] (not a training variable)", 50, 1.6, 2.1}
    };

    const JetPtBin jetBins[] = {
        {"1-20 GeV/c", 1.0, 20.0, "jetPt1_20"},
        {"20-30 GeV/c", 20.0, 30.0, "jetPt20_30"},
        {"30+ GeV/c", 30.0, 1.0e9, "jetPt30plus"}
    };
    const int nJetBins = sizeof(jetBins) / sizeof(jetBins[0]);

    TFile input(inputName.c_str(), "READ");
    if (input.IsZombie()) {
        std::cerr << "Could not open " << inputName << "\n";
        return 1;
    }

    TTree* tree = nullptr;
    input.GetObject("Candidates", tree);
    if (!tree) {
        std::cerr << "Could not find TTree 'Candidates' in " << inputName << "\n";
        return 1;
    }

    int isSignal = 0;
    float weight = 1.0f;
    float sameJet = 0.0f;
    tree->SetBranchAddress("isSignal", &isSignal);
    tree->SetBranchAddress("weight", &weight);
    if (tree->SetBranchAddress("sameJet", &sameJet) < 0) {
        std::cerr << "Required jet branches missing. Regenerate candidates with buildtree_jetbins.\n";
        return 1;
    }

    // Selection-only branch: never added to the ML feature list.
    bool sameSign = false;
    if (tree->SetBranchAddress("sameSign", &sameSign) < 0) {
        std::cerr << "Branch 'sameSign' not found. Regenerate the candidate file.\n";
        return 1;
    }

    std::vector<float> values(features.size(), 0.0f);
    for (size_t k = 0; k < features.size(); ++k) {
        if (tree->SetBranchAddress(features[k].name, &values[k]) < 0) {
            std::cerr << "Branch '" << features[k].name << "' not found. "
                      << "Regenerate candidates with the updated buildtree.\n";
            return 1;
        }
    }

    // jetPt is read through the feature vector below. Do not bind the
    // same ROOT branch to a second address: SetBranchAddress would
    // replace the first address and leave the scalar jetPt unchanged.
    size_t jetPtIndex = features.size();
    for (size_t k = 0; k < features.size(); ++k) {
        if (features[k].name == std::string("jetPt")) {
            jetPtIndex = k;
            break;
        }
    }
    if (jetPtIndex == features.size()) {
        std::cerr << "Internal error: jetPt is not in the feature list.\\n";
        return 1;
    }

    std::vector<std::vector<TH1F*> > hSig(nJetBins), hBkg(nJetBins);
    for (int ib = 0; ib < nJetBins; ++ib) {
        for (const auto& f : features) {
            const std::string suffix = std::string("_") + jetBins[ib].fileTag;
            TH1F* s = new TH1F((std::string("sig_") + f.name + suffix).c_str(),
                (std::string(";") + f.label + ";Fraction of candidates").c_str(),
                f.nBins, f.lo, f.hi);
            TH1F* b = new TH1F((std::string("bkg_") + f.name + suffix).c_str(),
                (std::string(";") + f.label + ";Fraction of candidates").c_str(),
                f.nBins, f.lo, f.hi);
            s->SetDirectory(nullptr);
            b->SetDirectory(nullptr);
            s->Sumw2();
            b->Sumw2();
            hSig[ib].push_back(s);
            hBkg[ib].push_back(b);
        }
    }

    long long nSig[nJetBins] = {0, 0, 0};
    long long nBkg[nJetBins] = {0, 0, 0};
    const Long64_t nEntries = tree->GetEntries();

    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        if (sameJet < 0.5f) continue; // enforce same-jet selection defensively
        if (sameSign) continue;       // enforce OS selection defensively

        int ib = -1;
        for (int b = 0; b < nJetBins; ++b) {
            if (values[jetPtIndex] >= jetBins[b].lo &&
                values[jetPtIndex] < jetBins[b].hi) {
                ib = b;
                break;
            }
        }
        if (ib < 0) continue;

        if (isSignal == 1) ++nSig[ib];
        else ++nBkg[ib];

        for (size_t k = 0; k < features.size(); ++k) {
            if (isSignal == 1) hSig[ib][k]->Fill(values[k]);
            else hBkg[ib][k]->Fill(values[k], weight);
        }
    }
    input.Close();

    gStyle->SetOptStat(0);
    const int nCols = 6;
    const int nPads = static_cast<int>(features.size()) + 1; // + legend
    const int nRows = (nPads + nCols - 1) / nCols;

    for (int ib = 0; ib < nJetBins; ++ib) {
        std::cout << "\n===== shared-jet pT " << jetBins[ib].label << " =====\n"
                  << "Signal rows: " << nSig[ib] << "\n"
                  << "Background rows: " << nBkg[ib] << "\n"
                  << "Separation (0 = identical shapes, 1 = no overlap):\n";

        struct Score { std::string name; double sep; };
        std::vector<Score> scores;
        for (size_t k = 0; k < features.size(); ++k) {
            TH1F* s = hSig[ib][k];
            TH1F* b = hBkg[ib][k];
            const double iS = s->Integral(0, s->GetNbinsX() + 1);
            const double iB = b->Integral(0, b->GetNbinsX() + 1);
            if (iS > 0) s->Scale(1.0 / iS);
            if (iB > 0) b->Scale(1.0 / iB);

            double sep = 0.0;
            for (int bin = 0; bin <= s->GetNbinsX() + 1; ++bin)
                sep += std::abs(s->GetBinContent(bin) - b->GetBinContent(bin));
            sep *= 0.5;
            scores.push_back({features[k].name, sep});
        }
        std::sort(scores.begin(), scores.end(),
                  [](const Score& a, const Score& b) { return a.sep > b.sep; });
        for (const auto& score : scores)
            std::cout << "  " << score.name << ": " << score.sep << "\n";

        TCanvas canvas((std::string("canvas_") + jetBins[ib].fileTag).c_str(),
                       (std::string("OS-only features, shared-jet pT ") + jetBins[ib].label).c_str(),
                       400 * nCols, 330 * nRows);
        canvas.Divide(nCols, nRows);

        for (size_t k = 0; k < features.size(); ++k) {
            canvas.cd(k + 1);
            TH1F* s = hSig[ib][k];
            TH1F* b = hBkg[ib][k];
            s->SetLineColor(kBlue + 1);
            s->SetLineWidth(2);
            b->SetLineColor(kRed + 1);
            b->SetLineWidth(2);
            s->SetMaximum(1.3 * std::max(s->GetMaximum(), b->GetMaximum()));
            s->SetMinimum(0);
            s->Draw("HIST");
            b->Draw("HIST SAME");
        }

        canvas.cd(features.size() + 1);
        TLegend leg(0.08, 0.35, 0.92, 0.65);
        leg.AddEntry(hSig[ib][0], "Signal: truth D^{0}#rightarrowK#pi", "l");
        leg.AddEntry(hBkg[ib][0], "Background: other opposite-sign same-jet K#pi pairs", "l");
        leg.Draw();

        const std::string outName =
            "plots/features_os" + tag + "_" + jetBins[ib].fileTag + ".png";
        canvas.SaveAs(outName.c_str());
        std::cout << "Saved " << outName << "\n";
    }

    for (int ib = 0; ib < nJetBins; ++ib) {
        for (TH1F* h : hSig[ib]) delete h;
        for (TH1F* h : hBkg[ib]) delete h;
    }
    return 0;
}
