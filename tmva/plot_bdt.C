
#include "TFile.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TList.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

namespace {

TH1* getHistogram(TDirectory* dir, const char* name) {
    TH1* original = nullptr;
    dir->GetObject(name, original);

    if (!original) {
        std::cerr << "Missing histogram: " << name << std::endl;
        return nullptr;
    }

    TH1* copy = dynamic_cast<TH1*>(
        original->Clone(Form("%s_plot", name))
    );

    copy->SetDirectory(nullptr);
    return copy;
}

void normalize(TH1* h) {
    if (h && h->Integral() > 0)
        h->Scale(1.0 / h->Integral());
}

void plotOne(const std::string& inputName) {

    // Derive model name from ROOT filename.
    std::string model = inputName.substr(0, inputName.size() - 5);

    std::string outputDir = "tmva/plots/" + model;

    std::cout << "\n====================================\n"
              << "Processing: " << inputName << "\n"
              << "Model: " << model << "\n"
              << "====================================\n";

    TFile f(("tmva/" + inputName).c_str(), "READ");

    if (f.IsZombie()) {
        std::cerr << "Cannot open file." << std::endl;
        return;
    }

    // Dataset directory is derived from the ROOT filename.
    std::string dataset =
        "data/dataset_" + model + "/Method_BDT/BDT";

    TDirectory* dir = f.GetDirectory(dataset.c_str());

    if (!dir) {
        std::cerr << "Cannot find BDT directory: "
                  << dataset << std::endl;
        return;
    }

    gSystem->mkdir(outputDir.c_str(), true);

    // 1. Classifier output
    {
        TH1* s = getHistogram(dir, "MVA_BDT_S");
        TH1* b = getHistogram(dir, "MVA_BDT_B");

        if (s && b) {
            normalize(s);
            normalize(b);

            TCanvas c("classifier", "Classifier output", 900, 700);

            s->SetLineColor(kBlue);
            b->SetLineColor(kRed);
            s->SetLineWidth(2);
            b->SetLineWidth(2);

            s->SetTitle(
                Form("%s: Classifier output;"
                     "BDT response;Normalized events",
                     model.c_str())
            );

            s->Draw("HIST");
            b->Draw("HIST SAME");

            TLegend legend(0.65, 0.75, 0.88, 0.88);
            legend.AddEntry(s, "Signal (test)", "l");
            legend.AddEntry(b, "Background (test)", "l");
            legend.Draw();

            c.SaveAs(Form("%s/classifier_output.png",
                          outputDir.c_str()));
        }

        delete s;
        delete b;
    }

    // 2. ROC curve
    {
        TH1* roc = getHistogram(dir, "MVA_BDT_rejBvsS");

        if (roc) {
            TCanvas c("roc", "ROC curve", 900, 700);

            roc->SetTitle(
                Form("%s: ROC curve;"
                     "Signal efficiency;Background rejection",
                     model.c_str())
            );

            roc->SetLineColor(kBlue);
            roc->SetLineWidth(3);
            roc->Draw("HIST");

            c.SaveAs(Form("%s/roc_curve.png",
                          outputDir.c_str()));
        }

        delete roc;
    }

    // 3. Efficiencies
    {
        TH1* s = getHistogram(dir, "MVA_BDT_effS");
        TH1* b = getHistogram(dir, "MVA_BDT_effB");

        if (s && b) {
            TCanvas c("efficiency", "Efficiencies", 900, 700);

            s->SetLineColor(kBlue);
            b->SetLineColor(kRed);
            s->SetLineWidth(2);
            b->SetLineWidth(2);

            s->SetTitle(
                Form("%s: Classification efficiencies;"
                     "BDT cut;Efficiency",
                     model.c_str())
            );

            s->Draw("HIST");
            b->Draw("HIST SAME");

            TLegend legend(0.65, 0.75, 0.88, 0.88);
            legend.AddEntry(s, "Signal", "l");
            legend.AddEntry(b, "Background", "l");
            legend.Draw();

            c.SaveAs(Form("%s/efficiencies.png",
                          outputDir.c_str()));
        }

        delete s;
        delete b;
    }

    // 4. Overtraining check
    {
        TH1* st = getHistogram(dir, "MVA_BDT_Train_S");
        TH1* bt = getHistogram(dir, "MVA_BDT_Train_B");
        TH1* sv = getHistogram(dir, "MVA_BDT_S");
        TH1* bv = getHistogram(dir, "MVA_BDT_B");

        if (st && bt && sv && bv) {
            normalize(st);
            normalize(bt);
            normalize(sv);
            normalize(bv);

            TCanvas c("overtraining", "Overtraining", 900, 700);

            st->SetLineColor(kBlue);
            bt->SetLineColor(kRed);
            sv->SetMarkerColor(kBlue);
            bv->SetMarkerColor(kRed);

            sv->SetMarkerStyle(20);
            bv->SetMarkerStyle(24);

            st->SetLineWidth(2);
            bt->SetLineWidth(2);

            st->SetTitle(
                Form("%s: Overtraining check;"
                     "BDT response;Normalized events",
                     model.c_str())
            );

            st->Draw("HIST");
            bt->Draw("HIST SAME");
            sv->Draw("E SAME");
            bv->Draw("E SAME");

            TLegend legend(0.58, 0.68, 0.88, 0.88);
            legend.AddEntry(st, "Signal (train)", "l");
            legend.AddEntry(sv, "Signal (test)", "p");
            legend.AddEntry(bt, "Background (train)", "l");
            legend.AddEntry(bv, "Background (test)", "p");
            legend.Draw();

            c.SaveAs(Form("%s/overtraining.png",
                          outputDir.c_str()));
        }

        delete st;
        delete bt;
        delete sv;
        delete bv;
    }

    std::cout << "Finished: " << model
              << "\nPlots saved to: " << outputDir << std::endl;
}

} // namespace

void plot_bdt() {

    gStyle->SetOptStat(0);

    TSystemDirectory directory("tmva", "tmva");
    TList* files = directory.GetListOfFiles();

    if (!files) {
        std::cerr << "Cannot read tmva/ directory." << std::endl;
        return;
    }

    std::vector<std::string> rootFiles;

    TIter next(files);
    while (TSystemFile* file =
           dynamic_cast<TSystemFile*>(next())) {

        std::string name = file->GetName();

        if (!file->IsDirectory() &&
            name.size() > 5 &&
            name.substr(name.size() - 5) == ".root") {
            rootFiles.push_back(name);
        }
    }

    std::sort(rootFiles.begin(), rootFiles.end());

    std::cout << "Found " << rootFiles.size()
              << " ROOT files in tmva/\n";

    for (const auto& file : rootFiles)
        plotOne(file);

    std::cout << "\nAll plotting attempts completed.\n";
}
