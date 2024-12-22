#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>

void rootmacro5_DirectPhoton() {

    // ROOT 파일 열기
    TFile *file = TFile::Open("stage_output.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open the file!" << std::endl;
        return;
    }

    // 트리 가져오기
    TTree *tree_alp = (TTree*)file->Get("ALPrealphotonEnergy");
    TTree *tree_stage2 = (TTree*)file->Get("Stage2");

    // 히스토그램 생성
    TH1F *h_alpPhotonEnergy = new TH1F("h_alpPhotonEnergy", "ALPrealphotonEnergy;Energy (MeV);Number of Particles", 100, 0, 350);
    TH1F *h_stage2GammaEnergy = new TH1F("h_stage2GammaEnergy", "Stage2 Gamma Energy;Energy (MeV);Number of Particles", 100, 0, 350);

    // 데이터 채우기
    tree_alp->Draw("ALPrealphotonEnergy >> h_alpPhotonEnergy", "ALPrealphotonEnergy != 0");
    tree_stage2->Draw("gamma_Energy >> h_stage2GammaEnergy", "gamma_Energy != 0");

    // 캔버스 생성
    TCanvas *canvas = new TCanvas("canvas", "Direct Photon Energy", 1800, 600);
    canvas->cd();
    canvas->SetLogy(); 
    // 최대값 설정하여 잘리지 않도록 조정
    Double_t maxY = std::max(h_alpPhotonEnergy->GetMaximum(), h_stage2GammaEnergy->GetMaximum());
    h_alpPhotonEnergy->SetMaximum(maxY * 1.2);  // 최대 높이의 20% 여유 설정

    // 히스토그램 그리기
    h_alpPhotonEnergy->SetLineColor(kBlue);
    h_alpPhotonEnergy->Draw();

    h_stage2GammaEnergy->SetLineColor(kRed);
    h_stage2GammaEnergy->Draw("SAME");

    // 범례 추가
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);  // 범례 위치 설정
    legend->AddEntry(h_alpPhotonEnergy, "ALP Real Photon Energy", "l");
    legend->AddEntry(h_stage2GammaEnergy, "Stage2 Gamma Energy", "l");
    legend->Draw();

    // 이미지로 저장
    canvas->SaveAs("output/DirectPhotonEnergy.png");

    // 파일 닫기
    file->Close();
    
}