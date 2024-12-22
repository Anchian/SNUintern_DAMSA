#include <cmath> // 수학 함수 사용을 위해 포함
#include <iostream> // std::cerr 사용을 위해 포함

int CalculateOptimalBins(int nEntries) {
    if (nEntries <= 0) return 10; // 데이터가 없는 경우 기본 값 반환
    return std::max(5, static_cast<int>(sqrt(nEntries))); // 최소 5개의 빈을 유지, sqrt(N) 반환
}


void rootmacro3_signal() {
    // ROOT 파일 열기
    TFile *file = TFile::Open("stage_output.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open the file!" << std::endl;
        return;
    }

    // 필요한 트리 가져오기
    TTree *treeALP = (TTree*)file->Get("ALPpair");
    TTree *treeInvMass = (TTree*)file->Get("InvariantMass");
    TTree *tree6 = (TTree*)file->Get("AnalyzeGamma");

    // 히스토그램 생성
    TH1F *histALPTime = new TH1F("histALPTime", "ALP Time Distribution;Time (ns);Number of ALPs", 100, 0, 2);
    TH1F *histRealInvMass = new TH1F("histRealInvMass", "Real Invariant Mass;Invariant Mass (MeV);Number of Events", 70, 10, 20);
    TH1F *histDetectorInvMass = new TH1F("histDetectorInvMass", "Detector Invariant Mass;Invariant Mass (MeV);Number of Events", 70, 10, 20);
    

    // 데이터 채우기
    treeALP->Draw("ALPpair_time >> histALPTime", "ALPpair_time != 0");
    treeInvMass->Draw("invariant_mass >> histRealInvMass", "invariant_mass != 0");
    tree6->Draw("Vertex_0 >> histDetectorInvMass", "Vertex_0 != 0");
    

    // 캔버스 생성 및 그래프 그리기
    TCanvas *canvas = new TCanvas("canvas", "ALP Analysis", 1800, 600);
    canvas->Divide(3, 1);

    // ALP Time Distribution
    canvas->cd(1);
    histALPTime->Draw();
    
    // 0.2보다 작은 값들의 갯수 계산
    int count_below_02 = histALPTime->Integral(1, histALPTime->FindBin(0.2) - 1);
    int total_count = histALPTime->GetEntries();
    double percentage = (double)count_below_02 / total_count * 100;
    
    TLatex *lat1 = new TLatex();
    lat1->SetNDC();
    lat1->SetTextSize(0.04);
    lat1->DrawLatex(0.15, 0.85, Form("Events < 0.2 ns: %d (%.2f%%)", count_below_02, percentage));

    // Real Invariant Mass
    canvas->cd(2);
    histRealInvMass->Draw();
    
    // 정확히 50 MeV인 이벤트 수 계산
    int count_50mev_real = treeInvMass->GetEntries("invariant_mass == 50");
    
    TLatex *lat2 = new TLatex();
    lat2->SetNDC();
    lat2->SetTextSize(0.04);
    lat2->DrawLatex(0.15, 0.85, Form("Events at 50 MeV: %d", count_50mev_real));
    lat2->DrawLatex(0.15, 0.80, "Real Invariant Mass");

    // Detector Invariant Mass
    canvas->cd(3);
    histDetectorInvMass->Draw();
    
    // 정확히 50 MeV인 이벤트 수 계산
    int count_50mev_detector = tree6->GetEntries("Vertex_0 == 50");
    
    TLatex *lat3 = new TLatex();
    lat3->SetNDC();
    lat3->SetTextSize(0.04);
    lat3->DrawLatex(0.15, 0.85, Form("Events at 50 MeV: %d", count_50mev_detector));
    lat3->DrawLatex(0.15, 0.80, "Detector Invariant Mass");

    // 이미지로 저장
    canvas->SaveAs("output/alp_analysis.png");

    // 파일 닫기
    file->Close();
}
