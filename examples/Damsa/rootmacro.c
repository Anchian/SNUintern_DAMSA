#include <cmath> // 수학 함수 사용을 위해 포함






// 빈의 개수를 자동으로 계산해주는 함수
int CalculateOptimalBins(int nEntries) {
    if (nEntries <= 0) return 10; // 데이터가 없는 경우 기본 값 반환
    return std::max(5, static_cast<int>(sqrt(nEntries))); // 최소 5개의 빈을 유지, sqrt(N) 반환
}

void rootmacro() {
    // ROOT 파일 열기
    TFile *file = TFile::Open("stage_output.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open the file!" << std::endl;
        return;
    }

    // Stage1 및 Stage2 트리 가져오기
    TTree *tree1 = (TTree*)file->Get("Stage1"); // Stage1 트리 이름 확인 및 수정 필요
    TTree *tree2 = (TTree*)file->Get("Stage2"); // Stage2 트리 이름 확인 및 수정 필요
    TTree *tree0 = (TTree*)file->Get("Stage0"); // Stage0 
    TTree *tree3 = (TTree*)file->Get("GammaTimeDistribution");
    TTree *tree4 = (TTree*)file->Get("GammaTimeDistributionCut_Pair");
    TTree *tree5 = (TTree*)file->Get("InvariantMass");
    TTree *tree6 = (TTree*)file->Get("AnalyzeGamma");
    TTree *tree7 = (TTree*)file->Get("Vertex_distance_distribution");
    TTree *tree8 = (TTree*)file->Get("ALPpair");

    int nEntries_gamma_stage1 = tree1->GetEntries("gamma_Energy != 0");
    int nEntries_electron_stage1 = tree1->GetEntries("electron_Energy != 0");
    int nEntries_neutron_stage1 = tree1->GetEntries("neutron_Energy != 0");

    int nEntries_gamma_stage2 = tree2->GetEntries("gamma_Energy != 0");
    int nEntries_electron_stage2 = tree2->GetEntries("electron_Energy != 0");
    int nEntries_neutron_stage2 = tree2->GetEntries("neutron_Energy != 0");

    int nEntries_gamma_stage0 = tree0->GetEntries("gamma_Energy != 0");

    int nEntries_gamma_stage2_timedistribution = tree3->GetEntries("gamma_time!=0");
    int nEntries_gamma_stage2_timedistribution_cut = tree4->GetEntries("gamma_time_cut!=0");

    int nEntries_invariant_mass = tree5->GetEntries("invariant_mass!=0");
    int nEntries_invariant_mass_cut = tree5->GetEntries("VertexCut_invariant_mass!=0");
    int nEntries_invariant_mass_inside_chamber = tree5->GetEntries("Vertex_Inside_Chamber!=0");
    int nEntries_invariant_mass_cut_chamber_steel = tree5->GetEntries("Vertex_Chamber_Steel!=0");
    int nEntries_invariant_mass_cut_target = tree5->GetEntries("Vertex_target!=0");
    int nEntries_invariant_mass_cut_outside_chamber = tree5->GetEntries("Vertex_Outside_Chamber>=10");
    
    int nEntries_analyze_gamma_invariant_totalmass = tree6->GetEntries("Vertex_0!=0");
    int nEntries_analyze_gamma_invariant_inside_chamber = tree6->GetEntries("Vertex_Inside_Chamber!=0");
    int nEntries_analyze_gamma_invariant_chamber_steel = tree6->GetEntries("Vertex_Chamber_Steel!=0");
    int nEntries_analyze_gamma_invariant_target = tree6->GetEntries("Vertex_target!=0");
    int nEntries_analyze_gamma_invariant_outside_chamber = tree6->GetEntries("Vertex_Outside_Chamber!=0");

    int nEntries_pair_distance_distribution = tree7->GetEntries("distance_distribution!=0");
    int nEntries_pair_transverse_distance_distribution = tree7->GetEntries("transverse_distance_distribution!=0");
    int nEntries_pair_z_distance_distribution = tree7->GetEntries("z_distance_distribution!=0");
    int nEntries_pair_distance_distribution_cut = tree7->GetEntries("cut_distance_distribution!=0");
    int nEntries_pair_transverse_distance_distribution_cut = tree7->GetEntries("cut_transverse_distance_distribution!=0");
    int nEntries_pair_z_distance_distribution_cut = tree7->GetEntries("cut_z_distance_distribution!=0");

    int nEntries_ALPpair_time = tree8->GetEntries("ALPpair_time!=0");
    
    


    // 빈의 개수 자동 계산
    int bins_gamma_stage1 = CalculateOptimalBins(nEntries_gamma_stage1);
    int bins_electron_stage1 = CalculateOptimalBins(nEntries_electron_stage1);
    int bins_neutron_stage1 = CalculateOptimalBins(nEntries_neutron_stage1);

    int bins_gamma_stage2 = CalculateOptimalBins(nEntries_gamma_stage2);
    int bins_electron_stage2 = CalculateOptimalBins(nEntries_electron_stage2);
    int bins_neutron_stage2 = CalculateOptimalBins(nEntries_neutron_stage2);

    int bins_gamma_stage0 = CalculateOptimalBins(nEntries_gamma_stage0);

    int bins_gamma_stage2_timedistribution = CalculateOptimalBins(nEntries_gamma_stage2_timedistribution);
    int bins_gamma_stage2_timedistribution_cut = CalculateOptimalBins(nEntries_gamma_stage2_timedistribution_cut);
    
    int bins_invariant_mass = CalculateOptimalBins(nEntries_invariant_mass);
    int bins_invariant_mass_cut = CalculateOptimalBins(nEntries_invariant_mass_cut);
    int bins_invariant_mass_cut_pair = CalculateOptimalBins(nEntries_invariant_mass_inside_chamber);
    int bins_invariant_mass_cut_pair2 = CalculateOptimalBins(nEntries_invariant_mass_cut_chamber_steel);
    int bins_invariant_mass_cut_pair3 = CalculateOptimalBins(nEntries_invariant_mass_cut_target);
    int bins_invariant_mass_cut_pair4 = CalculateOptimalBins(nEntries_invariant_mass_cut_outside_chamber);

    int bins_analyze_gamma_invariant_totalmass = CalculateOptimalBins(nEntries_analyze_gamma_invariant_totalmass);
    int bins_analyze_gamma_invariant_inside_chamber = CalculateOptimalBins(nEntries_analyze_gamma_invariant_inside_chamber);
    int bins_analyze_gamma_invariant_chamber_steel = CalculateOptimalBins(nEntries_analyze_gamma_invariant_chamber_steel);
    int bins_analyze_gamma_invariant_target = CalculateOptimalBins(nEntries_analyze_gamma_invariant_target);
    int bins_analyze_gamma_invariant_outside_chamber = CalculateOptimalBins(nEntries_analyze_gamma_invariant_outside_chamber);
    
    int bins_pair_distance_distribution = CalculateOptimalBins(nEntries_pair_distance_distribution);
    int bins_pair_transverse_distance_distribution = CalculateOptimalBins(nEntries_pair_transverse_distance_distribution);
    int bins_pair_z_distance_distribution = CalculateOptimalBins(nEntries_pair_z_distance_distribution);
    int bins_pair_distance_distribution_cut = CalculateOptimalBins(nEntries_pair_distance_distribution_cut);
    int bins_pair_transverse_distance_distribution_cut = CalculateOptimalBins(nEntries_pair_transverse_distance_distribution_cut);
    int bins_pair_z_distance_distribution_cut = CalculateOptimalBins(nEntries_pair_z_distance_distribution_cut);

    int bins_ALPpair_time = CalculateOptimalBins(nEntries_ALPpair_time);

    TH1F *hist0_gamma = new TH1F("hist0_gamma", "Stage0 Gamma Energy;Energy (MeV);Number of Particles", bins_gamma_stage0, 0, 5);

    // Stage1 브랜치에 대한 히스토그램 생성 및 데이터 필터링
    TH1F *hist1_gamma = new TH1F("hist1_gamma", "Stage1 Gamma Energy;Energy (MeV);Number of Particles", bins_gamma_stage1, 0, 20);
    TH1F *hist1_electron = new TH1F("hist1_electron", "Stage1 Electron Energy;Energy (MeV);Number of Particles", 50, 0, 20);
    TH1F *hist1_neutron = new TH1F("hist1_neutron", "Stage1 Neutron Energy;Energy (MeV);Number of Particles", 200, 0, 20);
    
    // Stage2 브랜치에 대한 히스토그램 생성 및 데이터 필터링
    TH1F *hist2_gamma = new TH1F("hist2_gamma", "Stage2 Gamma Energy;Energy (MeV);Number of Particles", 30, 0, 20);
    TH1F *hist2_electron = new TH1F("hist2_electron", "Stage2 Electron Energy;Energy (MeV);Number of Particles", 10, 0, 20);
    TH1F *hist2_neutron = new TH1F("hist2_neutron", "Stage2 Neutron Energy;Energy (MeV);Number of Particles", 5, 0, 20);
    

    TH1F *hist2_gamma_timedistribution = new TH1F("hist2_gamma_timedistribution", "Stage2 Gamma Time Distribution;Time (ns);Number of Particles", bins_gamma_stage2_timedistribution, 1.0,3.5 );
    TH1F *hist2_gamma_timedistribution_cut = new TH1F("hist2_gamma_timedistribution_cut", "Stage2 Gamma Time Distribution (Cut);Time ", bins_gamma_stage2_timedistribution_cut, 0, 0.15);
    
    TH1F *hist2_invariant_mass = new TH1F("hist2_invariant_mass", "Invariant Mass;Invariant Mass (MeV);Number of Particles", bins_invariant_mass, 0, 5);
    TH1F *hist2_invariant_mass_cut = new TH1F("hist2_invariant_mass_cut", "Invariant Mass (Cut);Invariant Mass (MeV);Number of Particles", bins_invariant_mass_cut, 0, 5);
    TH1F *hist2_invariant_mass_cut_inside_chamber = new TH1F("hist2_invariant_mass_cut_inside_chamber", "Invariant Mass (chamber);Invariant Mass (MeV);Number of Particles", bins_invariant_mass_cut_pair, 0, 5);
    TH1F *hist2_invariant_mass_cut_chamber_steel = new TH1F("hist2_invariant_mass_cut_chamber_steel", "Invariant Mass (chamber_steel);Invariant Mass (MeV);Number of Particles", bins_invariant_mass_cut_pair2, 0, 5);
    TH1F *hist2_invariant_mass_cut_target = new TH1F("hist2_invariant_mass_cut_target", "Invariant Mass (target);Invariant Mass (MeV);Number of Particles", bins_invariant_mass_cut_pair3, 0, 5);
    TH1F *hist2_invariant_mass_cut_outside_chamber = new TH1F("hist2_invariant_mass_cut_outside_chamber", "Invariant Mass (outside_chamber);Invariant Mass (MeV);Number of Particles", bins_invariant_mass_cut_pair4, 0, 5);
    TH1F *hist2_invariant_mass_cut_highE = new TH1F("hist2_invariant_mass_cut_highE", "Invariant Mass (target_highE);Invariant Mass (MeV);Number of Particles",bins_invariant_mass_cut_pair3 , 10, 30);
    
    TH1F *hist2_analyze_gamma_invariant_totalmass = new TH1F("hist2_analyze_gamma_invariant_totalmass", "Invariant Mass (totalmass);Invariant Mass (MeV);Number of Particles", bins_analyze_gamma_invariant_totalmass, 0, 1);
    TH1F *hist2_analyze_gamma_invariant_inside_chamber = new TH1F("hist2_analyze_gamma_invariant_inside_chamber", "Invariant Mass (inside_chamber);Invariant Mass (MeV);Number of Particles", bins_analyze_gamma_invariant_inside_chamber, 0, 1);
    TH1F *hist2_analyze_gamma_invariant_chamber_steel = new TH1F("hist2_analyze_gamma_invariant_chamber_steel", "Invariant Mass (chamber_steel);Invariant Mass (MeV);Number of Particles", bins_analyze_gamma_invariant_chamber_steel, 0, 1);
    TH1F *hist2_analyze_gamma_invariant_target = new TH1F("hist2_analyze_gamma_invariant_target", "Invariant Mass (target);Invariant Mass (MeV);Number of Particles", bins_analyze_gamma_invariant_target, 0, 1);
    TH1F *hist2_analyze_gamma_invariant_outside_chamber = new TH1F("hist2_analyze_gamma_invariant_outside_chamber", "Invariant Mass (outside_chamber);Invariant Mass (MeV);Number of Particles", bins_analyze_gamma_invariant_outside_chamber, 0, 1);
    
    TH1F *hist2_pair_distance_distribution = new TH1F("hist2_pair_distance_distribution", "Pair Distance Distribution;Distance (mm);Number of Particles", bins_pair_distance_distribution, 0, 30);
    TH1F *hist2_pair_transverse_distance_distribution = new TH1F("hist2_pair_transverse_distance_distribution", "Pair Transverse Distance Distribution;Transverse Distance (mm);Number of Particles", bins_pair_transverse_distance_distribution, 0, 30);
    TH1F *hist2_pair_z_distance_distribution = new TH1F("hist2_pair_z_distance_distribution", "Pair Z Distance Distribution;Z Distance (mm);Number of Particles", bins_pair_z_distance_distribution, 0, 0.1);
    TH1F *hist2_pair_distance_distribution_cut = new TH1F("hist2_pair_distance_distribution_cut", "Pair Distance Distribution (Cut);Distance (mm);Number of Particles", bins_pair_distance_distribution_cut, 0, 0.1);
    TH1F *hist2_pair_transverse_distance_distribution_cut = new TH1F("hist2_pair_transverse_distance_distribution_cut", "Pair Transverse Distance Distribution (Cut);Transverse Distance (mm);Number of Particles", bins_pair_transverse_distance_distribution_cut, 0, 0.1);
    TH1F *hist2_pair_z_distance_distribution_cut = new TH1F("hist2_pair_z_distance_distribution_cut", "Pair Z Distance Distribution (Cut);Z Distance (mm);Number of Particles", bins_pair_z_distance_distribution_cut, 0, 0.001);
    
    TH1F *hist2_ALPpair_time = new TH1F("hist2_ALPpair_time", "ALPpair Time Distribution;Time (ns);Number of Particles", 100, 0, 2);
    // 데이터 필터링 및 히스토그램 채우기
    tree0->Draw("gamma_Energy >> hist0_gamma", "gamma_Energy != 0");
    // Stage1 트리
    tree1->Draw("gamma_Energy >> hist1_gamma", "gamma_Energy != 0");
    tree1->Draw("electron_Energy >> hist1_electron", "electron_Energy != 0");
    tree1->Draw("neutron_Energy >> hist1_neutron", "neutron_Energy != 0");
    
    // Stage2 트리
    tree2->Draw("gamma_Energy >> hist2_gamma", "gamma_Energy != 0");
    tree2->Draw("electron_Energy >> hist2_electron", "electron_Energy != 0");
    tree2->Draw("neutron_Energy >> hist2_neutron", "neutron_Energy != 0");
    
    tree3->Draw("gamma_time >> hist2_gamma_timedistribution", "gamma_time!=0");
    tree4->Draw(" gamma_time_cut >> hist2_gamma_timedistribution_cut", "gamma_time_cut!=0");

    tree5->Draw("invariant_mass >> hist2_invariant_mass", "invariant_mass!=0");
    tree5->Draw("VertexCut_invariant_mass >> hist2_invariant_mass_cut", "VertexCut_invariant_mass!=0");
    tree5->Draw("Vertex_Inside_Chamber >> hist2_invariant_mass_cut_inside_chamber", "Vertex_Inside_Chamber!=0");
    tree5->Draw("Vertex_Chamber_Steel >> hist2_invariant_mass_cut_chamber_steel", "Vertex_Chamber_Steel!=0");
    tree5->Draw("Vertex_target >> hist2_invariant_mass_cut_target", "Vertex_target!=0");
    tree5->Draw("Vertex_Outside_Chamber >> hist2_invariant_mass_outside_chamber", "Vertex_Outside_Chamber!=0");
    tree5->Draw("Vertex_target >> hist2_invariant_mass_cut_highE", "Vertex_target>=10");

    tree6->Draw("Vertex_0 >> hist2_analyze_gamma_invariant_totalmass", "Vertex_0!=0");
    tree6->Draw("Vertex_Inside_Chamber >> hist2_analyze_gamma_invariant_inside_chamber", "Vertex_Inside_Chamber!=0");
    tree6->Draw("Vertex_Chamber_Steel >> hist2_analyze_gamma_invariant_chamber_steel", "Vertex_Chamber_Steel!=0");
    tree6->Draw("Vertex_target >> hist2_analyze_gamma_invariant_target", "Vertex_target!=0");
    tree6->Draw("Vertex_Outside_Chamber >> hist2_analyze_gamma_invariant_outside_chamber", "Vertex_Outside_Chamber!=0");

    tree7->Draw("distance_distribution >> hist2_pair_distance_distribution", "distance_distribution!=0");
    tree7->Draw("transverse_distance_distribution >> hist2_pair_transverse_distance_distribution", "transverse_distance_distribution!=0");
    tree7->Draw("z_distance_distribution >> hist2_pair_z_distance_distribution", "z_distance_distribution!=0");
    tree7->Draw("cut_distance_distribution >> hist2_pair_distance_distribution_cut", "cut_distance_distribution!=0");
    tree7->Draw("cut_transverse_distance_distribution >> hist2_pair_transverse_distance_distribution_cut", "cut_transverse_distance_distribution!=0");
    tree7->Draw("cut_z_distance_distribution >> hist2_pair_z_distance_distribution_cut", "cut_z_distance_distribution!=0");

    tree8->Draw("ALPpair_time >> hist2_ALPpair_time", "ALPpair_time!=0");
    // 캔버스 생성 및 그래프 그리기
    TCanvas *canvas1 = new TCanvas("canvas1", "Stage1 Branches", 1600, 1200);
    canvas1->Divide(2, 2); // 캔버스를 2x2로 분할
    canvas1->cd(1);
    hist1_gamma->SetLineColor(kBlue);
    hist1_gamma->Draw();
    TLatex *latex1 = new TLatex();
    latex1->SetNDC(); // Normalized device coordinates
    latex1->SetTextSize(0.04);
    latex1->DrawLatex(0.2, 0.8, Form("Above 20: %f", hist1_gamma->Integral(hist1_gamma->FindBin(20), hist1_gamma->GetNbinsX())));

    canvas1->cd(2);
    hist1_electron->SetLineColor(kRed);
    hist1_electron->Draw();
    TLatex *latex2 = new TLatex();
    latex2->SetNDC();
    latex2->SetTextSize(0.04);
    latex2->DrawLatex(0.2, 0.8, Form("Above 20: %f", hist1_electron->Integral(hist1_electron->FindBin(20), hist1_electron->GetNbinsX())));

    canvas1->cd(3);
    hist1_neutron->SetLineColor(kGreen);
    hist1_neutron->Draw();
    TLatex *latex3 = new TLatex();
    latex3->SetNDC();
    latex3->SetTextSize(0.04);
    latex3->DrawLatex(0.2, 0.8, Form("Above 20: %f", hist1_neutron->Integral(hist1_neutron->FindBin(20), hist1_neutron->GetNbinsX())));

    canvas1->cd(4);
    hist1_gamma->SetLineColor(kBlue);
    hist1_gamma->Draw();
    hist1_electron->SetLineColor(kRed);
    hist1_electron->Draw("same");
    hist1_neutron->SetLineColor(kGreen);
    hist1_neutron->Draw("same");

    
    // Stage2의 그래프를 다른 캔버스에 그리기
    TCanvas *canvas2 = new TCanvas("canvas2", "Stage2 Branches", 1600, 1200);
    canvas2->Divide(2, 2);

    canvas2->cd(1);
    hist2_gamma->SetLineColor(kBlue);
    hist2_gamma->Draw();
    TLatex *latex4 = new TLatex();
    latex4->SetNDC();
    latex4->SetTextSize(0.04);
    latex4->DrawLatex(0.2, 0.8, Form("Above 20: %f", hist2_gamma->Integral(hist2_gamma->FindBin(20), hist2_gamma->GetNbinsX())));

    canvas2->cd(2);
    hist2_electron->SetLineColor(kRed);
    hist2_electron->Draw();
    TLatex *latex5 = new TLatex();
    latex5->SetNDC();
    latex5->SetTextSize(0.04);
    latex5->DrawLatex(0.2, 0.8, Form("Above 20: %f", hist2_electron->Integral(hist2_electron->FindBin(20), hist2_electron->GetNbinsX())));

    canvas2->cd(3);
    hist2_neutron->SetLineColor(kGreen);
    hist2_neutron->Draw();
    TLatex *latex6 = new TLatex();
    latex6->SetNDC();
    latex6->SetTextSize(0.04);
    latex6->DrawLatex(0.2, 0.8, Form("Above 20: %f", hist2_neutron->Integral(hist2_neutron->FindBin(20), hist2_neutron->GetNbinsX())));

    canvas2->cd(4);
    hist2_gamma->SetLineColor(kBlue);
    hist2_gamma->Draw();
    hist2_electron->SetLineColor(kRed);
    hist2_electron->Draw("same");
    hist2_neutron->SetLineColor(kGreen);
    hist2_neutron->Draw("same");


    TCanvas *canvas0 = new TCanvas("canvas0", "Stage0 Branches", 1600, 1200);
    canvas0->Divide(1, 1);

    canvas0->cd(1);
    gPad->SetLogy();
    hist0_gamma->SetLineColor(kBlue);
    hist0_gamma->Draw();
    TLatex *latex7 = new TLatex();
    latex7->SetNDC();
    latex7->SetTextSize(0.04);
    latex7->DrawLatex(0.2, 0.8, Form("Above 200: %f", hist0_gamma->Integral(hist0_gamma->FindBin(200), hist0_gamma->GetNbinsX())));

    TCanvas *canvas3 = new TCanvas("canvas3", "Stage2 Time Distribution", 1600, 1200);
    canvas3->Divide(2, 2);

    canvas3->cd(1);
    hist2_gamma_timedistribution->SetLineColor(kBlue);
    hist2_gamma_timedistribution->Draw();
    TLatex *latex10 = new TLatex();
    latex10->SetNDC();
    latex10->SetTextSize(0.04);
    latex10->DrawLatex(0.2, 0.8, Form("Above 3: %f", hist2_gamma_timedistribution->Integral(hist2_gamma_timedistribution->FindBin(3), hist2_gamma_timedistribution->GetNbinsX())));

    canvas3->cd(2);
    hist2_gamma_timedistribution_cut->SetLineColor(kRed);
    hist2_gamma_timedistribution_cut->Draw();
    TLatex *latex11 = new TLatex();
    latex11->SetNDC();
    latex11->SetTextSize(0.04);
    latex11->DrawLatex(0.2, 0.8, Form("Above 0.15: %f", hist2_gamma_timedistribution_cut->Integral(hist2_gamma_timedistribution_cut->FindBin(0.15), hist2_gamma_timedistribution_cut->GetNbinsX())));

    canvas3->cd(3);
    hist2_invariant_mass->SetLineColor(kGreen);
    hist2_invariant_mass->Draw();
    TLatex *latex8 = new TLatex();
    latex8->SetNDC();
    latex8->SetTextSize(0.04);
    latex8->DrawLatex(0.2, 0.8, Form("Above 50: %f", hist2_invariant_mass->Integral(hist2_invariant_mass->FindBin(50), hist2_invariant_mass->GetNbinsX())));

    canvas3->cd(4);
    hist2_invariant_mass_cut->SetLineColor(kBlack);
    hist2_invariant_mass_cut->Draw();
    TLatex *latex9 = new TLatex();
    latex9->SetNDC();
    latex9->SetTextSize(0.04);
    latex9->DrawLatex(0.2, 0.8, Form("Above 50: %f", hist2_invariant_mass_cut->Integral(hist2_invariant_mass_cut->FindBin(50), hist2_invariant_mass_cut->GetNbinsX())));
    
    TCanvas *canvas4 = new TCanvas("canvas4", "Stage2 Vertex Cut", 1600, 1200);
    canvas4->Divide(2, 2);
    canvas4->cd(1);
    hist2_invariant_mass_cut->SetLineColor(kBlack);
    hist2_invariant_mass_cut->Draw();

    //canvas4->cd(2);
    //hist2_invariant_mass_cut_inside_chamber->SetLineColor(kBlack);
    //hist2_invariant_mass_cut_inside_chamber->Draw();

    canvas4->cd(2);
    hist2_invariant_mass_cut_chamber_steel->SetLineColor(kBlack);
    hist2_invariant_mass_cut_chamber_steel->Draw();
    TLatex *latex12 = new TLatex();
    latex12->SetNDC();
    latex12->SetTextSize(0.04);
    latex12->DrawLatex(0.2, 0.8, Form("Above 20: %f", hist2_invariant_mass_cut_chamber_steel->Integral(hist2_invariant_mass_cut_chamber_steel->FindBin(20), hist2_invariant_mass_cut_chamber_steel->GetNbinsX())));

    canvas4->cd(3);
    hist2_invariant_mass_cut_target->SetLineColor(kBlack);
    hist2_invariant_mass_cut_target->Draw();
    TLatex *latex13 = new TLatex();
    latex13->SetNDC();
    latex13->SetTextSize(0.04);
    latex13->DrawLatex(0.2, 0.8, Form("Above 40: %f", hist2_invariant_mass_cut_target->Integral(hist2_invariant_mass_cut_target->FindBin(40), hist2_invariant_mass_cut_target->GetNbinsX())));

    canvas4->cd(4);
    hist2_invariant_mass_cut_highE->SetLineColor(kBlack);
    hist2_invariant_mass_cut_highE->Draw();
    TLatex *latex14 = new TLatex();
    latex14->SetNDC();
    latex14->SetTextSize(0.04);
    latex14->DrawLatex(0.2, 0.8, Form("Above 10: %f", hist2_invariant_mass_cut_highE->Integral(hist2_invariant_mass_cut_highE->FindBin(10), hist2_invariant_mass_cut_highE->GetNbinsX())));

    
    //canvas4->cd(5);
    //hist2_invariant_mass_cut_outside_chamber->SetLineColor(kBlack);
    //hist2_invariant_mass_cut_outside_chamber->Draw();

    TCanvas *canvas5 = new TCanvas("canvas5", "Stage2 Analyze Gamma Invariant Mass", 1600, 1200);
    canvas5->Divide(2, 3);
    
    canvas5->cd(1);
    hist2_analyze_gamma_invariant_totalmass->SetLineColor(kBlack);
    hist2_analyze_gamma_invariant_totalmass->Draw();
    TLatex *latex15 = new TLatex();
    latex15->SetNDC();
    latex15->SetTextSize(0.04);
    latex15->DrawLatex(0.2, 0.8, Form("Above 10: %f", hist2_analyze_gamma_invariant_totalmass->Integral(hist2_analyze_gamma_invariant_totalmass->FindBin(10), hist2_analyze_gamma_invariant_totalmass->GetNbinsX())));

    canvas5->cd(2);
    hist2_analyze_gamma_invariant_inside_chamber->SetLineColor(kBlack);
    hist2_analyze_gamma_invariant_inside_chamber->Draw();
    
    canvas5->cd(3);
    hist2_analyze_gamma_invariant_outside_chamber->SetLineColor(kBlack);
    hist2_analyze_gamma_invariant_outside_chamber->Draw();

    canvas5->cd(4);
    hist2_analyze_gamma_invariant_chamber_steel->SetLineColor(kBlack);
    hist2_analyze_gamma_invariant_chamber_steel->Draw();
    TLatex *latex16 = new TLatex();
    latex16->SetNDC();
    latex16->SetTextSize(0.04);
    latex16->DrawLatex(0.2, 0.8, Form("Above 10: %f", hist2_analyze_gamma_invariant_chamber_steel->Integral(hist2_analyze_gamma_invariant_chamber_steel->FindBin(10), hist2_analyze_gamma_invariant_chamber_steel->GetNbinsX())));

    canvas5->cd(5);
    hist2_analyze_gamma_invariant_target->SetLineColor(kBlack);
    hist2_analyze_gamma_invariant_target->Draw();
    TLatex *latex17 = new TLatex();
    latex17->SetNDC();
    latex17->SetTextSize(0.04);
    latex17->DrawLatex(0.2, 0.8, Form("Above 10: %f", hist2_analyze_gamma_invariant_target->Integral(hist2_analyze_gamma_invariant_target->FindBin(10), hist2_analyze_gamma_invariant_target->GetNbinsX())));

    TCanvas *canvas6 = new TCanvas("canvas6", "Stage2 Pair Distribution", 1600, 1200);
    canvas6->Divide(2, 3);

    canvas6->cd(1);
    hist2_pair_distance_distribution->SetLineColor(kBlack);
    hist2_pair_distance_distribution->Draw();

    canvas6->cd(2);
    hist2_pair_distance_distribution_cut->SetLineColor(kBlack);
    hist2_pair_distance_distribution_cut->Draw();

    canvas6->cd(3);
    hist2_pair_transverse_distance_distribution->SetLineColor(kBlack);
    hist2_pair_transverse_distance_distribution->Draw();

    canvas6->cd(4);
    hist2_pair_transverse_distance_distribution_cut->SetLineColor(kBlack);
    hist2_pair_transverse_distance_distribution_cut->Draw();

    canvas6->cd(5);
    hist2_pair_z_distance_distribution->SetLineColor(kBlack);
    hist2_pair_z_distance_distribution->Draw();

    canvas6->cd(6);
    hist2_pair_z_distance_distribution_cut->SetLineColor(kBlack);
    hist2_pair_z_distance_distribution_cut->Draw();



    TCanvas *canvas7 = new TCanvas("canvas7", "Stage2 ALPpair Time Distribution", 1600, 1200);
    canvas7->Divide(1, 1);
    canvas7->cd(1);
    hist2_ALPpair_time->SetLineColor(kBlack);
    hist2_ALPpair_time->Draw();
    





    // 이미지로 저장
    canvas1->SaveAs("output/stage1_branch_histograms.png");
    canvas2->SaveAs("output/stage2_branch_histograms.png");
    canvas0->SaveAs("output/stage0_branch_histogram.png");
    canvas3->SaveAs("output/stage2_time_distribution.png");
    canvas4->SaveAs("output/stage2_real_vertex.png");
    canvas5->SaveAs("output/stage2_detector_vertex.png");
    canvas6->SaveAs("output/stage2_distance_distribution.png");
    canvas7->SaveAs("output/stage2_ALPpair_time_distribution.png");
    
    // 파일 닫기
    file->Close();
}

void extract_alp_time() {
    TFile *file = TFile::Open("stage_output.root");
    TTree *tree = (TTree*)file->Get("ALPpair");

    float alp_time;
    tree->SetBranchAddress("ALPpair_time", &alp_time);

    std::ofstream outfile("alp_time_data.txt");

    for (int i = 0; i < tree->GetEntries(); ++i) {
        tree->GetEntry(i);
        outfile << alp_time << std::endl;
    }

    outfile.close();
    file->Close();
}

