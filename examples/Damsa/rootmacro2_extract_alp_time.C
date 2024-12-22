#include <iostream>
#include <fstream>

void rootmacro2_extract_alp_time() {
    TFile *file = TFile::Open("stage_output.root");
    TTree *tree = (TTree*)file->Get("ALPpair");

    float alp_time;
    tree->SetBranchAddress("ALPpair_time", &alp_time);

    // 특정 폴더에 파일 저장
    std::ofstream outfile("output/alp_time_data.txt");

    for (int i = 0; i < tree->GetEntries(); ++i) {
        tree->GetEntry(i);
        outfile << alp_time << std::endl;
    }

    outfile.close();
    file->Close();
}

