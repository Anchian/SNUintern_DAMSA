#include "RunAction.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "G4TrackingManager.hh"
#include "G4SteppingManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G4AnalysisManager.hh"
#include "EventAction.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4ThreeVector.hh"
#include "DetectorConstruction.hh"
#include "G4AccumulableManager.hh"

#include <iostream>
#include <iomanip>
#include <numeric>



RunAction::RunAction()
: G4UserRunAction(){
    
    // SteppingAction을 가져옴
    fSteppingAction = (SteppingActionDMG4*)G4RunManager::GetRunManager()->GetUserSteppingAction();
    // stepping Action 형태의 변수인데 , 그 변수가 지금 실제로 사용되는 steppingaction과 연결해줌 
     // stage2의 root file 
    auto analysisManager = G4AnalysisManager::Instance();
    
    analysisManager->SetDefaultFileType("root");
    analysisManager->SetVerboseLevel(3);
    //analysisManager->SetFileName("stage_output.root");
    // Create a ROOT file
    analysisManager->SetNtupleMerging(true);
    analysisManager->OpenFile("stage_output.root");
    
    //stage0 ntuple 은 Ntuple 2로 설정함


    // Stage 1 Ntuple 생성
    analysisManager->CreateNtuple("Stage1", "Stage 1 Particle Data"); 
    analysisManager->CreateNtupleFColumn("gamma_Energy"); // 에너지
    analysisManager->CreateNtupleFColumn("electron_Energy"); 
    analysisManager->CreateNtupleFColumn("neutron_Energy"); 
    analysisManager->FinishNtuple(0);

    // Stage 2 Ntuple 생성
    analysisManager->CreateNtuple("Stage2", "Stage 2 Particle Data");
    analysisManager->CreateNtupleFColumn("gamma_Energy"); // 에너지
    analysisManager->CreateNtupleFColumn("electron_Energy"); 
    analysisManager->CreateNtupleFColumn("neutron_Energy"); 
    analysisManager->FinishNtuple(1);

    //stage0 ntuple 은 Ntuple 2로 설정함
    analysisManager->CreateNtuple("Stage0","Target inside gamma production");
    analysisManager->CreateNtupleFColumn("gamma_Energy");
    analysisManager->FinishNtuple(2);   

    analysisManager->CreateNtuple("GammaTimeDistributionCut_Pair","Gamma Time Distribution Cut Pairs");
    analysisManager->CreateNtupleFColumn("gamma_time_cut");
    analysisManager->FinishNtuple(3);

    analysisManager->CreateNtuple("GammaTimeDistribution","Gamma Time Distribution");
    analysisManager->CreateNtupleFColumn("gamma_time");
    analysisManager->FinishNtuple(4);

    analysisManager->CreateNtuple("InvariantMass","Invariant Mass Distribution");
    analysisManager->CreateNtupleFColumn("invariant_mass");
    analysisManager->CreateNtupleFColumn("VertexCut_invariant_mass");

    analysisManager->CreateNtupleFColumn("Vertex_Inside_Chamber");//
    analysisManager->CreateNtupleFColumn("Vertex_Chamber_Steel"); 
    analysisManager->CreateNtupleFColumn("Vertex_target");
    analysisManager->CreateNtupleFColumn("Vertex_Outside_Chamber");
    analysisManager->CreateNtupleFColumn("invariant_mass_0");

    analysisManager->CreateNtupleSColumn("Vertex_inside_Chamber_process"); //7
    analysisManager->CreateNtupleSColumn("Vertex_Chamber_Steel_process");
    analysisManager->CreateNtupleSColumn("Vertex_target_process");
    analysisManager->CreateNtupleSColumn("Vertex_Outside_Chamber_process");
    
    
    analysisManager->FinishNtuple(5);

    analysisManager->CreateNtuple("AnalyzeGamma", "Aspect_to_Detector");
    analysisManager->CreateNtupleFColumn("Vertex_0");
    analysisManager->CreateNtupleFColumn("Vertex_Inside_Chamber");//
    analysisManager->CreateNtupleFColumn("Vertex_Chamber_Steel"); 
    analysisManager->CreateNtupleFColumn("Vertex_target");
    analysisManager->CreateNtupleFColumn("Vertex_Outside_Chamber");

    analysisManager->CreateNtupleSColumn("Vertex_inside_Chamber_process");//5
    analysisManager->CreateNtupleSColumn("Vertex_Chamber_Steel_process");
    analysisManager->CreateNtupleSColumn("Vertex_target_process");
    analysisManager->CreateNtupleSColumn("Vertex_Outside_Chamber_process");


    analysisManager->FinishNtuple(6);

    analysisManager->CreateNtuple("Vertex_distance_distribution", "Aspect_Detector");
    analysisManager->CreateNtupleFColumn("distance_distribution" ); //0 
    analysisManager->CreateNtupleFColumn("transverse_distance_distribution");
    analysisManager->CreateNtupleFColumn("z_distance_distribution");
    analysisManager->CreateNtupleFColumn("cut_distance_distribution");//3
    analysisManager->CreateNtupleFColumn("cut_transverse_distance_distribution");
    analysisManager->CreateNtupleFColumn("cut_z_distance_distribution");
    analysisManager->FinishNtuple(7);

    analysisManager->CreateNtuple("ALPpair", "ALPpairtimedistribution"); 
    analysisManager->CreateNtupleFColumn("ALPpair_time"); // 에너지
    analysisManager->CreateNtupleFColumn("ALPpair_invariant_mass");
    analysisManager->CreateNtupleFColumn("E1");
    analysisManager->CreateNtupleFColumn("E2");
    analysisManager->CreateNtupleFColumn("p1");
    analysisManager->CreateNtupleFColumn("p2");
    analysisManager->FinishNtuple(8);

    analysisManager->CreateNtuple("ALPrealphotonEnergy","ALPrealphotonEnergy");
    analysisManager->CreateNtupleFColumn("ALPrealphotonEnergy");
    analysisManager->CreateNtupleFColumn("ALPrealphotonTime");
    analysisManager->CreateNtupleFColumn("ALPrealphotonArrivedEnergy");
    analysisManager->FinishNtuple(9);
    


}

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
    
    // RunManager verbose level 설정
    G4RunManager::GetRunManager()->SetVerboseLevel(0);

    // EventManager verbose level 설정
    G4EventManager::GetEventManager()->SetVerboseLevel(0);

    // TrackingManager verbose level 설정
    G4TrackingManager* trackingManager = G4EventManager::GetEventManager()->GetTrackingManager();
    trackingManager->SetVerboseLevel(0);

    // SteppingManager verbose level 설정
    G4SteppingManager* steppingManager = trackingManager->GetSteppingManager();
    steppingManager->SetVerboseLevel(0);
}


////////진행률 나타내기 ////////////




/////////////////////////////////



RunAction::~RunAction()
{
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->Write();
    analysisManager->CloseFile();
    
}



//void RunAction::BeginOfRunAction(const G4Run* /*run*/)
//{ 
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  
  // Get analysis manager
//  auto analysisManager = G4AnalysisManager::Instance();
//  analysisManager->OpenFile();
//}





std::string GetVolumeNameAtPosition(const G4ThreeVector& position) 
{
    // Navigator 인스턴스를 가져옴
    G4Navigator* navigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();

    // 주어진 좌표에서 물리적 볼륨 찾기
    G4VPhysicalVolume* physicalVolume = navigator->LocateGlobalPointAndSetup(position);

    if (physicalVolume) {
        // 물리적 볼륨에 대응하는 논리적 볼륨의 이름을 반환
        G4LogicalVolume* logicalVolume = physicalVolume->GetLogicalVolume();
        return logicalVolume->GetName();
    } else {
        // 해당 좌표가 어떤 볼륨에도 속하지 않는 경우
        return "The point is outside of any volume!";
    }
}


G4double CalculateInvariantMass(G4double E1, G4double E2, const G4ThreeVector& p1, const G4ThreeVector& p2) {
    G4double totalEnergy = E1 + E2;
    G4ThreeVector totalMomentum = E1*p1 + E2*p2;
    auto analysisManager = G4AnalysisManager::Instance();
    // 불변 질량 계산: M^2 = (E1 + E2)^2 - |p1 + p2|^2
    G4double invariantMassSquared = totalEnergy * totalEnergy - totalMomentum.mag2();
    if (invariantMassSquared < 0) {
        analysisManager->FillNtupleFColumn(5, 6, invariantMassSquared/ MeV);
    }
    return (invariantMassSquared > 0) ? std::sqrt(invariantMassSquared) : 0.0;
}


std::vector<std::vector<G4double>> MinimumDistanceBetweenLines(const G4ThreeVector& p1, const G4ThreeVector& v1, const G4ThreeVector& p2, const G4ThreeVector& v2) 
{
    // 두 방향 벡터의 외적 (cross product)
    G4ThreeVector crossProduct = v1.cross(v2);

    // 직선 사이의 최소 거리 계산 (일반적인 경우)
    G4ThreeVector delta = p2 - p1;

    // 파라미터 t1과 t2 계산 (직선 1과 직선 2에서 각각 가장 가까운 점)
    G4double t1 = (delta.dot(v2) * v1.dot(v2) - delta.dot(v1) * v2.dot(v2)) / (v1.dot(v1) * v2.dot(v2) - pow(v1.dot(v2), 2));
    G4double t2 = (t1 * v1.dot(v2) + delta.dot(v2)) / v2.dot(v2);

    // 가장 가까운 점 계산
    G4ThreeVector closestPoint1 = p1 + t1 * v1;  // 직선 1에서 가장 가까운 점
    G4ThreeVector closestPoint2 = p2 + t2 * v2;  // 직선 2에서 가장 가까운 점
    G4double distance = (closestPoint1 - closestPoint2).mag();

    std::vector<std::vector<G4double>> result = {
        {closestPoint1.x(), closestPoint1.y(), closestPoint1.z()},  // 직선 1에서 가장 가까운 점
        {closestPoint2.x(), closestPoint2.y(), closestPoint2.z()},  // 직선 2에서 가장 가까운 점
        {distance}  // 두 점 사이의 거리
    };

    return result;
}

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{


    
    
    /*
    // 저장된 gamma 에너지 출력
    const std::vector<G4double>& gammaEnergies = fSteppingAction->GetGammaEnergies();
    G4cout << "Gamma energies at the end of the run:" << G4endl;

    for (G4double energy : gammaEnergies) {
        G4cout << energy / MeV << " MeV" << G4endl;
    }

    G4cout << "Total number of gamma particles: " << gammaEnergies.size() << G4endl;
    */
    auto analysisManager = G4AnalysisManager::Instance();
    // stage 1을 지난 모든 입자들의 이름과 에너지를 출력함
    const std::vector<G4double>& stage1particleEnergy = fSteppingAction->GetStage1ParticlesEnergy();
    const std::vector<G4String>& stage1particleName = fSteppingAction->GetStage1ParticlesName();
    
    
//////////////////////////////////////////////



    
    
///////////////////////////////////

    




    /*
    G4cout << "---- Summary of Stage 1 particles recorded in this run ----" << G4endl;
    for (size_t i = 0; i < stage1particleName.size(); ++i) {
        G4cout << "Particle: " << stage1particleName[i] << ", Energy: " 
                << stage1particleEnergy[i] / MeV << " MeV" << G4endl;
    }*/
    G4cout << "-------------------------------------------------" << stage1particleName.size()<< G4endl;

     // stage 2을 지난 모든 입자들의 이름과 에너지를 출력함
    const std::vector<G4double>& stage2particleEnergy = fSteppingAction->GetStage2ParticlesEnergy();
    const std::vector<G4String>& stage2particleName = fSteppingAction->GetStage2ParticlesName();
    /*
    G4cout << "---- Summary of Stage 2 particles recorded in this run ----" << G4endl;
    for (size_t i = 0; i < stage2particleName.size(); ++i) {
        G4cout << "Particle: " << stage2particleName[i] << ", Energy: " 
                << stage2particleEnergy[i] / MeV << " MeV" << G4endl;
    }*/
    G4cout << "-------------------------------------------------"<<stage2particleName.size() << G4endl;
    
    const std::vector<G4double>& stage0particleEnergy = fSteppingAction->GetStage0ParticlesEnergy();
    /*
    G4cout << "---- Summary of Stage 0 particles recorded in this run ----" << G4endl;
    for (size_t i = 0; i < stage0particleEnergy.size(); ++i) {
        G4cout << "Particle: " << "Target" << ", Energy: "
                << stage0particleEnergy[i] /MeV << " MeV"  << G4endl;
    }*/
    G4cout << "-------------------------------------------------"<<stage0particleEnergy.size() << G4endl;
    

    //////// Stage2의 time cut & distribution mass cut &  mass distribution
    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////

    std::vector<std::vector<G4double>> St2TiEnMo = fSteppingAction->GetStage2TimeCut4Vec();
    std::vector<std::vector<G4double>> VertexInfo = fSteppingAction->GetVertexPosition();
    // 시간을 기준으로 정렬 (첫 번째 열이 time, 두 번째가 energy, 세 번째가 momentum)
    if (St2TiEnMo.size() < 6) {
        G4cout << "No gamma pairs found." << G4endl;
        return;} 
    else {
        // 시간을 기준으로 오름차순 정렬
        std::sort(St2TiEnMo.begin(), St2TiEnMo.end(),
        [](const std::vector<G4double>& a, const std::vector<G4double>& b) {
            return a[0] < b[0];  // 첫 번째 열이 time
        });

    int pairCount = 0;
    int vertexCount = 0;
    G4double timeDifference = 0.0;
    G4double E1 = 0.0;
    G4double E2 = 0.0;
    G4ThreeVector p1 = {0.0,0.0,0.0};
    G4ThreeVector p2 = {0.0,0.0,0.0};
    
    G4double p1x = 0.0;
    G4double p1y = 0.0;
    G4double p1z = 0.0;

    G4double p2x = 0.0;
    G4double p2y = 0.0;
    G4double p2z = 0.0;

    G4double invariantMass = 0.0;

    G4ThreeVector VertexPosition(0, 0, 0);
    G4ThreeVector VertexPosition2(0, 0, 0);

    G4ThreeVector x1 = {0.0,0.0,0.0};
    G4ThreeVector x2 = {0.0,0.0,0.0};



    // 시간 차이가 0.1ns 이하인 쌍을 찾기
    size_t totalPairs = St2TiEnMo.size() * (St2TiEnMo.size() - 1) / 2;
    size_t processedPairs = 0;

    for (size_t i = 0; i < St2TiEnMo.size() - 1; ++i) {
        for (size_t j = i + 1; j < St2TiEnMo.size(); ++j) {
            processedPairs++;
            if (processedPairs % 1000 == 0) { // 1000번째 쌍마다 진행률 출력
                G4cout << "Progress: " << (100.0 * processedPairs / totalPairs) << "%\r" << std::flush;
            }

            timeDifference = std::abs(St2TiEnMo[j][0] - St2TiEnMo[i][0]);

            // 시간 차이가 0.1ns 이하인 경우
            if (timeDifference <= 0.1 * ns) {
                pairCount++;

                // Ntuple에 시간 차이 기록
                analysisManager->FillNtupleFColumn(3, 0, timeDifference / ns);
                analysisManager->AddNtupleRow(3);

                // 운동량과 에너지를 이용해 불변 질량 계산
                E1 = St2TiEnMo[i][1];  // 첫 번째 감마 입자의 에너지
                E2 = St2TiEnMo[j][1];  // 두 번째 감마 입자의 에너지

                p1x = St2TiEnMo[i][2];  // 첫 번째 감마 입자의 운동량 (예시)
                p1y = St2TiEnMo[i][3];
                p1z = St2TiEnMo[i][4];
                p1 = {p1x, p1y, p1z};

                p2x = St2TiEnMo[j][2];
                p2y = St2TiEnMo[j][3];
                p2z = St2TiEnMo[j][4];
                p2 = {p2x,p2y,p2z}; // 두 번째 감마 입자의 운동량 (예시)

                x1 = G4ThreeVector(St2TiEnMo[i][5], St2TiEnMo[i][6], St2TiEnMo[i][7]);
                x2 = G4ThreeVector(St2TiEnMo[j][5], St2TiEnMo[j][6], St2TiEnMo[j][7]);

                std::vector<G4String>  process = fSteppingAction->GetPhysicalProcess();
                
                // 불변 질량 계산
                invariantMass = CalculateInvariantMass(E1, E2, p1, p2);
                // 0.1ns 안에 들어온것 중에 invariantmass 가 0보다 크다는 보장이 없기 때문에 
                // timedifference갯수 보다 적음 
                
                // 불변 질량을 Ntuple에 저장
                analysisManager->FillNtupleFColumn(5, 0, invariantMass / MeV);
                analysisManager->AddNtupleRow(5);  // Ntuple에 행 추가
                
                //////////////////////////////////vertex cut 설정 //////////////////////////////////
                

                G4double vertexoptionX = VertexInfo[i][0] - VertexInfo[j][0];
                G4double vertexoptionY = VertexInfo[i][1] - VertexInfo[j][1];
                G4double vertexoptionZ = VertexInfo[i][2] - VertexInfo[j][2];
                G4double vertexoptionXY = std::sqrt(vertexoptionX*vertexoptionX + vertexoptionY*vertexoptionY);
                G4double vertexoptionZZ = std::sqrt(vertexoptionZ*vertexoptionZ);

                

                // 쌍의 vertex point 가 일치하는거 찾아야 하지만 systmatic uncertainty 때문이라고 한다면 
                if (vertexoptionXY < 0.1*cm && vertexoptionZZ < 0.1*cm){
                    analysisManager->FillNtupleFColumn(5, 1, invariantMass / MeV);

                    analysisManager->AddNtupleRow(5);
                    vertexCount++;
                    // 어디에서 생성 된건지 확인 
                    VertexPosition = G4ThreeVector(VertexInfo[i][0], VertexInfo[i][1], VertexInfo[i][2]);
                    VertexPosition2 = G4ThreeVector(VertexInfo[j][0], VertexInfo[j][1], VertexInfo[j][2]);
                    std::string volumeName = GetVolumeNameAtPosition(VertexPosition);
                    std::string volumeName2 = GetVolumeNameAtPosition(VertexPosition2);

                    
                    if ( volumeName == "decaychamberLogical"|| volumeName2 == "decaychamberLogical" )
                    {
                        analysisManager->FillNtupleFColumn(5, 2, invariantMass / MeV);
                        analysisManager->FillNtupleSColumn( 5,7, process[i]+process[j]);
                        analysisManager->AddNtupleRow(5); // e+ e- 예상
                        
                    }
                    else if( volumeName == "decaychamberstLogical" || volumeName2 == "decaychamberstLogical" )
                    {
                        analysisManager->FillNtupleFColumn(5, 3, invariantMass / MeV);
                        analysisManager->FillNtupleSColumn(5,8, process[i]+process[j]);
                        analysisManager->AddNtupleRow(5); // neutron 예상
                        
                    }
                    else if( volumeName == "targetLogical" || volumeName2 == "targetLogical" )
                    {
                        analysisManager->FillNtupleFColumn(5, 4, invariantMass / MeV);
                        analysisManager->FillNtupleSColumn(5,9, process[i]+process[j]);
                        analysisManager->AddNtupleRow(5); // 초기 감마 예상 
                        
                    }
                    else if( volumeName == "worldLogical" || volumeName2 == "worldLogical" )
                    {
                        analysisManager->FillNtupleFColumn(5, 5, invariantMass / MeV);
                        analysisManager->FillNtupleSColumn(5,10, process[i]+process[j]);
                        analysisManager->AddNtupleRow(5); // 몰루 
                        
                    }
                    else{}
                }
                // 쌍의 vertex point가 완전히 일치하는거 
                //그리고 나서 케이스 분리 , vertex 위치가 챔버 내인지 ( stage1에서 온건지 아니면 steel의 neutron에서 온건지)
                // 다른 케이스는 챔버내에서 온건데 이건 나올 수 있는게 전자와 양전자로 부터 나온 감마 
                //그리고 나서 완전히 일치하는 경우의 invmass를 그리고 그다음 각각 케이스에 해당하는 invmass 합해서 비교 
                // invmass를 정확히 분석했는지 확인하면 될듯?

                std::vector<std::vector<G4double>> distance = MinimumDistanceBetweenLines(x1,p1,x2,p2);
                G4double delta_X = std::abs(distance[0][0] - distance[1][0]);
                G4double delta_Y = std::abs(distance[0][1] - distance[1][1]);
                G4double distance_transverse = std::sqrt(delta_X * delta_X + delta_Y * delta_Y);
                G4double distance_z = std::abs(distance[0][2] - distance[1][2]);

                //if (distance[2][0] < 3*cm ) {
                //G4cout << "dist<5cm-------"<<distance[2][0] /cm <<"distance_transverse" <<distance_transverse/cm
                //<<" distance_z"<< distance_z /cm<< G4endl;
                //}
                analysisManager->FillNtupleFColumn(7,0,distance[2][0]/cm);
                analysisManager->FillNtupleFColumn(7, 1, distance_transverse/cm);
                analysisManager->FillNtupleFColumn(7, 2, distance_z/cm);
                analysisManager->AddNtupleRow(7);



                if ( distance_transverse < 0.1*cm && distance_z < 0.1*cm ){
                    analysisManager->FillNtupleFColumn(6, 0, invariantMass / MeV);
                    analysisManager->AddNtupleRow(6);
                    
                    analysisManager->FillNtupleFColumn(7, 3, distance[2][0]/cm);
                    analysisManager->FillNtupleFColumn(7, 4, distance_transverse/cm);
                    analysisManager->FillNtupleFColumn(7, 5, distance_z/cm);

                
                    G4ThreeVector signallike_position1 = {distance[0][0],distance[0][1],distance[0][2]};
                    G4ThreeVector signallike_position2 = {distance[1][0],distance[1][1],distance[1][2]};
                    std::string volumeName1 = GetVolumeNameAtPosition(signallike_position1);
                    std::string volumeName2 = GetVolumeNameAtPosition(signallike_position2);
                    if ( volumeName1  == "decaychamberLogical" && volumeName2 == "decaychamberLogical")
                    {
                        analysisManager->FillNtupleFColumn(6, 1, invariantMass / MeV);
                        analysisManager->FillNtupleSColumn(6,5, process[i]+process[j]);
                        analysisManager->AddNtupleRow(6); // e+ e- 예상
                        
                    }
                    else if( volumeName1 == "decaychamberstLogical" && volumeName2 == "decaychamberstLogical"  )
                    {
                        analysisManager->FillNtupleFColumn(6, 2, invariantMass / MeV);
                        analysisManager->FillNtupleSColumn(6,6, process[i]+process[j]);
                        analysisManager->AddNtupleRow(6); // neutron 예상
                        
                    }
                    else if( volumeName1 == "targetLogical" && volumeName2 == "targetLogical" )
                    {
                        analysisManager->FillNtupleFColumn(6, 3, invariantMass / MeV);
                        analysisManager->FillNtupleSColumn(6,7, process[i]+process[j]);
                        analysisManager->AddNtupleRow(6); // 초기 감마 예상 
                        
                    }
                    else if( volumeName1 == "worldLogical" && volumeName2== "worldLogical"  )
                    {
                        analysisManager->FillNtupleFColumn(6, 4, invariantMass / MeV);
                        analysisManager->FillNtupleSColumn(6,8, process[i]+process[j]);
                        analysisManager->AddNtupleRow(6); // 몰루 
                        
                    }
                    else{}
                }

        
                
                // vertex cut ok


                // pair 중에  
            }
        }
    }
    // EventAction 객체에 접근하기 위해 G4RunManager 사용
    auto eventAction = (EventAction*)G4RunManager::GetRunManager()->GetUserEventAction();
    // EventAction에서 PhotonPair 관련 데이터를 가져옴
    //G4int ALPpair = eventAction->GetPhotonPairCount();


    std::vector<G4double> alpdecayposition = fSteppingAction->GetAlpDecayPosition();
    G4double sum = std::accumulate(alpdecayposition.begin(), alpdecayposition.end(), 0.0);
    G4double average = (alpdecayposition.size() > 0) ? sum / alpdecayposition.size() : 0.0;
    G4cout << "Average of alpdecayposition: " << average << " cm" << G4endl;

    G4double maxDecayPosition = *std::max_element(alpdecayposition.begin(), alpdecayposition.end());
    // 최대값 출력
    G4cout << "Maximum ALP decay position: " << maxDecayPosition << " cm" << G4endl;

    G4cout << "Total pairs found: " << pairCount << G4endl;
    G4cout << "Vertex cut pairs found: " << vertexCount << G4endl;
    G4cout << "ALP pair found " << ALPpairnumber << G4endl;
    //G4cout << "ALP pairs detected: " << ALPpair << G4endl;
    }

    }
















