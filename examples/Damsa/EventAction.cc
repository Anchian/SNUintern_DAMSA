#include "globals.hh"

#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "SteppingActionDMG4.hh"

#include "DarkMatter.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4UImanager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G4AnalysisManager.hh"
#include <iostream>
G4int ALPpairnumber = 0; 


void ShowProgressBar(G4int currentEventID, G4int totalEvents) 
{
    const int barWidth = 50;  // 진행률 바의 너비
    G4double progress = double(currentEventID + 1) / totalEvents;
    int pos = barWidth * progress;

    std::cout << "[";
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();

    if (currentEventID + 1 == totalEvents) {
        std::cout << std::endl;  // 모든 이벤트가 끝났을 때 줄바꿈
    }
}



EventAction::EventAction(DetectorConstruction* myDC, DarkMatter* DMPointer)
: myDetector(myDC), myDarkMatter(DMPointer), NEmissions(0) 
{;}


EventAction::~EventAction()
{
  G4cout << "Total number of DM emissions = " << NEmissions << G4endl;
  ofstream outFile("Report.txt");
  if(NEmissions >= 3) outFile << "Total number of DM emissions = " << NEmissions << G4endl;
  outFile.close();
}


void EventAction::BeginOfEventAction(const G4Event* event)
{
  theSteppingAction->Reset();
  ALPphotonData.clear();
  ALPphotonList.clear();
  CheckIN.clear();
  pairsNumber = 0;
  myDetector->SetAEmission(0);
  
  //초기화 함수 
  ALPphotonData.clear();
  ALPphotonList.clear();
  CheckIN.clear();
  
    // 매 이벤트 시작 시, 변수들 초기화
  
    // photonPairCount 초기화 제거
}


void EventAction::EndOfEventAction(const G4Event* event)
{
    
  theSteppingAction->Finalize();
   // 총 이벤트 수와 현재 이벤트 ID 가져오기
  G4int totalEvents = G4RunManager::GetRunManager()->GetCurrentRun()->GetNumberOfEventToBeProcessed();
  G4int currentEventID = event->GetEventID();

    // 진행률 바를 터미널에 표시
  ShowProgressBar(currentEventID, totalEvents);
  //G4cout << "number of ALP pairs" << ALPpairnumber << G4endl;
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////



  


  /////////////////////////////// 

  //////이거 다시 짜야 함

  
  /*if (HasPhotonPair()) {
    G4double time1 = GetFirstPhotonTime();
    G4double time2 = GetSecondPhotonTime();
    G4double timeDifference = std::abs(time1 - time2);
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillNtupleFColumn(8, 0, timeDifference/ns);
    analysisManager->AddNtupleRow(8);
    

    auto firstPhotonData = alpPhotonData[photonTrackID1];
    auto secondPhotonData = alpPhotonData[photonTrackID2];

    G4double E1 = std::get<0>(firstPhotonData);  // 첫 번째 광자의 에너지
    G4ThreeVector p1(std::get<1>(firstPhotonData), std::get<2>(firstPhotonData), std::get<3>(firstPhotonData)); // 첫 번째 광자의 운동량

    G4double E2 = std::get<0>(secondPhotonData); // 두 번째 광자의 에너지
    G4ThreeVector p2(std::get<1>(secondPhotonData), std::get<2>(secondPhotonData), std::get<3>(secondPhotonData)); // 두 번째 광자의 운동량


    // invariant mass 계산
    G4double invariantMassSquared = (E1 + E2) * (E1 + E2) - (p1*E1 + p2*E2).mag2();
    G4double invariantMass = (invariantMassSquared > 0) ? std::sqrt(invariantMassSquared) : 0.0;

    analysisManager->FillNtupleFColumn(8, 1, invariantMass/MeV);
    analysisManager->FillNtupleFColumn(8, 2, E1/MeV);
    analysisManager->FillNtupleFColumn(8, 3, E2/MeV);
    analysisManager->FillNtupleFColumn(8, 4, p1.mag()/MeV);
    analysisManager->FillNtupleFColumn(8, 5, p2.mag()/MeV);
    analysisManager->AddNtupleRow(8);

    

    //G4cout << "Time difference between two photons: " << timeDifference << " ns" << G4endl;
    //G4cout << "//////////////////////////////////////////////////////////////// //" << G4endl;
    
    

}
    //G4cout << "check track id  is not empty" << alpTrackIDs.size() << G4endl;*/
}


void EventAction::AddALPPhotonData(G4int trackID, G4int pairedTrackID, G4int pairNumber, 
                                   G4double energy, G4double px, G4double py, G4double pz) {
    ALPphotonData[trackID] = std::make_tuple(pairedTrackID, pairNumber, energy, px, py, pz);
    ALPphotonList[trackID] = pairNumber;
}

void EventAction::SetPairNumber(G4int pairNumber) {
    pairsNumber = pairNumber;
}

G4int EventAction::GetPairNumber() const {
    return pairsNumber;
}

void EventAction::AddCheckIN(G4int pairNumber) {
    if (pairNumber >= CheckIN.size()) {
        CheckIN.resize(pairNumber + 1, 0); // 새로운 쌍에 대해 CheckIN 초기화
    }
    CheckIN[pairNumber] +=1;
}
