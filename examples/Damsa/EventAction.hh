#pragma once

#include "G4UserEventAction.hh"

class DetectorConstruction;
class SteppingActionDMG4;

class DarkMatter;

class G4Event;
#include <vector>
#include <map>
#include "G4ios.hh"

#ifndef EVENTACTION_HH
#define EVENTACTION_HH
extern G4int ALPpairnumber;  // 전역 변수의 선언
#endif


class EventAction : public G4UserEventAction
{
  public:
    EventAction(DetectorConstruction* myDC, DarkMatter* DMPointer);
    ~EventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event* event) override;
    void SetSteppingAction(SteppingActionDMG4* action) {theSteppingAction = action;}
    DarkMatter* GetDarkMatterPointer() {return myDarkMatter;}
    void CountEmission() {NEmissions++;}

    ////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////ALP photon data 저장용 멤버 함수
    // 데이터 저장용 멤버 함수
    void AddALPPhotonData(G4int trackID, G4int pairedTrackID, G4int pairNumber, G4double energy, 
                          G4double px, G4double py, G4double pz);

    void SetPairNumber(G4int pairNumber);
    G4int GetPairNumber() const;
    
    void AddCheckIN(G4int pairNumber);

    // 데이터 구조
    std::map<G4int, std::tuple<G4int, G4int, G4double, G4double, G4double, G4double>> ALPphotonData;
    std::map<G4int, G4int> ALPphotonList;
    std::vector<G4int> CheckIN;

    
    ////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////


  private:
    DetectorConstruction* myDetector;
    SteppingActionDMG4* theSteppingAction;
    DarkMatter* myDarkMatter;
    
    
    G4int pairsNumber;

    G4int NEmissions;
    
};
