#pragma once

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include <vector>
#include "G4LogicalVolume.hh"
#include "G4Step.hh"
#include "DetectorConstruction.hh"



class DetectorConstruction;
class EventAction;


class SteppingActionDMG4 : public G4UserSteppingAction
{
  public:
    SteppingActionDMG4(DetectorConstruction* myDC, EventAction* myEA);
    virtual ~SteppingActionDMG4(){};

    virtual void UserSteppingAction(const G4Step*);
    
    virtual void Reset();

    virtual void Finalize();

        std::vector<G4double> fGammaEnergies;
        std::vector<G4String> fStage1ParticlesName;
        std::vector<G4double> fStage1ParticlesEnergy;
        std::vector<G4String> fStage2ParticlesName;
        std::vector<G4double> fStage2ParticlesEnergy;
        std::vector<G4double> fStage2GammaTimeDistribution;

        std::vector<G4double> fStage0ParticlesEnergy;

        std::vector<G4double> fStage2GammaInvariantMassDistribution;

        std::vector<std::vector<G4double>> fStage2TimeCut4Vec;
        std::vector<std::vector<G4double>> fVertexPosition;

        std::vector<G4String> fphysicalprocess;

        std::vector<G4double> fAlpDecayPosition;
        


        const std::vector<G4double>& GetGammaEnergies() const{return fGammaEnergies;}

        //stage 1 에너지와 입자의 이름 저장하는 vector
        const std::vector<G4double>& GetStage1ParticlesEnergy() const{return fStage1ParticlesEnergy; }
        const std::vector<G4String>& GetStage1ParticlesName() const{return fStage1ParticlesName; } 
        // stages2 에너지와 이름 저장하는 vector
        const std::vector<G4double>& GetStage2ParticlesEnergy() const{return fStage2ParticlesEnergy; }
        const std::vector<G4String>& GetStage2ParticlesName() const{return fStage2ParticlesName; }
        const std::vector<G4double>& GetStage2GammaTimeDistribution() const{return fStage2GammaTimeDistribution; }

        //stage 0 에��지 저장하는 vector

        const std::vector<G4double>& GetStage0ParticlesEnergy() const{return fStage0ParticlesEnergy; }

        const std::vector<std::vector<G4double>>& GetStage2TimeCut4Vec() const{return fStage2TimeCut4Vec; }

        const std::vector<std::vector<G4double>>& GetVertexPosition() const{return fVertexPosition; }

        const std::vector<G4String>& GetPhysicalProcess() const{return fphysicalprocess; }
        
        const std::vector<G4double>& GetAlpDecayPosition() const{return fAlpDecayPosition; }
        

  private:

    EventAction* eventAction;
};
