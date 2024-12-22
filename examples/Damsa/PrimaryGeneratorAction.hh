#pragma once
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

#include "G4ios.hh"

class G4ParticleGun;
class G4ParticleDefinition;
class DetectorConstruction;
class G4GeneralParticleSource;
class G4VPrimaryGenerator;
class G4Event;


class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction(DetectorConstruction* myDC);    
    ~PrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);

  private:
    G4GeneralParticleSource* particleGun;
    DetectorConstruction* myDetector;

    void DefineCommands();

    G4ParticleGun* fParticleGun;
    G4ParticleDefinition* fPositron;
    G4ParticleDefinition* fMuon;
    G4ParticleDefinition* fPion;
    G4ParticleDefinition* fKaon;
    G4ParticleDefinition* fProton;
};
