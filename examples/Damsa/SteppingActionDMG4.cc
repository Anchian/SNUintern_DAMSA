#include "globals.hh"
#include "G4ios.hh"

#include "SteppingActionDMG4.hh"
//#include "RunAction.hh"
//#include "RunActionDMG4.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "DarkMatter.hh"
#include "DarkPhotons.hh"
#include "DarkScalars.hh"
#include "ALP.hh"

#include "DMParticles/DMParticle.hh"

#include "G4SteppingManager.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"

#include "G4ProcessManager.hh"
#include "G4ProcessType.hh"
#include "Randomize.hh"

#include "G4ParticleTypes.hh"
#include "G4DynamicParticle.hh"
#include "G4EventManager.hh"
#include "G4TrackVector.hh"
#include "G4SystemOfUnits.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"

#include "G4ParticleTypes.hh"
#include "G4DynamicParticle.hh"
#include "G4EventManager.hh"
#include "G4TrackVector.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"   // 단위 정의
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"
#include "RunAction.hh"
#include <algorithm>
#include <cmath>
#include <string>


SteppingActionDMG4::SteppingActionDMG4(DetectorConstruction* myDC, EventAction* myEA)
: eventAction(myEA)
{
  eventAction->SetSteppingAction(this);
}


void SteppingActionDMG4::UserSteppingAction(const G4Step* step)
{
    auto analysisManager = G4AnalysisManager::Instance();
    G4StepPoint* SPointPreStep = step->GetPreStepPoint();


    // DM process인 경우하는 경우 찾기////////////////////////////////////////////////////////////////////////////////////////////
    if(step->GetPostStepPoint()->GetProcessDefinedStep() != 0) {
    if((step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()).find(string("DMProcess")) != string::npos) {

    eventAction->CountEmission();
    
    G4cout  <<"Dark Matter production at E = " << SPointPreStep->GetKineticEnergy()/MeV << G4endl;

    auto secondaries = step->GetSecondaryInCurrentStep();

    for (const auto* secondary : *secondaries) {
        G4String particleName = secondary->GetDefinition()->GetParticleName();

    // DMParticleALP가 생성된 경우에만 에너지를 출력
        if (particleName == "DMParticleALP") {
        G4double alpEnergy = secondary->GetKineticEnergy() / MeV;
        G4cout << "DMParticleALP generated with energy = " << alpEnergy << " MeV" << G4endl;
        }
    }
    ////////////////////////////////////////////////////////////////

    }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////ALP 생성 과정 추적 

    G4StepPoint* preStepPoint = step->GetPreStepPoint();
    G4StepPoint* postStepPoint = step->GetPostStepPoint();
    G4ThreeVector preStepPos = preStepPoint->GetPosition();
    G4ThreeVector postStepPos = postStepPoint->GetPosition();


    const G4Track* track = step->GetTrack();
    
    // 최초입자가 아닌경우에서 입자가 DM paricle인경우 찾고 거기서 나온 2차입자인 감마 찾기 
    G4int trackID1global ;
    G4int trackID2global ;

    G4int parentID = track->GetParentID();
    G4int alpTrackID = 0;
        if (parentID != 0) {
        // We can now use GetParentID to know the parent ID but not the parent track object itself.
        G4String particleName = track->GetDefinition()->GetParticleName();
        G4double particleE = preStepPoint->GetKineticEnergy();
        // Check if the particle's parent is DMParticleALP
        if (particleName == "DMParticleALP") {
            alpTrackID = track->GetTrackID();
            //eventAction->StoreALPTrackID(alpTrackID); 

            //G4cout << "alpTrackID: " << alpTrackID << G4endl;
            //G4String volumeName = postStepPoint->GetPhysicalVolume()->GetName(); 
            //G4cout << "ALP decay volume: " << volumeName << G4endl;
            //G4StepStatus stepStatus = postStepPoint->GetStepStatus();
            //G4cout << "Step Status: " << stepStatus << G4endl;
            G4ThreeVector decayPosition = postStepPoint->GetPosition();
        
            //G4ThreeVector appearPosition = preStepPoint->GetPosition();
            //G4cout << "ALP appear at:"<< appearPosition.z()<< G4endl;
            double decayPositionZ = (decayPosition.z() + 50)/10;
            G4cout << "ALP decay at:"<<  decayPositionZ  << " cm"<< G4endl;
            fAlpDecayPosition.push_back(decayPositionZ); // ALP붕괴위치가 어디서 인지 확인 
            
            //if (postStepPoint->GetStepStatus() == fStopAndKill) {
            //    G4ThreeVector decayPosition = postStepPoint->GetPosition();
            //    G4String volumeName = postStepPoint->GetPhysicalVolume()->GetName(); 
            //    G4cout << "ALP decay volume: " << volumeName << G4endl;
            //}

            ///// alp particle 생성후 감마입자로 invmass check 
            const std::vector<const G4Track*>* secondaries = step->GetSecondaryInCurrentStep();
            if (secondaries && secondaries->size() == 2) {  // Check for exactly two secondaries
                const G4Track* gamma1 = (*secondaries)[0];
                const G4Track* gamma2 = (*secondaries)[1];
                
                // Ensure both secondaries are gamma particles
                if (gamma1->GetDefinition() == G4Gamma::GammaDefinition() && 
                    gamma2->GetDefinition() == G4Gamma::GammaDefinition()) {

                    G4double E1 = gamma1->GetKineticEnergy();
                    G4double E2 = gamma2->GetKineticEnergy();
                    G4ThreeVector p1 = gamma1->GetMomentum();
                    G4ThreeVector p2 = gamma2->GetMomentum();

                    // Calculate invariant mass
                    G4double totalEnergy = E1 + E2;
                    G4ThreeVector totalMomentum = p1 + p2;
                    G4double invariantMassSquared = totalEnergy * totalEnergy - totalMomentum.mag2();
                    G4double invariantMass = (invariantMassSquared > 0) ? std::sqrt(invariantMassSquared) : 0.0;

                    G4double trackID1 = gamma1->GetTrackID();
                    G4double trackID2 = gamma2->GetTrackID();
                    
                    trackID1global = trackID1;
                    trackID2global = trackID2;
                    
                    //eventAction->AddALPPhotonData(trackID1,E1,p1.x(),p1.y(),p1.z());
                    //eventAction->AddALPPhotonData(trackID2,E2,p2.x(),p2.y(),p2.z());
                    
                    //if (eventAction->IsALPTrack(parentTrackID)) {
                //    G4double photonTime = preStepPoint->GetGlobalTime();
                //    eventAction->RecordPhotonTime(photonTrackID, photonTime);
                //    eventAction->RecordALPphoton(photonTrackID,gammaEnergy,momentum.x(),momentum.y(),momentum.z());

                    


                    // Output ALP particle and secondaries information
                    //G4cout << "ALP particle: " << particleName << ", Energy: " << particleE / MeV << " MeV" << G4endl;
                    //G4cout << "Secondary Particle 1: " << gamma1->GetDefinition()->GetParticleName() 
                    //       << ", Energy: " << E1 / MeV << " MeV" << G4endl;
                    //G4cout << "Secondary Particle 2: " << gamma2->GetDefinition()->GetParticleName() 
                    //       << ", Energy: " << E2 / MeV << " MeV" << G4endl;
                    //G4cout << "Invariant Mass of the two photons: " << invariantMass / MeV << " MeV" << G4endl;
                    
                    
                    G4String volumeName1 = gamma1->GetTouchableHandle()->GetVolume()->GetName();
                    //G4String volumeName2 = gamma2->GetTouchableHandle()->GetVolume()->GetName();

                    G4cout << "Gamma 1 created in volume: " << volumeName1 << G4endl;
                    //G4cout << "Gamma 2 created in volume: " << volumeName2 << G4endl;
                    
                    analysisManager->FillNtupleFColumn(9, 0, E1);
                    analysisManager->AddNtupleRow(9);
                    analysisManager->FillNtupleFColumn(9, 0, E2);
                    analysisManager->AddNtupleRow(9);
                }
            }


        }
        }
    G4int maybe_alpTrackID = track->GetTrackID();
    // Track ID 를 가진 감마입자의 post 위치가 Ecal 이면 출력하게 만들기 
    if ( maybe_alpTrackID == trackID1global  ||  maybe_alpTrackID == trackID2global ) {
        G4cout<< "trackID matches"<< G4endl;
        G4String ALPvolumename = preStepPoint->GetTouchableHandle()->GetVolume()->GetName();
        if ( ALPvolumename == "EcalPhysical"){
            G4cout << "There is a signal " << G4endl;
        }
    }
    


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// track ID 이용하여 감마입자 저장
    

    G4int pairsnumber = 0; // 이번호는 eventaction에서 관리 해야겠다. 
    //G4vector<G4int> CheckIN;
    //G4vector<std::vector<G4double>> ALPphotonData; //이 안에 자기 trackID, 같이 생긴 TrackID,  에너지 , 운동량 저장
    //G4vector<G4double> ALPphotonlist; // 이 안에 자기 track IDm, pair number 저장
    G4String particleName = track->GetDefinition()->GetParticleName();

    if (particleName == "DMParticleALP") {
        alpTrackID = track->GetTrackID();
        const std::vector<const G4Track*>* secondaries = step->GetSecondaryInCurrentStep();
            if (secondaries && secondaries->size() == 2) {  // Check for exactly two secondaries
                
                
                const G4Track* gamma1 = (*secondaries)[0];
                const G4Track* gamma2 = (*secondaries)[1];
                // Ensure both secondaries are gamma particles
                if (gamma1->GetDefinition() == G4Gamma::GammaDefinition() && 
                    gamma2->GetDefinition() == G4Gamma::GammaDefinition()) {

                    G4double E1 = gamma1->GetKineticEnergy();
                    G4double E2 = gamma2->GetKineticEnergy();
                    G4ThreeVector p1 = gamma1->GetMomentum();
                    G4ThreeVector p2 = gamma2->GetMomentum();
                    
                    G4double p1x = p1.x();
                    G4double p1y = p1.y();
                    G4double p1z = p1.z();
                    G4double p2x = p2.x();
                    G4double p2y = p2.y();
                    G4double p2z = p2.z();

                    G4int trackID1 = gamma1->GetTrackID();
                    G4int trackID2 = gamma2->GetTrackID();


                    pairsnumber = eventAction->GetPairNumber();
                    pairsnumber++; // 이전에 쌍 번호가 현재 뭔지 불러오고 1 더해서 인덱스 증가

                    eventAction-> SetPairNumber(pairsnumber);
                    //CheckIN[pairsnumber] = 0; // 쌍 이름으로 0으로 생성
                    //ALPphotonData[trackID1] = (trackID2,pairsnumber,E1,p1.x(),p1.y(),p1.z());
                    //ALPphotonData[trackID2] = (trackID1,pairsnumber,E2,p2.x(),p2.y(),p2.z());
                    //ALPphotonlist[trackID1] = pairsnumber;
                    //ALPphotonlist[trackID2] = pairsnumber;
                    eventAction->AddCheckIN(pairsnumber);

                    eventAction->AddALPPhotonData(trackID1, trackID2, pairsnumber, E1,p1x,p1y,p1z );
                    eventAction->AddALPPhotonData(trackID2, trackID1, pairsnumber, E2, p2x,p2y,p2z);
                    //이걸 전부다 eventaction에 저장해야함 
                    }

            }        

        

    }


    /////////////////////////// 감마입자가 ecal에 들어간 경우 track ID 찾아보고 체크인 하기 
    if(step->IsFirstStepInVolume() && step->GetPreStepPoint()->GetTouchable()->GetVolume()->GetName() == "EcalPhysical")
    {
        const G4Track* track = step->GetTrack();
        if (particleName == "gamma" ){
            G4int trackID = track->GetTrackID();
            auto it = eventAction->ALPphotonList.find(trackID);
            if ( it != eventAction->ALPphotonList.end() ) {
                G4int pairnumber = eventAction->ALPphotonList[trackID];
                eventAction->CheckIN[pairnumber] += 1;
                if (eventAction-> CheckIN[pairnumber] == 2){
                    ALPpairnumber += 1; // eventaction 에서 초기값 0으로 설정함 
                }
        }
        }
    }












        /*
        G4ParticleDefinition* particleType = track->GetDefinition();
        G4int particleID = particleType->GetPDGEncoding();   
        G4String volumeName = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName();
        if (volumeName == "EcalPhysical" && particleID == 22) {
            G4int trackID = track->GetTrackID();
            G4cout << "trackID: " << trackID << G4endl;
            
            if (trackID==alpTrackID ){
                G4cout << " //////////////////////////////////" << G4endl;
            }
            if (eventAction->IsALPPhotonTrack(trackID)) {
                G4cout << "////////////////////////////////////" << G4endl;
                G4double photonTime = preStepPoint->GetGlobalTime();
                eventAction->RecordPhotonTime(trackID, photonTime);
                G4cout << "photon detected at Ecal" << G4endl;
                G4double arrivedEnergy = preStepPoint->GetKineticEnergy();
                analysisManager->FillNtupleFColumn(9, 1, photonTime);
                analysisManager->FillNtupleFColumn(9, 2, arrivedEnergy);
                analysisManager->AddNtupleRow(9);
            }
        
        }

    // ALP 에서 붕괴한 감마입자에 대한 정보     
    if (particleID == 22) {
    G4int currentTrackID = track->GetTrackID();
    // Check if the current track ID is in the ALP photon data
    if (eventAction->HasALPPhotonData(currentTrackID)) {
        G4cout << "///////////////////////////////" << G4endl;
        G4StepPoint* postStepPoint = step->GetPostStepPoint();
        G4String volumeName = postStepPoint->GetPhysicalVolume()->GetName();
        
        G4ThreeVector position = postStepPoint->GetPosition();
        G4cout << "Gamma with Track ID " << currentTrackID
               << " in volume: " << volumeName << G4endl;
        }
    }
    

    
    */
    
    // 두 감마 입자의 트랙 ID를 확인하여 맞는 경우에만 출력
    
    
    
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// stage  조건///////////////////////////////////////////////////////////////////////////////


    if (step->IsFirstStepInVolume()) {
        // 현재 볼륨 이름 가져오기
        G4String volumeName = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName();
        const G4Track* track = step->GetTrack();
        // 타겟 볼륨 내에서 첫 번째 스텝인 경우
        if (volumeName =="targetPhysical" && track->GetParticleDefinition() == G4Gamma::GammaDefinition()) {
            // 감마 입자의 운동 에너지 가져오기
            G4double gammaEnergy = preStepPoint->GetKineticEnergy()/MeV;

            // 감지된 입자의 에너지를 저장 (switch 문 활용)

            fStage0ParticlesEnergy.push_back(preStepPoint->GetKineticEnergy());
            // 에너지 저장 (히스토그램 또는 Ntuple 등으로 저장)
            // 예시: Ntuple의 첫 번째 컬럼에 에너지를 저장
            analysisManager->FillNtupleFColumn(2, 0, gammaEnergy);
            analysisManager->AddNtupleRow(2);

            ////////////////////////////////////////////////////////////////
        }    



    // prestepPos의 z 좌표가 5cm 보다작고 postStepPos의 z좌표가 5cm 보다 큰 경우
    if (preStepPos.z() < 5.0*cm && postStepPos.z() > 5.0*cm) {
        const G4Track* track = step->GetTrack();
        G4ParticleDefinition* particleType = track->GetDefinition();
        // 입자의 이름과 에너지 출력
        //G4cout << " stage 1 Particle : " << particleType->GetParticleName() << "Energy" << preStepPoint->GetKineticEnergy() / MeV << "Mev" << G4endl;
        fStage1ParticlesEnergy.push_back(preStepPoint->GetKineticEnergy());
        fStage1ParticlesName.push_back(particleType->GetParticleName());




        //입자의 이름과 에너지를 rootfile에 저장
        G4int particleID = particleType->GetPDGEncoding();
        
        //////////////////////////////////////////////////////////////////////////// Stage 1 데이터를 저장하는 경우
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        G4double gammaEnergy = 0.0;
        G4double electronEnergy = 0.0;
        G4double neutronEnergy = 0.0;

        // 감지된 입자의 에너지를 저장 (switch 문 활용)
        switch (particleID) {
            case 22:  // gamma
                gammaEnergy = postStepPoint->GetKineticEnergy()/MeV;
                break;
            case 11:  // electron
                electronEnergy = postStepPoint->GetKineticEnergy()/MeV;
                break;
            case 2112:  // neutron
                neutronEnergy = postStepPoint->GetKineticEnergy()/MeV;
                break;
        }

        // 모든 Column에 값을 설정하고 나서 한 번에 Row를 추가
        analysisManager->FillNtupleFColumn(0, 0, gammaEnergy);
        analysisManager->FillNtupleFColumn(0, 1, electronEnergy);
        analysisManager->FillNtupleFColumn(0, 2, neutronEnergy);
        analysisManager->AddNtupleRow(0);
        
        

    }
    
    

    // stage 2 조건///////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // PreStepPoint가 Stage 2에 해당하지 않고, PostStepPoint가 Stage 2에 해당하는 경우에만 실행
    if(step->IsFirstStepInVolume() && step->GetPreStepPoint()->GetTouchable()->GetVolume()->GetName() == "EcalPhysical")
    {
        const G4Track* track = step->GetTrack();
        G4ParticleDefinition* particleType = track->GetDefinition();

        //G4cout << " stage 2 Particle : " << particleType->GetParticleName() << "Energy" << preStepPoint->GetKineticEnergy() / MeV << "Mev" << G4endl;
        fStage2ParticlesEnergy.push_back(preStepPoint->GetKineticEnergy());
        fStage2ParticlesName.push_back(particleType->GetParticleName()); 

         

        //입자의 이름과 에너지를 rootfile에 저장
        G4int particleID = particleType->GetPDGEncoding();
        // Stage 1 데이터를 저장하는 경우
        G4double gammaEnergy = 0.0;
        G4double electronEnergy = 0.0;
        G4double neutronEnergy = 0.0;
        std::vector<G4double> row = {0.0, 0.0, 0.0,0.0,0.0};
        std::vector<G4double> vertexrow = {0.0, 0.0, 0.0};
        G4double globalTime = 0.0;
        G4ThreeVector momentum = {0.0, 0.0, 0.0};
        
        G4ThreeVector vertexposition = {0.0, 0.0, 0.0};

        G4ThreeVector arriveposition = {0.0, 0.0, 0.0};
        
        G4String process =" ";
        const G4VProcess* creatorProcess;

        G4int photonTrackID = 0;
        G4int parentTrackID = 0;

        // 감지된 입자의 에너지를 저장 (switch 문 활용)
        switch (particleID) {
            case 22:  // gamma
                gammaEnergy = preStepPoint->GetKineticEnergy() / MeV;
                globalTime = preStepPoint->GetGlobalTime() / ns;  // 나노초(ns)
                momentum = preStepPoint->GetMomentumDirection() / MeV;
                vertexposition = track->GetVertexPosition();
                arriveposition = preStepPoint->GetPosition();   
                
                creatorProcess = track->GetCreatorProcess();
                process = creatorProcess->GetProcessName(); 

                photonTrackID = track->GetTrackID();
                parentTrackID = track->GetParentID();

                row = {globalTime, gammaEnergy, momentum.x(), momentum.y(), momentum.z() ,arriveposition.x(), arriveposition.y(),arriveposition.z() };
                fStage2TimeCut4Vec.push_back(row);

                vertexrow = {vertexposition.x(), vertexposition.y(), vertexposition.z()};
                fVertexPosition.push_back(vertexrow);
                
                fphysicalprocess.push_back(process);
                
                //if (eventAction->IsALPTrack(parentTrackID)) {
                //    G4double photonTime = preStepPoint->GetGlobalTime();
                //    eventAction->RecordPhotonTime(photonTrackID, photonTime);
                //    eventAction->RecordALPphoton(photonTrackID,gammaEnergy,momentum.x(),momentum.y(),momentum.z());
    



                break;
            case 11:  // electron
                electronEnergy = preStepPoint->GetKineticEnergy() / MeV;
                break;
            case 2112:  // neutron
                neutronEnergy = preStepPoint->GetKineticEnergy() / MeV;
                break;
        }
        


        
        // gamma invariant mass distribution
        
            
        



        // 모든 Column에 값을 설정하고 나서 한 번에 Row를 추가
        analysisManager->FillNtupleFColumn(1, 0, gammaEnergy);
        analysisManager->FillNtupleFColumn(1, 1, electronEnergy);
        analysisManager->FillNtupleFColumn(1, 2, neutronEnergy);
        analysisManager->AddNtupleRow(1);

        analysisManager->FillNtupleFColumn(4, 0, globalTime);
        analysisManager->AddNtupleRow(4);




  
}
}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 초기화 함수 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SteppingActionDMG4::Reset()
{
  eventAction->GetDarkMatterPointer()->ResetNEmissions();
}


void SteppingActionDMG4::Finalize()
{
}



/*
// 딸입자 판별기 
    G4String parentName = step->GetTrack()->GetDefinition()->GetParticleName();

    if (parentName=="DMParticleALP") {
        const std::vector<const G4Track*>* secondaries = step->GetSecondaryInCurrentStep();

    // 2. 자식 입자들이 존재하는지 확인
    if (secondaries->size() > 0) {
        G4int nSecondaryParticles = secondaries->size(); // 생성된 자식 입자 수
        G4cout << "Number of secondary particles generated in this step: " 
               << nSecondaryParticles << G4endl;
        
        // 3. 자식 입자들의 이름 출력
        for (size_t i = 0; i < secondaries->size(); i++) {
            G4String secondaryName = (*secondaries)[i]->GetDefinition()->GetParticleName();
            G4cout << "Secondary particle " << i + 1 << ": " << secondaryName << G4endl;
        }
    } else {
        G4cout << "No secondary particles generated in this step." << G4endl;
    }
        }

*/