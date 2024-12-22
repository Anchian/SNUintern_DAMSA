#include "globals.hh"
#include "Randomize.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "EventAction.hh"
#include "SteppingActionDMG4.hh"
#include "ActionInitialization.hh"
#include "RunAction.hh"

#include "QGSP_BERT.hh"
#include "FTFP_BERT.hh"

#include "DarkMatterPhysics.hh"
#include "DarkMatter.hh"
#include "G4StepLimiterPhysics.hh"

#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4RunManagerFactory.hh"
#include "G4RunManager.hh"
#include "G4PhysListFactory.hh"
#include "G4MTRunManager.hh" // 멀티스레드 런 매니저


#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#include "G4ios.hh"
#include <chrono> // 시간 측정 라이브러리 포함
#include <iostream>

int main(int argc,char** argv) {
  auto start = std::chrono::high_resolution_clock::now();
  
  // 난수 엔진 초기화 및 시드 설정
    CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
    G4long seed = time(nullptr);  // 현재 시간 기반 시드
    CLHEP::HepRandom::setTheSeed(seed);
  // Run manager
  G4RunManager * runManager = new G4RunManager;
  //auto* runManager = G4RunManagerFactory::CreateRunManager();
  // UserInitialization classes
  DetectorConstruction* mkexp = new DetectorConstruction;
  runManager->SetUserInitialization(mkexp);
  runManager->SetNumberOfThreads(4);
  // ___ Here the "extension" part starts ___
  G4PhysListFactory factory;
  G4VModularPhysicsList * phys = factory.GetReferencePhysList("FTFP_BERT");
  // ^^^ most of the standard physics lists are available by this interface

//  G4PhysicsListHelper * phLHelper = G4PhysicsListHelper::GetPhysicsListHelper();
//  phLHelper->DumpOrdingParameterTable();

  DarkMatterPhysics* myPhysics = new DarkMatterPhysics();
  phys->RegisterPhysics(myPhysics);
  // ^^^ Here the "extension" part ends ^^^

  runManager->SetUserInitialization(phys);  // init phys


#ifdef G4VIS_USE
  // Visualization, if you choose to have it!
  G4VisManager* visManager = new G4VisExecutive;

  visManager->Initialize();
#endif
   
  // UserAction classes
  //runManager->SetUserAction(runAction);
  runManager->SetUserAction(new PrimaryGeneratorAction(mkexp));

  EventAction* myEA = new EventAction(mkexp, myPhysics->GetDarkMatterPointer());
  runManager->SetUserAction(myEA);
  runManager->SetUserAction(new SteppingActionDMG4(mkexp, myEA));
  runManager->SetUserInitialization(new ActionInitialization());
  // User interactions
  G4UImanager * UI = G4UImanager::GetUIpointer();  

  if(argc==1)
  // Define (G)UI terminal for interactive mode  
  { 
    // G4UIterminal is a (dumb) terminal.
    G4UIsession * session = new G4UIterminal;
    UI->ApplyCommand("/control/execute prerun.g4mac");    
    session->SessionStart();
    delete session;
  }
  else
  // Batch mode
  { 
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UI->ApplyCommand(command+fileName);
  }

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;
    auto end = std::chrono::high_resolution_clock::now();

    // 소요된 시간 계산 (밀리초 단위)
    auto duration_sec = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    std::cout << "Simulation run time: " << duration_sec.count() << " seconds " << std::endl;

  return 0;
}
