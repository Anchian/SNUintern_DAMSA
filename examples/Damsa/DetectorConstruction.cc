#include "DetectorConstruction.hh"

#include "G4Element.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GeometryManager.hh"
#include "G4TransportationManager.hh"
#include "G4UserLimits.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4GenericMessenger.hh"
#include "G4NistManager.hh"
#include "G4ios.hh"


DetectorConstruction::DetectorConstruction()
{
}


DetectorConstruction::~DetectorConstruction()
{
}


G4VPhysicalVolume* DetectorConstruction::Construct()
{

  G4double a, iz, z, density;
  G4String name, symbol;
  G4int nel;

// Elements: ----------------

  a = 14.01*g/mole;
  G4Element* elN = new G4Element(name="Nitrogen", symbol="N", iz=7., a);

  a = 16.00*g/mole;
  G4Element* elO = new G4Element(name="Oxygen", symbol="O", iz=8., a);

// Materials: ----------------

    ConstructMaterials();
    G4Material* vaccum = G4Material::GetMaterial("G4_Galactic");
    G4Material* tungsten = G4Material::GetMaterial("G4_W");
    G4Material* air = G4Material::GetMaterial("G4_AIR");
    G4Material* steel = G4Material::GetMaterial("Steel");
    G4bool checkOverlaps = true;
  // Vacuum
    a = 1.*g/mole;
    density = 1.e-15*g/cm3;
    GeneralVacuum = new G4Material(name="GeneralVacuum",z=1.,a,density);

//--------- Sizes of the principal geometrical components (solids)  ---------

//--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------

  //------------------------------
  // World
  //------------------------------
    G4VSolid* worldSolid 
      = new G4Box("worldBox",5.*m,5.*m,5.*m);
    G4LogicalVolume* worldLogical
      = new G4LogicalVolume(worldSolid,air,"worldLogical");
    G4VPhysicalVolume* worldPhysical
      = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"worldPhysical",0,
                          false,0,checkOverlaps);

    G4VSolid* target 
      = new G4Box("targetBox",5.*1/2*cm,5.*1/2*cm,10.*1/2*cm);
    G4LogicalVolume* targetlogical
      = new G4LogicalVolume(target,tungsten,"targetLogical");
    G4VPhysicalVolume* targetphysical =new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),targetlogical,
    "targetPhysical",worldLogical,false,0,checkOverlaps);

    G4VSolid* decaychamber 
      = new G4Tubs("decaychamber",0.,10.*cm,30*1/2*cm,0.*deg,360.*deg);
    G4LogicalVolume* decaychamberLogical
      = new G4LogicalVolume(decaychamber,GeneralVacuum,"decaychamberLogical");
    G4VPhysicalVolume*decaychambervolume = new G4PVPlacement(0,G4ThreeVector(0.,0.,20.*cm),decaychamberLogical,
    "chamber1Physical",worldLogical,false,0,checkOverlaps);
    
    
    G4VSolid* decaychamberst 
    = new G4Tubs("decaychamberst",10.*cm,10.6*cm,30*1/2*cm,0.*deg,360.*deg);
    G4LogicalVolume* decaychamberstLogical
    = new G4LogicalVolume(decaychamberst,steel,"decaychamberstLogical"); 
    G4VPhysicalVolume*decaychamberstvolume = new G4PVPlacement(0,G4ThreeVector(0.,0.,20.*cm),decaychamberstLogical,
    "chamberst1Physical",worldLogical,false,0,checkOverlaps);
    


    G4VSolid* Ecal
      = new G4Box("Ecal",12.*1/2*cm,12.*1/2*cm,44.*1/2*cm);
    G4LogicalVolume * Ecallogical
      = new G4LogicalVolume(Ecal,vaccum,"Ecallogical");
    G4VPhysicalVolume* Ecalphysical =new G4PVPlacement(0,G4ThreeVector(0.,0.,57.*cm),Ecallogical,
    "EcalPhysical",worldLogical,false,0,checkOverlaps);
    

  
  return worldPhysical;
}

void DetectorConstruction::ConstructMaterials()
{
    G4NistManager* nistManager = G4NistManager::Instance();

    // Air 
    nistManager->FindOrBuildMaterial("G4_AIR");
  
    // Argon gas
    nistManager->FindOrBuildMaterial("G4_W");

    // Scintillator
    // (PolyVinylToluene, C_9H_10)
    nistManager->FindOrBuildMaterial("G4_Galactic");

    // 철(Fe), 크롬(Cr), 니켈(Ni) 원소 정의
    G4Material* Fe = nistManager->FindOrBuildMaterial("G4_Fe");  // 철
    G4Material* Cr = nistManager->FindOrBuildMaterial("G4_Cr");  // 크롬
    G4Material* Ni = nistManager->FindOrBuildMaterial("G4_Ni");  // 니켈

// 강철 (SS304) 비율로 합금 정의
    G4double density = 8.00*g/cm3;  // SS304 밀도
    G4Material* Steel = new G4Material("Steel", density, 3);
    Steel->AddMaterial(Fe, 0.70);  // 철 70%
    Steel->AddMaterial(Cr, 0.18);  // 크롬 18%
    Steel->AddMaterial(Ni, 0.12);  // 니켈 12%
    


    G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}