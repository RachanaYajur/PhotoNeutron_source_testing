//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4RotationMatrix.hh"
#include "G4tgbRotationMatrix.hh"

#include "G4SubtractionSolid.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),
 fAbsorMaterial(0), fLAbsor(0), fContainMaterial(0), fLContain(0),
 fWorldMaterial(0), fPWorld(0), fDetectorMessenger(0)
{
  //fAbsorRadius = 15*mm;
  fAbsorRadius = 5*mm;
  fAbsorLength = 4.5*mm;
  fContainThickness = 150*mm;
  DefineMaterials();
  SetAbsorMaterial  ("BeO");
  SetContainMaterial("Vacuum");
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ delete fDetectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  G4int ncomponents, natoms;
  
  G4Element* Be = new G4Element("Beryllium","Be" ,  4.,  9.01*g/mole);
  G4Element* N  = new G4Element("Nitrogen" ,"N"  ,  7., 14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen"   ,"O"  ,  8., 16.00*g/mole);
  G4Element* Cr = new G4Element("Chromium" ,"Cr" , 24., 51.99*g/mole);
  G4Element* Fe = new G4Element("Iron"     ,"Fe" , 26., 55.84*g/mole);
  G4Element* Ni = new G4Element("Nickel"   ,"Ni" , 28., 58.69*g/mole);
    


  G4Material* BeO = 
  new G4Material("BeO", 3.05*g/cm3, ncomponents=2);
  BeO->AddElement(Be, natoms=1);
  BeO->AddElement( O, natoms=1);
  
 // G4Material* Pb =
 // new G4Material("Lead",


  G4Material* inox = 
  new G4Material("Stainless-Steel", 8*g/cm3, ncomponents=3);
  inox->AddElement(Fe, 74*perCent);
  inox->AddElement(Cr, 18*perCent);
  inox->AddElement(Ni,  8*perCent);

  G4Material* Air = 
  new G4Material("Air", 1.290*mg/cm3, ncomponents=2);
  Air->AddElement(N, 70.*perCent);
  Air->AddElement(O, 30.*perCent);

 
  

  G4double atomicNumber = 1.;   
  G4double massOfMole = 1.008*g/mole;   
  G4double density = 1.e-25*g/cm3; 
  G4double temperature = 2.73*kelvin; 
  G4double pressure = 3.e-18*pascal; 
  G4Material* Vacuum =   new G4Material("Vacuum", atomicNumber, massOfMole, density, kStateGas,temperature, pressure);
 

  fWorldMaterial = Vacuum;

  G4NistManager* nist = G4NistManager::Instance();
  Sb  = nist->FindOrBuildMaterial("G4_Sb");
  Pb  = nist->FindOrBuildMaterial("G4_Pb");
  Epoxy = new G4Material("Epoxy", 0.9828*g/cm3, 4);
  G4Isotope* H1 = new G4Isotope("H1", 1, 1);
  G4Element* elH = new G4Element("Hydrogen", "H", 1);
  elH->AddIsotope(H1, 100. * perCent);
  G4Isotope* N14 = new G4Isotope("N14", 7, 14);
  G4Element* elN = new G4Element("Nitrogen", "N", 1);
  elN->AddIsotope(N14, 100. * perCent);
  G4Isotope* O16 = new G4Isotope("O16", 8, 16);
  G4Element* elO = new G4Element("Oxygen", "O", 1);
  elO->AddIsotope(O16, 100. * perCent);



  G4Element* elC = nist->FindOrBuildElement("C", true);

  Epoxy->AddElement(elH, 0.0797);
  Epoxy->AddElement(elC, 0.7016);
  Epoxy->AddElement(elN, 0.0142);
  Epoxy->AddElement(elO, 0.2045);




 
  ///G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::MaterialWithSingleIsotope( G4String name,
                           G4String symbol, G4double density, G4int Z, G4int A)
{
 // define a material from an isotope
 //
 G4int ncomponents;
 G4double abundance, massfraction;

 G4Isotope* isotope = new G4Isotope(symbol, Z, A);
 
 G4Element* element  = new G4Element(name, symbol, ncomponents=1);
 element->AddIsotope(isotope, abundance= 100.*perCent);
 
 G4Material* material = new G4Material(name, density, ncomponents=1);
 material->AddElement(element, massfraction=100.*perCent);

 return material;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  
  // compute dimensions
  G4double ContainRadius = fAbsorRadius + fContainThickness;
  G4double ContainLength = fContainThickness;
  
  G4double WorldSizeXY = 2.4*ContainRadius;
  G4double WorldSizeZ  = 1.2*ContainLength;


  // Rotation matrices
  //
 

  std::vector<G4double> rID{ 0    ,  0    , 0};
  std::vector<G4double> rX{ 90*deg,  0    , 0};
  std::vector<G4double> rY{  0    , 90*deg, 0};
  G4tgbRotationMatrix rot;
  G4RotationMatrix* rotID = rot.BuildG4RotMatrixFrom3(rID);
  G4RotationMatrix* rotX  = rot.BuildG4RotMatrixFrom3(rX );
  G4RotationMatrix* rotY  = rot.BuildG4RotMatrixFrom3(rY );



  
  // World
  //
  G4Box*
  sWorld = new G4Box("World",                                    //name
              5*m,5*m,5*m);   //dimensions
                   
  G4LogicalVolume*
  lWorld = new G4LogicalVolume(sWorld,                  //shape
                             fWorldMaterial,            //material
                             "World");                  //name

  fPWorld = new G4PVPlacement(0,                        //no rotation
                            G4ThreeVector(),            //at (0,0,0)
                            lWorld,                     //logical volume
                            "World",                    //name
                            0,                          //mother volume
                            false,                      //no boolean operation
                            0);                         //copy number

  // Container
  //
 // G4Tubs* 
 // sContain = new G4Tubs("Container",                            //name
  //           0., ContainRadius, 0.5*ContainLength, 0., twopi);  //dimensions


  G4Box* sContain = new G4Box("Container",ContainLength,ContainLength,ContainLength);
  fLContain = new G4LogicalVolume(sContain,            //shape
                       fContainMaterial,               //material
                       fContainMaterial->GetName());   //name

           new G4PVPlacement(0,                        //no rotation
                       G4ThreeVector(),                //at (0,0,0)
                       fLContain,                      //logical volume
                       fContainMaterial->GetName(),    //name
                       lWorld,                         //mother  volume
                       false,                          //no boolean operation
                       0);                             //copy number

  // Absorber
  //
  G4Tubs* 
  sAbsor = new G4Tubs("Absorber",                                //name
               0., fAbsorRadius, 0.5*fAbsorLength, 0., twopi);    //dimensions

  fLAbsor = new G4LogicalVolume(sAbsor,                //shape
                       fAbsorMaterial,                 //material
                       fAbsorMaterial->GetName());     //name

           new G4PVPlacement(0,                        //no rotation
                       G4ThreeVector(),                //at (0,0,0)
                       fLAbsor,                        //logical volume
                       fAbsorMaterial->GetName(),      //name
                       fLContain,                      //mother  volume
                       false,                          //no boolean operation
                       0);                             //copy number

/*

//  Lead Shielding


 G4Box* solidLeadFull         = new G4Box("LeadFull"        , 127. *mm, 101.6*mm, 177.8*mm);
  G4Box* solidLeadInHole       = new G4Box("LeadInHole"      ,  50.8*mm,  50.8*mm,  25.4*mm);
  G4Box* solidLeadOutFrontHole = new G4Box("LeadOutFrontHole",  25.5*mm, 102. *mm,  89. *mm);
  G4Box* solidLeadOutBackHole  = new G4Box("LeadOutBackHole" ,  25.5*mm, 102. *mm,  38.2*mm);
  G4SubtractionSolid* solidLead = new G4SubtractionSolid("Lead", solidLeadFull, solidLeadInHole      , 0, G4ThreeVector( -25.4*mm, 0,  50.8*mm));
  solidLead                     = new G4SubtractionSolid("Lead", solidLead    , solidLeadOutFrontHole, 0, G4ThreeVector(-101.6*mm, 0, -89. *mm));
  solidLead                     = new G4SubtractionSolid("Lead", solidLead    , solidLeadOutBackHole , 0, G4ThreeVector(-101.6*mm, 0, 139.8*mm));

 // G4LogicalVolume* logicLead = new G4LogicalVolume(solidLead, Pb, "Lead");

  fLContain = new G4LogicalVolume(solidLead,fContainMaterial ,fContainMaterial->GetName());
  leadPhys.push_back(new G4PVPlacement(0, G4ThreeVector(-17.8*mm, 0, 307.*mm), fLContain, fContainMaterial->GetName(), logicWorld, false, 0,0));

// bricks



  G4Box* solidLeadBrick1 = new G4Box("LeadBrick1", 25.4*mm, 70. *mm,  50.8*mm);
  G4Box* solidLeadBrick2 = new G4Box("LeadBrick2", 50.8*mm, 50.8*mm,  25.4*mm);
  G4Box* solidLeadBrick3 = new G4Box("LeadBrick3", 50.8*mm, 26.3*mm,  25.4*mm);

  G4LogicalVolume* logicLeadBrick1 = new G4LogicalVolume(solidLeadBrick1, Pb, "LeadBrick1");
  G4LogicalVolume* logicLeadBrick2 = new G4LogicalVolume(solidLeadBrick2, Pb, "LeadBrick2");
  G4LogicalVolume* logicLeadBrick3 = new G4LogicalVolume(solidLeadBrick3, Pb, "LeadBrick3");

  leadPhys.push_back(new G4PVPlacement(0, G4ThreeVector(-170.2*mm,   0     , 357.8*mm), logicLeadBrick1, "LeadBrick", logicWorld, false, 0, 0));
  leadPhys.push_back(new G4PVPlacement(0, G4ThreeVector(-144.8*mm,   0     , 434. *mm), logicLeadBrick2, "LeadBrick", logicWorld, false, 1, 0));
  leadPhys.push_back(new G4PVPlacement(0, G4ThreeVector(-144.8*mm,   0     , 281.6*mm), logicLeadBrick2, "LeadBrick", logicWorld, false, 2, 0));
  leadPhys.push_back(new G4PVPlacement(0, G4ThreeVector(   7.6*mm, 127.9*mm, 357.8*mm), logicLeadBrick3, "LeadBrick", logicWorld, false, 3, 0));
  leadPhys.push_back(new G4PVPlacement(0, G4ThreeVector(   7.6*mm,-127.9*mm, 357.8*mm), logicLeadBrick3, "LeadBrick", logicWorld, false, 4, 0));



// BeO


G4VSolid* solidBeOTop        = new G4Tubs("BeOTop"       , 0     , 13.75*mm,  4.5 *mm, 0, 2*M_PI);
  G4VSolid* solidBeOTube       = new G4Tubs("BeOTube"      , 6.5*mm, 13.  *mm, 12.75*mm, 0, 2*M_PI);
  G4VSolid* solidBeOBottomFull = new G4Tubs("BeOBottomFull", 0     , 22.6 *mm, 12.5 *mm, 0, 2*M_PI);
  G4VSolid* solidBeOBottomHole = new G4Tubs("BeOBottomHole", 0     , 13.75*mm,  8.1 *mm, 0, 2*M_PI);
  G4SubtractionSolid* solidBeOBottom = new G4SubtractionSolid("BeOBottom", solidBeOBottomFull, solidBeOBottomHole, 0, G4ThreeVector(0, 0, 4.6*mm));

G4LogicalVolume* logicBeOTop    = new G4LogicalVolume(solidBeOTop   , fAbsorMaterial, "BeOTop"   );
  G4LogicalVolume* logicBeOTube   = new G4LogicalVolume(solidBeOTube  , fAbsorMaterial, "BeOTube"  );
  G4LogicalVolume* logicBeOBottom = new G4LogicalVolume(solidBeOBottom, fAbsorMaterial, "BeOBottom");

  BeOPhys.push_back(new G4PVPlacement(rotY, G4ThreeVector(-31.4 *mm, 0, 357.8*mm), logicBeOTop   , "BeOTop"   , logicWorld, false, 0, 0));
  BeOPhys.push_back(new G4PVPlacement(rotY, G4ThreeVector(-14.15*mm, 0, 357.8*mm), logicBeOTube  , "BeOTube"  , logicWorld, false, 0, 0));
  BeOPhys.push_back(new G4PVPlacement(rotY, G4ThreeVector( -4.9 *mm, 0, 357.8*mm), logicBeOBottom, "BeOBottom", logicWorld, false, 0, 0));


// Antimony source 
//

 G4VSolid* solidSbSource  = new G4Tubs("SbSource" , 0, 1.5875*mm,  7.5 *mm, 0, 2*M_PI);
  G4VSolid* solidEpoxyFull = new G4Tubs("EpoxyFull", 0, 6.5   *mm, 12.75*mm, 0, 2*M_PI);
  G4SubtractionSolid* solidEpoxy = new G4SubtractionSolid("Epoxy", solidEpoxyFull, solidSbSource, 0, G4ThreeVector(0, 0,-5.25*mm));

  G4LogicalVolume* logicSbSource = new G4LogicalVolume(solidSbSource, Sb   , "SbSource");
  G4LogicalVolume* logicEpoxy    = new G4LogicalVolume(solidEpoxy   , Epoxy, "Epoxy"   );

  new G4PVPlacement(rotY, G4ThreeVector( -8.9 *mm, 0, 357.8*mm), logicSbSource, "SbSource", logicWorld, false, 0, 0);
  new G4PVPlacement(rotY, G4ThreeVector(-14.15*mm, 0, 357.8*mm), logicEpoxy   , "Epoxy"   , logicWorld, false, 0, 0);

*/



  PrintParameters();
  
  //always return the root volume
  //
  return fPWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
 G4cout << "\n The Absorber  is a cylinder of " << fAbsorMaterial->GetName()
        << "  radius = " << G4BestUnit(fAbsorRadius,"Length")
        << "  length = " << G4BestUnit(fAbsorLength,"Length") << G4endl;
 G4cout << " The Container is a cylinder of " << fContainMaterial->GetName()
        << "  thickness = " << G4BestUnit(fContainThickness,"Length") << G4endl;
 
 G4cout << "\n" << fAbsorMaterial << G4endl;
 G4cout << "\n" << fContainMaterial << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  
  if (pttoMaterial) { 
    fAbsorMaterial = pttoMaterial;
    if(fLAbsor) { fLAbsor->SetMaterial(fAbsorMaterial); }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetAbsorMaterial : "
           << materialChoice << " not found" << G4endl;
  }              
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetContainMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  
  if (pttoMaterial) { 
    fContainMaterial = pttoMaterial;
    if(fLContain) { fLContain->SetMaterial(fContainMaterial); }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetContainMaterial : "
           << materialChoice << " not found" << G4endl;
  }              
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorRadius(G4double value)
{
  fAbsorRadius = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorLength(G4double value)
{
  fAbsorLength = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetContainThickness(G4double value)
{
  fContainThickness = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

