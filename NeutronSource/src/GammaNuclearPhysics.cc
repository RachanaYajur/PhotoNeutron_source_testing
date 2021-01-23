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
/// \file GammaNuclearPhysics.cc
/// \brief Implementation of the GammaNuclearPhysics class
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "GammaNuclearPhysics.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

// Processes

//#include "G4PhotoNuclearProcess.hh"
//#include "G4CascadeInterface.hh"

#include "G4SystemOfUnits.hh"

#include "G4EmExtraPhysics.hh"
#include "G4BertiniElectroNuclearBuilder.hh"
#include "G4LENDModel.hh"
#include "G4LENDInelastic.hh"
#include "G4LENDorBERTModel.hh"
#include "G4LENDCombinedCrossSection.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GammaNuclearPhysics::GammaNuclearPhysics(const G4String& name)
:  G4VPhysicsConstructor(name)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GammaNuclearPhysics::~GammaNuclearPhysics()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GammaNuclearPhysics::ConstructProcess()
{
   G4ProcessManager* pManager = G4Gamma::Gamma()->GetProcessManager();
   //
   G4PhotoNuclearProcess* process = new G4PhotoNuclearProcess();
   //
   //
   //G4LENDModel* lend = new G4LENDModel();
  // process->RegisterMe(lendmodel);  

  G4LENDInelastic* lend = new G4LENDInelastic(G4Gamma::Gamma());
  //G4LENDorBERTModel* lend = new G4LENDorBERTModel(G4Gamma::Gamma());
  lend->SetMaxEnergy(10*MeV);
  process->RegisterMe(lend);
  G4LENDCombinedCrossSection* lendXS = new G4LENDCombinedCrossSection(G4Gamma::Gamma());
   process->AddDataSet(lendXS);
   //G4LENDInelastic* lend = new G4LENDInelastic();
   

   /*
   G4CascadeInterface* bertini = new G4CascadeInterface();
   bertini->SetMaxEnergy(10*GeV);
   process->RegisterMe(bertini);
   //
  
   G4ElectroVDNuclearModel* model = new G4ElectroVDNuclearModel();
   process->RegisterMe(model);
   //
   G4TheoFSGenerator* generator = new G4TheoFSGenerator();
   process->RegisterMe(generator);
   //
   G4GeneratorPrecompoundInterface* interface = new G4GeneratorPrecompoundInterface();  
   process->RegisterMe(interface);
   //
   G4QGSModel< G4GammaParticipants >* smodel = new G4QGSModel< G4GammaParticipants >();
   process->RegisterMe(smodel);   
   //
   //G4QGSMFragmentation* frag = new G4QGSMFragmentation();
   
   //
  // G4ExcitedStringDecay* string = new G4ExcitedStringDecay();
  // frag->RegisterMe(string);   
  // process->RegisterMe(frag);
  */
   pManager->AddDiscreteProcess(process);
   


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
