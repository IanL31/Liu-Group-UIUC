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
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "PrimaryGeneratorAction.hh"

/*#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Geantino.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Gamma.hh"
#include "G4Neutron.hh"
#include "G4Positron.hh"
#include "Randomize.hh"
#include "DetectorConstruction.hh"*/

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4Geantino.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  fParticleGun->SetParticleDefinition(
               particleTable->FindParticle(particleName="neutron"));

  fParticleGun->SetParticleEnergy(0.003325*eV);
  fParticleGun->SetParticlePosition(G4ThreeVector(0, 0, 17.52*cm));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4double ShellRadius = (18.52 - 1); //Static cast from DetConst class
  G4double zLength     = 3;
  G4bool   halfShell   = true;        //Static cast from DetConst class
  G4bool   useTarget   = true;        //Static cast from DetConst class
  G4int    evtNb       = anEvent->GetEventID();

  if(useTarget && evtNb%2==1) //Events for boron carbide target
  {
	//Split odd events to impact target instead of isotope shell
    fParticleGun->SetParticlePosition(G4ThreeVector(0,
                                                    0,
                                                    -15*cm));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,
                                                             0,
                                                             -1));
  }

  else //Events for background spectrum
  {
 
      //Determine the position & direction of the neutrons
    G4double randN1 = 2*G4UniformRand() - 1;
    G4double randN2 = 2*G4UniformRand() - 1;
    G4double randN3 = 2*G4UniformRand() - 1;

      //Normalize the vectors to just inside the Isotope Shell inner radius
    G4double Norm = pow((pow(randN1,2) + pow(randN2,2) + pow(randN3,2)),0.5)/ShellRadius;

    G4double NormRN1 = randN1/Norm;
    G4double NormRN2 = randN2/Norm;
    G4double NormRN3 = randN3/Norm + zLength/2; //To adjust for the offset center

    if(halfShell && NormRN3<0)
    {
  	//to increase simulation efficiency, find a way to restrict neutron location to inside of the shell  
      fParticleGun->SetParticlePosition(G4ThreeVector(1*NormRN1*cm,
			  			    1*NormRN2*cm,
						    -1*NormRN3*cm));
      fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1*NormRN1,
			  				     1*NormRN2,
							     -1*NormRN3));
      //G4cout << "Neutron z_init: " << NormRN3 << G4endl;
      //G4cout << "Neutron loc: " << (NormRN1, NormRN2, -1*NormRN3) << G4endl;
							     
    }

    else
    {
      fParticleGun->SetParticlePosition(G4ThreeVector(NormRN1*cm,
                                                    NormRN2*cm,
                                                    NormRN3*cm));
      fParticleGun->SetParticleMomentumDirection(G4ThreeVector(NormRN1,
                                                             NormRN2,
                                                             NormRN3));
							     
    }
  }

  fParticleGun->GeneratePrimaryVertex(anEvent);
   
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
