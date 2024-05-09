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
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4GenericMessenger.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    DetectorConstruction();
   ~DetectorConstruction();

    //Scoring Volume
    G4LogicalVolume *GetScoringVolume() const { return fScoringVolume; }

    virtual     
    G4VPhysicalVolume* Construct();
                        
    G4double GetWorldSize() {return fWorldSize;}; 
    

  private:
 
    //DetectorConstruction variables 
    G4GenericMessenger *fMessenger;

    G4double fWorldSize, rThickness, zLength, dDiameter, outerShieldDiameter, ShellSweep;
     
    G4bool isShielded, isHydrogen, halfShell, useTarget;

//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~

    //DefineMaterials variables
    G4NistManager *man;

    G4Material *Galactic, *Ge, *Cu, *Pb, *Al, *IsotopeShellMat, *BC;

    void DefineMaterials();

    G4Isotope *H1, *Cl35, *Cr50, *Cr52, *Cr53, *Fe56, *Ni58, *Ni60, *Ni62, *Ni63, *Ni64, *Na24, *Mn56;

    G4Element *elH, *elCl, *elCr, *elFe, *elNi, *elNa, *elMn;

//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~

    //Construct variables
    G4double zDistanceFromOrigin, Sidel, Capl, innerShellRadius, ShellThickness, Di1, Do1, zL1, DistOrigin1, Di2, Do2, zL2, DistOrigin2, Di3, Do3, zL3, DistOrigin3, TarLength, TarWidth, TarDistance;

    G4Box *solidWorld, *solidTarget;
    G4Tubs *solidDet, *solidAl1, *solidAl2, *solidShield1, *solidShield2, *solidShield3;
    G4Sphere *solidIsoSphere;
    
    G4LogicalVolume *logicWorld, *logicDet, *fScoringVolume, *logicAl1, *logicAl2, *logicIsoSphere, *logicShield1, *logicShield2, *logicShield3, *logicTarget;

    G4VPhysicalVolume *physiWorld, *physDet, *physAl1, *physAl2_1, *physAl2_2, *physShield1, *physShield2, *physShield3, *physIsoSphere, *physTarget;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

