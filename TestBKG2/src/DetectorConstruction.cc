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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction()
{
  fWorldSize = 1*m;
  
  //Initialize custom messenger to adjust parameters
  fMessenger = new G4GenericMessenger(this, "/detcon/", "Detector Construction");

  fMessenger->DeclareProperty("isShielded", isShielded, "Load in Lead & Copper Shielding");
  fMessenger->DeclareProperty("rThickness", rThickness, "If isShielded is true, set Thickness of Lead Shield. Unit = cm.");
  fMessenger->DeclareProperty("dDiameter", dDiameter, "Set Diameter of HPGe. Unit = cm.");
  fMessenger->DeclareProperty("zLength", zLength, "Set Length of HPGe. Unit = cm.");
  fMessenger->DeclareProperty("outerShieldDiameter", outerShieldDiameter, "Set Outer Diameter of the Outer Shield (S3). Unit = in.");
  fMessenger->DeclareProperty("halfShell", halfShell, "Set True to study only half of the isotope shell");
  fMessenger->DeclareProperty("isHydrogen", isHydrogen, "Load Isotope Shell material as hydrogen instead of isotope mixture");


  //Set default values here (these can be changed here for debugging)
  isShielded = true;
  isHydrogen = false;
  rThickness = 5;
  dDiameter = 8;
  zLength = 3;
  outerShieldDiameter = 8; //This will be changed //Default: 8
  halfShell = true;
  useTarget = true;

  DefineMaterials();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  man = G4NistManager::Instance();

    //Vacuum using NIST Manager //From Jin's Simulation
  Galactic = new G4Material("Galactic", universe_mean_density, 1);
  Galactic->AddElement(man->FindOrBuildElement("H"),1);

//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~

    //Imported Materials
      //Detector
  Ge = man->FindOrBuildMaterial("G4_Ge");

      //Shielding
  Cu = man->FindOrBuildMaterial("G4_Cu");
  Pb = man->FindOrBuildMaterial("G4_Pb");
  Al = man->FindOrBuildMaterial("G4_Al");

      //Target
  BC = man->FindOrBuildMaterial("G4_BORON_CARBIDE");

//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~

    //Isotope Shell - Source of background spectrum
      //Hydrogen for initial test
    H1 = new G4Isotope("H1", 1, 1);//name, iz, n
   elH = new G4Element("isopure hydrogen", "H", 1);//name, symbol, # of isotopes
  	elH->AddIsotope(H1, 100*perCent); //isotope, abundance

      //Jin's List
  Cl35 = new G4Isotope("Cl35", 17, 35);
  Cr50 = new G4Isotope("Cr50", 24, 50);
  Cr52 = new G4Isotope("Cr50", 24, 52);
  Cr53 = new G4Isotope("Cr50", 24, 53);
  Fe56 = new G4Isotope("Cr50", 26, 56);
  Ni58 = new G4Isotope("Cr50", 28, 58);
  Ni60 = new G4Isotope("Cr50", 28, 60);
  Ni62 = new G4Isotope("Cr50", 28, 62);
  Ni63 = new G4Isotope("Cr50", 28, 63);
  Ni64 = new G4Isotope("Cr50", 28, 64);
  Na24 = new G4Isotope("Cr50", 11, 24);
  Mn56 = new G4Isotope("Cr50", 25, 56);

  elCl = new G4Element("MatCl", "Cl", 1);
  	elCl->AddIsotope(Cl35, 100*perCent);
  elCr = new G4Element("MatCr", "Cr", 3);
  	elCr->AddIsotope(Cr50, 33.33*perCent);
	elCr->AddIsotope(Cr52, 33.33*perCent);
	elCr->AddIsotope(Cr53, 33.34*perCent);
  elFe = new G4Element("MatFe", "Fe", 1);
  	elFe->AddIsotope(Fe56, 100*perCent);
  elNi = new G4Element("MatNi", "Ni", 5);
  	elNi->AddIsotope(Ni58, 20*perCent);
	elNi->AddIsotope(Ni60, 20*perCent);
	elNi->AddIsotope(Ni62, 20*perCent);
	elNi->AddIsotope(Ni63, 20*perCent);
	elNi->AddIsotope(Ni64, 20*perCent);
  elNa = new G4Element("MatNa", "Na", 1);
  	elNa->AddIsotope(Na24, 100*perCent);
  elMn = new G4Element("MatMn", "Mn", 1);
  	elMn->AddIsotope(Mn56, 100*perCent);
  
      //Construct Isotope Shell Material
  if(isHydrogen) //Hydrogen Testing
  {
    IsotopeShellMat = new G4Material("IsotopeShellMat",
		             	     0.07085 * g/cm3,
				     1);
    	IsotopeShellMat->AddElement(elH, 100*perCent);
  }

  else //From Jin's List
  {  
    IsotopeShellMat = new G4Material("IsotopeShellMat", //name
		  		     7.85 * g/cm3, //density of material
				     6); //number of elements
  	  IsotopeShellMat->AddElement(elCl, 8.34*perCent); //I gave them all roughly the same abundance & sum to 100 percent
  	  IsotopeShellMat->AddElement(elCr, 24.99*perCent);
  	  IsotopeShellMat->AddElement(elFe, 8.34*perCent);
  	  IsotopeShellMat->AddElement(elNi, 41.65*perCent);
  	  IsotopeShellMat->AddElement(elNa, 8.34*perCent);
  	  IsotopeShellMat->AddElement(elMn, 8.34*perCent);
  }

}

//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~
//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~
//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    //Define Constants
  
  zDistanceFromOrigin = 0; //In this simulation, the detector is centered at the origin
  G4cout << "zLength = " << zLength << G4endl;

//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~

    //World Volume
  solidWorld = new G4Box("World",
		  	 fWorldSize/2,
			 fWorldSize/2,
			 fWorldSize/2);
  logicWorld = new G4LogicalVolume(solidWorld,
		  		   Galactic,
				   "World");
  physiWorld = new G4PVPlacement(0,
		  		 G4ThreeVector(),
				 logicWorld,
				 "World",
				 0,
				 false,
				 0);

//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~

    //HPGe Detector
  solidDet = new G4Tubs("Det",
		 	0,
			(dDiameter/2)*cm,
			(zLength/2)*cm,
			0*deg,
			360*deg);
  logicDet = new G4LogicalVolume(solidDet,
		  		 Ge,
				 "Det");
  fScoringVolume = logicDet;
  physDet = new G4PVPlacement(0,
		  	      G4ThreeVector(0,
				      	    0,
					    (zDistanceFromOrigin+(zLength/2))*cm),
			      logicDet,
			      "Det",
			      logicWorld,
			      false,
			      0,
			      false);

//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~

    //Al Shield (for e+e- fix)
    //Shell
  Sidel = (zLength/2 + 1)*cm;
  solidAl1 = new G4Tubs("solidAl1",
		  	(dDiameter*cm/2 + 1*cm),
			(dDiameter*cm/2 + 1.1*cm),
			Sidel,
			0*deg,
			360*deg);
  logicAl1 = new G4LogicalVolume(solidAl1,
		  		 Al,
				 "logicAl1");
  physAl1 = new G4PVPlacement(0,
		  	      G4ThreeVector(0,
				      	    0,
					    (zDistanceFromOrigin+(zLength/2))*cm),
			      logicAl1,
			      "physAl1",
			      logicWorld,
			      false,
			      0,
			      false);

    //Caps
  Capl = 1*mm;
  solidAl2 = new G4Tubs("solidAl1",
		        0,
			(dDiameter*cm/2 + 1.1*cm),
			Capl/2,
			0*deg,
			360*deg);
  logicAl2 = new G4LogicalVolume(solidAl2,
		  		 Al,
				 "logicAl2");
  physAl2_1 = new G4PVPlacement(0,
		  		G4ThreeVector(0,
					      0,
					      (zDistanceFromOrigin+(zLength/2))*cm+(Sidel+Capl/2)),
				logicAl2,
				"physAl2_1",
				logicWorld,
				false,
				0,
				false);
  physAl2_2 = new G4PVPlacement(0,
		  		G4ThreeVector(0,
					      0,
					      (zDistanceFromOrigin+(zLength/2))*cm-(Sidel+Capl/2)),
				logicAl2,
				"physAl2_2",
				logicWorld,
				false,
				0,
				false);

//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~

    //Isotope Shell
    innerShellRadius = 18.52*cm;
    ShellThickness   = 1*cm;
    if(halfShell){ShellSweep = 90*deg;}
    else {ShellSweep = 180*deg;}


    solidIsoSphere = new G4Sphere("solidIsoSphere",
		                  innerShellRadius,
				  innerShellRadius+ShellThickness,
				  0*deg,
				  360*deg,
				  0*deg,
				  ShellSweep);
    logicIsoSphere = new G4LogicalVolume(solidIsoSphere,
		    			 IsotopeShellMat,
					 "logicIsoSphere");
    physIsoSphere = new G4PVPlacement(0,
		    		      G4ThreeVector(0,
					      	    0,
						    (zDistanceFromOrigin+zLength/2)*cm),
				      logicIsoSphere,
				      "physIsoSphere",
				      logicWorld,
				      false,
				      0,
				      false);

//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~

    //Shielding
  if(isShielded)
  {
	//Shield 1 (Lead) ~~~~|~~~~|~~~~|~~~~|~~~~|~~~~|~~~~|~~~~
	Di1 = 1; //Inches
	Do1 = 5;
	zL1 = 2;
	DistOrigin1 = zDistanceFromOrigin/2.54 + (zL1/2) + zLength /2.54 + 1; // +1 will be changed
  	solidShield1 = new G4Tubs("solidShield1",
			  	  Di1*2.54/2*cm,
				  Do1*2.54/2*cm,
				  zL1*2.54/2*cm,
				  0*deg,
				  360*deg);
	logicShield1 = new G4LogicalVolume(solidShield1,
					   Pb,
					   "logicShield1");
	physShield1 = new G4PVPlacement(0,
					G4ThreeVector(0,
						      0,
						      DistOrigin1*2.54*cm),
					logicShield1,
					"physShield1",
					logicWorld,
					false,
					0,
					false);
    
	//Shield 2 (Copper) ~~~~|~~~~|~~~~|~~~~|~~~~|~~~~|~~~~|~~~~
        Di2 = 4.125; //Inches
        Do2 = 4.25;
        zL2 = zLength/2.54 + 0.5; //This will be changed
        DistOrigin2 = zDistanceFromOrigin/2.54 + (zL2/2);
	solidShield2 = new G4Tubs("solidShield2",
				  Di2*2.54/2*cm,
				  Do2*2.54/2*cm,
				  zL2*2.54/2*cm,
				  0*deg,
				  360*deg);
        logicShield2 = new G4LogicalVolume(solidShield2,
					   Cu,
					   "logicShield2");
        physShield2 = new G4PVPlacement(0,
					G4ThreeVector(0,
						      0,
						      DistOrigin2*2.54*cm),
					logicShield2,
					"physShield2",
					logicWorld,
					false,
					0,
					false);

	//Shield 3 (Lead) ~~~~|~~~~|~~~~|~~~~|~~~~|~~~~|~~~~|~~~~
        Di3 = 5.125; //Inches This will be changed
 	Do3 = outerShieldDiameter; //Inches
        zL3 = zL1 + zLength/2.54 + 1; //This will be changed
        DistOrigin3 = zDistanceFromOrigin/2.54+(zL3/2); //This will be changed
        solidShield3 = new G4Tubs("solidShield3",
				  Di3*2.54/2*cm,
				  Do3*2.54/2*cm,
				  zL3*2.54/2*cm,
				  0*deg,
				  360*deg);
        logicShield3 = new G4LogicalVolume(solidShield3,
					   Pb,
					   "logicShield3");
        physShield3 = new G4PVPlacement(0,
					G4ThreeVector(0,
						      0,
						      DistOrigin3*2.54*cm),
					logicShield3,
					"physShield3",
					logicWorld,
					false,
					0,
					false);
  }

//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~//~~~~

    //Boron Carbide Target
    TarLength = 10*cm;
    TarWidth = 1*mm;
    TarDistance = 20*cm;
    if(useTarget)
    {
    	solidTarget = new G4Box("solidTarget",
                                TarLength/2,
                                TarLength/2,
                                TarWidth/2);
    	logicTarget = new G4LogicalVolume(solidTarget,
                                          BC,
                                          "logicTarget");
    	physTarget = new G4PVPlacement(0,
                                       G4ThreeVector(0,
					   	     0,
						     zDistanceFromOrigin-TarDistance),
                                       logicTarget,
                                       "physTarget",
                                       logicWorld,
                                       false,
                                       0,
				       false);
    }

  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
