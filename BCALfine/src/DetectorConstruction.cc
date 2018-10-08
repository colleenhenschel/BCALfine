/*
 	Detector Construction file has definitions for variables and classes in DetectorConstruction.hh
	The .hh file is located under the include directory within the simulation directory
  ie.  /BCALfine/include/DetectorConstruction.hh 

	The file is used to define the geometry of the detector as well as the materials
*/

#include "DetectorConstruction.hh"

#include "MCEvent.hh"
#include "DataManager.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
#include "G4NistManager.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"


#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
using namespace std;

//-----------------------------------------------------------------------------
 
DetectorConstruction::DetectorConstruction(G4double worldSize)
 
{
  solidWorld = 0;
  logicWorld = 0; 
  physiWorld = 0;

	solidScintFiber = 0;
	logicScintFiber = 0;
	

  stepLimit = 0;

  NISTManager         = G4NistManager::Instance();

  World_VisAtt   = 0;
  VisAttScintFiber = 0; 
 
  fWorldLength        = worldSize;
  //field size devided by (fiber_diameter+fiber_gap) has to be even number to obtain symetric arangement of fibers
  field_size          = 40.455*cm;//10.44*cm;  //40.455*cm for 931 fibers per row; 
                                   //4.35*cm for 101 fibers per row; 
                                   //2.61*cm for 61 fibers per row; 
                                   //2.61*mm for 7 fibers per row; 
                                   //100.*cm;

 	fiber_length = 390.0*cm;
	fiber_radius = 0.5*mm;
	module_outter = 114.27*mm;
	module_inner = 85.24*mm;
	module_length = 390.0*cm;
	module_height = 221.94*mm;
	fibervertspacing = 1.18*mm;
	halfglueboxheight = 0.1165*mm;
	midradius = 77.75; //76.26; // distance from center of module to center of BCAL (in centimeters), set to avoid overlaps with baseplates
	bp_height = module_height + 7.93*mm + 31.75*mm; // height from bottom of inner base plate to top of outter base plate
  baseplate_inner = 84.75*mm;
  baseplate_outter = 119.0*mm;
  //-----------------------------------------------------------------;
  //DataManager* dataManager = DataManager::GetInstance();
}

//-----------------------------------------------------------------------------
 
DetectorConstruction::~DetectorConstruction()
{
  // if(World_VisAtt!=0)  {delete World_VisAtt;   World_VisAtt   = 0;}
  // if(VisAttWrapping!=0){delete VisAttWrapping; VisAttWrapping = 0;}
  // if(VisAttScint!=0)   {delete VisAttScint;    VisAttScint    = 0;}
  // if(VisAttSensor!=0)  {delete VisAttSensor;   VisAttSensor   = 0;}  

  // if(solidWrapping!=0){delete solidWrapping; solidWrapping = 0;}
  // if(logicWrapping!=0){delete logicWrapping; logicWrapping = 0;}
  // if(physiWrapping!=0){delete physiWrapping; physiWrapping = 0;}

  // if(solidScint!=0){delete solidScint; solidScint = 0;}
  // if(logicScint!=0){delete logicScint; logicScint = 0;}
  // if(physiScint!=0){delete physiScint; physiScint = 0;}

  // if(solidSensor!=0){delete solidSensor; solidSensor = 0;}
  // if(logicSensor!=0){delete logicSensor; logicSensor = 0;}
  // if(physiSensor!=0){delete physiSensor; physiSensor = 0;}
  
  // if(solidWorld!=0){delete solidWorld; solidWorld = 0;}
  // if(logicWorld!=0){delete logicWorld; logicWorld = 0;}
  // if(physiWorld!=0){delete physiWorld; physiWorld = 0;}

  if(stepLimit!=0){delete stepLimit; stepLimit =0;}
}

//-----------------------------------------------------------------------------
 
G4VPhysicalVolume* DetectorConstruction::Construct()
{  
  //------------------------------
  // Construct Materials
  //------------------------------
  ConstructMaterials();

  //------------------------------ 
  // World
  //------------------------------ 
  G4GeometryManager::GetInstance()->SetWorldMaximumExtent(fWorldLength);
  G4cout << "Computed tolerance = "
         << G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/mm
         << " mm" << G4endl;

  G4double HalfWorldLength = 0.5*fWorldLength;

  solidWorld = new G4Box("sWorld",HalfWorldLength,HalfWorldLength,HalfWorldLength);
  logicWorld = new G4LogicalVolume(solidWorld, NISTManager->FindOrBuildMaterial("G4_AIR"), "lWorld", 0, 0, 0);
  
  //  Must place the World Physical volume unrotated at (0,0,0).
  physiWorld = new G4PVPlacement(0, G4ThreeVector(0.*mm,0.*mm,0.*mm), logicWorld, "pvWorld", 0, false, 0);   


  //Start building the detector
  //------------------------------------------------------------------

	G4double modx = 0.;
	G4double mody = midradius;
	G4double modz = 0.;
	
	

	//Defining module:
	
	G4RotationMatrix *Rot = new G4RotationMatrix;
	Rot->rotateX(90.0*deg); 
	Rot->rotateZ(90.0*deg);
	//Rot->rotateY(90.0*deg);


	solidModule = new G4Trd("sModule", module_inner/2, module_outter/2, module_length/2, module_length/2, module_height/2);
	logicModule = new G4LogicalVolume(solidModule, NISTManager->FindOrBuildMaterial("G4_Pb"), "lModule", 0, 0, 0);

	solidBasePlate = new G4Trd("sBasePlate", baseplate_inner/2, baseplate_outter/2, module_length/2, module_length/2, bp_height/2);
 	logicBasePlate = new G4LogicalVolume(solidBasePlate, NISTManager->FindOrBuildMaterial("G4_Al"), "lBasePlate", 0, 0, 0);

	physiModule = new G4PVPlacement(0, G4ThreeVector(0.*cm, 0.*cm, -(3.175-0.793)/2.*cm), logicModule, "pvModule", logicBasePlate, false, 0);

		
 	// Definition of G4Tubs -> cylinder (name, inner radius, outer radius, half length, starting phi, segment angle)

	G4RotationMatrix *rot2 = new G4RotationMatrix;
	rot2->rotateY(90.0*deg);
	rot2->rotateX(90.0*deg);

	solidScintFiber = new G4Tubs("sScintFiber", 0.0, 0.96*mm/2, fiber_length/2, 0.0, 360.0*deg);
	logicScintFiber = new G4LogicalVolume(solidScintFiber, NISTManager->FindOrBuildMaterial("plasticScint"), "lScintFiber", 0, 0, 0);

	solidCladding1 = new G4Tubs("sCladding1", 0, 0.99*mm/2, fiber_length/2, 0.0, 360.0*deg);
	logicCladding1 = new G4LogicalVolume(solidCladding1, NISTManager->FindOrBuildMaterial("PMMA"), "lCladding1", 0, 0, 0);

	solidCladding2 = new G4Tubs("sCladding2", 0, 0.5*mm, fiber_length/2, 0.0, 360.0*deg);
	logicCladding2 = new G4LogicalVolume(solidCladding2, NISTManager->FindOrBuildMaterial("FPethylene"), "lCladding2", 0, 0, 0);

	solidGlueRing = new G4Tubs("sGlueRing", 0, 0.55*mm, fiber_length/2, 0.0, 360.0*deg);
	logicGlueRing = new G4LogicalVolume(solidGlueRing, NISTManager->FindOrBuildMaterial("epoxy"), "lGlueRing", 0, 0, 0);

	solidGlueBox = new G4Box("sGlueBox", 0.125*mm, module_length/2, halfglueboxheight);
	logicGlueBox = new G4LogicalVolume(solidGlueBox, NISTManager->FindOrBuildMaterial("epoxy"), "lGlueBox", 0, 0, 0);

  new G4PVPlacement(0, G4ThreeVector(0.*mm, 0.*mm, 0.*mm), logicScintFiber, "pvScintFiber", logicCladding1, false, 0);
	new G4PVPlacement(0, G4ThreeVector(0.*mm, 0.*mm, 0.*mm), logicCladding1, "pvCladding1", logicCladding2, false, 0);
	new G4PVPlacement(0, G4ThreeVector(0.*mm, 0.*mm, 0.*mm), logicCladding2, "pvCladding2", logicGlueRing, false, 0);
					
					

	G4double currentwidthb;
	G4double currentwidtht;

	G4double currentxpos;
	G4int copyNofiber = 0;
	for(int j = 1; j< 188; j++)
	{
		currentwidthb = (module_inner) + 2*(14.5*mm/module_height)*(j*fibervertspacing/* +1.145*mm */ - halfglueboxheight);
		currentwidtht = (module_inner) + 2*(14.5*mm/module_height)*(j*fibervertspacing/* +1.145*mm */ + halfglueboxheight);
		
		
		VisualizationAttributes();
		if(j%2 == 0)
		{
			currentxpos = 0;
			while(currentxpos > -currentwidthb/2 + 0.55)
			{
				new G4PVPlacement(rot2, G4ThreeVector(currentxpos*mm, 0.*mm, -module_height/2+j*fibervertspacing/* +1.145*mm */), logicGlueRing, "pvGlueRing", logicModule, false, copyNofiber);
				copyNofiber++;
				if(currentxpos > -currentwidthb/2 + 0.55 + 0.25)					
					new G4PVPlacement(0, G4ThreeVector(currentxpos*mm - 0.675, 0.*mm, -module_height/2+j*fibervertspacing/* +1.145*mm */), logicGlueBox, "pvGlueBox", logicModule, false, 0);			
				currentxpos -= 1.35;
			}
			currentxpos = 1.35;
			while(currentxpos < currentwidthb/2 - 0.55)
			{
				new G4PVPlacement(rot2, G4ThreeVector(currentxpos*mm, 0.*mm, -module_height/2+j*fibervertspacing/* +1.145*mm */), logicGlueRing, "pvGlueRing", logicModule, false, copyNofiber);
				copyNofiber++;
				new G4PVPlacement(0, G4ThreeVector(currentxpos*mm - 0.675, 0.*mm, -module_height/2+j*fibervertspacing/* +1.145*mm */), logicGlueBox, "pvGlueBox", logicModule, false, 0);	
				if(currentxpos < currentwidthb/2 - 0.55 - 0.25)
					new G4PVPlacement(0, G4ThreeVector(currentxpos*mm + 0.675, 0.*mm, -module_height/2+j*fibervertspacing/* +1.145*mm */), logicGlueBox, "pvGlueBox", logicModule, false, 0);					
				currentxpos += 1.35;
			}
		}
		else
		{
			currentxpos = -0.675;
			while(currentxpos > -currentwidthb/2 + 0.55)
			{
				new G4PVPlacement(rot2, G4ThreeVector(currentxpos*mm, 0.*mm, -module_height/2+j*fibervertspacing/* +1.145*mm */), logicGlueRing, "pvGlueRing", logicModule, false, copyNofiber);
				copyNofiber++;
				if(currentxpos > -currentwidthb/2 + 0.55 + 0.25)					
					new G4PVPlacement(0, G4ThreeVector(currentxpos*mm - 0.675, 0.*mm, -module_height/2+j*fibervertspacing/* +1.145*mm */), logicGlueBox, "pvGlueBox", logicModule, false, 0);				
				currentxpos -= 1.35;
			}
			currentxpos = 0.675;
			while(currentxpos < currentwidthb/2 - 0.55)
			{
				new G4PVPlacement(rot2, G4ThreeVector(currentxpos*mm, 0.*mm, -module_height/2+j*fibervertspacing/* +1.145*mm */), logicGlueRing, "pvGlueRing", logicModule, false, copyNofiber);
				copyNofiber++;
				new G4PVPlacement(0, G4ThreeVector(currentxpos*mm - 0.675, 0.*mm, -module_height/2+j*fibervertspacing/* +1.145*mm */), logicGlueBox, "pvGlueBox", logicModule, false, 0);	
				if(currentxpos < currentwidthb/2 - 0.55 - 0.25)
					new G4PVPlacement(0, G4ThreeVector(currentxpos*mm + 0.675, 0.*mm, -module_height/2+j*fibervertspacing/* +1.145*mm */), logicGlueBox, "pvGlueBox", logicModule, false, 0);			
				currentxpos += 1.35;
			}
		}			
	}
 
	G4int nummods = 48;

	paramBasePlate = new BasePlateParameterisation(Rot, NISTManager->FindOrBuildMaterial("G4_Al"), nummods, midradius);
  physiBasePlate = new G4PVParameterised("pvBasePlate",            // Its name
                                         logicBasePlate,             // Its logical volume
                                         logicWorld,                // Mother logical volume
                                         kZAxis,               // Allow default voxelising -- no axis
                                         nummods,                  // Number of modules
                                         paramBasePlate, false);    // The parameterisation

  	
		

  //--------- example of User Limits -------------------------------

  // below is an example of how to set tracking constraints in a given
  // logical volume(see also in PhysicsList how to setup the processes
  // G4StepLimiter or G4UserSpecialCuts).
    
  // Sets a max Step length in the tracker region, with G4StepLimiter
  //
  // G4double maxStep = 0.5*cm;
  // stepLimit = new G4UserLimits(maxStep);
  // logicWorld->SetUserLimits(stepLimit);
  
  // Set additional contraints on the track, with G4UserSpecialCuts
  //
  // G4double maxLength = 2*fTrackerLength, maxTime = 0.1*ns, minEkin = 10*MeV;
  // logicTracker->SetUserLimits(new G4UserLimits(maxStep,maxLength,maxTime,
  //                                               minEkin));
  
  return physiWorld;
}
 
//======================================================================
void DetectorConstruction::ConstructMaterials()
{
  NISTManager=G4NistManager::Instance();

  //-----------------------------
  G4Material* scintMat = new G4Material ("plasticScint", 1.06*g/cm3, 1, kStateSolid);
  scintMat->AddMaterial(G4NistManager::Instance()->FindOrBuildMaterial("G4_POLYSTYRENE"), 1.0);

  //------------------- Material Properties Table for Scint -----------------------------

  //========================================================================//
  //======               Material properties for the Scintillation      ====//
  //====== Refractive index values are also used for cherenkov generation ==//
  //========================================================================//

  const G4int ScintENTRIES = 13;
  G4double ppscint[ScintENTRIES] = {2.48*eV, 2.43*eV, 2.41*eV, 2.40*eV, 2.38*eV, 2.34*eV, 2.26*eV, 
				    2.21*eV, 2.18*eV, 2.10*eV, 2.07*eV, 1.98*eV, 1.91*eV};
    
  G4double rindex_scint[ScintENTRIES]={1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6};
  G4double absorption[ScintENTRIES]={3.5*m, 3.5*m, 3.5*m, 3.5*m, 3.5*m, 3.5*m, 3.5*m, 3.5*m, 3.5*m, 3.5*m, 3.5*m, 3.5*m, 3.5*m};
  G4double ScintilFast[ScintENTRIES]={0.019, 0.037, 0.074, 0.111, 0.148, 0.185, 0.148, 
				      0.110, 0.074, 0.037, 0.033, 0.018, 0.006};

  G4MaterialPropertiesTable* MPT_scint = new G4MaterialPropertiesTable();
  MPT_scint->AddProperty("RINDEX", ppscint, rindex_scint, ScintENTRIES);
  MPT_scint->AddProperty("ABSLENGTH",ppscint,absorption, ScintENTRIES);
  MPT_scint->AddProperty("FASTCOMPONENT",ppscint, ScintilFast, ScintENTRIES);
  MPT_scint->AddConstProperty("SCINTILLATIONYIELD",7100./MeV);
  MPT_scint->AddConstProperty("RESOLUTIONSCALE",1.0);
  MPT_scint->AddConstProperty("FASTTIMECONSTANT", 7.0*ns);
  scintMat->SetMaterialPropertiesTable(MPT_scint);
  scintMat->GetIonisation()->SetBirksConstant(0.07943*mm/MeV);
  G4cout<<"Polystyrene radLength = "<<scintMat->GetRadlen()/cm<<G4endl;


  //----====================
  // build Epoxy
  //----====================
	G4double a = 12.011*g/mole;
  G4Element* elC = new G4Element("Carbon","C",6., a);

	a = 14.007*g/mole;
  G4Element* elN = new G4Element("Nitrogen","N",7., a);

	a = 15.999*g/mole;
  G4Element* elO = new G4Element("Oxygen","O", 8., a);

	a = 1.008*g/mole;
  G4Element* elH = new G4Element("Hydrogen", "H", 1., a);

	G4Material* EPOXY = new G4Material("epoxy", 1.18*g/cm3, 4);
	EPOXY->AddElement(elC, 82.2*perCent);
	EPOXY->AddElement(elH, 9.1*perCent);
	EPOXY->AddElement(elN, 3.2*perCent);
	EPOXY->AddElement(elO, 5.5*perCent);

	//----====================
	// build PMMA
	//----====================
	G4NistManager* fNistMan;
	fNistMan = G4NistManager::Instance();
	const G4int nEntries = 50;

	
   G4double photonEnergy[nEntries] =
   {2.00*eV,2.03*eV,2.06*eV,2.09*eV,2.12*eV,
    2.15*eV,2.18*eV,2.21*eV,2.24*eV,2.27*eV,
    2.30*eV,2.33*eV,2.36*eV,2.39*eV,2.42*eV,
    2.45*eV,2.48*eV,2.51*eV,2.54*eV,2.57*eV,
    2.60*eV,2.63*eV,2.66*eV,2.69*eV,2.72*eV,
    2.75*eV,2.78*eV,2.81*eV,2.84*eV,2.87*eV,
    2.90*eV,2.93*eV,2.96*eV,2.99*eV,3.02*eV,
    3.05*eV,3.08*eV,3.11*eV,3.14*eV,3.17*eV,
    3.20*eV,3.23*eV,3.26*eV,3.29*eV,3.32*eV,
    3.35*eV,3.38*eV,3.41*eV,3.44*eV,3.47*eV};

	std::vector<G4String> elements;
	std::vector<G4int> natoms;
	G4double density;
	G4Material* fPMMA;

	elements.push_back("C");     natoms.push_back(5);
  elements.push_back("H");     natoms.push_back(8);
  elements.push_back("O");     natoms.push_back(2);
  density = 1.190*g/cm3;
  fPMMA = fNistMan->ConstructNewMaterial("PMMA", elements, natoms, density);
  elements.clear();
  natoms.clear();
	
	 G4double refractiveIndexWLSfiber[nEntries] =
   { 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
     1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
     1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
     1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
     1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49};
 
   G4double absWLSfiber[nEntries] =
   {5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,
    5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,
    5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,1.10*m,
    1.10*m,1.10*m,1.10*m,1.10*m,1.10*m,1.10*m, 1.*mm, 1.*mm, 1.*mm, 1.*mm,
     1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm};
 
   G4double emissionFib[nEntries] =
   {0.05, 0.10, 0.30, 0.50, 0.75, 1.00, 1.50, 1.85, 2.30, 2.75,
    3.25, 3.80, 4.50, 5.20, 6.00, 7.00, 8.50, 9.50, 11.1, 12.4,
    12.9, 13.0, 12.8, 12.3, 11.1, 11.0, 12.0, 11.0, 17.0, 16.9,
    15.0, 9.00, 2.50, 1.00, 0.05, 0.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00};
 
   // Add entries into properties table
   G4MaterialPropertiesTable* mptWLSfiber = new G4MaterialPropertiesTable();
   mptWLSfiber->
            AddProperty("RINDEX",photonEnergy,refractiveIndexWLSfiber,nEntries);
   // mptWLSfiber->AddProperty("ABSLENGTH",photonEnergy,absWLSfiber,nEntries);
   mptWLSfiber->AddProperty("WLSABSLENGTH",photonEnergy,absWLSfiber,nEntries);
   mptWLSfiber->AddProperty("WLSCOMPONENT",photonEnergy,emissionFib,nEntries);
   mptWLSfiber->AddConstProperty("WLSTIMECONSTANT", 0.5*ns);
 
   fPMMA->SetMaterialPropertiesTable(mptWLSfiber);

	// CLADDING 1?

	G4Material* fPethylene;

	elements.push_back("C");     natoms.push_back(2);
  elements.push_back("H");     natoms.push_back(4);

  density = 1.190*g/cm3;

  fPethylene = fNistMan->
          ConstructNewMaterial("Pethylene", elements, natoms, density);

  elements.clear();
  natoms.clear();

	G4double refractiveIndexClad1[nEntries] =
  { 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
    1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
    1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
    1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
    1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49};
 
  G4double absClad[nEntries] =
  {20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,
   20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,
   20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,
   20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,
   20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m};

   // Add entries into properties table
  G4MaterialPropertiesTable* mptClad1 = new G4MaterialPropertiesTable();
  mptClad1->AddProperty("RINDEX",photonEnergy,refractiveIndexClad1,nEntries);
  mptClad1->AddProperty("ABSLENGTH",photonEnergy,absClad,nEntries);

  fPethylene->SetMaterialPropertiesTable(mptClad1);

	//----===============
	// build cladding 2
	//----===============
	
	G4Material* fFPethylene;

	elements.push_back("C");     natoms.push_back(2);
  elements.push_back("H");     natoms.push_back(4);
 
  density = 1.430*g/cm3;

  fFPethylene = fNistMan->
          ConstructNewMaterial("FPethylene", elements, natoms, density);

  elements.clear();
  natoms.clear();

	G4double refractiveIndexClad2[nEntries] =
  { 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
    1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
    1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
    1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
    1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42};

 // Add entries into properties table
  G4MaterialPropertiesTable* mptClad2 = new G4MaterialPropertiesTable();
  mptClad2->AddProperty("RINDEX",photonEnergy,refractiveIndexClad2,nEntries);
  mptClad2->AddProperty("ABSLENGTH",photonEnergy,absClad,nEntries);

  fFPethylene->SetMaterialPropertiesTable(mptClad2);


}

//======================================================================

void DetectorConstruction::VisualizationAttributes()
{
  World_VisAtt = new G4VisAttributes();
  World_VisAtt->SetForceWireframe(true);  
  // logicWorld->SetVisAttributes(G4VisAttributes::Invisible);
  logicWorld->SetVisAttributes(World_VisAtt);

  VisAttScintFiber = new G4VisAttributes(true, G4Colour::Cyan());
  VisAttScintFiber->SetForceWireframe(false);
  logicScintFiber->SetVisAttributes(VisAttScintFiber);

	VisAttModule = new G4VisAttributes(true, G4Colour::Grey());
  VisAttModule->SetForceWireframe(false);
  logicModule->SetVisAttributes(VisAttModule);

	VisAttGlueRing = new G4VisAttributes(true, G4Colour::Red());
  VisAttGlueRing->SetForceWireframe(false);
  logicGlueRing->SetVisAttributes(VisAttGlueRing);

	VisAttCladding1 = new G4VisAttributes(true, G4Colour::Blue());
  VisAttCladding1->SetForceWireframe(false);
  logicCladding1->SetVisAttributes(VisAttCladding1);
	
	VisAttCladding2 = new G4VisAttributes(true, G4Colour::White());
  VisAttCladding2->SetForceWireframe(false);
  logicCladding2->SetVisAttributes(VisAttCladding2);

	VisAttGlueBox = new G4VisAttributes(true, G4Colour::Red());
	VisAttGlueBox->SetForceWireframe(false);
	logicGlueBox->SetVisAttributes(VisAttGlueBox);

  // VisAttSensor = new G4VisAttributes(true, G4Colour::White());  
  // VisAttSensor->SetForceWireframe(true);  
  // logicSensor->SetVisAttributes(VisAttSensor);
}
