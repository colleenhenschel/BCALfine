/*
	DetectorConstruction.hh
	This file defines the DetectorConstruction class
	The functions can be found in DetectorConstruction.cc located in the source directory
	ie. /BCALfine/src/DetectorConstruction.cc

	The simulation will create a finegrained BCAL with the following volume types:
			World
			Module			
			ScintFiber
			GlueBox
			GlueRing
			
*/

#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include <math.h> 
#include "G4PVParameterised.hh"
#include "ModuleParameterisation.hh"
#include "BasePlateParameterisation.hh"

#define PI 3.14159265358979

class Fiber;
class LightAbsorber;
class CPID_info;
class G4Box;
//class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4NistManager;
class G4UserLimits;
class G4VisAttributes;

//---------------------------------------------------------------------------
class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  DetectorConstruction(G4double worldSize=0);
  ~DetectorConstruction();
  
public:
  virtual G4VPhysicalVolume* Construct();
  G4double GetWorldFullLength()      {return fWorldLength;}
  
  const G4LogicalVolume* GetWorldLVolume()  const { return logicWorld;}

  G4double GetFieldSize()            { return field_size;}

  G4double GetSourceToCPIDdistance()   { return source_to_cpid_dist;}  

  void SetFieldDimenssion(const G4double &size) { field_size = size;}
  
private:    
  void ConstructMaterials();
  void VisualizationAttributes();

// World Volume:

  G4Box*             	solidWorld;        // pointer to the solid envelope 
  G4LogicalVolume*   	logicWorld;        // pointer to the logical envelope
  G4VPhysicalVolume* 	physiWorld;        // pointer to the physical envelope
 
// Module Volume:

	G4Trd*							solidModule;
	G4LogicalVolume*		logicModule;
	G4VPhysicalVolume*	physiModule;	

	G4Trd*							solidBasePlate;
	G4LogicalVolume*    logicBasePlate;

// Glue Box volumes:

	G4Box*							solidGlueBox;
	G4LogicalVolume*		logicGlueBox;
	G4VPhysicalVolume*	physiGlueBox;	

// Scintillating Fiber Volume:

	G4Tubs*							solidScintFiber;
	G4LogicalVolume*		logicScintFiber;
	G4VPhysicalVolume*	physiScintFiber;


	G4Tubs*							solidCladding1;
	G4LogicalVolume*		logicCladding1;
	G4VPhysicalVolume*	physiCladding1;

	G4Tubs*							solidCladding2;
	G4LogicalVolume*		logicCladding2;
	G4VPhysicalVolume*	physiCladding2;

	G4Tubs*							solidGlueRing;
	G4LogicalVolume*		logicGlueRing;
	G4VPhysicalVolume*	physiGlueRing;

//ModuleParameterisation*		paramModule; //These are needed if you want to remove the base plates
	//G4PVParameterised*				physiModule;

	BasePlateParameterisation*		paramBasePlate;
	G4PVParameterised*						physiBasePlate;


  G4UserLimits*      stepLimit;         // adjusting max. step size for 
                                        // more detailed simulation of 
                                        // the energy deposition
  
  G4NistManager*     NISTManager;

  G4VisAttributes* World_VisAtt;
  G4VisAttributes* VisAttScintFiber;
	G4VisAttributes* VisAttModule;
	G4VisAttributes* VisAttGlueRing;
	G4VisAttributes* VisAttCladding1;
	G4VisAttributes* VisAttCladding2;
	G4VisAttributes* VisAttGlueBox;

  G4double           	fWorldLength;      // Full length of the world volume 
  G4double           	field_size;
  G4double           	source_to_cpid_dist;
  G4double						fiber_length;
	G4double						fiber_radius;
	G4double						module_outter;
	G4double						module_inner;
	G4double						module_length;
	G4double						module_height;
	G4double 						fibervertspacing;
	G4double						halfglueboxheight;
	G4double 						midradius;
	G4double						bp_height; // distance between outter edge of inner baseplate to outter edge of outter baseplate
	G4double						baseplate_inner; // small side width of inner baseplate
	G4double						baseplate_outter; // wide side width of outter baseplate
};

//-----------------------------------------------------------------------------
#endif
