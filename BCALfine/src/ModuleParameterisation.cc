#include "ModuleParameterisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4VTouchable.hh"
#include "G4NistManager.hh"
#include "MCEvent.hh"
#include "G4PVPlacement.hh"
#include <vector>
#include <algorithm>
#include "G4SystemOfUnits.hh"

#include <iostream>
#include <fstream>
#include <math.h> 
using namespace std;

#define PI 3.14159265358979

//______________________________________________________________________________

ModuleParameterisation::ModuleParameterisation(G4RotationMatrix* r, G4Material* m, G4int n, G4double midrad)
{   
  rot = r;
	mat = m;
	nmodules = n;
	midradius = midrad;
}

//________________________________________________________________________________
ModuleParameterisation::~ModuleParameterisation()
{
 
}

//________________________________________________________________________________
G4int ModuleParameterisation::GetMaterialIndex(const G4int parentCopyNo) const
{
  return 1;    
}

//________________________________________________________________________________
void ModuleParameterisation::ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const
{
	G4double anglebetweenmodules = copyNo*7.5;
	G4double modz;
	G4double mody;
	if(anglebetweenmodules != 0 && anglebetweenmodules != 180)
	{
		modz = 0. + (midradius*sin(anglebetweenmodules*PI/180)/(sin(((180-anglebetweenmodules)/2)*PI/180)))*cos((anglebetweenmodules/2)*PI/180);
		mody = (midradius*sin((90-anglebetweenmodules)*PI/180));
	}
	else if (anglebetweenmodules == 0)
	{ 
		modz = 0;
		mody = midradius;
	}
	else
	{
		modz = 0;
		mody = -midradius;
	}
	//cout << "x, y, z = 0, " << mody << ", " << modz << endl;
	G4RotationMatrix *Rot = new G4RotationMatrix;
 	Rot->rotateX(90.0*deg - copyNo*7.5*deg); 
	Rot->rotateZ(90.0*deg);
  physVol->SetTranslation(G4ThreeVector(0.*cm, mody*cm, modz*cm));
  physVol->SetRotation(Rot);
}


//________________________________________________________________________________
G4Material* ModuleParameterisation::ComputeMaterial(G4VPhysicalVolume*, const G4int, const G4VTouchable *parentTouch)
{
  return mat;
}
//_________________________________________________________________________________
G4Material* ModuleParameterisation::GetMaterial(G4int idx) const 
{
  return mat;
}
//_________________________________________________________________________________
G4int ModuleParameterisation::GetNumberOfMaterials() const 
{
  return 1;
}
