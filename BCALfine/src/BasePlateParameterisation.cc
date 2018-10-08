#include "BasePlateParameterisation.hh"

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

BasePlateParameterisation::BasePlateParameterisation(G4RotationMatrix* r, G4Material* m, G4int n, G4double midrad)
{   
  rot = r;
	mat = m;
	nBasePlates = n;
	midradius = midrad;
}

//________________________________________________________________________________
BasePlateParameterisation::~BasePlateParameterisation()
{
 
}

//________________________________________________________________________________
G4int BasePlateParameterisation::GetMaterialIndex(const G4int parentCopyNo) const
{
  return 1;    
}

//________________________________________________________________________________
void BasePlateParameterisation::ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const
{
	G4double anglebetweenBasePlates = copyNo*7.5;
	G4double modz;
	G4double mody;
	if(anglebetweenBasePlates != 0 && anglebetweenBasePlates != 180)
	{
		modz = 0. + (midradius*sin(anglebetweenBasePlates*PI/180)/(sin(((180-anglebetweenBasePlates)/2)*PI/180)))*cos((anglebetweenBasePlates/2)*PI/180);
		mody = (midradius*sin((90-anglebetweenBasePlates)*PI/180));
	}
	else if (anglebetweenBasePlates == 0)
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
 	Rot->rotateX(90.0*deg); 
	Rot->rotateZ(90.0*deg);
	Rot->rotateZ(90.0*deg);
	Rot->rotateY(copyNo*7.5*deg);
  physVol->SetTranslation(G4ThreeVector(modz*cm, mody*cm, 0.*cm));
  physVol->SetRotation(Rot);
}


//________________________________________________________________________________
G4Material* BasePlateParameterisation::ComputeMaterial(G4VPhysicalVolume*, const G4int, const G4VTouchable *parentTouch)
{
  return mat;
}
//_________________________________________________________________________________
G4Material* BasePlateParameterisation::GetMaterial(G4int idx) const 
{
  return mat;
}
//_________________________________________________________________________________
G4int BasePlateParameterisation::GetNumberOfMaterials() const 
{
  return 1;
}
