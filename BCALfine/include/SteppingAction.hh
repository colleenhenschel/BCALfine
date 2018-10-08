            
#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class DataManager;
class fiberHit;
class opticalHit;
class G4Track;
class G4StepPoint;
class G4VPhysicalVolume;
class G4LogicalVolume;


//-----------------------------------------------------------------------------

class SteppingAction : public G4UserSteppingAction
{
public:
  SteppingAction();
  ~SteppingAction();

  void UserSteppingAction(const G4Step*);

  void Clear();

private:
  DataManager *dataManager;
  G4Track*     aTrack;

  G4StepPoint* preStepPoint;
  G4VPhysicalVolume* preStepPV;
  G4String preStepPVname;

	//G4LogicalVolume* mothervol;
  
	G4ThreeVector hitPos;
  G4ThreeVector vertPos;

  G4String proc_name;
  G4String particle_name;

  G4StepPoint* postStepPoint;
  G4VPhysicalVolume* postStepPV;
  G4String postStepPVname;
  G4int trackID;

	G4int curcopynum = 0;
	G4int curmodule = 0;
	
	G4int eventid = 0;

  fiberHit   *particleHit;
};

//-----------------------------------------------------------------------------

#endif
