

#include "SteppingAction.hh"
#include "G4SteppingManager.hh"

#include "G4SystemOfUnits.hh"

#include "G4ThreeVector.hh"
#include "G4ios.hh"
#include "G4GeometryTolerance.hh"
#include "G4UnitsTable.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"

#include "DataManager.hh"
#include "MCEvent.hh"

#include "G4OpticalPhoton.hh"
#include "G4TransportationManager.hh"

//------------------------------------------------------------------------------

SteppingAction::SteppingAction()
{ 
  dataManager = DataManager::GetInstance();

  aTrack = 0;
  preStepPoint  = 0;
  preStepPV     = 0;
  preStepPVname = "";
  hitPos  = G4ThreeVector(0.*mm, 0.*mm, 0.*mm);
  vertPos = G4ThreeVector(0.*mm, 0.*mm, 0.*mm);

  proc_name = "";
  particle_name = "";

  postStepPoint = 0;
  postStepPV    = 0;
  preStepPVname ="";

  trackID=0;
	curcopynum = 0;
	curmodule = 0;
  particleHit = new fiberHit();
}

//-------------------------------------------------------------------------------
SteppingAction::~SteppingAction()
{
  if(particleHit) {delete particleHit; particleHit=0;}
}
//-------------------------------------------------------------------------------
void SteppingAction::Clear()
{
  aTrack = 0;
  preStepPoint = 0;
  preStepPV = 0;
  preStepPVname = "";
  hitPos = G4ThreeVector(0.*mm, 0.*mm, 0.*mm);
  vertPos = G4ThreeVector(0.*mm, 0.*mm, 0.*mm);

  proc_name = "";
  particle_name = "";

  postStepPoint = 0;
  postStepPV = 0;
  postStepPVname = "";

  trackID = 0;
	//curcopynum = 0;
	//curmodule = 0;
 	particleHit->Clear("C");
}


void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  Clear();
  aTrack = aStep->GetTrack();
  preStepPoint     = aStep->GetPreStepPoint();
  preStepPV  = preStepPoint->GetTouchableHandle()->GetVolume();  
  preStepPVname = preStepPV->GetName();

	postStepPoint = aStep->GetPostStepPoint();
  postStepPV    = postStepPoint->GetTouchableHandle()->GetVolume();
  if(postStepPV!=0)
    postStepPVname = postStepPV->GetName();
  else
    postStepPVname ="";

	if(preStepPVname.contains("pvScintFiber") )
	{
		particleHit->SetTrackID(aTrack->GetTrackID());
	  particleHit->SetParentID(aTrack->GetParentID());
	  particleHit->SetGlobalTime(aTrack->GetGlobalTime()/ns);
	  particleHit->SetLocalTime(aTrack->GetLocalTime()/ns);
		particleHit->SetEnergyDep(aStep->GetTotalEnergyDeposit()/MeV);
		vertPos = preStepPoint->GetPosition();      
	  particleHit->SetVertexPos(TVector3(vertPos.x()/cm, 
						     vertPos.y()/cm, 
						     vertPos.z()/cm));
		if(aTrack->GetTrackID()==1) 
	    particleHit->SetCreatorProc("generator");
	 	else 
	  {
	    proc_name = aTrack->GetCreatorProcess()->GetProcessName();
	    particleHit->SetCreatorProc(proc_name.data());
		}
	 	particle_name = aTrack->GetParticleDefinition()->GetParticleName();
	 	particleHit->SetTrackName(particle_name.data());
		dataManager->AddFiberHit(particleHit);
		
	}

}
