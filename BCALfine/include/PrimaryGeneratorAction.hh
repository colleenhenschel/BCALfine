
 
#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

#include <vector>

class G4ParticleGun;
class DetectorConstruction;
class G4Event;
class G4ParticleTable;
class G4ParticleDefinition;

class PrimaryGeneratorMessenger;
class DataManager;
class EventPrimary;

//-------------------------------------------------------------------------------
 
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction(DetectorConstruction* myDC);    
  ~PrimaryGeneratorAction();
  
public:
  virtual void GeneratePrimaries(G4Event* anEvent);
  
  void SetParticleEnergy(G4double e)   {particleEnergy = e;}
  void SetFieldSize(G4double fs)       
  {
    field_size = fs;
    minCosTheta = 1./sqrt(1+field_size*field_size/(4.*source_to_cpid_dist*source_to_cpid_dist));
    maxCoordinate = field_size/2.;
  }
  void SetFixedEnergyFlag(G4bool flag) {fixedEnergyFlag = flag;}
  void SetParallelBeamFlad(G4bool flag) {parallelBeamFlag = flag;}
  void SetBeamPositionOffset(G4ThreeVector offset) {beamPositionOffset = offset;}   
private:
  void SampleEnergy(G4double probability); 
  void SampleDirection();
  void SampleBeamStartXY();

public:  
  void InitializeSpectrum(G4String name);

  G4double GetSpectrumIntegral() {return integral;}

private:
  DetectorConstruction  *myDetector;
  G4ParticleTable       *particleTable;
  G4ParticleDefinition  *particle;
  G4ParticleGun         *particleGun;

  DataManager           *dataManager;
  EventPrimary          *eventPrimary;

  // particle direction sampler
  G4double               source_to_cpid_dist;
  G4double               field_size;
  G4double               minCosTheta;
  G4double               maxCoordinate;
	G4double							 theta;
	G4double 							 phi;
  G4ThreeVector          particleDirection;
	G4ThreeVector					 particlePosition;
  G4ThreeVector          beamPositionOffset;
  
  G4bool                 parallelBeamFlag;


  // particle energy sampler
  G4double               particleEnergy;
  G4bool                 fixedEnergyFlag;
  G4double               integral;    //used by generator internally

  std::vector<G4double> energyV;       //used by generator internally
  std::vector<G4double> intensityV;    //used by generator internally
  std::vector<G4double> intensityVCI;  //used by generator internally
  std::pair<std::vector<G4double>::iterator, std::vector<G4double>::iterator> position;

  // Messenger for communication with UI:
  PrimaryGeneratorMessenger* messenger;
};

//--------------------------------------------------------------------------------

#endif


