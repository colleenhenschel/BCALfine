
#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

#include "G4ParticleGun.hh"
#include "Randomize.hh"
#include "G4Poisson.hh"

#include "DataManager.hh"
#include "MCEvent.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <math.h> 

//----------------------------------------------------------------------------

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* myDC)
  :G4VUserPrimaryGeneratorAction(), myDetector(myDC)
{
  // default particle
  particleTable = G4ParticleTable::GetParticleTable();
  particle = particleTable->FindParticle("mu-");
  
  particleGun = new G4ParticleGun(particle, 1);
  particleGun->SetParticlePosition(G4ThreeVector(0.*mm, 0.*mm, 0.*mm));
  
  //initialize private memebers
  source_to_cpid_dist = myDC->GetSourceToCPIDdistance();
  field_size          = myDC->GetFieldSize();
  minCosTheta = 1./sqrt(1+field_size*field_size/(4.*source_to_cpid_dist*source_to_cpid_dist));
  maxCoordinate = field_size/2.;

  G4cout<< "\n\n MIN_THETA = "<< acos(minCosTheta)/deg <<" deg \n\n"<<G4endl;

  particleDirection = G4ThreeVector(0., 0., 1.);
	particlePosition = G4ThreeVector(0.,0.,0.);
  beamPositionOffset = G4ThreeVector(0.*um, 0.*um, 0.*um);

  particleEnergy = 0.*MeV;
  integral       = 0.; 
  
  //Instantiation of messenger
  messenger = new PrimaryGeneratorMessenger(this);

  fixedEnergyFlag  = true;
  parallelBeamFlag = false;
  if(fixedEnergyFlag==false)InitializeSpectrum("./energy_spectra//pinacle_6MV.txt");

  //Data Management
  dataManager  = DataManager::GetInstance();
  eventPrimary = new EventPrimary;
}

//-----------------------------------------------------------------------------

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete eventPrimary;
  delete particleGun;
  delete messenger;
}

//-----------------------------------------------------------------------------
void PrimaryGeneratorAction::InitializeSpectrum(G4String name)
{
  if(fixedEnergyFlag == true) return;
  G4double energy=0., intensity=0.;

  energyV.clear();
  intensityV.clear();
  intensityVCI.clear();

  std::ifstream spectrum_file;
  
  spectrum_file.open(name.data());
  
  if(!spectrum_file.good())    
    { 
      G4cerr<<G4endl<<__FILE__<<": "<<__LINE__<<G4endl
  	    <<"failed to open spectrum file : \""<<name<<"\""<<G4endl<<G4endl;
      exit(-1);
    }
  else
    {
      G4cout<<"Initializing X-ray spectrum to: "<<name<<G4endl;
      while(!spectrum_file.eof())
  	{
  	  spectrum_file>>energy>>intensity;
	  if(spectrum_file.eof())break;
  	  energyV.push_back(energy*MeV);
	  intensityV.push_back(intensity);
  	}
    }

  spectrum_file.close();

  for(unsigned int i=0; i<intensityV.size(); i++)
    {
      if(i>0)
	integral += 0.5*(intensityV[i] + intensityV[i-1])*(energyV[i]-energyV[i-1]);
      else if(i==0)
	integral += 0.5*(intensityV[i])*(energyV[i]);
    }

  G4double CIvalue = 0.;

  for(unsigned int i=0; i<intensityV.size(); i++)
    {
      if(i>0)
	CIvalue += 0.5*(intensityV[i] + intensityV[i-1])*(energyV[i]-energyV[i-1])/integral;
      else if (i==0)
	CIvalue += 0.5*(intensityV[i])*(energyV[i])/integral;
      intensityVCI.push_back(CIvalue);
    }
}
//-----------------------------------------------------------------------------
void PrimaryGeneratorAction::SampleEnergy(G4double probability)
{
	particleEnergy = 1*GeV + probability*6*GeV; // uniformly generates muons with energy 1 to 7 GeV
}

double cdf(double x) //integral of pdf
{
    return (2/pi)*(x/2 + sin(2*x)/4) + 0.5; //from Wolfram Alpha
}

double inverse_cdf(double u)
{   //bisection, not 100% accurate
    double low  = -halfpi;
    double high = halfpi;
    double epsilon = 1e-10; //any small number, e.g. 1e-15
    while (high - low > epsilon)
    {
        double mid = (low + high) / 2;
        if (cdf(mid) == u) return mid;
        if (cdf(mid) < u) low = mid; else high = mid;
    }
    return (low + high) / 2;
}

//-----------------------------------------------------------------------------
void PrimaryGeneratorAction::SampleDirection()
{
	phi = twopi*G4UniformRand();
	
	// Using bisection method to sample theta
	G4double u = G4UniformRand();

	theta = inverse_cdf(u);
 
  particleDirection = G4ThreeVector(sin(theta)*cos(phi),-cos(theta), sin(theta)*sin(phi));
}

//----------------------------------------------------------------------------
void PrimaryGeneratorAction::SampleBeamStartXY()
{
  G4double x = (G4UniformRand() - 0.5)*600*cm; 
  G4double z = (G4UniformRand() - 0.5)*600*cm; 
	
	particlePosition = G4ThreeVector(x, 118.*cm, z); 

  particleGun->SetParticlePosition(G4ThreeVector(x, 118.*cm, z));   
}


//-----------------------------------------------------------------------------

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{ 
  if(particle==0) return;

  dataManager->CleanUpEvent();
  dataManager->SetEventNumber(anEvent->GetEventID());
  dataManager->SetEventName("singlePhotonEvent");

  if(fixedEnergyFlag==false)
    { SampleEnergy(G4UniformRand()); }
  
	SampleBeamStartXY();
	SampleDirection();
	bool checkevent = true; 
	while(checkevent)
	{
		G4double y0 = particlePosition.getY();
		G4double x0 = particlePosition.getX();
		G4double z0 = particlePosition.getZ();
		G4double vx = particleDirection.x();
		G4double vy = particleDirection.y();
		G4double vz = particleDirection.z();
		if(vy == 0)
		{	
			SampleBeamStartXY();
			SampleDirection();
			continue;
		}			
		G4double t1 = (0*cm - y0)/vy;
		G4double t2 = (45*cm - y0)/vy;
		G4double t3 = (-45*cm - y0)/vy;
		G4double t4 = (70*cm - y0)/vy;
		G4double t5 = (-70*cm - y0)/vy;
		G4double x1 = x0 + t1*vx; 
		G4double z1 = z0 + t1*vz;
		G4double x2 = x0 + t2*vx;
		G4double z2 = z0 + t2*vz;
		G4double x3 = x0 + t3*vx;
		G4double z3 = z0 + t3*vz;
		G4double x4 = x0 + t4*vx;
		G4double z4 = z0 + t4*vz;
		G4double x5 = x0 + t5*vx;
		G4double z5 = z0 + t5*vz;
		if((fabs(x1) < 92.*cm && fabs(z1) < 200.*cm)||(fabs(x2) < 80*cm && fabs(z2) < 200.*cm)|| (fabs(x3) < 80*cm && fabs(z3) < 200.*cm) || (fabs(x4) < 60*cm && fabs(z4) < 200.*cm) || (fabs(x5) < 60*cm && fabs(z5) < 200.*cm)) // checking if event passes through BCAL at different planes
		{
			checkevent = false;
		}
		else
		{
			SampleBeamStartXY();
			SampleDirection();
		}
				
	}	
  
  particleGun->SetParticleMomentumDirection(particleDirection);
  particleGun->SetParticleEnergy(particleEnergy);

  eventPrimary->Clear("");
  
  eventPrimary->SetTrackID(0);
  eventPrimary->SetName(particle->GetParticleName());
  eventPrimary->SetMass(particle->GetPDGMass());
  eventPrimary->SetCharge(particle->GetPDGCharge());
  eventPrimary->SetTotalEnergy(particleEnergy/MeV);
  //   eventPrimary->SetLocalTime();
  //   eventPrimary->SetGlobalTime();
  //eventPrimary->SetPolarization(TVector3(fBeamPolarization.x(),
  //   					 fBeamPolarization.y(),
  //					   fBeamPolarization.z())); 
  eventPrimary->SetMomentum(TVector3(particleEnergy*particleDirection.x()/MeV,
  				     particleEnergy*particleDirection.y()/MeV,
  				     particleEnergy*particleDirection.z()/MeV));
	eventPrimary->SetVertexPosition(TVector3(particlePosition.getX(), particlePosition.getY(), particlePosition.getZ()));
	eventPrimary->SetTheta(theta);
  
  dataManager->AddPrimary(eventPrimary);


  particleGun->GeneratePrimaryVertex(anEvent);
}

