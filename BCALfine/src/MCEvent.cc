
#include "MCEvent.hh"

ClassImp(fiberHit)
ClassImp(opticalHit)
ClassImp(EventPrimary)
ClassImp(MCEvent)
ClassImp(MODULE_info)

//-------------------------------------------------------------------------------
fiberHit::fiberHit():TObject()
{
  trackID         = -1;
  parentID        = -1;
  local_time      = 0.;
  global_time     = 0.;
  edep            = 0.;
  vertexPos       = TVector3();
  creator_process = "";
  track_name      = "";
}

//--------------------------------------------------------------------------------
fiberHit::fiberHit(const fiberHit & orig):TObject(orig)
{
  trackID         = orig.trackID;
  parentID        = orig.parentID;
  local_time      = orig.local_time;
  global_time     = orig.global_time;
  edep            = orig.edep;
  vertexPos       = orig.vertexPos;
  creator_process = orig.creator_process;
  track_name      = orig.track_name;
}
//--------------------------------------------------------------------------------

void fiberHit::Clear(Option_t *)
{
  trackID         = -1;
  parentID        = -1;
  local_time      = 0.;
  global_time     = 0.;
  edep            = 0.;
  vertexPos       = TVector3();
  creator_process = "";
  track_name      = "";
}

//-------------------------------------------------------------------------------
opticalHit::opticalHit():TObject()
{
  trackID         = -1;
  parentID        = -1;
  local_time      = 0.;
  global_time     = 0.;
  vertexPos       = TVector3();
  hitPos          = TVector3();
  enenergy        = 0.;
}

//--------------------------------------------------------------------------------
opticalHit::opticalHit(const opticalHit & orig):TObject(orig)
{
  trackID         = orig.trackID;
  parentID        = orig.parentID;
  local_time      = orig.local_time;
  global_time     = orig.global_time;
  vertexPos       = orig.vertexPos;
  hitPos          = orig.hitPos;
  enenergy        = orig.enenergy;
}
//--------------------------------------------------------------------------------

void opticalHit::Clear(Option_t *)
{
  trackID         = -1;
  parentID        = -1;
  local_time      = 0.;
  global_time     = 0.;
  vertexPos       = TVector3();
  hitPos          = TVector3();
  enenergy        = 0.;
}

//________________________________________________________________________________
EventPrimary::EventPrimary(): TObject()
{
  fTrackID       = 0;
  fName          = ""; 
  fMass          = 0.;
  fCharge        = 0.;
  fTotalEnergy   = 0.; 
  fMomentum.SetXYZ(0.,0.,0.);
  fVertextPos.SetXYZ(0., 0., 0.);
	fTheta         = 0.;
}

//______________________________________________________________________________
EventPrimary::EventPrimary(const EventPrimary& orig):TObject(orig)
{
  fTrackID       = orig.fTrackID;
  fName          = orig.fName; 
  fMass          = orig.fMass;
  fCharge        = orig.fCharge;
  fTotalEnergy   = orig.fTotalEnergy; 
  fMomentum      = orig.fMomentum;
  fVertextPos    = orig.fVertextPos;
	fTheta         = orig.fTheta;
}

//________________________________________________________________________________
void EventPrimary::Clear(Option_t*)
{
  fTrackID       = 0; 
  fName          = "";
  fMass          = 0; 
  fCharge        = 0.;
  fTotalEnergy   = 0.;
  fMomentum.SetXYZ(0.,0.,0.);
  fVertextPos.SetXYZ(0., 0., 0.);
	fTheta         = 0.;
}

//______________________________________________________________________________
MCEvent::MCEvent() : TObject()
{
  fEventName        = "";
  fEventNumber      = 0;
  fN_OptCreated     = 0;
  fN_Primaries      = 0;
  fN_FiberHits      = 0;
  fN_OptTransmitted = 0;
  fEventPrimaries   = new TClonesArray("EventPrimary", 20);
  fFiberHits        = new TClonesArray("fiberHit", 100);
  //fOpticalHits      = new TClonesArray("opticalHit", 100);
}
//______________________________________________________________________________

MCEvent::MCEvent(const MCEvent& orig) : TObject(orig)
{
  fEventName          = orig.fEventName;
  fEventNumber        = orig.fEventNumber;
  fN_OptCreated       = orig.fN_OptCreated;
  fN_Primaries        = orig.fN_Primaries;
  fN_FiberHits        = orig.fN_FiberHits;
  fN_OptTransmitted   = orig.fN_OptTransmitted;
  fEventPrimaries     = orig.fEventPrimaries;
  fFiberHits          = orig.fFiberHits;
  //fOpticalHits        = orig.fOpticalHits;
}

//______________________________________________________________________________
MCEvent::~MCEvent()
{
  Clear("C");
  if(fEventPrimaries)    {fEventPrimaries->Clear(); delete fEventPrimaries; fEventPrimaries=0;}  
  if(fFiberHits)       {fFiberHits->Clear(); delete fFiberHits; fFiberHits=0;}
  //if(fOpticalHits)       {delete fOpticalHits; fOpticalHits=0;}
}

//______________________________________________________________________________
EventPrimary* MCEvent::AddPrimary(EventPrimary* primaryPtr)
{
   TClonesArray &evPrimaries = *fEventPrimaries;
   EventPrimary *primary = new(evPrimaries[fN_Primaries++]) (EventPrimary)(*primaryPtr);
  
   return primary;
}

//______________________________________________________________________________
fiberHit *MCEvent::AddFiberHit(fiberHit *hitPtr)
{
   TClonesArray &tgtHits = *fFiberHits;
   fiberHit* hit = new(tgtHits[fN_FiberHits++]) (fiberHit)(*hitPtr);

   return hit;
}

//______________________________________________________________________________
//opticalHit *MCEvent::AddOptTransmitted(opticalHit *hitPtr)
//{
//   TClonesArray &tgtHits = *fOpticalHits;
//   opticalHit *hit = new(tgtHits[fN_OptTransmitted++]) (opticalHit)(*hitPtr);

//   return hit;
//}

//______________________________________________________________________________
void MCEvent::Clear(Option_t *opt)
{
  fEventName          = "";
  fEventNumber        = 0;
  fN_OptCreated       = 0;
  fN_Primaries        = 0;
  fN_FiberHits      = 0;
  fN_OptTransmitted   = 0;
  if(fEventPrimaries)fEventPrimaries->Clear(opt);
  if(fFiberHits)fFiberHits->Clear(opt);
  //if(fOpticalHits)fOpticalHits->Clear(opt);
}

//_______________________________________________________________________________
MODULE_info::MODULE_info(const Int_t& num)    
{
  SetNumFibers(num);
  fNFibers1D =0;
  
  fFiber_size_x=0.;
  fFiber_size_y=0.;
  fFiber_size_z=0.;

  reflector_thickness = 0.; 
  side_wall_thickness = 0.;   
  light_sensor_thickness=0.;
  light_guide_thickness=0.;
  light_guide_assembly_gap=0.;
  fiber_to_fiber_gap=0.;
  distance_from_centre=0.;

  fAssembly_size_x =0.;
  fAssembly_size_y =0.;
  fAssembly_size_z =0.;

  fSourceRadius = 0;
  vSourcePos.SetXYZ(0.,0.,0.);
}
//-------------------------------------------------------------------------------
MODULE_info::MODULE_info(const MODULE_info &orig):TNamed(orig)
{
  fFiberPosVector      = orig.fFiberPosVector;
  fNFibers1D = orig.fNFibers1D;
  
  fFiber_size_x=orig.fFiber_size_x;
  fFiber_size_y=orig.fFiber_size_y;
  fFiber_size_z=orig.fFiber_size_z;

  reflector_thickness      = orig.reflector_thickness; 
  side_wall_thickness      = orig.side_wall_thickness;   
  light_sensor_thickness   = orig.light_sensor_thickness;
  light_guide_thickness    = orig.light_guide_thickness;
  light_guide_assembly_gap = orig.light_guide_assembly_gap;
  fiber_to_fiber_gap   = orig.fiber_to_fiber_gap;
  distance_from_centre     = orig.distance_from_centre;

  fAssembly_size_x = orig.fAssembly_size_x;
  fAssembly_size_y = orig.fAssembly_size_y;
  fAssembly_size_z = orig.fAssembly_size_z;

  fSourceRadius = orig.fSourceRadius;
  vSourcePos = orig.vSourcePos;
}
