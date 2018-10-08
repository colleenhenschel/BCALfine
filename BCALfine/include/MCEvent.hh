#ifndef MCEvent_H
#define MCEvent_H
                                            
#include "TObject.h"
#include "TClonesArray.h"
#include "TRefArray.h"
#include "TVector3.h"


//--------------------------------------------------------------------
class fiberHit : public TObject {

private:  
  Int_t    trackID;
  Int_t    parentID;
  Double_t local_time;
  Double_t global_time;
  Double_t edep;
  TVector3 vertexPos;
  TString  creator_process;
  TString  track_name;

public:
  fiberHit();
  fiberHit(const fiberHit& orig);
  virtual ~fiberHit() {}

  virtual void Clear(Option_t*);

  void AddEnergyDep(const Double_t& e) {edep += e;}
  void SetTrackID(const Int_t &id)      {trackID = id;}
  void SetParentID(const Int_t &id)     {parentID = id;}
  void SetGlobalTime(const Double_t &t) {global_time = t;}
  void SetLocalTime(const Double_t &t)  {local_time = t;}
  void SetEnergyDep(const Double_t &e)  {edep =e;}
  void SetVertexPos(const TVector3 &v)  {vertexPos = v;}
  void SetCreatorProc(const TString &s) {creator_process = s;}
  void SetTrackName(const TString &s)   {track_name = s;}

  Int_t    GetTrackID()     {return trackID;}
  Int_t    GetParentID()    {return parentID;}
  Double_t GetGlobalTime()  {return global_time;}
  Double_t GetLocalTime()   {return local_time;}
  Double_t GetEnergyDep()   {return edep;}
  TVector3 GetVertexPos()   {return vertexPos;}
  TString  GetCreatorProc() {return creator_process;}
  TString  GetTrackName()   {return track_name;}

  ClassDef(fiberHit,1)  //[Analyze] Fiber Event /*do not remove this comment*/
};

//------------------------------------------------------------
class opticalHit : public TObject {

private:  
  Int_t    trackID;
  Int_t    parentID;
  Double_t local_time;
  Double_t global_time;
  Double_t enenergy;
  TVector3 vertexPos;
  TVector3 hitPos;

public:
  opticalHit();
  opticalHit(const opticalHit& orig);
  virtual ~opticalHit() {}

  virtual void Clear(Option_t*);

  void SetTrackID(const Int_t &id)      {trackID = id;}
  void SetParentID(const Int_t &id)     {parentID = id;}
  void SetGlobalTime(const Double_t &t) {global_time = t;}
  void SetLocalTime(const Double_t &t)  {local_time = t;}
  void SetVertexPos(const TVector3 &v)  {vertexPos = v;}
  void SetHitPos(const TVector3 &v)     {hitPos = v;}
  void SetEnergy(const Double_t &e)     {enenergy = e;}

  Int_t    GetTrackID()    {return trackID;}
  Int_t    GetParentID()   {return parentID;}
  Double_t GetGlobalTime() {return global_time;}
  Double_t GetLocalTime()  {return local_time;}
  TVector3 GetVertexPos()  {return vertexPos;}
  TVector3 GetHitPos()     {return hitPos;}
  Double_t GetWaveLnegth() {return 1240./enenergy;}
  ClassDef(opticalHit,1)  //[Analyze] Fiber Event /*do not remove this comment*/
};


//------------------------------------------------------------------------------
class EventPrimary : public TObject 
{
public:
  EventPrimary();
  EventPrimary(const EventPrimary& orig);
  virtual ~EventPrimary() {};

  virtual void Clear(Option_t*);

  void SetTrackID(const Int_t &id)                { fTrackID       = id;       }
  void SetName(const TString &name)               { fName          = name;     }
  void SetMass(const Double_t &mass)              { fMass          = mass;     } 
  void SetCharge(const Double_t &charge)          { fCharge        = charge;   }
  void SetTotalEnergy(const Double_t &energy)     { fTotalEnergy   = energy;   }
  void SetMomentum(const TVector3 &momentum)      { fMomentum      = momentum; }
  void SetVertexPosition(const TVector3 &pos)     { fVertextPos    = pos;      }
	void SetTheta(const Double_t &theta)						{ fTheta         = theta;    }

  Int_t    GetPrimaryTrackID() { return fTrackID;      }
  TString  GetPrimaryName()    { return fName;         }
  Double_t GetMass()           { return fMass;         }
  Double_t GetCharge()         { return fCharge;       }
  Double_t GetTotalEnergy()    { return fTotalEnergy;  }
  TVector3 GetMomentum()       { return fMomentum;     }
  TVector3 GetVertexPosition() { return fVertextPos;   } 

private:
  Int_t    fTrackID; 
  TString  fName;
  Double_t fMass;
  Double_t fCharge;
  Double_t fTotalEnergy;
  TVector3 fVertextPos;
  TVector3 fMomentum;
	Double_t fTheta;

  ClassDef(EventPrimary,1)  //[Analyze] event primary /*do not remove this comment*/
};


class MCEvent : public TObject {

private:  
  TString        fEventName;       //name in character format
  Int_t          fEventNumber;  
  Int_t          fN_OptCreated;

  Int_t          fN_Primaries;
  Int_t          fN_FiberHits;  
  Int_t          fN_OptTransmitted;

  TClonesArray*  fEventPrimaries;    //->array with primaries
  TClonesArray*  fFiberHits;       //->array with fiber hits /*do not remove comment*/
  //TClonesArray*  fOpticalHits;      //->array of transmitted Cerenkov photons/*do not remove this comment*/

public:
  MCEvent();
  MCEvent(const MCEvent& orig);
  virtual ~MCEvent();

  void Clear(Option_t *option ="");
  
  void SetEventNumber(Int_t e_num )     { fEventNumber = e_num; }
  void SetEventName(TString name)       { fEventName   = name; }  

  void          AddNumOptCreated(const Int_t &n){fN_OptCreated += n;}
  EventPrimary *AddPrimary(EventPrimary *primaryPtr);
  fiberHit   *AddFiberHit(fiberHit *hitPtr);
  //opticalHit   *AddOptTransmitted(opticalHit *hitPtr);

  TString       GetEventName()      {return fEventName;}
  Int_t         GetEventNumber()    {return fEventNumber;}
  Int_t         GetN_OptCreated()   {return fN_OptCreated;}
  Int_t         GetN_Primaries()    {return fN_Primaries;}
  Int_t         GetN_FiberHits()  {return fN_FiberHits;}
  Int_t         GetN_OptTransm()    {return fN_OptTransmitted;}

  TClonesArray *GetPrimaries()       { return fEventPrimaries; }
  TClonesArray *GetFiberHits()     { return fFiberHits;}
  //TClonesArray *GetOptTransmitted()  { return fOpticalHits;}

  ClassDef(MCEvent,1)  //[Analyze] MCEvent structure /*do not remove this comment*/
};


#include <utility>
#include <vector>

typedef std::pair<Double_t, Double_t> Coordintate_t;
typedef std::vector<Coordintate_t>  FiberXY_t;

class MODULE_info:public TNamed{
public:
  MODULE_info(){}
  MODULE_info(const Int_t& num);
  MODULE_info(const MODULE_info &orig);
  ~MODULE_info(){};

public:
  Int_t GetNumFibers() const {return fFiberPosVector.size();}

  Double_t GetFiberX(const int &id) const { return (fFiberPosVector.at(id)).first;}
  Double_t GetFiberY(const int &id) const { return (fFiberPosVector.at(id)).second;}

  Double_t GetFiberSizeX(){return fFiber_size_x;}
  Double_t GetFiberSizeY(){return fFiber_size_y;}
  Double_t GetFiberSizeZ(){return fFiber_size_z;}

  Double_t GetReflectorThickness()     {return reflector_thickness;}
  Double_t GetSideWallThickness()      {return side_wall_thickness;}
  Double_t GetLightSensorThickness()   {return light_sensor_thickness;}
  Double_t GetLightGuideThickness ()   {return light_guide_thickness;}
  Double_t GetLightGuideAssemblyGap()  {return light_guide_assembly_gap;}
  Double_t GetFiberToFiberGap()    {return fiber_to_fiber_gap;}
  Double_t GetDistanceToCentere()      {return distance_from_centre;}

  Int_t GetNumFibers1D(){return fNFibers1D;}

  Double_t GetAssemblySizeX(){return fAssembly_size_x;}
  Double_t GetAssemblySizeY(){return fAssembly_size_y;}
  Double_t GetAssemblySizeZ(){return fAssembly_size_z;}

  Double_t GetSourceRadius()      {return fSourceRadius;}
  TVector3 GetSourcePos()         {return vSourcePos;}
  Bool_t   GetSourcePhantomFlag() {return fSourcePhantomFlag;}

public:
  void SetCenterXY(const Int_t &id, const Double_t &X, const Double_t &Y)
  { Coordintate_t xy = std::make_pair(X, Y); fFiberPosVector[id]=xy;}
  void SetFiberSizeX(const Double_t& x)  {fFiber_size_x = x;}
  void SetFiberSizeY(const Double_t& y)  {fFiber_size_y = y;}
  void SetFiberSizeZ(const Double_t& z)  {fFiber_size_z = z;}

  void SetReflectorThickness(const Double_t &t)   {reflector_thickness=t;}
  void SetSideWallThickness(const Double_t &t)    {side_wall_thickness =t;}
  void SetLightSensorThickness(const Double_t &t) {light_sensor_thickness = t;}
  void SetLightGuideThickness (const Double_t &t) {light_guide_thickness = t;}
  void SetLightGuideAssemblyGap(const Double_t &t){light_guide_assembly_gap =t;}
  void SetFiberToFiberGap(const Double_t &t)  {fiber_to_fiber_gap = t;}
  void SetDistanceToCentere(const Double_t &d)    {distance_from_centre = d;}

  void SetNumFibers1D(const Int_t &n){fNFibers1D = n;}

  void SetAssemblySizeX(const Double_t &x)   {fAssembly_size_x = x;}
  void SetAssemblySizeY(const Double_t &y)   {fAssembly_size_y = y;}
  void SetAssemblySizeZ(const Double_t &z)   {fAssembly_size_z = z;}
  void SetSourceRadius(const Double_t &r)    {fSourceRadius = r;}
  void SetSourcePos(const TVector3 &v)       {vSourcePos = v;}
  void SetSourcePhantomFlag(const Bool_t &f) {fSourcePhantomFlag = f;}

private:
  void SetNumFibers(const Int_t& num_fibers) {fFiberPosVector.resize(num_fibers);}

private:
  
  FiberXY_t fFiberPosVector;

  Double_t fFiber_size_x;
  Double_t fFiber_size_y;
  Double_t fFiber_size_z;

  Double_t reflector_thickness; 
  Double_t side_wall_thickness;   
  Double_t light_sensor_thickness;
  Double_t light_guide_thickness;
  Double_t light_guide_assembly_gap;
  Double_t fiber_to_fiber_gap;
  Double_t distance_from_centre;

  Int_t fNFibers1D;

  Double_t fAssembly_size_x;
  Double_t fAssembly_size_y;
  Double_t fAssembly_size_z;

  Double_t fSourceRadius;
  TVector3 vSourcePos;
  Bool_t   fSourcePhantomFlag;
  ClassDef(MODULE_info,1)
};
#endif


