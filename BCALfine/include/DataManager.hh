#ifndef DataManager_h
#define DataManager_h 1

#include "globals.hh"
#include "MCEvent.hh"

class TTree;
class TFile;

class DataManager //a singleton class
{
public:
    static DataManager* GetInstance();

private:
  DataManager():fFile(NULL), fFileName(0), fEventTree(NULL), fEvent(NULL) {}
  DataManager(DataManager const& );             //not defined, not copyable
  DataManager& operator= (DataManager const& ); //not defined, not assignable
  ~DataManager();

private:
  TFile*       fFile; 
  TString      fFileName;
  TTree*       fEventTree;
  // static Int_t CerenTrackID;
  // static Int_t ScintTrackID;
  static Int_t TrackID;

  MCEvent*     fEvent;
  MODULE_info*  ModuleInfo;

public:
  void open(Int_t RN);
  void open(const char* name);

  void close();
  void SetFileName(TString name)      { fFileName=name; }

  void CleanUpEvent() { 
    if(fEvent) fEvent->Clear("C"); 
    // CerenTrackID =-1;
    // ScintTrackID =-1;
    TrackID =-1;
  }

  // void SetCerenkovTrackID(const Int_t &id) {CerenTrackID = id;}
  // void SetScintTrackID(const Int_t &id)    {ScintTrackID = id;}
  void SetTrackID(const Int_t &id)  {TrackID = id;}

  void SetEventName(TString name) { fEvent->SetEventName(name); }
  void SetEventNumber(Int_t numb) { fEvent->SetEventNumber(numb);}

  void AddNumOptCreated(const Int_t &n){fEvent->AddNumOptCreated(n);}
  void AddPrimary(EventPrimary* primary); 
  void AddFiberHit(fiberHit *hitPtr);   
  //void AddOptTransmitted(opticalHit *hitPtr);   

  void SetEvent(MCEvent* ev) { fEvent = ev; }

  void SetModuleInfo(MODULE_info* xy) { ModuleInfo = xy; }
  MODULE_info *GetModuleInfo()        { return ModuleInfo; }
  
  // Int_t GetCerenkovTrackID() {return CerenTrackID;}
  // Int_t GetScintTrackID()    {return ScintTrackID;}
  Int_t GetTrackID()   {return TrackID;}

  void FillEvent();
//  protected:
  MCEvent* GetEvent() { return fEvent; }
};

namespace
{ struct ForceSingletonInitialization
  { ForceSingletonInitialization() { DataManager::GetInstance(); } } GetInstance;
}

#endif




