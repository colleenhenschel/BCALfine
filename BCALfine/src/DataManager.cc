
#include "DataManager.hh"

#include "TFile.h"
#include "TTree.h"
#include "TString.h"


DataManager* DataManager::GetInstance()
{ static DataManager *pointerToTheSingletonInstance = new DataManager();
  return pointerToTheSingletonInstance;
}

// Int_t DataManager::CerenTrackID = -1;
// Int_t DataManager::ScintTrackID = -1;
Int_t DataManager::TrackID = -1;
//------------------------------------------------------------------------

DataManager::~DataManager()
{
  if(fEvent)     delete fEvent;
  if(fEventTree) delete fEventTree;
  if(fFile)      delete fFile;
}

//-----------------------------------------------------------------------
void DataManager::open(G4int RN)
{
  fFileName = "mc_data_";
  std::ostringstream tmp_str;
  tmp_str << RN;
  fFileName+= tmp_str.str();
  fFileName+= ".root";
  
  fFile = new TFile(fFileName, "RECREATE");

  G4String tree_comment="tree, Run # ";
  tree_comment+=tmp_str.str();

  fEventTree = new TTree("event_tree", tree_comment);
  fEventTree->SetDirectory(fFile);
  fEventTree->SetAutoSave(5000000); // autosave when every 5Mb
  fEventTree->SetCacheSize(10000000);  //set a 10 MBytes cache (useless when writing local files)
  fEventTree->Branch("EventBranch", "MCEvent", &fEvent, 512000, 99); 
  //--------------------------//

  // displacement_histo = new TH2F("displacement_histo", "; Lateral displacement/mm; Vertex energy (keV)", 
  // 				100, 0., 10., 20000, 0., 20.);
  // displacement_histo->SetDirectory(fFile);
}

//-----------------------------------------------------------------------
void DataManager::open(const char* name)
{
  fFileName = name;
  
  fFile = new TFile(fFileName, "RECREATE");

  G4String tree_comment="tree, Run # ";
  tree_comment+= fFileName;

  fEventTree = new TTree("event_tree", tree_comment);
  fEventTree->SetDirectory(fFile);
  fEventTree->SetAutoSave(5000000); // autosave when every 5Mb
  fEventTree->SetCacheSize(10000000);  //set a 10 MBytes cache (useless when writing local files)
  fEventTree->Branch("EventBranch", "MCEvent", &fEvent, 512000, 99); 
  //--------------------------//

  // displacement_histo = new TH2F("displacement_histo", "; Lateral displacement/mm; Vertex energy (keV)", 
  // 				100, 0., 10., 20000, 0., 20.);
  // displacement_histo->SetDirectory(fFile);
}


//---------------------------------------------------
void DataManager::close()
{
  if(fFile)
    {
      //ModuleInfo->Write("ModuleInfo");
      fFile->Write(); 
      fFile->Close();
    }
}

//-------------------------------------------
void DataManager::FillEvent()
{
  fEventTree->Fill();
}

//------------------------------------------------
void DataManager::AddPrimary(EventPrimary* primary)  
{
  fEvent->AddPrimary(primary);
}

//--------------------------------------------------
void DataManager::AddFiberHit(fiberHit *hitPtr)    
{ 
	fEvent->AddFiberHit(hitPtr); 
}



//------------------------------------------------
//void DataManager::AddOptTransmitted(opticalHit *hitPtr)  
//{
//  fEvent->AddOptTransmitted(hitPtr);
//}
