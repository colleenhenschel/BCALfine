#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TF2.h>
#include <TProfile.h>
#include <TSpectrum.h>
#include <TStyle.h>
#include <TChain.h>
#include <TPostScript.h>
#include <TLatex.h>
#include <TMath.h>
#include <TPaveLabel.h>
#include <TMinuit.h>

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <utility>
#include <vector>
#include <map>

#include <string>
#include <cstring>
#include <sstream>
#include <inttypes.h>

#include "TPaveStats.h"
#include "TGraphPainter.h"
#include "TString.h"
#include "TCollection.h"

#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TKey.h"
#include "TDatime.h"
#include "TAxis.h"
#include "TLine.h"
#include "TTree.h"
#include "TBranch.h"
#include "TROOT.h"
#include "TSystem.h"

using namespace std;

//********************************************************************************

void cosmics_2018_datafilter_points(int runnum=4101,int nfiles=23) {


  float Ethre = 0.010 ; //Energy threshold of 10keV (for v2 it was 1eV)

  int imo = 12; //dummy
  int ico = 1;  //dummy 
  
//Histograms Definition

  
  TH1F *ehis1 = new TH1F("ehis1"," ",100,0.,20.);
   ehis1->SetTitle("");
   ehis1->GetXaxis()->SetTitle("Number of Fired Modules");
   ehis1->GetYaxis()->SetTitle("Number of Events");
  
  TH1F *ehis2 = new TH1F("ehis2"," ",120,0.,60.);
   ehis2->SetTitle("");
   ehis2->GetXaxis()->SetTitle("Number of Fired Cells");
   ehis2->GetYaxis()->SetTitle("Number of Events");

//************************************************************************************
//************************************************************************************
  
//************************************************************************************
//Output Vector File Init
  char *prnname = new char[90];
  sprintf(prnname, "cosmics_2018points_%i.dat",runnum);
  ofstream prnfile(prnname);

//************************************************************************************
  
  int nmodu; //Number of the fired modules
  int lastmodu; //The last fired module
  int ncell; //Number of the fired cells
  int lastcell; //The last fired cell
  
  int cellcode;
  
  int ncell2;
  int Dmod[1536],Dsec[1536],Dlay[1536];
  double Dup[1536],Ddo[1536],Denergy[1536],Dz[1536];      
      
  int icurrentevent;
//************************************************************************************


//Root file open
  char *rootfilename = new char[90];
 
  TChain *bcaltree = new TChain("BCALpoint");
//  sprintf(rootfilename, "/work/halld/beattite/pointcosmics/hists/hd_root.root");
//  bcaltree->Add(rootfilename);
  for  (int iru=0; iru<nfiles; iru++){
//  for  (int iru=0; iru<1; iru++){
    if (iru < 10) sprintf(rootfilename, "/work/halld/beattite/pointcosmics/hists/Run0%d/hd_root_Run0%d_00%d.root",runnum,runnum,iru);
    else if (iru < 100) sprintf(rootfilename, "/work/halld/beattite/pointcosmics/hists/Run0%d/hd_root_Run0%d_0%d.root",runnum,runnum,iru);
    else sprintf(rootfilename, "/work/halld/beattite/pointcosmics/hists/Run0%d/hd_root_Run0%d_%d.root",runnum,runnum,iru);
    bcaltree->Add(rootfilename); 
  }
  cout<<"Run #"<<runnum<<" is connected: "<<rootfilename<<endl;

   int eventnum;
   bcaltree->SetBranchAddress("eventnum",&eventnum); bcaltree->SetBranchStatus("eventnum",1);
   double point_E;
   bcaltree->SetBranchAddress("point_E",&point_E); bcaltree->SetBranchStatus("point_E",1);
   double point_z;
   bcaltree->SetBranchAddress("point_z",&point_z); bcaltree->SetBranchStatus("point_z",1);
   double hitU_E;
   bcaltree->SetBranchAddress("hitU_E",&hitU_E); bcaltree->SetBranchStatus("hitU_E",1);
   double hitD_E;
   bcaltree->SetBranchAddress("hitD_E",&hitD_E); bcaltree->SetBranchStatus("hitD_E",1);
   int module;
   bcaltree->SetBranchAddress("module",&module); bcaltree->SetBranchStatus("module",1);
   int layer;
   bcaltree->SetBranchAddress("layer",&layer); bcaltree->SetBranchStatus("layer",1);
   int sector;
   bcaltree->SetBranchAddress("sector",&sector); bcaltree->SetBranchStatus("sector",1);

  unsigned long nevent = bcaltree->GetEntries();
  cout<<"Rootfiles contain "<<nevent<<" records."<<endl;
//************************************************************************************
  icurrentevent = -999;
  
  nmodu = 0;
  lastmodu = -999;
  
  ncell = 0;
  lastcell = -999;
  

  
//Data Reading

  for(Int_t i = 0; i < nevent; i++){          //start of records loop

//    trig_time = 0.;
    if (1000000*(i/1000000)==i) cout<<i<<" records processed"<<endl; 

    bcaltree->GetEntry(i);

    cellcode = 100*module + 10*sector +layer;

    if (eventnum != icurrentevent) {   // End of event: Making of decision

      if (ncell2>0) {
        ehis1->Fill(1.*nmodu);
        ehis2->Fill(1.*ncell);

        prnfile<<runnum<<" "<<icurrentevent<<" "<<endl;
        prnfile<<"  "<<nmodu<<" "<<ncell<<" "<<ncell2<<endl;
        for (int jj=0; jj<ncell2; jj++){
          prnfile<<"         "<<Dmod[jj]<<" "<<Dsec[jj]<<" "<<Dlay[jj]<<" "<<Dup[jj]<<" "<<Ddo[jj]<<" "<<Denergy[jj]<<" "<<Dz[jj]<<endl;
        }
      } 	 

      icurrentevent = eventnum ;      // Next event init
      
      for (int jj=0; jj<1536; jj++){
        Dmod[jj]=-999;
        Dsec[jj]=-999;
        Dlay[jj]=-999;
        Dup[jj]=-999.;
        Ddo[jj]=-999.;
        Denergy[jj]=-999.;
        Dz[jj]=-999.;
        ncell2 = 0;
      }
     

    } //End of making decision/ End of event


    if (lastmodu != module && module>0&&module<49) {
      lastmodu = module;
      nmodu += 1;
    }

    if (lastcell != cellcode && module>0&&module<49) {
      lastcell = cellcode;
      ncell += 1;
    }

    Dmod[ncell2] = module;
    Dsec[ncell2] = sector;
    Dlay[ncell2] = layer;
    Dup[ncell2] = hitU_E;
    Ddo[ncell2] = hitD_E;
    Denergy[ncell2] = point_E;
    Dz[ncell2] = point_z;
    ncell2++;    
  }
  
  
  cout<<"We have "<<icurrentevent<<" cosmic events."<<endl;
  
  prnfile.close();

//************************************************************************************

}
