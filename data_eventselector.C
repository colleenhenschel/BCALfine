#if 1
#include <iostream>
#include <vector>
#include <math.h>

#include <iomanip>

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>

#include "TTree.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TFile.h"
#include "TPolyLine.h"
#include "TImage.h"

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
#include <TColor.h>

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

#endif

double fline(Double_t *x, Double_t *par)
{
  // Put the formula for your PDF's value here. Use the pre-computed
  // value of _norm to normalize the result.

  double qa=0,qb=0,qc=0,qx=0,qy=0;

  double xx = x[0];
  double pr1 = par[1];
  if (pr1==0) pr1 = 0.000000000001;
  return par[0]+xx/pr1;

}
//********************************************************************************
//********************************************************************************
double fline2(Double_t *x, Double_t *par)
{
  // Put the formula for your PDF's value here. Use the pre-computed
  // value of _norm to normalize the result.

  double qa=0,qb=0,qc=0,qx=0,qy=0;

  double xx = x[0];
  double pr1 = par[1];
  double pr2 = par[1];
  if (pr1==0) pr1 = 0.000000000001;
  if (pr2==0) pr2 = 0.000000000001;
  pr1 = 1./pr1;
  pr2 = 1./pr2;
  if ((pr1>0.&&xx<par[2])||(pr1<0.&&xx>par[2])) pr2=pr2*(1.+par[3]);
  double pc = par[0] + (pr1-pr2)*par[2];

  return pc+xx*pr2;

}

//********************************************************************************
// main function:
//********************************************************************************
void data_eventselector(int runnum=1002, int iplot=1) {

  int c;
  int imo,ise,ila;

  int cellenergy;

  Int_t  bd1   = 5;
  Int_t  bd2   = 6;

  //  Int_t  bd1   = 15;
  //  Int_t  bd2   = 16;

  TTree *outtree = new TTree("cosmics","cosmics");
  
  double mevpercm;
  double crossdist;
  double mev;
  double outz;
  double upstreame;
  double downstreame;
  int outlayer;
  int outsector;
  int outmodule;
  int cellid;
  int outrunnumber;
  
  outtree->Branch("mevpercm", &mevpercm, "mevpercm/D");
  outtree->Branch("crossdist", &crossdist, "crossdist/D");
  outtree->Branch("mev", &mev, "mev/D");
  outtree->Branch("z", &outz, "z/D"); 
  outtree->Branch("upE", &upstreame, "upE/D");
  outtree->Branch("downE", &downstreame, "downE/D");
  outtree->Branch("layer", &outlayer, "layer/I");
  outtree->Branch("sector", &outsector, "sector/I");
  outtree->Branch("module", &outmodule, "module/I");
  outtree->Branch("cellid", &cellid, "cellid/I");
  outtree->Branch("runnum", &outrunnumber, "runnum/I");

  double xx[768];
  double exx[768];
  double yy[768];
  double eyy[768];
  double yyy[768]; // y values for hits that z can be calculated for
  double eyyy[768];
  double zzz[768];
  double ezzz[768];
  int ipp, izz;

//********************************************************************************
 

//************************************************************************************

//************************************************************************************
  // BCAL

  Float_t pi = 3.1415926;

  TPolyLine *cell_cont[49];
  TPolyLine *cell[49][5][5];

  double xce[49][5][5],yce[49][5][5];
  double xR[5][49][5][5],yR[5][49][5][5];

  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

//************************************************************************************
  TCanvas *c1 = new TCanvas("BCAL Upstream","BCAL Upstream",700,700);
//  gPad->Range(-120,-120,120,120);
//  gPad->Range(-150,-150,150,150);

  gPad->Range(-100,-100,100,100);
  gStyle->SetOptStat(0);
//  gStyle->SetOptFit(0);
  gStyle->SetOptFit(1110);
  gStyle->SetStatX(0.68);
  gStyle->SetStatY(0.6);

  TCanvas *c2 = new TCanvas ("Y vs. Z","Y vs. Z",700,366);

  if (iplot != 0) {
    delete c1;
    c1 = new TCanvas("c1", " c1",12,29,499,699);
    c1->Divide(2,2);
    c1->Range(0.1875,-81.7688,3.3125,735.919);
    c1->SetBorderSize(2);
    c1->SetFrameFillColor(0);
    c1->SetGridx();
    c1->SetGridy();
    c1->SetTitle("Datafile: ");
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetOptStat(1110);
    gStyle->SetOptFit(1110);
    gStyle->SetPalette(1);
    gStyle->SetStatX(0.95);
    gStyle->SetStatY(0.95);

  }

//************************************************************************************
//Color Re-Definition

  double dcol = 1./47.00000001 ;
  double xcol = 1.;
  double xcol2 = 0.;
  for (int icol=1; icol<=49; icol++) {
    TColor *color=(TColor*)(gROOT->GetListOfColors()->At(icol));
    if (icol<49&&icol>1) {
//    color->SetRGB(1,xcol,xcol); //Red Gradient
//    color->SetRGB(xcol2,0,xcol); //Blue-to-Red Gradient
      color->SetRGB(1,xcol,0); //Yellow-to-Red Gradient
    } else {
      color->SetRGB(0,0,0);
    }
    xcol = xcol - dcol;
    xcol2 = xcol2 + dcol;
  }
  TColor *color50=(TColor*)(gROOT->GetListOfColors()->At(50));
  color50->SetRGB(0,0,1);
  TColor *color51=(TColor*)(gROOT->GetListOfColors()->At(51));
  color51->SetRGB(0,1,0);


  Float_t rad_min  = 65.462;
  Float_t rad_max  = rad_min + 22.294;

  Float_t  dphi  = 7.5*pi/180;

  Float_t x_cell[5];
  Float_t y_cell[5];

  for(Int_t module = 1; module <= 48; module++) {

    Float_t phi_min =  pi - (module - 1)*dphi  + dphi*0.5;
    Float_t phi_max =  phi_min - dphi;


    x_cell[0] = rad_min*cos(phi_min);
    y_cell[0] = rad_min*sin(phi_min);

    x_cell[1] = rad_max*cos(phi_min);
    y_cell[1] = rad_max*sin(phi_min);

    x_cell[3] = rad_min*cos(phi_max);
    y_cell[3] = rad_min*sin(phi_max);

    x_cell[2] = rad_max*cos(phi_max);
    y_cell[2] = rad_max*sin(phi_max);

    x_cell[4] = rad_min*cos(phi_min);
    y_cell[4] = rad_min*sin(phi_min);

    cell_cont[module] = new TPolyLine(5, x_cell, y_cell);
    cell_cont[module]->SetLineColor(49);
    cell_cont[module]->SetLineWidth(3);
  }


  Float_t dr = (rad_max - rad_min)/10.;
  Float_t dr_step[5];
//  dr_step[0] = 0;
//  dr_step[1] = dr;
//  dr_step[2] = 3*dr;
//  dr_step[3] = 6*dr;
//  dr_step[4] = 10*dr;

  dr_step[0] = 0.;
  dr_step[1] = 2.0574;
  dr_step[2] = 3.*2.0574;
  dr_step[3] = 6.*2.0574;
  dr_step[4] = 22.294;


  Float_t X_IN1, Y_IN1, X_IN2, Y_IN2, X_OUT1, Y_OUT1, X_OUT2, Y_OUT2;

  for(Int_t module = 1; module <= 48; module++) {

    //  Float_t phi_min =  pi - (module - 1)*dphi + dphi/2.;
    //  Float_t phi_max =  phi_min - dphi;

    Float_t phi_min =  pi - (module - 1)*dphi  + dphi*0.5;
    Float_t phi_max =  phi_min - dphi;

    for(Int_t layer = 1; layer <= 4; layer++) {

      X_IN1 = (rad_min + dr_step[layer-1])*cos(phi_min);
      Y_IN1 = (rad_min + dr_step[layer-1])*sin(phi_min);

      X_IN2 = (rad_min + dr_step[layer-1])*cos(phi_max);
      Y_IN2 = (rad_min + dr_step[layer-1])*sin(phi_max);

      X_OUT1 = (rad_min + dr_step[layer])*cos(phi_min);
      Y_OUT1 = (rad_min + dr_step[layer])*sin(phi_min);

      X_OUT2 = (rad_min + dr_step[layer])*cos(phi_max);
      Y_OUT2 = (rad_min + dr_step[layer])*sin(phi_max);


      Float_t L_IN = sqrt( pow(X_IN2 - X_IN1,2) + pow(Y_IN2 - Y_IN1,2));
      Float_t L_OUT = sqrt( pow(X_OUT2 - X_OUT1,2) + pow(Y_OUT2 - Y_OUT1,2));

      //    cout << L_IN << " " << L_OUT << endl;

      Float_t sin1 = (Y_IN2 - Y_IN1)/L_IN;
      Float_t cos1 = (X_IN2 - X_IN1)/L_IN;

      for(Int_t column = 1; column <= 4; column++) {

        x_cell[0] =  X_IN1 + (column - 1)*L_IN*cos1/4;
        y_cell[0] =  Y_IN1 + (column - 1)*L_IN*sin1/4;

        x_cell[1] =  X_OUT1 + (column - 1)*L_OUT*cos1/4;
        y_cell[1] =  Y_OUT1 + (column - 1)*L_OUT*sin1/4;

        x_cell[2] =  X_OUT1 + column*L_OUT*cos1/4;
        y_cell[2] =  Y_OUT1 + column*L_OUT*sin1/4;

        x_cell[3] =  X_IN1 + column*L_IN*cos1/4;
        y_cell[3] =  Y_IN1 + column*L_IN*sin1/4;

        x_cell[4] =  X_IN1 + (column - 1)*L_IN*cos1/4;
        y_cell[4] =  Y_IN1 + (column - 1)*L_IN*sin1/4;

        cell[module][layer][column] = new TPolyLine(5, x_cell, y_cell);
        cell[module][layer][column]->SetLineColor(49);
        cell[module][layer][column]->SetLineWidth(1);

        cell[module][layer][column]->SetFillStyle(1001);

        xce[module][layer][column]=(x_cell[0]+x_cell[1]+x_cell[2]+x_cell[3])/4.; //Coordinates of cell centers
        yce[module][layer][column]=(y_cell[0]+y_cell[1]+y_cell[2]+y_cell[3])/4.;

        xR[0][module][layer][column] = x_cell[0]; //Coordinates of cell corners
        yR[0][module][layer][column] = y_cell[0];
        xR[1][module][layer][column] = x_cell[1];
        yR[1][module][layer][column] = y_cell[1];
        xR[2][module][layer][column] = x_cell[2];
        yR[2][module][layer][column] = y_cell[2];
        xR[3][module][layer][column] = x_cell[3];
        yR[3][module][layer][column] = y_cell[3];
        xR[4][module][layer][column] = x_cell[4];
        yR[4][module][layer][column] = y_cell[4];
      }
    }

  }

//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
  int icurrentevent, imod,isec,ilay, xrunnum;
  TH1F *hncell = new TH1F("ncellspassed","ncellspassed", 800, 0, 800);
  
  int ncell2;
  
  double eneU[49][5][5], eneCore[49][5][5];
  double Z[49][5][5];
  int cellcolorU[49][5][5];
  double eneD[49][5][5];
  double cross[49][5][5], eneUnorm[49][5][5], eneDnorm[49][5][5], eneCorenorm[49][5][5];
  double dist[49][5][5];
  double enenorm[49][5][5];

  double errX[49][5][5];
  double errY[49][5][5];

  double ene[49][5][5];

  int out_flag[49][5][5];

  int cellcolorD[49][5][5];
  int nmodu; //Number of the fired modules
  int ncell; //Number of the fired cells

  double xmax,xmin,ymax,ymin,zmax,zmin,p1est,p0est;
  double xx1,xx2, yy1,yy2;

  double A1,A2,B1,B2,C1,C2, pp, xCR1,yCR1, xCR2,yCR2, xCR,yCR;
  int cross_flag;

  double errcell[5]= {0.5, 0.5, 1., 1.5, 2.};
  double errcell1[5]= {0.5, 0.25, 1., 1.5, 2.};
  double errcell0[5]= {2., 2., 2., 2., 2.};

  double distcut[5]= {2., 3., 3.5, 4.5, 5.5};

  double menergy[5]= {20., 20., 40., 40., 40.};

  int kill_flag, itop, ibot;

  int ipage = 0;

  double weightfun, enr;

  double x_int1,x_int2,y_int1,y_int2;
  bool good_int1,good_int2,kill_track;

  double XX1,XX2,XX3,XX4,YY1,YY2,YY3,YY4,D1,D2,D3,D4,DD0;
  double A1a,B1a,C1a, A1b,B1b,C1b, Xkink,Ykink;
  double CRX1,CRY1,CRX2,CRY2; //Begin(up) and end(down) region for MultScatt

  double R1 = rad_min;
  double R2 = rad_max;

//************************************************************************************
//Input Data File Init

  char *prnname = new char[90];

  sprintf(prnname, "/run/media/beattite/04c81856-39a1-4c05-837d-5f479ec8a84b/Cosmics/cosmics_2018points_%i.root",runnum);
  //sprintf(prnname, "cosmics_2018points_%i.root",runnum);
  TFile *fin = TFile::Open(prnname);
  TTree *t;
	fin->GetObject("Events", t);
	
	vector<int>* Dmod;
	vector<int>* Dsec;
	vector<int>* Dlay;
	vector<double>* Dup;
	vector<double>* Ddo;
	vector<double>* Dcore;
	vector<double>* Dz;
	t->SetBranchAddress("runnum", &xrunnum);
	t->SetBranchAddress("icurrentevent", &icurrentevent);
	t->SetBranchAddress("nmodu", &nmodu);
	t->SetBranchAddress("ncell", &ncell);
	t->SetBranchAddress("ncell2", &ncell2);
	t->SetBranchAddress("Dmod", &Dmod);
	t->SetBranchAddress("Dlay", &Dlay);
	t->SetBranchAddress("Dsec", &Dsec);
	t->SetBranchAddress("Dz", &Dz);
	t->SetBranchAddress("Denergy", &Dcore);	
	t->SetBranchAddress("Dup", &Dup);
	t->SetBranchAddress("Ddo", &Ddo); 

//************************************************************************************
//Output Vector File Init

  char *prnname2 = new char[90];
  sprintf(prnname2, "cosmics_2018_%i.vec",runnum);
  ofstream outf2(prnname2);

//************************************************************************************
//Postscript Drawing
    char *pdfname = new char[90];
    sprintf(pdfname, "cosmics_2018_event7h_%i.pdf",runnum);
    char *pdfname1 = new char[90];
    sprintf(pdfname1, "cosmics_2018_event7h_%i.pdf[",runnum);
    char *pdfname2 = new char[90];
    sprintf(pdfname2, "cosmics_2018_event7h_%i.pdf]",runnum);
  if(iplot == 0)
  {
    c1->Update();
    c1->Print(pdfname1);
  }
  
  int nevents = t->GetEntries();
  
  int ii=0;
//  while (true) {
  for (int ijk=0; ijk<nevents; ijk++) {
//  for (int ijk=0; ijk<100; ijk++) {
    t->GetEntry(ijk);
    if(ncell2 > 100) // skip events where more than 100 cells fire
      continue;
    if (1000*(ijk/1000)==ijk) cout<<"Event #"<<ijk<<endl;
    ipp = 0; //Number of Data points
    izz = 0;

  
    for(Int_t module = 1; module <= 48; module++) { // Init
      for(Int_t layer = 1; layer <= 4; layer++) {
        for(Int_t column = 1; column <= 4; column++) {
          eneU[module][layer][column] = 0.;
          eneCore[module][layer][column] = 0.;
          cellcolorU[module][layer][column] = 0;
          eneD[module][layer][column] = 0.;
          cellcolorD[module][layer][column] = 0;
          ene[module][layer][column] = 0.;

          cross[module][layer][column] = -999.;
          eneUnorm[module][layer][column] = -999.;
          eneDnorm[module][layer][column] = -999.;
          eneCorenorm[module][layer][column] = -999.;
          enenorm[module][layer][column] = -999.;
          Z[module][layer][column] = -999.;
          dist[module][layer][column] = 999.;
        }
      }
    }
    for (int jj=0; jj<Dup->size(); jj++) {
      Dup->at(jj) = Dup->at(jj)*1000; // converting energies from GeV to MeV
      Ddo->at(jj) = Ddo->at(jj)*1000;
      
      Dcore->at(jj) = Dcore->at(jj)*1000;
      //Dz->at(jj) = Dz->at(jj) - 147.0; // pointz 0 at center of target      
      
      imod = Dmod->at(jj);
      isec = Dsec->at(jj);
      ilay = Dlay->at(jj);
      if (Dcore->at(jj) > 0.) {
        eneCore[imod][ilay][isec] = Dcore->at(jj);
        
        Z[imod][ilay][isec] = Dz->at(jj) - 147.0;
      }
      if (Dup->at(jj)>0.) {
        eneU[imod][ilay][isec] = Dup->at(jj);
        ene[imod][ilay][isec] += Dup->at(jj);
      }

      if (Ddo->at(jj)>0.) {
        eneD[imod][ilay][isec] = Ddo->at(jj);
        ene[imod][ilay][isec] += Ddo->at(jj);
      }

      if (ene[imod][ilay][isec]>0.) {
        cellcolorU[imod][ilay][isec] = ene[imod][ilay][isec]*(48./1.2);
        if (cellcolorU[imod][ilay][isec]>48) cellcolorU[imod][ilay][isec]=48;
        xx[ipp] = xce[imod][ilay][isec];
        yy[ipp] = yce[imod][ilay][isec];
        exx[ipp] = errcell0[ilay]*menergy[ilay]*menergy[ilay]/(ene[imod][ilay][isec]*ene[imod][ilay][isec]);
        eyy[ipp] = errcell0[ilay]*menergy[ilay]*menergy[ilay]/(ene[imod][ilay][isec]*ene[imod][ilay][isec]);
        if(eneCore[imod][ilay][isec] > 0.) {
          zzz[izz] = Z[imod][ilay][isec];
          yyy[izz] = yce[imod][ilay][isec];
          ezzz[izz] = errcell0[ilay]*menergy[ilay]*menergy[ilay]/(ene[imod][ilay][isec]*ene[imod][ilay][isec]);
          eyyy[izz] = errcell0[ilay]*menergy[ilay]*menergy[ilay]/(ene[imod][ilay][isec]*ene[imod][ilay][isec]);
          izz++;
        }
        
        ipp++;
      }
    }
    //cout<<"Event #"<<ijk<<endl;
    ii++;

//************************************************************************************
//*************************************************************************************************************** Fit #1
//************************************************************************************


    xmin = 120;
    ymin = 120;
    xmax = -120;
    ymax = -120;

    for(int jj = 0; jj <ipp; jj++) {
      if(yy[jj] > ymax) {
        ymax = yy[jj];
        xmax = xx[jj];
      }
      if(yy[jj] < ymin) {
        ymin = yy[jj];
        xmin = xx[jj];
      }
    }

    p1est = (xmax-xmin)/(ymax-ymin);


    p0est = ymax - xmax/p1est;



//************************************************************************************

    TGraphErrors *gryz0 = new TGraphErrors(ipp,zzz,yyy,ezzz,eyyy);
    gryz0->GetXaxis()->SetLimits(-230.,230.);
    gryz0->SetMinimum(-120.);
    gryz0->SetMaximum(120.);
    gryz0->SetMarkerStyle(20);
    gryz0->SetMarkerSize(0.9);
    gryz0->SetMarkerColor(48);

    TGraphErrors *gr11 = new TGraphErrors(ipp,xx,yy,exx,eyy);
    gr11->GetXaxis()->SetLimits(-120.,120.);
    gr11->SetMinimum(-120.);
    gr11->SetMaximum(120.);
    gr11->SetMarkerStyle(20);
    gr11->SetMarkerSize(0.6);
    gr11->SetMarkerColor(1);

    TF1 *fu0 = new TF1("fline",fline,-120.,120.,2);
    fu0->SetLineColor(50);
    fu0->SetParameter(0,p0est);
    fu0->SetParameter(1,p1est);
//  fu0->SetParLimits(1,-0.59,0.44);
//  fu0->SetParLimits(1,-0.9,0.9);
    fu0->SetNpx(1000);
    gr11->Fit("fline","QMER");
//*************************************************************
//************************************************************************************* Fit #2
//*************************************************************
//Parameters of the track line
    A1 = 1./(fu0->GetParameter(1));
    B1 = -1.;
    C1 = (fu0->GetParameter(0));

//*************************************************************
//Clean Data Set

    ipp = 0;

    itop = 0;
    ibot = 0;
    kill_flag = 0;

    for(Int_t module = 1; module <= 48; module++) {
      for(Int_t layer = 1; layer <= 4; layer++) {
        for(Int_t column = 1; column <= 4; column++) {

          xx1=xce[module][layer][column]; //Distance from the fit line
          yy1=yce[module][layer][column];
          dist[module][layer][column]=abs(A1*xx1+B1*yy1+C1)/sqrt(A1*A1+B1*B1);

          if (dist[module][layer][column]<2.*distcut[layer]&&ene[module][layer][column]>0.) {
            xx[ipp] = xce[module][layer][column];
            yy[ipp] = yce[module][layer][column];
            exx[ipp] = errcell[layer]*menergy[layer]/(ene[module][layer][column]);
            eyy[ipp] = errcell[layer]*menergy[layer]/(ene[module][layer][column]);
            if (yy[ipp]>0) itop++;
            if (yy[ipp]<=0) ibot++;
            ipp++;
          }


        }
      }
    }

//************************************************************************************
    p0est = fu0->GetParameter(0);
    p1est = fu0->GetParameter(1);

//************************************************************************************


    TGraphErrors *gr12 = new TGraphErrors(ipp,xx,yy,exx,eyy);
    gr12->GetXaxis()->SetLimits(-120.,120.);
    gr12->SetMinimum(-120.);
    gr12->SetMaximum(120.);
    gr12->SetMarkerStyle(20);
    gr12->SetMarkerSize(0.6);
    gr12->SetMarkerColor(1);

    TF1 *fu2 = new TF1("fline2",fline,-120.,120.,2);
    if (kill_flag==0) {
      fu2->SetLineColor(50);
    } else {
      fu2->SetLineColor(51);
    }
    fu2->SetParameter(0,p0est);
    fu2->SetParameter(1,p1est);
//  fu2->SetParLimits(1,-0.9,0.9);
    fu2->SetNpx(1000);
    gr12->Fit("fline2","QMER");

//*************************************************************
//************************************************************************************** Fit #3
//*************************************************************
//Parameters of the track line
    A1 = 1./(fu2->GetParameter(1));
    B1 = -1.;
    C1 = (fu2->GetParameter(0));

//************************
//Checking Crossing with Outer Radius for Through-Tracks
    x_int1 = (-2*A1*C1 + sqrt((2*A1*C1)*(2*A1*C1) - 4*(C1*C1 + A1*A1*C1*C1 - A1*A1*87*87 - 87*87))) / (2*(1 + A1*A1));
    x_int2 = (-2*A1*C1 - sqrt((2*A1*C1)*(2*A1*C1) - 4*(C1*C1 + A1*A1*C1*C1 - A1*A1*87*87 - 87*87))) / (2*(1 + A1*A1));
    y_int1 = A1*x_int1 + C1;
    y_int2 = A1*x_int2 + C1;

//************************
//Clean Data Set Again

    ipp = 0;

    itop = 0;
    ibot = 0;
    kill_flag = 0;
    kill_track = false;
    good_int1 = false;
    good_int2 = false;

    for(Int_t module = 1; module <= 48; module++) {
      for(Int_t layer = 1; layer <= 4; layer++) {
        for(Int_t column = 1; column <= 4; column++) {

          xx1=xce[module][layer][column]; //Distance from the fit line
          yy1=yce[module][layer][column];
          dist[module][layer][column]=abs(A1*xx1+B1*yy1+C1)/sqrt(A1*A1+B1*B1);

          if (!good_int1&&ene[module][layer][column]>0.) if (sqrt((xx1 - x_int1)*(xx1 - x_int1) + (yy1 - y_int1)*(yy1 - y_int1)) < 6.0) good_int1 = true;
          if (!good_int2&&ene[module][layer][column]>0.) if (sqrt((xx1 - x_int2)*(xx1 - x_int2) + (yy1 - y_int2)*(yy1 - y_int2)) < 6.0) good_int2 = true;

          if (dist[module][layer][column]<distcut[layer]&&ene[module][layer][column]>0.) {
            xx[ipp] = xce[module][layer][column];
            yy[ipp] = yce[module][layer][column];
            exx[ipp] = errcell1[layer]*menergy[layer]/(ene[module][layer][column]);
            eyy[ipp] = errcell1[layer]*menergy[layer]/(ene[module][layer][column]);
            if (yy[ipp]>0) itop++;
            if (yy[ipp]<=0) ibot++;
            ipp++;
          }


        }
      }
    }

    if (!good_int1 || !good_int2) kill_track = true; // The track does not go completely through the BCAL

//************************************************************************************
    p0est = fu2->GetParameter(0);
    p1est = fu2->GetParameter(1);

//************************************************************************************
    TGraphErrors *gr14 = new TGraphErrors(ipp,xx,yy,exx,eyy);
    gr14->GetXaxis()->SetLimits(-120.,120.);
    gr14->SetMinimum(-120.);
    gr14->SetMaximum(120.);
    gr14->SetMarkerStyle(20);
    gr14->SetMarkerSize(0.6);
    gr14->SetMarkerColor(1);

    TF1 *fu3 = new TF1("fline3",fline,-120.,120.,2);
    if (!kill_track&&kill_flag==0) {
      fu3->SetLineColor(50);
    } else {
      fu3->SetLineColor(51);
    }
    fu3->SetParameter(0,p0est);
    fu3->SetParameter(1,p1est);
    fu3->SetParLimits(0,0.8*p0est,1.1*p0est);
    fu3->SetParLimits(1,0.9*p1est,1.1*p1est);
    fu3->SetNpx(1000);
    gr14->Fit("fline3","QMER");


//*************************************************************
//Parameters of the track line
    A1 = 1./(fu3->GetParameter(1));
    B1 = -1.;
    C1 = (fu3->GetParameter(0));
//*************************************************************
//Search for the crossing with cells
    for(Int_t module = 1; module <= 48; module++) {
      for(Int_t layer = 1; layer <= 4; layer++) {
        for(Int_t column = 1; column <= 4; column++) {
          cross_flag = 0;
          for (int jj=1; jj<5; jj++) {
            xx1 = xR[jj-1][module][layer][column];
            yy1 = yR[jj-1][module][layer][column];
            xx2 = xR[jj][module][layer][column];
            yy2 = yR[jj][module][layer][column];
            if (xx2!=xx1) {
              pp = (yy2-yy1)/(xx2-xx1);
            } else {
              pp = 1000000000000.;
            }
            A2 = pp;
            B2 = -1.;
            C2 = yy1-pp*xx1;


            xCR = (-C1*B2+C2*B1)/(A1*B2-A2*B1);
            yCR = (-A1*C2+A2*C1)/(A1*B2-A2*B1);

            xmin = min(xx1,xx2);
            xmax = max(xx1,xx2);
            ymin = min(yy1,yy2);
            ymax = max(yy1,yy2);
            if (cross_flag==1&&xCR>=xmin&&xCR<=xmax&&yCR>=ymin&&yCR<=ymax) {
              xCR2 = xCR;
              yCR2 = yCR;
              cross_flag++;
            }  else if (cross_flag==0&&xCR>=xmin&&xCR<=xmax&&yCR>=ymin&&yCR<=ymax) {
              xCR1 = xCR;
              yCR1 = yCR;
              cross_flag++;
            }

          }

          if (cross_flag==2) {
            cross[module][layer][column]=sqrt((xCR1-xCR2)*(xCR1-xCR2)+(yCR1-yCR2)*(yCR1-yCR2));
            if (cross[module][layer][column]>0.2) {
              eneUnorm[module][layer][column] = eneU[module][layer][column]/cross[module][layer][column];
              eneDnorm[module][layer][column] = eneD[module][layer][column]/cross[module][layer][column];
              if(eneCore[module][layer][column] > 0)
                eneCorenorm[module][layer][column] = eneCore[module][layer][column]/cross[module][layer][column];
              enenorm[module][layer][column]=eneUnorm[module][layer][column]+eneDnorm[module][layer][column];
            }
          }

        }
      }
    }

//*************************************************************
//************************************************************************************** Fit #4
//*************************************************************
//Clean Data Set Once Again

    ipp = 0;
    itop = 0;
    ibot = 0;
    kill_flag = 0;

    for(Int_t module = 1; module <= 48; module++) {
      for(Int_t layer = 1; layer <= 4; layer++) {
        for(Int_t column = 1; column <= 4; column++) {

          xx1=xce[module][layer][column]; //Distance from the fit line
          yy1=yce[module][layer][column];
          dist[module][layer][column]=abs(A1*xx1+B1*yy1+C1)/sqrt(A1*A1+B1*B1);

          if (dist[module][layer][column]<distcut[layer]&&ene[module][layer][column]>0.) {
            xx[ipp] = xce[module][layer][column];
            yy[ipp] = yce[module][layer][column];

            enr=(eneCorenorm[module][layer][column]/200.)-1.;
            weightfun = 1.;
            if (cross[module][layer][column]>0.2) {
              weightfun = 0.5 + 0.5*exp(-enr/0.5);
            }
            exx[ipp] = 2.*weightfun*errcell1[layer]*menergy[layer]/(ene[module][layer][column]);
            eyy[ipp] = 2.*weightfun*errcell1[layer]*menergy[layer]/(ene[module][layer][column]);
            errX[module][layer][column] = exx[ipp];
            errY[module][layer][column] = eyy[ipp];

            if (yy[ipp]>0) itop++;
            if (yy[ipp]<=0) ibot++;
            ipp++;
          }


        }
      }
    }


//************************************************************************************
    p0est = fu3->GetParameter(0);
    p1est = fu3->GetParameter(1);
//if (p1est>0.45) p1est = 0.43;
//if (p1est<-0.60) p1est = -0.58;

//************************************************************************************
    TGraphErrors *gr15 = new TGraphErrors(ipp,xx,yy,exx,eyy);
    gr15->GetXaxis()->SetLimits(-120.,120.);
    gr15->SetMinimum(-120.);
    gr15->SetMaximum(120.);
    gr15->SetMarkerStyle(20);
    gr15->SetMarkerSize(0.6);
    gr15->SetMarkerColor(1);

    TF1 *fu4 = new TF1("fline4",fline,-120.,120.,2);
    if(!kill_track && kill_flag == 0) {
      fu4->SetLineColor(50);
    } else {
      fu4->SetLineColor(51);
    }
    fu4->SetParameter(0,p0est);
    fu4->SetParameter(1,p1est);
    fu4->SetParLimits(0,0.8*p0est,1.2*p0est);
    //fu4->SetParLimits(1,-0.60,0.45);
    fu4->SetNpx(1000);
    gr15->Fit("fline4","QMER");

//*************************************************************
//*************************************************************
//************************************************************************************** Fit #5
//*************************************************************
//Parameters of the track line
    A1 = 1./(fu4->GetParameter(1));
    B1 = -1.;
    C1 = (fu4->GetParameter(0));

//Crossing with the BCAL inner radius
    D1 = A1*A1+B1*B1;
    D2 = 2.*A1*C1;
    D3 = C1*C1 - B1*B1*R1*R1;

    DD0 = D2*D2 - 4.*D1*D3;
    XX1 = -999.;
    XX2 = -999.;
    if (DD0>0) {
      XX1 = (-D2+sqrt(DD0))/(2.*D1) ;
      YY1 = (-C1-A1*XX1)/B1;
      XX2 = (-D2-sqrt(DD0))/(2.*D1) ;
      YY2 = (-C1-A1*XX2)/B1;
      if (YY1>YY2) {
        if (XX1>XX2) {
          CRX1 = XX1;
          CRX2 = XX1+(XX1-XX2)/20.;
        } else {
          CRX2 = XX1;
          CRX1 = XX1+(XX1-XX2)/20.;
        }
      } else {
        if (XX2>XX1) {
          CRX1 = XX2;
          CRX2 = XX2+(XX2-XX1)/20.;
        } else {
          CRX2 = XX2;
          CRX1 = XX2+(XX2-XX1)/20.;
        }
      }
    } else {
//Crossing with the BCAL outer radius
      D1 = A1*A1+B1*B1;
      D2 = 2.*A1*C1;
      D3 = C1*C1 - B1*B1*R2*R2;

      DD0 = D2*D2 - 4.*D1*D3;
      XX3 = -999.;
      XX4 = -999.;
      if (DD0>0) {
        XX3 = (-D2+sqrt(DD0))/(2.*D1) ;
        YY3 = (-C1-A1*XX3)/B1;
        XX4 = (-D2-sqrt(DD0))/(2.*D1) ;
        YY4 = (-C1-A1*XX4)/B1;
        if (XX3>XX4) {
          CRX1 = XX4 + 0.4*(XX3-XX4);
          CRX2 = XX4 + 0.6*(XX3-XX4);
        } else {
          CRX1 = XX3 + 0.4*(XX4-XX3);
          CRX2 = XX3 + 0.6*(XX4-XX3);
        }
      } else {
        CRX1 = -999. ;
        CRX2 = -998. ;
      }

    }

//************************
//Clean Data Set Again

    ipp = 0;
    izz = 0;

    itop = 0;
    ibot = 0;
    kill_flag = 0;

    for(Int_t module = 1; module <= 48; module++) {
      for(Int_t layer = 1; layer <= 4; layer++) {
        for(Int_t column = 1; column <= 4; column++) {

          xx1=xce[module][layer][column]; //Distance from the fit line
          yy1=yce[module][layer][column];
          dist[module][layer][column]=abs(A1*xx1+B1*yy1+C1)/sqrt(A1*A1+B1*B1);

          if (dist[module][layer][column]<distcut[layer]&&ene[module][layer][column]>0.) {
            xx[ipp] = xce[module][layer][column];
            yy[ipp] = yce[module][layer][column];


            enr=(eneCorenorm[module][layer][column]/200.)-1.;
            weightfun = 1.;
            if (cross[module][layer][column]>0.2) {
              weightfun = 0.5 + 0.5*exp(-enr/0.5);
            }

            exx[ipp] = 2.*weightfun*errcell1[layer]*menergy[layer]/(ene[module][layer][column]);
            eyy[ipp] = 2.*weightfun*errcell1[layer]*menergy[layer]/(ene[module][layer][column]);
            errX[module][layer][column] = exx[ipp];
            errY[module][layer][column] = eyy[ipp];

            if(eneCore[module][layer][column] > 0.) {
              zzz[izz] = Z[module][layer][column]; 
              yyy[izz] = yce[module][layer][column];
              ezzz[izz] = 2.*weightfun*errcell1[layer]*menergy[layer]/(ene[module][layer][column]);
              eyyy[izz] = 2.*weightfun*errcell1[layer]*menergy[layer]/(ene[module][layer][column]);
              izz++;
            }
            if (yy[ipp]>0) itop++;
            if (yy[ipp]<=0) ibot++;
            ipp++;
          }


        }
      }
    }

// if (itop<3||ibot<3) kill_flag = 100;

//************************
//************************
// Y vs. Z plot
    TGraphErrors *gryz = new TGraphErrors(izz,zzz,yyy,ezzz,eyyy);
    gryz->GetXaxis()->SetLimits(-230., 230.);
    gryz->SetMinimum(-120.);
    gryz->SetMaximum(120.);
    gryz->SetMarkerStyle(20);
    gryz->SetMarkerSize(0.6);
    gryz->SetMarkerColor(1);

    zmin = 230;
    ymin = 120;
    zmax = -230;
    ymax = -120;

    for(int jj =0; jj <izz; jj++)
    {
      if(yyy[jj] > ymax)
      {
        ymax = yyy[jj];
        zmax = zzz[jj];
      }
      if(yyy[jj] < ymin)
      {
        ymin = yyy[jj];
        zmin = zzz[jj];
      }
    }

    if(zmax == zmin || ymax == ymin) p1est = 0.00000001;
    else p1est = (zmax-zmin)/(ymax -ymin);
    p0est = ymax - zmax/p1est;

    TF1 *fuyz = new TF1("flineyz", fline, -230.0, 230.0, 2);
    fuyz->SetParameter(0, p0est);
    fuyz->SetParameter(1, p1est);
    fuyz->SetNpx(400);
    gryz->Fit("flineyz", "QW");
    gryz->Fit("flineyz", "QMER");

    double A1z, B1z, C1z, distz, zCR1, zCR2;
    int kill_flagz = 0;

    A1z = 1./(fuyz->GetParameter(1));
    B1z = -1.;
    C1z = (fuyz->GetParameter(0));

    double chi_fuyz = (fuyz->GetChisquare())/(fuyz->GetNDF());
    if(chi_fuyz > 30.) kill_flagz += 1000;

    if (fabs(A1z) < 0.25) kill_flagz += 100; // y vs. z angle cut, not needed if using 3D crossing

    for(int i = 0; i < izz; i++)
    {
      distz = abs(A1z*zzz[i]+B1z*yyy[i]+C1z)/sqrt(A1z*A1z + B1z*B1z);
      if(distz > 6.0) kill_flagz++;
    }

//************************************************************************************
    p0est = fu4->GetParameter(0);
    p1est = fu4->GetParameter(1);
//if (p1est>0.45) p1est = 0.43;
//if (p1est<-0.60) p1est = -0.58;

//************************************************************************************
    TGraphErrors *gr16 = new TGraphErrors(ipp,xx,yy,exx,eyy);
    gr16->GetXaxis()->SetLimits(-120.,120.);
    gr16->SetMinimum(-120.);
    gr16->SetMaximum(120.);
    gr16->SetMarkerStyle(20);
    gr16->SetMarkerSize(0.6);
    gr16->SetMarkerColor(1);

    TF1 *fu5 = new TF1("fline3",fline2,-120.,120.,4);
    /*  if (kill_flag==0) {
      fu5->SetLineColor(50);
      } else {
      fu5->SetLineColor(51);
      } */
    fu5->SetParameter(0,p0est);
    fu5->SetParameter(1,p1est);
    fu5->SetParLimits(0,0.8*p0est,1.1*p0est);
//  fu5->SetParLimits(1,-0.65,0.50);
    fu5->SetParameter(2,(CRX1+CRX2)*0.5);
    fu5->SetParLimits(2,CRX1,CRX2);
    fu5->SetParameter(3,0.);
//  fu5->SetParLimits(3,-0.9,0.9);
    fu5->SetNpx(10000);
    gr16->Fit("fline3","QMER");

//*************************************************************
//Parameters of the track lines
    A1a = 1./(fu5->GetParameter(1));
    B1a = -1.;
    C1a = (fu5->GetParameter(0));

    Xkink = (fu5->GetParameter(2));
    Ykink = (-C1a-A1a*Xkink)/B1a;

    A1b = (1./(fu5->GetParameter(1)))*(1.+(fu5->GetParameter(3)));
    B1b = -1.;
    C1b = C1a + (A1a-A1b)*Xkink;

//*************************************************************
//Search for the crossing with cells
    for(Int_t module = 1; module <= 48; module++) {
      for(Int_t layer = 1; layer <= 4; layer++) {
        for(Int_t column = 1; column <= 4; column++) {
          cross_flag = 0;
          for (int jj=1; jj<5; jj++) {
            xx1 = xR[jj-1][module][layer][column];
            yy1 = yR[jj-1][module][layer][column];
            xx2 = xR[jj][module][layer][column];
            yy2 = yR[jj][module][layer][column];

            if ((yy1+yy2)*0.5>Ykink) {
              A1 = A1a;
              B1 = B1a;
              C1 = C1a;
            } else {
              A1 = A1b;
              B1 = B1b;
              C1 = C1b;
            }

            if (xx2!=xx1) {
              pp = (yy2-yy1)/(xx2-xx1);
            } else {
              pp = 1000000000000.;
            }
            A2 = pp;
            B2 = -1.;
            C2 = yy1-pp*xx1;


            xCR = (-C1*B2+C2*B1)/(A1*B2-A2*B1);
            yCR = (-A1*C2+A2*C1)/(A1*B2-A2*B1);

            xmin = min(xx1,xx2);
            xmax = max(xx1,xx2);
            ymin = min(yy1,yy2);
            ymax = max(yy1,yy2);
            if (cross_flag==1&&xCR>=xmin&&xCR<=xmax&&yCR>=ymin&&yCR<=ymax) {
              xCR2 = xCR;
              yCR2 = yCR;
              cross_flag++;
            }  else if (cross_flag==0&&xCR>=xmin&&xCR<=xmax&&yCR>=ymin&&yCR<=ymax) {
              xCR1 = xCR;
              yCR1 = yCR;
              cross_flag++;
            }

          }

          if (cross_flag==2) {
            zCR1 = (yCR1 - C1z) / A1z;
            zCR2 = (yCR2 - C1z) / A1z;
            //cross[module][layer][column]=sqrt((xCR1-xCR2)*(xCR1-xCR2)+(yCR1-yCR2)*(yCR1-yCR2));
            cross[module][layer][column] = sqrt((xCR1-xCR2)*(xCR1-xCR2)+(yCR1-yCR2)*(yCR1-yCR2) + (zCR1 - zCR2)*(zCR1 - zCR2)); // 3D crossing
            if (cross[module][layer][column]>0.2) {
              eneUnorm[module][layer][column] = eneU[module][layer][column]/cross[module][layer][column];
              eneDnorm[module][layer][column] = eneD[module][layer][column]/cross[module][layer][column];
              if(eneCore[module][layer][column] > 0.)
                eneCorenorm[module][layer][column] = eneCore[module][layer][column]/cross[module][layer][column];
              enenorm[module][layer][column]=eneUnorm[module][layer][column]+eneDnorm[module][layer][column];
            }
          }
          if (cross[module][layer][column]<0.&&ene[module][layer][column]>0.) kill_flag++;

        }
      }
    }

    double chi_fu5 = (fu5->GetChisquare())/(fu5->GetNDF());
    if(chi_fu5>30.) kill_flag += 1000;



    if (kill_track || kill_flag > 4||kill_flagz > 4) {
      fuyz->SetLineColor(51);
      fu5->SetLineColor(51);
    } else {
      fuyz->SetLineColor(50);
      fu5->SetLineColor(50);
    }



//************************
//************************
//Draw the Picture
    if (iplot==0) c1->Clear();
    for(Int_t module = 1; module <= 48; module++) {
      for(Int_t layer = 1; layer <= 4; layer++) {
        for(Int_t column = 1; column <= 4; column++) {

          cellenergy = cellcolorU[module][layer][column];
          cell[module][layer][column]->SetFillColor(cellenergy);

          if (iplot==0) {
            cell[module][layer][column]->Draw("f");
            cell[module][layer][column]->Draw("Same");
          }

        }
      }
      if (iplot==0) cell_cont[module]->Draw();
    }
    if (iplot==0) {

      gStyle->SetStatX(0.68);
      gStyle->SetStatY(0.6);

      gr16->Draw("P");
      fu5->Draw("Same");
      if (ncell2>0)  {
        c1->Update();
        if(iplot == 0)c1->Print(pdfname);
        c1->Clear();
        ipage++;
        if (ipage==67) {
          //cout<<"A1="<<A1<<"   B1="<<B1<<"   C1="<<C1<<endl;
        }
      }
    }
    if(iplot == 0) {
      c2->cd();

      gStyle->SetStatX(0.9);
      gStyle->SetStatY(0.9);

      gryz0->Draw("AP");
      gryz->Draw("P");
      fuyz->Draw("Same");
      if(ncell2>0) {
        c2->Update();
        if(iplot == 0)c2->Print(pdfname);
        c2->Clear();
        if (ipage == 4) {
        }

      }
    }
//*************************************************************
    if (iplot!=0) {
      for(Int_t module = 1; module <= 48; module++) {
        for(Int_t layer = 1; layer <= 4; layer++) {
          for(Int_t column = 1; column <= 4; column++) {
            imo = module-1;
            ise = column-1;
            ila = layer -1;
            if (!kill_track&&kill_flag<4&&kill_flagz<4&&cross[module][layer][column]>0.4 && Z[module][layer][column] > -195 && Z[module][layer][column] < 195) {
              if (eneCore[module][layer][column]>0.)
              {
                mevpercm = eneCorenorm[module][layer][column];
                mev = eneCore[module][layer][column];
                crossdist = cross[module][layer][column];
                outz = Z[module][layer][column];
                upstreame = eneU[module][layer][column];
                downstreame = eneD[module][layer][column];
                outlayer = layer;
                outsector = column;
                outmodule = module;
                cellid = (module -1)*48 + (layer - 1)*4 + column;
                outrunnumber = runnum;
                outtree->Fill();
              }
            }
          }
        }
      }
    }

    if(gryz0) gryz0->Delete();
    if(gr11) gr11->Delete();
    if(gr12) gr12->Delete();
    if(gr14) gr14->Delete();
    if(gr15) gr15->Delete();
    if(gr16) gr16->Delete();
    if(gryz) gryz->Delete();
    
    if(fu0)  fu0->Delete();
    if(fu2) fu2->Delete();
    if(fu3) fu3->Delete();
    if(fu4) fu4->Delete();
    if(fu5) fu5->Delete();
    if(fuyz) fuyz->Delete();
//*************************************************************
//*************************************************************
//*************************************************************
  }//End of event loop
//*************************************************************
//*************************************************************
//*************************************************************
  

  
  
//*************************************************************
  c1->Update();
  if(iplot == 0) c1->Print(pdfname2);
  
  TFile *myoutfile = new TFile(Form("run%i_cosmics_spectra.root",runnum), "RECREATE");
  outtree->SetDirectory(myoutfile);
  outtree->Print();
  myoutfile->Write();
  myoutfile->Close();
//************************************
}

