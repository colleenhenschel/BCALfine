#include<iostream>
#include<fstream>

#include"TFile.h"
#include"TF1.h"
#include"TCanvas.h"
#include"TTree.h"
#include"TH1F.h"
#include"TStyle.h"

using namespace std;

void data_fits()
{
  TFile *file = TFile::Open("currentdata.root");
  TCanvas *c1 = new TCanvas("c1","canvas", 1300, 1000);
  c1->Divide(4,4);
  gStyle->SetOptFit(111);
  ofstream outfile;
  outfile.open("data_fitvaluesnew.txt");
  
  if (file == 0) 
  {
    // if we cannot open the file, print an error message and return immediatly
    printf("Error: cannot open root file\n");
    return;
  }

  TH1F* hist;

  TH1F* hists[48][4][4];

  TTree* t = (TTree*)file->Get("cosmics");

  TF1 *f = new TF1("f", "landau(0) + pol3(3)", 2, 18.0);
  TF1 *land = new TF1("land", "landau", 2, 18.0);
  land->SetLineColor(3);
  TF1 *bg = new TF1("bg", "pol3", 2, 18.0);
  bg->SetLineColor(8);
  double p0, p1, p2, p3, p4, p5, p6;

  p0 = 1000;
  p1 = 7.5;
  p2 = 1.1;
  p3 = 0.0;
  p4 = 0.0;
  p5 = 0.0;
  p6 = 0.0;

  for(int i = 0; i<48; i++)
  {
    for(int j = 0; j<4; j++)
    {
      for(int k = 0; k < 4; k++)
      {
        c1->cd((k*4+j+1));
        if(k == 3)
          f->Delete();
          f = new TF1("f", "landau(0) + pol3(3)", 4, 18.0);
          t->Draw("mevpercm>>h(100, 0.0, 20.0)",Form("module == %i && sector == %i && layer == %i", i+1, j+1, k+1));
        else
          t->Draw("mevpercm>>h(100, 0.0, 20.0)",Form("module == %i && sector == %i && layer == %i", i+1, j+1, k+1));
        


        hist = (TH1F*)gDirectory->Get("h");  

        hists[i][j][k] = (TH1F*)hist->Clone(Form("m%i_s%i_l%i", i+1, j+1, k+1));

        hists[i][j][k]->SetTitle(Form("Data point energy: Module %i, Sector %i, Layer %i; MeV/cm",i+1, j+1, k+1));

        f->SetParameter(0, p0);
        f->SetParameter(1, p1);
        f->SetParameter(2, p2);
        f->SetParameter(3, p3);
        f->SetParameter(4, p4);
        f->SetParameter(5, p5);
        f->SetParameter(6, p6);
    

        hists[i][j][k]->Fit("f", "R");
        bg->SetParameters(f->GetParameter(3), f->GetParameter(4), f->GetParameter(5), f->GetParameter(6));
        land->SetParameters(f->GetParameter(0), f->GetParameter(1), f->GetParameter(2));
        hists[i][j][k]->SetLineColor(4);

        c1->cd((k*4+j+1));
        hists[i][j][k]->Draw();
        bg->DrawCopy("Same");
        //land->Draw("Same");

        outfile << setw(7) << i+1 << setw(7) << j+1  << setw(7) << k+1 << setprecision(3) << setw(7) << f->GetMaximumX(2.0, 18.0) << setprecision(3) << setw(7)  << f->GetParameter(1) << endl;
      }
    } 
    c1->Print(Form("data_mod%i.pdf",i+1));   
  }
    
   outfile.close();
}
