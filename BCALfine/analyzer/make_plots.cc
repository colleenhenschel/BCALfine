

#include <iostream>
#include <cstdlib>

#include "MCEvent.hh"

#include "TFile.h"
#include "TString.h"
#include "TH1F.h"
#include "TVector3.h"
#include "TTree.h"
#include "TBenchmark.h"
#include "TChain.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
using namespace std;

#define PI 3.14159265358979

void Usage(void);
void ParseCommandLineArguments(int narg, char* argv[]);

TString outname = "";
TString inname = "";

int main(int narg, char* argv[])
{
  TBenchmark* bench = new TBenchmark();
  bench->Reset();
  bench->Start("benchmark");

  //parse command line arguements;
  ParseCommandLineArguments(narg, argv);

  //open the output root file
  TFile *outFile = new TFile(outname, "RECREATE");
  //TH1F* edep_total = new TH1F("total_edep", "; Edep/eV", 1000, -0.5, 999.5);
  //edep_total->SetDirectory(outFile);

	TTree *t = new TTree("Events", "Events");
	
	
	Int_t n_event = 0;
	Double_t mcEventTotalEdep = 0.;
	Double_t muinitE = 0.;
	Double_t mux = 0.;
	Double_t muz = 0.;

	t->Branch("event_num", &n_event, "Event Number/I");
	t->Branch("total_energy", &mcEventTotalEdep, "Total Energy/D");

	t->SetDirectory(outFile);

  MCEvent      *event    = 0;
  fiberHit   	 *fibHit   = 0;
  EventPrimary *primary  = 0;

  TClonesArray *primariesArray = 0;
  TClonesArray *fiberHitArray = 0;

  TChain * treeChain = new TChain("event_tree");
  treeChain->Add(inname);
  treeChain->SetBranchAddress("EventBranch", &event);  

	vector<int> module;
	vector<int> sector;
	vector<int> layer;
	vector<int> cellid;
	vector<double> cell_energy;
	vector<double> average_z;
	
	vector<int> fiber_cellid;
	vector<double> fiber_energy;
	vector<double> fiber_z;

	float E_core[4][4][48] = {0};
	float E_up[4][4][48] = {0};
	float E_down[4][4][48] = {0};
	float Z[4][4][48] = {0};

	t->Branch("module", &module);
	t->Branch("layer", &layer);
	t->Branch("sector", &sector);
	t->Branch("cellid", &cellid);
	t->Branch("cell_tot_energy", &cell_energy);
	t->Branch("average_z", &average_z);
	t->Branch("mu_initialenergy", &muinitE, "Muon Initial Energy/D");
	t->Branch("mu_initialx", &mux, "Muon Start X/D");
	t->Branch("mu_initialz", &muz, "Muon Start Z/D");
	t->Branch("E_core", &E_core, "E_core[4][4][48]/F");	
	t->Branch("E_up", &E_up, "E_up[4][4][48]/F");
	t->Branch("E_down", &E_down, "E_down[4][4][48]/F");
	t->Branch("Z", &Z, "Z[4][4][48]/F");


  for(Int_t i=0; i<treeChain->GetEntries(); i++)
  {
    treeChain->GetEntry(i);
      
    n_event = i;

    mcEventTotalEdep = 0.;
		module.clear();
		sector.clear();
		layer.clear();
		cellid.clear();
		cell_energy.clear();
		average_z.clear();

		fiber_cellid.clear();
		fiber_energy.clear();
		fiber_z.clear();
		for(int mm = 0; mm < 48; mm++)
		{
			for(int ss = 0; ss < 4; ss++)
			{
				E_core[ss][0][mm] = 0;
				E_up[ss][0][mm] = 0;
				E_down[ss][0][mm] = 0;
				Z[ss][0][mm] = 0;

				E_core[ss][1][mm] = 0;
				E_up[ss][1][mm] = 0;
				E_down[ss][1][mm] = 0;
				Z[ss][1][mm] = 0;

				E_core[ss][2][mm] = 0;
				E_up[ss][2][mm] = 0;
				E_down[ss][2][mm] = 0;
				Z[ss][2][mm] = 0;

				E_core[ss][3][mm] = 0;
				E_up[ss][3][mm] = 0;
				E_down[ss][3][mm] = 0;
				Z[ss][3][mm] = 0;
			}
		}
	
		primariesArray = event->GetPrimaries();
		for(Int_t nn = 0; nn<primariesArray->GetEntries(); nn++)
		{
			primary = (EventPrimary*) primariesArray->At(nn);
			muinitE = primary->GetTotalEnergy();
			mux = primary->GetVertexPosition().X();
			muz = primary->GetVertexPosition().Z();
		}


    fiberHitArray = event->GetFiberHits();
    for(Int_t kk=0; kk<fiberHitArray->GetEntries(); kk++)
    {
      fibHit = (fiberHit*) fiberHitArray->At(kk);
			
			Double_t x = fibHit->GetVertexPos().X();
			Double_t y = fibHit->GetVertexPos().Y();
			Double_t z = fibHit->GetVertexPos().Z();
			if(y > 100)
			{
				//cout << "Error: Initial y > 100, skipping event" << endl;
				continue;
			}
			Double_t phi;
			if(x > 0 && y >= 0)
				phi = atan(y/x) * 180 / PI; // angle phi in degrees
			else if(x == 0 && y > 0)
				phi = 90;
			else if(x < 0 && y >= 0)
				phi = 180 + atan(y/x) * 180 / PI;
			else if(x < 0 && y < 0)
				phi = 180 + atan(y/x) * 180 / PI;
			else if(x == 0 && y < 0)
				phi = 270;
			else
				phi = 360 + atan(y/x) * 180 / PI;

			Double_t r = sqrt(x*x + y*y);
			Int_t temp_mod;			
			if(phi >= 0 && phi < 356.25)
				temp_mod = trunc((phi + 3.75)/7.5) + 1;
			else if(phi >= 356.25 && phi <= 360)
				temp_mod = 1;
			else
			{
				cout << "Problem with module: " << endl;
				cout << "Phi = " << phi << endl;
				cout << "R = " << r << endl;
			}

			Double_t delta_phi = phi - (temp_mod - 1)*7.5;
			Int_t temp_sect = -1;
			if(phi >= 358.125 && phi <= 360)
			{
				temp_sect = 2;
			}
			else if(phi >= 356.25 && phi < 358.125)
			{
				temp_sect = 1;
			}
			if(delta_phi > 1.875 && delta_phi <= 3.75)
				temp_sect = 4;
			else if(delta_phi > 0 && delta_phi <= 1.875)
				temp_sect = 3;
			else if(delta_phi > -1.875 && delta_phi <= 0)
				temp_sect = 2;
			else if(delta_phi >= -3.75 && delta_phi <= -1.875)
				temp_sect = 1;

			Int_t temp_layer = -1;
			if( r >= 65.462 && r <= 67.5194) // bottom 3 layer light guides 0.81 inches tall, starting r comes from midradius - 3.175/2 (offset of module to baseplate center) - module_height/2
				temp_layer = 1;
			else if(r > 67.5194 && r <= 71.6342)
				temp_layer = 2;
			else if(r > 71.6342 && r <= 77.8064)
				temp_layer = 3;
			else if(r > 77.8064 && r <= 87.756) //top light guides 0.97 inches tall
				temp_layer = 4;
			else
			{
				cout << "Unphysical layer, ignore hit" << endl;
				cout << "Y = " << y << "   X = " << x << "    R = " << r << endl;
				continue;
			}
			fiber_energy.push_back(fibHit->GetEnergyDep());
			if(temp_mod < 1 || temp_mod > 48)
			{
				cout << "Calculated module problem: " << temp_mod << "    r = " << r << "   phi = " << phi << endl;
			}
			if(temp_layer < 1 || temp_layer > 4)
			{
				cout << "Calculated layer problem: " << temp_layer << "    r = " << r << "   phi = " << phi << endl;
			}
			if(temp_sect < 1 || temp_sect > 4)
			{
				cout << "Calculated sector problem: " << temp_sect << "    r = " << r << "   phi = " << phi << endl;
			}
			fiber_cellid.push_back((temp_mod - 1)*16 + (temp_layer - 1)*4 + temp_sect);
			fiber_z.push_back(z);

  	  mcEventTotalEdep += fibHit->GetEnergyDep();
    }

		std::vector<int>::iterator it;
		
		for(int jj = 0; jj < fiber_cellid.size(); jj++)
		{
			if(fiber_cellid.at(jj) != 0)
			{
				Double_t summed_E = 0;
				Double_t sum_z = 0;
				Double_t fib_count = 0;
				Int_t cur_cell = fiber_cellid.at(jj);
				for(int kk = 0; kk < fiber_cellid.size(); kk++)
				{
					if(fiber_cellid.at(kk) == cur_cell)
					{
						summed_E += fiber_energy.at(kk);
						sum_z += fiber_z.at(kk);
						fib_count++;
						fiber_cellid.at(kk) = 0;
					}
				}
				cellid.push_back(cur_cell);
				cell_energy.push_back(summed_E);
				Int_t mod = ceil(cur_cell/16.0);
				if(mod == 0 || mod == 49)
				{
					cout << "Weird module number, cell id = " << cur_cell << endl;
				}
				module.push_back(mod);
				Int_t lay = ceil((cur_cell - (mod - 1)*16.0)/4.0);
				layer.push_back(lay);
				Int_t sec = cur_cell - (mod-1)*16.0 - (lay-1)*4.0;
				sector.push_back(sec);
				average_z.push_back(sum_z/fib_count);
				E_core[sec - 1][lay - 1][mod - 1] += summed_E;
				E_up[sec - 1][lay - 1][mod - 1] += summed_E*exp(-(sum_z/fib_count + 390/2)/525.);
				E_down[sec - 1][lay - 1][mod - 1] += summed_E*exp(-(390/2 - sum_z/fib_count)/525.);
				Z[sec - 1][lay - 1][mod - 1] = sum_z/fib_count;
			}
		}

		t->Fill();    
		
		if(i%10==0)
  	{
          cout << "\r >>> "<< i <<" Events processed" << endl;
          //cout.flush();
  	}
  }

  cout<<endl;
	t->Print();
  outFile->Write();
  outFile->Close();	
  
  bench->Show("benchmark");
  delete bench;
  return 0;
}

//----------------------------
// ParseCommandLineArguments
//----------------------------
void ParseCommandLineArguments(int narg, char* argv[])
{
  if(narg<1)Usage();

  for(int i=1; i<narg; i++)
    {
      std::string arg = argv[i];
      if(arg=="-h" || arg=="--help")
        {
          Usage();
        }
      else if(arg=="-inF")
        {
          if(i==narg-1)
            { std::cout<<"-inF requires one argument!"<<std::endl;  Usage(); }

	  inname = argv[i+1];
        }
      else if(arg=="-outF")
        {
          if(i==narg-1)
            { std::cout<<"-outF requires one argument!"<<std::endl;  Usage(); }

	  outname = argv[i+1];
        }

    }

  if(outname==""||inname=="")
    {
      Usage();
    }

  std::cout<<" Input file name = "<<inname<<std::endl;
  std::cout<<" Output file name = "<<outname<<std::endl;

}


//----------------------------
// Usage
//----------------------------
void Usage(void)
{
  std::cout<<std::endl;
  std::cout<<"Usage:"<<std::endl;
  std::cout<<"      exec_name: -inF mc_data.root -outF output_name.root"<<std::endl;
  std::cout<<std::endl;
  std::cout<<"  options:"<<std::endl;
  std::cout<<"    -h                      print this help message"<<std::endl;
  std::cout<<"    -inF mc_data            set the name of MC data file (required)"<<std::endl;
  std::cout<<"    -outF mc_data           set the name of output file (required)"<<std::endl;

  exit(-2);
}
