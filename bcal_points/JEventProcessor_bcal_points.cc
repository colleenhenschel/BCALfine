// $Id$
//
//    File: JEventProcessor_bcal_points.cc
// Created: Mon May  5 15:20:49 EDT 2014
// Creator: davidl (on Darwin harriet.jlab.org 13.1.0 i386)
//

#include <iostream>
using namespace std;

#include "JEventProcessor_bcal_points.h"
using namespace jana;

#include <BCAL/DBCALHit.h>
#include <BCAL/DBCALPoint.h>

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_bcal_points());
}
} // "C"


//------------------
// JEventProcessor_bcal_points (Constructor)
//------------------
JEventProcessor_bcal_points::JEventProcessor_bcal_points()
{

}

//------------------
// ~JEventProcessor_bcal_points (Destructor)
//------------------
JEventProcessor_bcal_points::~JEventProcessor_bcal_points()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_bcal_points::init(void)
{
	BCALpoint = new TTree("BCALpoint","DBCALPoint objects for each channel and event");
	BCALpoint->Branch("channelnum",&channelnum,"channelnum/i");
	BCALpoint->Branch("eventnum",&eventnum,"eventnum/i");
	BCALpoint->Branch("point_E",&point_E,"point_E/d");
	BCALpoint->Branch("point_z",&point_z,"point_z/d");
	BCALpoint->Branch("hitU_E",&hitU_E,"hitU_E/d");
	BCALpoint->Branch("hitD_E",&hitD_E,"hitD_E/d");

	BCALpoint->Branch("module",&module,"module/i");
	BCALpoint->Branch("layer",&layer,"layer/i");
	BCALpoint->Branch("sector",&sector,"sector/i");

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_bcal_points::brun(JEventLoop *eventLoop, int32_t runnumber)
{
	// This is called whenever the run number changes
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_bcal_points::evnt(JEventLoop *loop, uint64_t eventnumber)
{
	channelnum = 0;
	eventnum = eventnumber;
	vector<int> cellids;

	// Get the DBCALPoint objects
	vector<const DBCALPoint*> bcalpoints;
	loop->Get(bcalpoints);
	// Get the DBCALHit objects
	vector<const DBCALHit*> bcalhits;
	loop->Get(bcalhits);

	japp->RootWriteLock();

	// Loop over DBCALPoint objects
	for(unsigned int i=0; i< bcalpoints.size(); i++){
		try {
			const DBCALPoint *bcalpoint = bcalpoints[i];

			module = bcalpoint->module();
			layer = bcalpoint->layer();
			sector = bcalpoint->sector();
			channelnum++;
			point_E = bcalpoint->E();
			point_z = bcalpoint->z();
			hitU_E = bcalpoint->E_US();
			hitD_E = bcalpoint->E_DS();

			cellids.push_back((module-1)*16 + (layer-1)*4 + (sector-1)); // cellid from 0 to 767

			// Fill tree
			BCALpoint->Fill();

		} catch (...) {}
	}
//cout << endl;
//for (int k = 0; k < cellids.size(); k++) cout << cellids[k] << "  ";
	// Loop over DBCALHit objects
	for(unsigned int i=0; i< bcalhits.size(); i++){
		try {
			const DBCALHit *bcalhit = bcalhits[i];

			module = bcalhit->module;
			layer = bcalhit->layer;
			sector = bcalhit->sector;
//cout << endl << (module-1)*16 + (layer-1)*4 + (sector-1);
			if (find(cellids.begin(), cellids.end(), (module-1)*16 + (layer-1)*4 + (sector-1)) == cellids.end()) {
//cout << " ! not in cellids vector";
				channelnum++;
				point_E = -999;
				point_z = -999;
				hitU_E = -999;
				hitD_E = -999;
				if (bcalhit->end == 0) hitU_E = bcalhit->E;
				else if (bcalhit->end == 1) hitD_E = bcalhit->E;

				// Fill tree
				BCALpoint->Fill();
			}

		} catch (...) {}
	}

	japp->RootUnLock();

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_bcal_points::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_bcal_points::fini(void)
{
	// Called before program exit after event processing is finished.
	return NOERROR;
}

