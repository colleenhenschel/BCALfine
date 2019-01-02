// $Id$
//
//    File: JEventProcessor_bcal_points.h
// Created: Mon May  5 15:20:49 EDT 2014
// Creator: davidl (on Darwin harriet.jlab.org 13.1.0 i386)
//

#ifndef _JEventProcessor_bcal_points_
#define _JEventProcessor_bcal_points_

#include <stdint.h>

#include <JANA/JEventProcessor.h>

#include <TTree.h>


class JEventProcessor_bcal_points:public jana::JEventProcessor{
	public:
		JEventProcessor_bcal_points();
		~JEventProcessor_bcal_points();
		const char* className(void){return "JEventProcessor_bcal_points";}

		TTree *BCALpoint;
		uint32_t channelnum;         ///< Arbitrary global channel number (sorted by crate, slot, channel)
		uint32_t eventnum;	         ///< Event number
		double point_E;
		double point_z;
		double hitU_E;
		double hitD_E;

		int module;
		int layer;
		int sector;

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_bcal_points_

