/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2013-2015 by
  o/.|.\o    E   nvironment    | LeMonADE Principal Developers (see AUTHORS)
    ooo                        |
----------------------------------------------------------------------------------

This file is part of LeMonADE.

LeMonADE is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

LeMonADE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with LeMonADE.  If not, see <http://www.gnu.org/licenses/>.

--------------------------------------------------------------------------------*/

#ifndef LEMONADE_ANALYZER_INITIAL_POSITIONS_H
#define LEMONADE_ANALYZER_INITIAL_POSITIONS_H

#include <string>

#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/analyzer/AbstractAnalyzer.h>
#include <LeMonADE/utility/ResultFormattingTools.h>
#include <LeMonADE/utility/MonomerGroup.h>

/*************************************************************************
 * definition of AnalyzerInitialPositions class
 * ***********************************************************************/

/**
 * @file
 *
 * @class AnalyzerInitialPositions
 *
 * @brief Analyzer for evaluating the squared radius of gyration IP
 *
 * @tparam IngredientsType Ingredients class storing all system information( e.g. monomers, bonds, etc).
 *
 * @details This analyzer calculates the squared radius of gyration of either the
 * IPplete system (if no subgroup of monomers was specified), or the average of
 * the IP of a set of subgroups (if subgroups were specified using the provided
 * function void setMonomerGroups(std::vector<MonomerGroup<molecules_type> >) )
 * The analyzer calculates IP in every timestep and saves the time series
 * of this data to disk in the format
 * mcs <IP_x> <IP_y> <IP_z>  <IP_tot>
 * The default output filename is IPTimeSeries.dat and it can be changed by argument
 * to the constructor or by using the setter function provided.
 * If a more sophisticated grouping of monomers into groups is required, one can
 * also write a new analyzer, inheriting from AnalyzerInitialPositions, and
 * overwriting the initialize function.
 */
template < class IngredientsType > class AnalyzerInitialPositions : public AbstractAnalyzer
{

private:
	//! typedef for the underlying container holding the monomers
	typedef typename IngredientsType::molecules_type molecules_type;
	//! reference to the IPplete system
	const IngredientsType& ingredients;
	//! IP is calculated for the groups in this vector
	std::vector<MonomerGroup<molecules_type> > groups;
	//! timeseries of the IP. Components: [0]-> IP_x, [1]->IP_y, [2]->IP_z, [3]:IP_tot
	std::vector< std::vector<double> > IPTimeSeries;
	//! vector of mcs times for writing the time series
	std::vector<double> MCSTimes;
	//! max length of internal buffer for each coordinate, before saving to disk
	uint32_t bufferSize;
	//! name of output files are outputFilePrefix_averages.dat and outputFilePrefix_timeseries.dat
	std::string outputFile;
	//! flag used in dumping time series output
	bool isFirstFileDump;
	//! save the current values in IPTimeSeriesX, etc., to disk
	uint32_t site;
	//! location of middle of macrodomain
	void dumpTimeSeries();
	//! calculate the Rg squared of the monomer group
	VectorDouble3 calculateIPComponents(const MonomerGroup<molecules_type>& group) const;
protected:
	//! Set the groups to be analyzed. This function is meant to be used in initialize() of derived classes.
	void setMonomerGroups(std::vector<MonomerGroup<molecules_type> > groupVector){groups=groupVector;}
public:

	//! constructor
	AnalyzerInitialPositions(const IngredientsType& ing,
			  std::string filename="IPTimeSeries.dat",uint32_t terC=0);

	//! destructor. does nothing
	virtual ~AnalyzerInitialPositions(){}
	//! Initializes data structures. Called by TaskManager::initialize()
	virtual void initialize();
	//! Calculates the IP for the current timestep. Called by TaskManager::execute()
	virtual bool execute();
	//! Writes the final results to file
	virtual void cleanup();
	//! Set the number of values, after which the time series is saved to disk
	void setBufferSize(uint32_t size){bufferSize=size;}
	//! Change the output file name
	void setOutputFile(std::string filename){outputFile=filename;isFirstFileDump=true;}

};

/*************************************************************************
 * implementation of members
 * ***********************************************************************/

/**
 * @param ing reference to the object holding all information of the system
 * @param filename output file name. defaults to "IPTimeSeries.dat".
 * */
template<class IngredientsType>
AnalyzerInitialPositions<IngredientsType>::AnalyzerInitialPositions(
	const IngredientsType& ing,
	std::string filename, uint32_t terC)
:ingredients(ing)
,site(terC)
,IPTimeSeries(ingredients.getMolecules().size()*3,std::vector<double>(0))
,bufferSize(100)
,outputFile(filename)
,isFirstFileDump(true)
{
}


/**
 * @details fills all monomers into the groups vector as a single group,
 * if the group size is zero. Normally this should always be the case,
 * unless the groups were already set explicitly by using the setter function
 * setMonomerGroups
 * */
template< class IngredientsType >
void AnalyzerInitialPositions<IngredientsType>::initialize()
{
	//! Select group
//    int min = (site-ingredients.getMolecules().size()/8)%ingredients.getMolecules().size();

//
//    if(groups.size()==0){
//        groups.push_back(MonomerGroup<molecules_type>(ingredients.getMolecules()));
//        for(size_t n=site;n<site+ingredients.getMolecules().size()/4;n++) {
//
//            groups[0].push_back(n%ingredients.getMolecules().size());
//
//        }
//    }
//! all monomers
    if(groups.size()==0){
        groups.push_back(MonomerGroup<molecules_type>(ingredients.getMolecules()));
        for(size_t n=site;n<ingredients.getMolecules().size();n++) {

            groups[0].push_back(n);

        }
    }
}

/**
 * @details Calculates the current IP, saves it in the
 * time series, and saves the time series to disk in regular intervals.
 * */
template< class IngredientsType >
bool AnalyzerInitialPositions<IngredientsType>::execute()
{
	VectorDouble3 IPComponents(0.0,0.0,0.0);

    for(size_t n=0;n<groups.size();n++)
    {
        //this vector will contain (IP_x, IP_y, IP_z), i.e. the squared Components!
        IPComponents+=calculateIPComponents(groups[n])/double(groups.size());
    }
    for(int i=0;i<ingredients.getMolecules().size())
	IPTimeSeries[0].push_back(ingredients.getMolecules()[i].getX());
	IPTimeSeries[1].push_back(IPComponents.getY());
	IPTimeSeries[2].push_back(IPComponents.getZ());
    IPTimeSeries[3].push_back(site);
	MCSTimes.push_back(ingredients.getMolecules().getAge());
	//save to disk in regular intervals
	if(MCSTimes.size()>=bufferSize)
		dumpTimeSeries();

	return true;
}


template<class IngredientsType>
void AnalyzerInitialPositions<IngredientsType>::cleanup()
{
	std::cout<<"AnalyzerInitialPositions::cleanup()...";
	//write the remaining data from the time series
	dumpTimeSeries();
	std::cout<<"done\n";
}

/**
 * @details Saves the current content of IPTimeSeries to the file
 * IPTimeSeries.dat. The output format is:
 * mcs IP_x IP_y IP_z IP_tot
 * */
template<class IngredientsType>
void AnalyzerInitialPositions<IngredientsType>::dumpTimeSeries()
{
	//fist make a single vector<vector<double> > for writing the results
	std::vector<std::vector<double> > resultsTimeseries=IPTimeSeries;
	resultsTimeseries.insert(resultsTimeseries.begin(),MCSTimes);

	//if it is written for the first time, include IPment in the output file
	if(isFirstFileDump){
		std::stringstream IPmentTimeSeries;
		IPmentTimeSeries<<"Created by AnalyzerInitialPositions\n";
		IPmentTimeSeries<<"file contains time series of average Rg_squared (IP) over all analyzed groups\n";
		IPmentTimeSeries<<"format: mcs\t IPX\t IPY\t IPZ\t IPTotal\n";

		ResultFormattingTools::writeResultFile(
			outputFile,
			ingredients,
			resultsTimeseries,
			IPmentTimeSeries.str()
		);

		isFirstFileDump=false;
	}
	//otherwise just append the new data
	else{
		ResultFormattingTools::appendToResultFile(outputFile,
							  resultsTimeseries);
	}
	//set all time series vectors back to zero size
	MCSTimes.resize(0);
	IPTimeSeries.resize(0);
	IPTimeSeries.resize(4,std::vector<double>(0));
}

/**
 * @details calculates the three Components IP_x, IP_y,IP_z and returns
 * them in a vector.
 * @return VectorDouble3 containing the Components IP_x, IP_y,IP_z, or (0.0,0.0,0.0) if group is empty)
 * @param group the monomer group of which the IP is calculated
 * */
template<class IngredientsType>
VectorDouble3 AnalyzerInitialPositions<IngredientsType>::calculateIPComponents(
	const MonomerGroup<molecules_type>& group) const
{
	//if group is empty, return zero vector
	if(group.size()==0){
		return VectorDouble3(0.0,0.0,0.0);
	}


	VectorDouble3 sum_sqr; VectorDouble3 IP_sum;
	//first calculate the center of mass
	for ( size_t n = 0; n < group.size(); ++n)
	{
		IP_sum.setX( IP_sum.getX() + group[n].getX() );
		IP_sum.setY( IP_sum.getY() + group[n].getY() );
		IP_sum.setZ( IP_sum.getZ() + group[n].getZ() );
	}
	double inv_N = 1.0 / double ( group.size() );

	VectorDouble3 IP (double ( IP_sum.getX() ) * inv_N,
			   double ( IP_sum.getY() ) * inv_N,
			   double ( IP_sum.getZ() ) * inv_N);
	
	return IP;

}
#endif


