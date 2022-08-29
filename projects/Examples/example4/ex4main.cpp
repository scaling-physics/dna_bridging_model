/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2013-2015 by
  o/.|.\o    E   nvironment    | LeMonADE Principal Developers (Hauke Rabbel)
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

/* *********************************************************************
* This example demonstrates the basic use of updaters and analyzers.
* A simple simulation program is put together and can be executed.
* In the example, a short polymer chain is placed in the system and
* simulated in a periodic box for 1E7 Monte Carlo steps. The
* configurations are saved to a file so they can be watched later, and
* the radius of gyration is analyzed.
* *********************************************************************/
#define LOG(x) std::cout<<x<<std::endl<<std::endl;

#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureBox.h>
#include <LeMonADE/feature/FeatureBondset.h>
#include <LeMonADE/feature/FeatureLattice.h>
#include <LeMonADE/utility/RandomNumberGenerators.h>
#include <LeMonADE/utility/TaskManager.h>
#include <LeMonADE/updater/UpdaterAddLinearChains.h>
#include <LeMonADE/updater/UpdaterSimpleSimulator.h>
#include <LeMonADE/updater/UpdaterSimpleSimulatorBridges.h>
#include <LeMonADE/updater/UpdaterSimpleSimulatorBridgesAll.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>
#include <LeMonADE/analyzer/AnalyzerPositionOfOri.h>
#include <LeMonADE/analyzer/AnalyzerPositionOfTer.h>
#include <LeMonADE/analyzer/AnalyzerRadiusOfGyration.h>
#include <LeMonADE/analyzer/AnalyzerCentreOfMass.h>
#include <LeMonADE/feature/FeatureNNInteractionSc.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include "LeMonADE/updater/Bridges.h"
#include <random>
#include "Flags.h"

int main(int argv, char **argc) {
    /* ****************************************************************
      * as described in the previous examples, we quickly seed the
      * random number generator, and set up a system in a box with the
      * standard BFM bond set. We add 10 monomers and connect them to a
      * linear polymer.
      * ***************************************************************/

    //first set up the random number generator
    RandomNumberGenerators randomNumbers;
    randomNumbers.seedAll();

    //now set up the system
    //here, we use one additional feature called FeatureMoleculesIO
    //strictly speaking this is not a feature in the sense that it
    //changes the simulation conditions. What it provides is the
    //basic functionalities for writing BFM files. So, in most cases
    //you are going to use this feature in simulations.
    typedef LOKI_TYPELIST_4(FeatureMoleculesIO, FeatureBox, FeatureExcludedVolumeSc<FeatureLattice<>>,
                            FeatureBoltzmann) Features;
    typedef ConfigureSystem<VectorInt3, Features, 2> Config;
    typedef Ingredients<Config> MyIngredients;
    MyIngredients mySystem;
    double nsteps=100;

    mySystem.setBoxX(88);
    mySystem.setBoxY(22);
    mySystem.setBoxZ(22);

    mySystem.setPeriodicX(false);
    mySystem.setPeriodicY(false);
    mySystem.setPeriodicZ(false);

    mySystem.modifyBondset().addBFMclassicBondset();
    mySystem.synchronize(mySystem);

    Bridges bridgeS;

    uint32_t simulation_number = atoi(argc[1]);

    bridgeS.r1 = atof(argc[2]);
    bridgeS.p1 = atof(argc[3])/atof(argc[2]);

    bridgeS.r2 = atof(argc[4]);
    bridgeS.p2 = atof(argc[5])/atof(argc[4]);

    bridgeS.theta = bridgeS.p2/bridgeS.p1;


    LOG(bridgeS.r1)
    LOG(bridgeS.p1)

//    uint32_t simulation_number = 1;
//    bridgeS.p1 = 0.1;
//    bridgeS.r1 = 10000;
    uint32_t polymerLength = 440;



    //! Initialize polymer

    mySystem.modifyMolecules().resize(polymerLength);

    std::ifstream file;
    int filenumber = randomNumbers.r250_rand32() % 400+1000;
    file.open("/scratch/Srikanth/ter_dynamics/2/IP_" + std::to_string(filenumber) + ".dat");

    //! ter and ori positions

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dis(0, polymerLength);

    uint32_t randomSite = dis(gen);



    int num;
    std::vector<int> vec;

    while (file >> num)
        vec.push_back(num);

    std::cout<<randomSite<<'\n';

//    std::cout<<vec[randomSite*3]<<'\t'<<vec[randomSite*3]<<'\t'<<vec[randomSite*3]<<'\n';
//
//    std::cout<<vec[0]<<'\t'<<vec[1]<<'\t'<<vec[2]<<'\n';
//
    std::rotate(vec.begin(),vec.begin()+3*(randomSite),vec.end()); //circular shift input vector

    std::cout<<vec[0]<<'\t'<<vec[1]<<'\t'<<vec[2];

    mySystem.modifyMolecules()[0].setX(vec[0]);
    mySystem.modifyMolecules()[0].setY(vec[1]);
    mySystem.modifyMolecules()[0].setZ(vec[2]);


    int j = 3;
    for (uint32_t i = 1; i < polymerLength; i++) {
        mySystem.modifyMolecules()[i].setX(vec[j]);
        mySystem.modifyMolecules()[i].setY(vec[j + 1]);
        mySystem.modifyMolecules()[i].setZ(vec[j + 2]);
        mySystem.modifyMolecules().connect(i - 1, i);
        j = j + 3;
    }

    mySystem.modifyMolecules().connect(0, polymerLength - 1);
    mySystem.synchronize(mySystem);


    int oriC = 0; //Position of Ori start
    uint32_t terC = polymerLength / 2;  // Position of ter start


//    mySystem.modifyMolecules()[0].setX(0);
//    mySystem.modifyMolecules()[0].setY(0);
//    mySystem.modifyMolecules()[0].setZ(0);
//
//
//    //intialize polymer
//    int n = 0, counter = 1;
//    for (uint32_t i = 1; i < mySystem.getBoxX() / 2; i++) {
//
//        mySystem.modifyMolecules()[counter].setX(i * 2);
//        mySystem.modifyMolecules()[counter].setY(n);
//        mySystem.modifyMolecules()[counter].setZ(0);
//
//        mySystem.modifyMolecules().connect(counter - 1, counter);
//        counter++;
//    }
//
//    while (n < 15) {
//        n += 2;
//        for (uint32_t i = 1; i < mySystem.getBoxX() / 2; i++) {
//            mySystem.modifyMolecules()[counter].setX((mySystem.getBoxX() / 2 - i) * 2);
//            mySystem.modifyMolecules()[counter].setY(n);
//            mySystem.modifyMolecules()[counter].setZ(0);
//
//            mySystem.modifyMolecules().connect(counter - 1, counter);
//            counter++;
//
//
//        }
//
//
//        n += 2;
//        for (uint32_t i = 1; i < mySystem.getBoxX() / 2; i++) {
//            mySystem.modifyMolecules()[counter].setX(i * 2);
//            mySystem.modifyMolecules()[counter].setY(n);
//            mySystem.modifyMolecules()[counter].setZ(0);
//            mySystem.modifyMolecules().connect(counter - 1, counter);
//            counter++;
//
//        }
//    }
//
//
//    n += 2;
//    for (uint32_t i = 1; i < mySystem.getBoxX() / 2; i++) {
//        mySystem.modifyMolecules()[counter].setX((mySystem.getBoxX() / 2 - i) * 2);
//        mySystem.modifyMolecules()[counter].setY(n);
//        mySystem.modifyMolecules()[counter].setZ(0);
//        mySystem.modifyMolecules().connect(counter - 1, counter);
//        counter++;
//
//    }
//
//
//    while (n > 0) {
//        mySystem.modifyMolecules()[counter].setX(0);
//        mySystem.modifyMolecules()[counter].setY(n);
//        mySystem.modifyMolecules()[counter].setZ(0);
//        mySystem.modifyMolecules().connect(counter - 1, counter);
//        counter++;
//        n -= 2;
//    }






    /* ****************************************************************
      * Now we can set up the task manager with the desired tasks.
      * We want to do three things:
      * 1) Simulate the system. For this we use the updater
      * UpdaterSimpleSimulator. The source code of this updater
      * can be found in src/updater/UpdaterSimpleSimulator.h
      * 2) Write the trajectory to a file. For this we use the analyzer
      * AnalyzerWriteBfmFile. The source code of this analyzer can be
      * found in src/analyzer/AnalyzerWriteBfmFile.h
      * 3) Calculate the radius of gyration. For this we use the
      * analyzer AnalyzerRadiusOfGyration, which  can be found in
      * src/analyzer/AnalyzerRadiusOfGyration.h
      *
      * For adding updaters and analyzers, the task manager provides
      * the functions addUpdater, addAnalyzer. These functions take two
      * arguments:
      * 1) a pointer to the updater/analyzer
      * 2) a number indicating how often the particular task should be
      * executed. For example, 1 means every cycle, 2 every second
      * cycle, etc. Setting the second argument to 0 means the task
      * will be executed only once in the first cycle. This is useful
      * if for example a polymer is added in the beginning by an updater.
      * Here, we have created the polymer by hand, so this is not
      * necessary.
      * ***************************************************************/

    //create the task manager
    TaskManager taskmanager;

    //add the simulators



    std::cout << bridgeS.p1 << '\t' << bridgeS.r1 << std::endl << std::endl;

    //! writing metadata
    std::ofstream myfile;

    myfile.open("Data/metadata_"+ std::to_string(simulation_number)  + ".dat");
    myfile<<polymerLength<<std::endl;
    myfile<<bridgeS.r1<<std::endl;
    myfile<<atof(argc[3])<<std::endl;
    myfile<<bridgeS.r2<<std::endl;
    myfile<<atof(argc[5])<<std::endl;


//
//    taskmanager.addUpdater(new UpdaterSimpleSimulator<MyIngredients, MoveLocalSc>(mySystem, nsteps,simulation_number), 1);


    taskmanager.addUpdater(new UpdaterSimpleSimulatorBridgesAll<MyIngredients, MoveLocalSc, Bridges>(mySystem, nsteps, oriC, terC,
                                                                                                     bridgeS, simulation_number), 0);


    taskmanager.addUpdater(new UpdaterSimpleSimulatorBridgesAll<MyIngredients, MoveLocalSc, Bridges>(mySystem, nsteps, oriC, terC,
            bridgeS, simulation_number), 1);

    //add the file output, the trajectory file name will be "polymer.bfm"
    taskmanager.addAnalyzer(new AnalyzerWriteBfmFile<MyIngredients>(
            "polymerCircular_" + std::to_string(simulation_number) + ".bfm", mySystem), 1);








    /* ****************************************************************
      * For running the desired tasks, the task manager provides four
      * functions:
      * - initialize() : should always be called at the beginning. It
      * initializes the updaters and analyzers. More on this topic in
      * the next example of the tutorial.
      * - run() : runs the circle of chosen tasks until all updaters
      * are finished
      * - run(nCircles) : runs the circle of chosen tasks nCircles times
      * - cleanup() : should be called at the end. finishes the tasks
      * of updaters and analyzers. More on this topic in the next
      * example of the tutorial.
      * ***************************************************************/

    //this will prepare and run the simulation. look in the directory
    //where you called the program for output files.

    taskmanager.initialize();
    taskmanager.run(1000);
    taskmanager.cleanup();

    /* **************************************************************
      * The program should produce among others a file called
      * polymer.bfm. This contains the trajectory. You can visualize
      * the trajectory using the provided LemonadeViewerFLTK from
      * the projects folder.
      * *************************************************************/

    return 0;
}
