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

#ifndef LEMONADE_UPDATER_UpdaterSimpleSimulatorBridgesAll_H
#define LEMONADE_UPDATER_UpdaterSimpleSimulatorBridgesAll_H

#include <LeMonADE/updater/AbstractUpdater.h>
#include <LeMonADE/updater/moves/MoveLocalBase.h>
#include "LeMonADE/updater/Bridges.h"
#include <iostream>
#include <fstream>
#include <math.h>
/**
 * @file
 *
 * @class UpdaterSimpleSimulatorBridges
 *
 * @brief Simple simulation updater for general purpose.
 *
 * @details It takes the type of move as template argument MoveType
 * and the number of mcs to be executed as argument for the constructor
 *
 * @tparam IngredientsType Ingredients class storing all system information( e.g. monomers, bonds, etc).
 * @tparam MoveType name of the specialized move.
 */
template<class IngredientsType, class MoveType, class BridgeType>
class UpdaterSimpleSimulatorBridgesAll : public AbstractUpdater {

public:
    /**
     * @brief Standard Constructor initialized with ref to Ingredients and MCS per cycle
     *
     * @param ing a reference to the IngredientsType - mainly the system
     * @param steps MCS per cycle to performed by execute()
     */
    UpdaterSimpleSimulatorBridgesAll(IngredientsType &ing, uint32_t steps, uint32_t oriC, uint32_t terC,
                                     Bridges &bridgeS, uint32_t simulation_number)
            : ingredients(ing), nsteps(steps), ori(oriC), ter(terC), bridge(bridgeS),
              simulation_number(simulation_number) {}

    /**
     * @brief This checks all used Feature and applies all Feature if all conditions are met.
     *
     * @details This function runs over \a steps MCS and performs the moves.
     * It setting the age of the system and prints a simple simple simulation speed
     * in the number of attempted monomer moves per s (tried and performed monomer moves/s).
     *
     * @return True if function are done.
     */
    bool execute() {
        time_t startTimer = time(NULL); //in seconds
        std::cout << "mcs " << ingredients.getMolecules().getAge() << " passed time "
                  << ((difftime(time(NULL), startTimer))) << std::endl;


        std::ofstream myfile;


        myfile.open("Polymer_coordinates/configs_" + std::to_string(simulation_number) + ".dat", std::ios_base::app);
        for (int i = 0; i < ingredients.getMolecules().size(); i++)
            myfile << ingredients.getMolecules()[i].getX() << '\t' << ingredients.getMolecules()[i].getY()
                   << '\t'
                   << ingredients.getMolecules()[i].getZ() << '\t';
        myfile << std::endl;
        myfile.close();
        bridge.count++;


         //! Coordination matrix for Hi-C
        myfile.open("Data/bridgeDistances_" + std::to_string(simulation_number) + ".dat",std::ios_base::app);
        for (int i = 0; i < ingredients.getMolecules().size(); i++)
        {
            int gen_dist = 0;
            if (bridge.getBridges(i)>0)
                gen_dist=abs(bridge.getBridges(i)-i);

            myfile << gen_dist<< '\t';}
        myfile << std::endl;
        myfile.close();


        int min = ter, max = ter + 1, range = max - min; // ter positions

        for (uint32_t n = 0; n < nsteps; n++) {

            for (size_t m = 0; m < ingredients.getMolecules().size(); m++) {


                // Local monte carlo moves
                move.init(ingredients, bridge);

                if (move.check(ingredients) == true) {

                    move.apply(ingredients);
                }


                // Bridge moves


                move.init(ingredients, bridge, ter);
//                for (size_t id = 0; id < ingredients.getMolecules().size(); id++) {
//                int id = randomNumbers.r250_rand32()%ingredients.getMolecules().size();
//                if (bridge.isBridged(id) == 1) {
//                    double p = randomNumbers.r250_drand();
//
//                    if (p < 1/(double)tauO) {
//                        {
////                            std::cout<<id<<std::endl;
//                            bridge.resetBridges(id);
//                        }
//                    }
////                }
//            }


            }
            for (size_t id = 0; id < ingredients.getMolecules().size(); id++) {
                if (bridge.isBridged(id) == 1) {
                    double p = randomNumbers.r250_drand();

                    if (p < 1/(double)tauO) {
                        {
//                            std::cout<<id<<std::endl;
                            bridge.resetBridges(id);
                        }
                    }
                }
            }

            /***Exponential bridge lifetime distribution***/


//            /*** Changed bridging lifetimes of ori and ter**/
//            for (size_t id = 0; id < ingredients.getMolecules().size(); id++) {
//
//                if (bridge.isBridged(id) == 1) {
//                    int partner = bridge.getBridges(id);
//
//                    bridge.lifetime[id]++;
//                    if ((unsigned) (id - min) >= range && (unsigned) (partner - min) >= range &&
//                        bridge.lifetime[id] % (tauO) == 0) {
//                        bridge.resetBridges(id);
//                        bridge.lifetime[id] = 0;
//                    }
//                    if ((unsigned) (id - min) >= range && (unsigned) (partner - min) <= range &&
//                        bridge.lifetime[id] % (int(sqrt(tauO * tauT))) == 0) {
//                        bridge.resetBridges(id);
//                        bridge.lifetime[id] = 0;
//                    }
//                    if ((unsigned) (id - min) <= range && (unsigned) (partner - min) >= range &&
//                        bridge.lifetime[id] % (int(sqrt(tauO * tauT)))== 0) {
//
//                        bridge.lifetime[id] = 0;
//                    }
//                    if ((unsigned) (id - min) <= range && (unsigned) (partner - min) <= range &&
//                        bridge.lifetime[id] % (tauT) == 0) {
//                        bridge.resetBridges(id);
//                        bridge.lifetime[id] = 0;
//                    }
//                }
//            }
        }


//                if (n % 1000 == 0) {

        //! Writing polymer coordinates



//                }




        int av_bridges = 0;

        for (size_t id = 0; id < ingredients.getMolecules().size(); id++) {
            if (bridge.isBridged(id) == 1) {

                bridge.totalbridges[id]++;

                av_bridges++;


                int partner = bridge.getBridges(id);
//                bridge.coordMatrix[id][partner]++;
            }


        }


        bridge.average_bridges.push_back(av_bridges / 2);


        ingredients.modifyMolecules().setAge(ingredients.modifyMolecules().getAge() + nsteps);


        std::cout << "mcs " << ingredients.getMolecules().getAge() << " with "
                  << (((1.0 * nsteps) * ingredients.getMolecules().size()) / (difftime(time(NULL), startTimer)))
                  << " [attempted moves/s]" << std::endl;
        std::cout << "mcs " << ingredients.getMolecules().getAge() << " passed time "
                  << ((difftime(time(NULL), startTimer))) << " with " << nsteps << " MCS " << std::endl;

        return true;
    }

    /**
     * @brief This function is called \a once in the beginning of the TaskManager.
     *
     * @details It´s a virtual function for inheritance.
     * Use this function for initializing tasks (e.g. init SDL)
     *
     **/
    virtual void initialize() {


    };

    /**
     * @brief This function is called \a once in the end of the TaskManager.
     *
     * @details It´s a virtual function for inheritance.
     * Use this function for cleaning tasks (e.g. destroying arrays, file outut)
     *
     **/
    virtual void cleanup() {

        std::ofstream myfile;

//        //! number of bridges made by a particular monomer
//        myfile.open("Data/totalbridges_" + std::to_string(simulation_number) + ".dat");
//        for (int i = 0; i < ingredients.getMolecules().size(); i++)
//            myfile << bridge.totalbridges[i] << std::endl;
//        myfile.close();

//      myfile.open ("Bridges/number_of_bridges_"+std::to_string(simulation_number)+".dat");
//      for(int i=0; i<ingredients.getMolecules().size();i++)
//        myfile <<bridge.number_of_bridges[i]<<std::endl;
//      myfile.close();
        //! average bridges made during the simulation on the entire polymer
        myfile.open("Data/average_bridges_" + std::to_string(simulation_number) + ".dat");
        myfile << simulation_number<<std::endl;
        myfile << ingredients.getMolecules().size() << std::endl;
        myfile << bridge.r1 << std::endl;
        myfile << bridge.p1 << std::endl;
        myfile << bridge.r2 << std::endl;
        myfile << bridge.p2 << std::endl;
        for (int i = 0; i < bridge.average_bridges.size(); i++)
            myfile << bridge.average_bridges[i] << '\t';
        myfile.close();

//        //! Coordination matrix for Hi-C
//        myfile.open("Data/coordMat_" + std::to_string(simulation_number) + ".dat");
//        for (int i = 0; i < ingredients.getMolecules().size(); i++) {
//            for (int j = 0; j < ingredients.getMolecules().size(); j++)
//                myfile << bridge.coordMatrix[i][j] << '\t';
//            myfile << std::endl;
//        }
//        myfile.close();
    };


private:
    //! A reference to the IngredientsType - mainly the system
    IngredientsType &ingredients;

    //! Specialized move to be used
    MoveType move;

    //! Number of mcs to be executed
    uint32_t nsteps;

    //! ori position
    uint32_t ori;

    //! ter position
    uint32_t ter;

    //! simulation number to be appended to file name
    uint32_t simulation_number;


    BridgeType &bridge;


    RandomNumberGenerators randomNumbers;

    int tauO = bridge.r1;

    int tauT = bridge.r2;

};

#endif
