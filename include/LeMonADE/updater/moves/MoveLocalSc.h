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

#ifndef LEMONADE_UPDATER_MOVES_MOVELOCALSC_H
#define LEMONADE_UPDATER_MOVES_MOVELOCALSC_H


#include <LeMonADE/updater/moves/MoveLocalBase.h>
#include "LeMonADE/updater/Bridges.h"
#include <math.h>


/*****************************************************************************/
/**
 * @file
 *
 * @class MoveLocalSc
 *
 * @brief Standard local bfm-move on simple cubic lattice for the scBFM.
 *
 * @details The class is a specialization of MoveLocalBase using the (CRTP) to avoid virtual functions.
 **/
/*****************************************************************************/
class MoveLocalSc : public MoveLocalBase<MoveLocalSc> {
public:
    MoveLocalSc() {
        steps[0] = VectorInt3(1, 0, 0);
        steps[1] = VectorInt3(-1, 0, 0);
        steps[2] = VectorInt3(0, 1, 0);
        steps[3] = VectorInt3(0, -1, 0);
        steps[4] = VectorInt3(0, 0, 1);
        steps[5] = VectorInt3(0, 0, -1);

        // Bridging vectors     // Including bond vectors upto length sqrt(6)
        bridgeSteps[0] = VectorInt3(2, 0, 0);
        bridgeSteps[1] = VectorInt3(-2, 0, 0);
        bridgeSteps[2] = VectorInt3(0, 2, 0);
        bridgeSteps[3] = VectorInt3(0, -2, 0);
        bridgeSteps[4] = VectorInt3(0, 0, 2);
        bridgeSteps[5] = VectorInt3(0, 0, -2);

        bridgeSteps[6] = VectorInt3(2, 1, 0);
        bridgeSteps[7] = VectorInt3(2, -1, 0);
        bridgeSteps[8] = VectorInt3(2, 1, 1);
        bridgeSteps[9] = VectorInt3(2, 1, -1);
        bridgeSteps[10] = VectorInt3(2, 0, 1);
        bridgeSteps[11] = VectorInt3(2, 0, -1);
        bridgeSteps[12] = VectorInt3(2, -1, 1);
        bridgeSteps[13] = VectorInt3(2, -1, -1);

        bridgeSteps[14] = VectorInt3(-2, 1, 0);
        bridgeSteps[15] = VectorInt3(-2, -1, 0);
        bridgeSteps[16] = VectorInt3(-2, 1, 1);
        bridgeSteps[17] = VectorInt3(-2, 1, -1);
        bridgeSteps[18] = VectorInt3(-2, 0, 1);
        bridgeSteps[19] = VectorInt3(-2, 0, -1);
        bridgeSteps[20] = VectorInt3(-2, -1, 1);
        bridgeSteps[21] = VectorInt3(-2, -1, -1);

        bridgeSteps[22] = VectorInt3(1, 2, 0);
        bridgeSteps[23] = VectorInt3(-1, 2, 0);
        bridgeSteps[24] = VectorInt3(0, 2, 1);
        bridgeSteps[25] = VectorInt3(0, 2, -1);
        bridgeSteps[26] = VectorInt3(1, 2, 1);
        bridgeSteps[27] = VectorInt3(-1, 0, -1);
        bridgeSteps[28] = VectorInt3(1, 2, -1);
        bridgeSteps[29] = VectorInt3(-1, 2, 1);

        bridgeSteps[30] = VectorInt3(1, -2, 0);
        bridgeSteps[31] = VectorInt3(-1, -2, 0);
        bridgeSteps[32] = VectorInt3(0, -2, 1);
        bridgeSteps[33] = VectorInt3(0, -2, -1);
        bridgeSteps[34] = VectorInt3(1, -2, 1);
        bridgeSteps[35] = VectorInt3(-1, -2, -1);
        bridgeSteps[36] = VectorInt3(1, -2, -1);
        bridgeSteps[37] = VectorInt3(-1, -2, 1);

        bridgeSteps[38] = VectorInt3(1, 0, 2);
        bridgeSteps[39] = VectorInt3(-1, 0, 2);
        bridgeSteps[40] = VectorInt3(0, 1, 2);
        bridgeSteps[41] = VectorInt3(0, -1, 2);
        bridgeSteps[42] = VectorInt3(1, 1, 2);
        bridgeSteps[43] = VectorInt3(-1, -1, 2);
        bridgeSteps[44] = VectorInt3(1, -1, 2);
        bridgeSteps[45] = VectorInt3(-1, 1, 2);

        bridgeSteps[46] = VectorInt3(1, 0, -2);
        bridgeSteps[47] = VectorInt3(-1, 0, -2);
        bridgeSteps[48] = VectorInt3(0, 1, -2);
        bridgeSteps[49] = VectorInt3(0, -1, -2);
        bridgeSteps[50] = VectorInt3(1, 1, -2);
        bridgeSteps[51] = VectorInt3(-1, -1, -2);
        bridgeSteps[52] = VectorInt3(1, -1, -2);
        bridgeSteps[53] = VectorInt3(-1, 1, -2);

    }

    // overload initialise function to be able to set the moves index and direction if neccessary
    template<class IngredientsType>
    void init(const IngredientsType &ing, Bridges &bridges);

    template<class IngredientsType>
    void init(const IngredientsType &ing, uint32_t index);

    template<class IngredientsType>
    void init(const IngredientsType &ing, Bridges &bridges, uint32_t ter);

    template<class IngredientsType>
    void init(const IngredientsType &ing, VectorInt3 dir);

    template<class IngredientsType>
    void init(const IngredientsType &ing, uint32_t index, VectorInt3 dir);

    template<class IngredientsType>
    bool check(IngredientsType &ing);

    template<class IngredientsType>
    void apply(IngredientsType &ing);

private:
    // holds the possible move directions
    /**
     * @brief Array that holds the 6 possible move directions
     *
     * @details In the scBFM the classic moves (dx,dy,dz) are along the lattice-axes as:
     * * steps   = (dx, dy, dz)
     * * steps[0]= ( 1,  0,  0);
     * * steps[1]= (-1,  0,  0);
     * * steps[2]= ( 0,  1,  0);
     * * steps[3]= ( 0, -1,  0);
     * * steps[4]= ( 0,  0,  1);
     * * steps[5]= ( 0,  0, -1);
     */
    VectorInt3 steps[6];
    VectorInt3 bridgeSteps[54];
};



/////////////////////////////////////////////////////////////////////////////
/////////// implementation of the members ///////////////////////////////////

/*****************************************************************************/
/**
 * @brief Initialize the move.
 *
 * @details Resets the move probability to unity. Dice a new random direction and
 * Vertex (monomer) index inside the graph.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 **/
template<class IngredientsType>
void MoveLocalSc::init(const IngredientsType &ing, Bridges &bridges) {
    this->resetProbability();

    //draw index
    uint32_t id = (this->randomNumbers.r250_rand32()) % (ing.getMolecules().size());


    //draw direction
    uint32_t randomDir = this->randomNumbers.r250_rand32() % 6;
    // check if the site is bridged
    if (bridges.isBridged(id) == 1) {

        VectorInt3 vertex(ing.getMolecules()[id]), partner(ing.getMolecules()[bridges.getBridges(id)]), newPos;
        newPos = vertex + steps[randomDir];


        int16_t d = pow((newPos[0] - partner[0]), 2) + pow((newPos[1] - partner[1]), 2) + pow((newPos[2] - partner[2]), 2);
        if (d > 6 ) {
//            std::cout<<d<<'\t';
            this->multiplyProbability(0);
//            std::cout<<id<<std::endl;
        }
    }
    this->setIndex(id);


    this->setDir(steps[randomDir]);
}


//
template<class IngredientsType>
void MoveLocalSc::init(const IngredientsType &ing, Bridges &bridges, uint32_t ter) {

    uint32_t id = (this->randomNumbers.r250_rand32()) %
                  (ing.getMolecules().size()); //random site where bridge should be attempted


    if (bridges.isBridged(id) == 0) {//bridging_everywhere
        std::vector<int> partnerPos;  // the coordinates of potential partner

        VectorInt3 vertexPos(ing.getMolecules()[id]); // coordinates of the selected monomer

//            VectorInt3 nearestPos[6]; // list that holds coordinates of all nearest neighbors

//            for (uint32_t i = 0; i < 6; i++)
//                nearestPos[i] = vertexPos + steps[i] + steps[i];


        uint32_t randomDir = this->randomNumbers.r250_rand32() % (54); // draw a random direction to check bridging

//            std::cout<<vertexPos+bridgeSteps[randomDir]<<std::endl;

        for (uint32_t m = 0; m < ing.getMolecules().size(); m++) {

            VectorInt3 a(ing.getMolecules()[m]);  //coordinate of every monomer

            if (m != (id - 1) && m != (id + 1) && a == (vertexPos + bridgeSteps[randomDir]) &&
                bridges.isBridged(m) == 0) { //bridging_everywhere

                double p = bridges.p1;
                double r = randomNumbers.r250_drand();

                if (r < p) {
                    bridges.setBridges(id, m);

                    break;
                }

            }
        }
//
    }
}



//template<class IngredientsType>
//void MoveLocalSc::init(const IngredientsType &ing, Bridges &bridges, uint32_t ter) {
//
//    uint32_t id = (this->randomNumbers.r250_rand32()) %
//                  (ing.getMolecules().size()); //random site where bridge should be attempted
//
//
//
//    if (bridges.isBridged(id) == 0 && id % 2 == 0) {//bridging_everywhere
////            std::vector<int> partnerPos;  // the coordinates of potential partner
//
//        VectorInt3 vertexPos(ing.getMolecules()[id]); // coordinates of the selected monomer
//
//        VectorInt3 nearestPos[6]; // list that holds coordinates of all nearest neighbors
//
//        for (uint32_t i = 0; i < 6; i++)
//            nearestPos[i] = vertexPos + steps[i] + steps[i];
//
//        uint32_t randomDir = this->randomNumbers.r250_rand32() % (6); // draw a random direction to check bridging
//
//
//        for (uint32_t m = 0; m < ing.getMolecules().size(); m++) {
//
//            VectorInt3 a(ing.getMolecules()[m]);  //coordinate of every monomer
//
//            if (m != (id - 1) && m != (id + 1) && a == nearestPos[randomDir] && bridges.isBridged(m) == 0 &&
//                m % 2 == 0) { //bridging_everywhere
//
//
//
//                double p = bridges.p1;
//                double r = randomNumbers.r250_drand();
//                if (r < p) {
//                    bridges.setBridges(id, m);
//
//                    break;
//                }
//            }
//
//        }
//    }
//
//}

/*****************************************************************************/
/**
 * @brief Initialize the move with a given monomer index.
 *
 * @details Resets the move probability to unity. Dice a new random direction.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 * @param index index of the monomer to be moved
 **/
template<class IngredientsType>
void MoveLocalSc::init(const IngredientsType &ing, uint32_t index) {
    this->resetProbability();

    //set index
    if ((index >= 0) && (index <= (ing.getMolecules().size() - 1)))
        this->setIndex(index);
    else
        throw std::runtime_error("MoveLocalSc::init(ing, index): index out of range!");


    //draw direction
    uint32_t randomDir = this->randomNumbers.r250_rand32() % 6;
    this->setDir(steps[randomDir]);

}

/*****************************************************************************/
/**
 * @brief Initialize the move with a given direction.
 *
 * @details Resets the move probability to unity. Dice a random monomer index and set move direction.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 * @param dir The direction of the move: must be one of the vectors P+-(1,0,0).
 **/
template<class IngredientsType>
void MoveLocalSc::init(const IngredientsType &ing, VectorInt3 dir) {
    this->resetProbability();

    //draw index
    this->setIndex((this->randomNumbers.r250_rand32()) % (ing.getMolecules().size()));


    //set direction
    if (dir == steps[0] ||
        dir == steps[1] ||
        dir == steps[2] ||
        dir == steps[3] ||
        dir == steps[4] ||
        dir == steps[5])
        this->setDir(dir);
    else
        throw std::runtime_error("MoveLocalSc::init(ing, dir): direction vector out of range!");

}

/*****************************************************************************/
/**
 * @brief Initialize the move with a given monomer index and a given direction.
 *
 * @details Resets the move probability to unity and set the move properties index and direction.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 * @param index index of the monomer to be moved
 * @param dir The direction of the move: must be one of the vectors P+-(1,0,0).
 **/
template<class IngredientsType>
void MoveLocalSc::init(const IngredientsType &ing, uint32_t index, VectorInt3 dir) {
    this->resetProbability();

    //set index
    if ((index >= 0) && (index <= (ing.getMolecules().size() - 1)))
        this->setIndex(index);
    else
        throw std::runtime_error("MoveLocalSc::init(ing, index, dir): index out of range!");

    //set direction
    if (dir == steps[0] ||
        dir == steps[1] ||
        dir == steps[2] ||
        dir == steps[3] ||
        dir == steps[4] ||
        dir == steps[5])
        this->setDir(dir);
    else
        throw std::runtime_error("MoveLocalSc::init(ing, index, dir): direction vector out of range!");
}

/*****************************************************************************/
/**
 * @brief Check if the move is accepted by the system.
 *
 * @details This function delegates the checking to the Feature.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 * @return True if move is valid. False, otherwise.
 **/
template<class IngredientsType>
bool MoveLocalSc::check(IngredientsType &ing) {
    //send the move to the Features to be checked
    return ing.checkMove(ing, *this);
}

/*****************************************************************************/
/**
 * @brief Apply the move to the system , e.g. add the displacement to Vertex (monomer) position.
 *
 * @details As first step: all Feature should apply the move using applyMove().\n
 * Second: Modify the positions etc. of the Vertex etc.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 **/
template<class IngredientsType>
void MoveLocalSc::apply(IngredientsType &ing) {
    ///@todo Think about the applying of move. Esp. make this independent of the order to avoid confusion!!

    //move must FIRST be applied to the features
    ing.applyMove(ing, *this);

    //THEN the position can be modified
    ing.modifyMolecules()[this->getIndex()] += this->getDir();

}

#endif
