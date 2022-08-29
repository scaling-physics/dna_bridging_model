//
// Created by subraman on 23.06.20.
//

#ifndef LEMONADE_BRIDGES_H
#define LEMONADE_BRIDGES_H


#include <vector>

class Bridges {

public:



    void setBridges (int vertex, int partner) {
        bridges[vertex] = partner;                      // stores the location of partner monomer at location of the current vertex
        bridges[partner] = vertex;

        isBridge[vertex]=1;                             // stores bridging information
        isBridge[partner]=1;
    }

    int getBridges (int vertex) {                       // get the bridge partner
        return bridges[vertex];

    }

    bool isBridged(int vertex){                         // checks if monomer at a particular location is bridged
        if(isBridge[vertex]==1)
            return true;
        else
            return false;
    }

    void resetBridges(int vertex)                   // delete bridges
    {

            int partner = bridges[vertex];
            isBridge[vertex] = 0;
            bridges[vertex] = -1;
            isBridge[partner] = 0;
            bridges[partner] = -1;

    }

    int lifetime[1000]={};                                 // bridge lifetime
    int number_of_bridges[1000] = {};
    int unbridgedtime = 0;                                 // unbridge lifetime
    int totalbridges[1000] = {};
    std::vector<int> average_bridges = {};
    int coordMatrix[1000][1000]={};
    int count =0;

    double p1;  // bridging probability of ori
    double p2;  // bridging probability of ter
    double theta;  // scaling ratio
    double r1; // bridge timescale change
    double r2; // bridge timescale change




private:

    int bridges[1000]={};       // makes a list for monomers and their partners
    int isBridge[1000]={};

};


#endif //LEMONADE_BRIDGES_H
