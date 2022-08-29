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
#include <fstream>
#include <boost/lexical_cast.hpp>

int main(int argv, char **argc) {
    int nChains(1), chainLength(400), type1(1), nMCS(500), nRuns(5000);

    RandomNumberGenerators randomNumbers;
    randomNumbers.seedAll();

    typedef LOKI_TYPELIST_4(
            FeatureMoleculesIO,
            FeatureAttributes<>,
            FeatureExcludedVolumeSc<FeatureLattice<>>, FeatureBoltzmann) Features;
    const uint max_bonds = 2;
    typedef ConfigureSystem<VectorInt3, Features, max_bonds> Config;
    typedef Ingredients<Config> IngredientsType;
    IngredientsType ingredients;

    RandomNumberGenerators rng;
    rng.seedAll();
    int L=95;
    ingredients.setBoxX(L);
    ingredients.setBoxY(L);
    ingredients.setBoxZ(L);
    ingredients.setPeriodicX(true);
    ingredients.setPeriodicY(true);
    ingredients.setPeriodicZ(true);
    ingredients.modifyBondset().addBFMclassicBondset();
    ingredients.synchronize();

    int oriC = 0; //Position of Ori start
    uint32_t terC = chainLength / 2;

//  u_int32_t simulation_number = 2;

    Bridges bridgeS;
    uint32_t simulation_number = atoi(argc[1]);
    bridgeS.r1 = atof(argc[2]);
    bridgeS.p1 = atof(argc[3]) / atof(argc[2]);
    bridgeS.r2 = bridgeS.r1;
    bridgeS.p2 = bridgeS.p1;
    TaskManager taskManager;

//    /// create polymer chain
//    int file_num = randomNumbers.r250_rand32()%8+1;
//    std::ifstream f("/scratch2/Srikanth/confinement_2/seed/Polymer_coordinates/configs_"+std::to_string(file_num)+".dat");
//    int lim=100;
//
//    int configs[300][1200];
//
//
//    for (int i = 0; i < lim; i++)
//    {
//        for (int j = 0; j < chainLength*3; j++) {
//            f >> configs[i][j];
////            std::cout<<configs[i][j]<<'\t';
//        }
////        std::cout<<std::endl;
//    }
//    f.close();
//
//
//    int num = randomNumbers.r250_rand32()%lim;
//
//    ingredients.modifyMolecules().resize(chainLength);
//    ingredients.modifyMolecules()[0].setX(configs[num][0]);
//    ingredients.modifyMolecules()[0].setY(configs[num][1]);
//    ingredients.modifyMolecules()[0].setZ(configs[num][2]);
//
//
//
//    int j=3;
//    for (int i = 1; i < chainLength; i++) {
//        ingredients.modifyMolecules()[i].setX(configs[num][j]);
//        ingredients.modifyMolecules()[i].setY(configs[num][j + 1]);
//        ingredients.modifyMolecules()[i].setZ(configs[num][j + 2]);
//        ingredients.modifyMolecules().connect(i - 1, i);
//        j = j + 3;
//    }
//
//    std::cout<<ingredients.getMolecules()[0].getX()<<'\t'<<ingredients.getMolecules()[1].getX();
//
//    ingredients.synchronize(ingredients);

    taskManager.addUpdater(new UpdaterAddLinearChains<IngredientsType>(ingredients, nChains, chainLength, type1, type1),
                           0);
    taskManager.addUpdater(
            new UpdaterSimpleSimulatorBridgesAll<IngredientsType, MoveLocalSc, Bridges>(ingredients, bridgeS.r1, oriC, terC,
                                                                                        bridgeS, simulation_number), 0);

    taskManager.addUpdater(
            new UpdaterSimpleSimulatorBridgesAll<IngredientsType, MoveLocalSc, Bridges>(ingredients, nMCS, oriC, terC,
                                                                                        bridgeS, simulation_number), 1);
    std::ofstream myfile;
    myfile.open("Polymer_coordinates/configs_" + std::to_string(simulation_number)+".dat");
    taskManager.addAnalyzer(new AnalyzerWriteBfmFile<IngredientsType>(
            "polymerLinear_" + boost::lexical_cast<std::string>(simulation_number) + ".bfm", ingredients, 1));


    taskManager.initialize();
    taskManager.run(nRuns);
    taskManager.cleanup();

    return 0;
} 
