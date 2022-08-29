#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/feature/FeatureNNInteractionScMod.h>

#include <LeMonADE/updater/UpdaterAddLinearChains.h>
#include <LeMonADE/updater/UpdaterSimpleSimulator.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>

#include <LeMonADE/utility/RandomNumberGenerators.h>
#include <LeMonADE/utility/TaskManager.h>

int main(int argc, char* argv[])
{
    int nChains(1),chainLength(48),type1(1),type2(2),nMCS(100),nRuns(1000), nSolvent(0);

    typedef LOKI_TYPELIST_4(
            FeatureMoleculesIO,
            FeatureAttributes<>,
            FeatureExcludedVolumeSc<FeatureLatticePowerOfTwo<uint8_t> >,
            FeatureNNInteractionSc<FeatureLatticePowerOfTwo>) Features;
    const uint max_bonds = 3;
    typedef ConfigureSystem<VectorInt3,Features,max_bonds> Config;
    typedef Ingredients<Config> IngredientsType;
    IngredientsType ingredients;

    RandomNumberGenerators rng;
    rng.seedAll();

 // Position of ter start

    ingredients.setBoxX(64);
    ingredients.setBoxY(64);
    ingredients.setBoxZ(32);
    ingredients.setPeriodicX(false);
    ingredients.setPeriodicY(false);
    ingredients.setPeriodicZ(false);

    ingredients.modifyBondset().addBFMclassicBondset();
    ingredients.setNNInteraction(1,2,1);

    ingredients.synchronize();

    TaskManager taskManager;
    taskManager.addUpdater(new UpdaterAddLinearChains<IngredientsType>(ingredients, nChains,chainLength,type1,type2),0);
//    taskManager.addUpdater(new UpdaterSimpleSimulator<IngredientsType,MoveLocalSc>(ingredients,nMCS));
    taskManager.addAnalyzer(new AnalyzerWriteBfmFile<IngredientsType>("config_nn_4.bfm",ingredients,AnalyzerWriteBfmFile<IngredientsType>::APPEND));

    taskManager.initialize();
    taskManager.run(nRuns);
    taskManager.cleanup();

    return 0;
}