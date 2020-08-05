#include <iostream>
#include <string>

#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/feature/FeatureNNInteractionSc.h>
#include <LeMonADE/feature/FeatureFixedMonomers.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>


#include <LeMonADE/updater/UpdaterSimpleSimulator.h>
#include <LeMonADE/updater/UpdaterReadBfmFile.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>

#include <LeMonADE/utility/RandomNumberGenerators.h>
#include <LeMonADE/utility/TaskManager.h>

#include "../../updater/UpdaterCreateCrossLink.h"

int main(int argc, char* argv[])
{
  int nRuns(1000), nSteps(100);
  std::string filenameI("crossLinkInit.bfm"); 
  std::string filenameR("crossLinkRun.bfm");

  typedef LOKI_TYPELIST_4(
    FeatureMoleculesIO, 
    FeatureAttributes<>,
    FeatureNNInteractionSc<FeatureLatticePowerOfTwo>,
    FeatureFixedMonomers) Features;
  const uint max_bonds=4;
  typedef ConfigureSystem<VectorInt3,Features,max_bonds> Config;
  typedef Ingredients<Config> IngredientsType;
  IngredientsType ingredients;

  RandomNumberGenerators rng;
  rng.seedAll();


  TaskManager taskManager;
  taskManager.addUpdater(new UpdaterReadBfmFile<IngredientsType>(filenameI, ingredients, UpdaterReadBfmFile<IngredientsType>::READ_STEPWISE));
  taskManager.addUpdater(new UpdaterSimpleSimulator<IngredientsType,MoveLocalSc>(ingredients, nSteps));
  taskManager.addAnalyzer(new AnalyzerWriteBfmFile<IngredientsType>(filenameR,ingredients,AnalyzerWriteBfmFile<IngredientsType>::APPEND));
  
  taskManager.initialize();
  taskManager.run(nRuns);
  taskManager.cleanup();
 
  return 0;
}
