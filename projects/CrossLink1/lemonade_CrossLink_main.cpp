#include <iostream>
#include <string>

#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/feature/FeatureNNInteractionSc.h>
#include <LeMonADE/feature/FeatureFixedMonomers.h>

#include <LeMonADE/updater/UpdaterSimpleSimulator.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>

#include <LeMonADE/utility/RandomNumberGenerators.h>
#include <LeMonADE/utility/TaskManager.h>

#include "UpdaterCreateCrossLink.h"

int main(int argc, char* argv[])
{
  int chainLength(10),nRuns(1);
  int boxX(64), boxY(64), boxZ(64);
  double interaction(1.);
  bool pX(0), pY(0),pZ(0);
  std::string filename("crossLink.bfm") 

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
  taskManager.addUpdater(new UpdaterCreateCrossLink<IngredientsType>(ingredients, boxX, boxY, boxZ, pX, pY, pZ, chainLenth, interaction),0);
  taskManager.addAnalyzer(new AnalyzerWriteBfmFile<IngredientsType>(filename,ingredients,AnalyzerWriteBfmFile<IngredientsType>::APPEND));
  
  taskManager.initialize();
  taskManager.run(nRuns);
  taskManager.cleanup();
 
  return 0;
}
