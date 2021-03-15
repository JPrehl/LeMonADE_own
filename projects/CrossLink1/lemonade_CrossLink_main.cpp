/*	LeMonADE Simulation 
 * 
 *  Test simulation to set an initial simualtion box
 */

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

#include "../../updater/UpdaterCreateCrossLink.h"

int main(int argc, char* argv[])
{
  int chainLength(10),nRuns(1);
  int boxX(32), boxY(32), boxZ(16);
  std::string modus("COHOMO");				// definition of polymer structure
  double nnInteraction(5.);					// WW (repulsive)
  bool pX(1), pY(1),pZ(0);					// periodicity
  std::string filename("crossLinkInit.bfm");// Output-File 

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
  taskManager.addUpdater(new UpdaterCreateCrossLink<IngredientsType>(ingredients, boxX, boxY, boxZ, pX, pY, pZ, nnInteraction, modus, chainLength),0);
  // -> initialize (periodicity) -> execute (set beads -> Walls, Chain) -> cleanup ()
  // -> (...,0) -> nPeriod=0 (default 1) wird eingerichtet aber nicht einmal ausgeführt
  taskManager.addAnalyzer(new AnalyzerWriteBfmFile<IngredientsType>(filename,ingredients,AnalyzerWriteBfmFile<IngredientsType>::OVERWRITE));
  // -> Bfm-File wird ausgegeben
  
  taskManager.initialize();
  taskManager.run(nRuns);
  taskManager.cleanup();
  
  // -> alles wird einmal durchlaufen um das System zu erstellen
 
  return 0;
}
