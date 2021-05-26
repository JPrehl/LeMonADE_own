/*	LeMonADE Simulation 
 * 
 *  Test simulation to set an initial simualtion box
 */

#include <iostream>
#include <string>
#include <stdlib.h> //for atoi
#include <unistd.h> //for getopt

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
  // parameters sorted by type
  std::string filename("crossLinkInit.bfm");// Output-File
  std::string modus("COHOMO");				// definition of polymer structure

  double nnInteraction(5.);					// WW (repulsive)

  int chainLength(10);
  int boxX(32), boxY(32), boxZ(16);

  // peridocity not taken from command line -> periodic in xy, non periodic in z
  bool pX(1), pY(1), pZ(0);					// periodicity

  // utility for getopt
  int option_char(0);
  //read in options by getopt
    while ((option_char = getopt (argc, argv, "f:m:i:c:x:y:z:h"))  != EOF){
        switch (option_char)
        {  
            case 'f':
                filename = std::string(optarg);
                break;
            case 'm':
                modus = std::string(optarg);
                break;
            case 'i':
                nnInteraction = atof(optarg);
                break;
            case 'c':
                chainLength = atoi(optarg);
                break;
            case 'x':
                boxX = atoi(optarg);
                break;
            case 'y':
                boxY = atoi(optarg);
                break;
            case 'z':
                boxZ = atoi(optarg);
                break;
            case 'h':
                std::cerr << "Usage:\n" << argv[0] << "\n[-f filename] [-c chainlength] [-x boxX] [-y boxY] [-z boxZ] [-t torque] [-w wallSwitch (1: with wall, 0: periodic box)]\n";
                return 0;
            case '?':
                printf("Unknown option: %c\n", optopt);
                break;
        }
    }
    // write out parameter list
    std::cout << 
    "filename = " << filename << 
    "\nnnInteraction = " << nnInteraction <<
    "\nchainlength = " << chainLength <<
    "\nboxX = " << boxX <<
    "\nboxY = " << boxY << 
    "\nboxZ = " << boxZ << 
    "\nmodus = " << modus <<std::endl;

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
  taskManager.run(1);
  taskManager.cleanup();
  
  // -> alles wird einmal durchlaufen um das System zu erstellen
 
  return 0;
}
