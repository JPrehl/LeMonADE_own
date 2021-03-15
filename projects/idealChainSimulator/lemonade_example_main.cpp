#include <iostream>

#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureAttributes.h>

#include <LeMonADE/updater/UpdaterAddLinearChains.h>
#include <LeMonADE/updater/UpdaterSimpleSimulator.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>

#include <LeMonADE/utility/RandomNumberGenerators.h>
#include <LeMonADE/utility/TaskManager.h>

int main(int argc, char* argv[])
{
  int nChains(3),chainLength(10),type1(1),nMCS(100),nRuns(10);
  
//  Molecules<VectorInt3,7,int>molecules1; 
  // Molecules with int-coord., max of seven bonds per particle + int-value per bonds
  
	
//  min=10;
//  max=13;
/*  for (int i=min, i<max-1, i++) {
  for (int j=min, j<max-1, j++) {
  for (int k=min, j<max-1, k++) {
	molecules1.connect(i,j,k);
  } } }
*/  

  typedef LOKI_TYPELIST_1(FeatureBox) Features1;
//  const uint max_bonds=4;
  typedef ConfigureSystem<VectorInt3,Features1> Config1;
  typedef Ingredients<Config1> IngredientsType1;
  IngredientsType1 ingredients1;
  
  
  int min=10;
  int max=13;
  int dn=max-min;
  ingredients1.modifyMolecules().resize(dn*dn*dn);

/* Idea: Set-Up a Grid of Particles representing a boundary material 
 * that may interact with linear chains
 */
  
  std::cout<<"ONW_TEST: Ingredients\n";
  int count=0;
  for (int i=min; i<max; i++) {
  for (int j=min; j<max; j++) {
  for (int k=min; k<max; k++) {
	ingredients1.modifyMolecules()[count].setX(i);
	ingredients1.modifyMolecules()[count].setY(j);
	ingredients1.modifyMolecules()[count].setZ(k);

	std::cout<< count << ". particle is at pos. "
			 << ingredients1.getMolecules()[count].getX() << " "
			 << ingredients1.getMolecules()[count].getY() << " "
			 << ingredients1.getMolecules()[count].getZ() << std::endl;

	count++;
  } } }  
  
  std::cout<< ingredients1.getMolecules().size() << std::endl;

  ingredients1.setBoxX(64);
  ingredients1.setBoxY(64);
  ingredients1.setBoxZ(64);
  ingredients1.setPeriodicX(true);
  ingredients1.setPeriodicY(true);
  ingredients1.setPeriodicZ(false);

  ingredients1.synchronize(ingredients1);



/*  typedef LOKI_TYPELIST_2(
    FeatureMoleculesIO, 
    FeatureAttributes<>) Features;
  const uint max_bonds=4;
  typedef ConfigureSystem<VectorInt3,Features,max_bonds> Config;
  typedef Ingredients<Config> IngredientsType;
  IngredientsType ingredients;

  RandomNumberGenerators rng;
  rng.seedAll();

  ingredients.setBoxX(64);
  ingredients.setBoxY(64);
  ingredients.setBoxZ(64);
  ingredients.setPeriodicX(true);
  ingredients.setPeriodicY(true);
  ingredients.setPeriodicZ(true);
  
  ingredients.modifyBondset().addBFMclassicBondset();
  ingredients.synchronize();

  TaskManager taskManager;
  taskManager.addUpdater(new UpdaterAddLinearChains<IngredientsType>(ingredients, nChains,chainLength,type1,type1),0);
  taskManager.addUpdater(new UpdaterSimpleSimulator<IngredientsType,MoveLocalSc>(ingredients,nMCS));
  taskManager.addAnalyzer(new AnalyzerWriteBfmFile<IngredientsType>("config.bfm",ingredients,AnalyzerWriteBfmFile<IngredientsType>::APPEND));
  
  taskManager.initialize();
  taskManager.run(nRuns);
  taskManager.cleanup();
 */ 
  return 0;
}