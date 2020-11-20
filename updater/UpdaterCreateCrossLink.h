/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2013-2015 by 
  o/.|.\o    E   nvironment    | LeMonADE Principal Developers
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

#ifndef LEMONADE_UPDATER_CREATE_CROSSLINK_H
#define LEMONADE_UPDATER_CREATE_CROSSLINK_H
/**
 * @file
 *
 * @class UpdaterCreateCrosslink
 *
 * @brief create Updater for systems of project Crosslink
 * 
 * @details ...
 *
 * @tparam IngredientsType
 *
 **/

#include <LeMonADE/updater/UpdaterAbstractCreate.h>
#include <LeMonADE/utility/Vector3D.h>

#include <LeMonADE/utility/RandomNumberGenerators.h>


template<class IngredientsType>
class UpdaterCreateCrossLink: public UpdaterAbstractCreate<IngredientsType>
{
  typedef UpdaterAbstractCreate<IngredientsType> BaseClass;
  
/**
 * @brief Member aus der UpdaterAbstractCreate-Klasse
 **/
  using BaseClass::ingredients;	

/**
 * @brief Methoden aus der UpdaterAbstractCreate-Klasse
 **/
/* ist das gleiche wie //! */
  using BaseClass::addMonomerToParent;
  using BaseClass::addSingleMonomer;
  using BaseClass::addMonomerAtPosition;
  using BaseClass::linearizeSystem;
  
  uint32_t boxX, boxY, boxZ;
  bool pX, pY, pZ;
  
  uint32_t chainLength; 
  std::string modus;
  double nnInteraction;
  
  /*	Number of different beads is set here 
   * -> if needed more flexible: CHANGE HERE
   */
  int typAlO;
  int typPA6;
  int typBlockA;
  int typBlockB;
  
  int lBlockA;
  int lBlockB;
  
  bool isExecuted;
  
  bool setWallMonomer();
  bool setLinearCoPolyChain();			// first verion
  bool setLinearBlockCoPolyChain();
										// several blocks of HomoPolymers
  bool setLinearCoHomoPolyChain();		// 2 HomoPolymers glued together
  bool setLinearStatCoPolyChain();		// statistical Polymer

protected:
	//! Random Number Generator (RNG)
	RandomNumberGenerators randomNumbers;

	
public:
  UpdaterCreateCrossLink(IngredientsType& ingredients_, uint32_t boxX_, uint32_t boxY_, uint32_t boxZ_, bool pX_, bool pY_, bool pZ_, double nnInteraction_, std::string modus_, uint32_t chainLength_, uint32_t lBlockA_=0, uint32_t lBlockB_=0);
  
  virtual void initialize();
  virtual bool execute();
  virtual void cleanup();

};

/*****************************************************************************
 *  Necessary (virtual) Methods + Constructor
 * ***************************************************************************
 */

template < class IngredientsType >
UpdaterCreateCrossLink<IngredientsType>::UpdaterCreateCrossLink(IngredientsType& ingredients_, uint32_t boxX_, uint32_t boxY_, uint32_t boxZ_, bool pX_, bool pY_, bool pZ_, double nnInteraction_, std::string modus_, uint32_t chainLength_, uint32_t lBlockA_, uint32_t lBlockB_):
  BaseClass(ingredients_), boxX(boxX_), boxY(boxY_), boxZ(boxZ_), pX(pX_), pY(pY_), pZ(pZ_), nnInteraction(nnInteraction_), modus(modus_), chainLength(chainLength_), lBlockA(lBlockA_),lBlockB(lBlockB_), isExecuted(false)
{}


template < class IngredientsType >
void UpdaterCreateCrossLink<IngredientsType>::initialize(){
//! Erzeugen von der Box, für Header von bfm-File

/* Version um die Standard-Output umzuleiten(!) fürs Rechenzentrum -> A
    // supress std::cout output of addBondset and synchronize
    std::streambuf *old = std::cout.rdbuf(); // <-- save
    std::stringstream ss;
    std::cout.rdbuf (ss.rdbuf());       // <-- redirect
*/
	
//	std::cout << "Updater CL: init" << std::endl;
	//set box size
    ingredients.setBoxX(boxX);
    ingredients.setBoxY(boxY);
    ingredients.setBoxZ(boxZ);

	//set periodicity
    ingredients.setPeriodicX(pX);
    ingredients.setPeriodicY(pY);
    ingredients.setPeriodicZ(pZ);

    //set Bondset
    ingredients.modifyBondset().addBFMclassicBondset();
	
	typPA6=1;
	typAlO=2;
	
	typBlockA=3;
	typBlockB=4;
	
    ingredients.setNNInteraction(typBlockA,typAlO,nnInteraction); // repulsive WW
    ingredients.setNNInteraction(typBlockB,typPA6,nnInteraction);
	// Note: attraction with negative nnInteraction
	
	// call synchronize
    ingredients.synchronize(ingredients);
/* Version um die Standard-Output umzuleiten(!) fürs Rechenzentrum -> B
    std::cout.rdbuf (old);
*/
	execute();
} 


template < class IngredientsType >
bool UpdaterCreateCrossLink<IngredientsType>::execute(){

//	std::cout << "Updater CL: execute -- rng" << randomNumbers.r250_drand() << std::endl;
	
	if(!isExecuted) {
//	std::cout << "Updater CL: execute" << std::endl;
		setWallMonomer();
		setLinearCoPolyChain();
		
		isExecuted=true;
	}
	return true;
}


template < class IngredientsType >
void UpdaterCreateCrossLink<IngredientsType>::cleanup(){

	std::cout << "Updater CCL: cleanup" << std::endl;
}

/*****************************************************************************
 *  Non-virtual Methods
 * ***************************************************************************
 */
template < class IngredientsType >
bool UpdaterCreateCrossLink<IngredientsType>::setWallMonomer(){

	int anzX(boxX/2);
	int anzY(boxY/2);
	int zmin(0), zmax(boxZ-2);
	// theortisch check auf Box-Größe%2==0
	
	for (int i=0; i<anzX; i++) {
		for (int j=0; j<anzY; j++) {
			int x = i*2, y = j*2;
			//! Add PA6-monomer
			addMonomerAtPosition(VectorInt3(x,y,zmin),typPA6);
			//! changes the latest added Monomer to fixed
			ingredients.modifyMolecules()[ingredients.getMolecules().size()-1].setMovableTag(false);
			
			//! Add first AlO-monomer
			addMonomerAtPosition(VectorInt3(x,y,zmax),typAlO);
			//! changes the latest added Monomer to fixed
			ingredients.modifyMolecules()[ingredients.getMolecules().size()-1].setMovableTag(false);
		}
	}
	return true;
}


template < class IngredientsType >
bool UpdaterCreateCrossLink<IngredientsType>::setLinearCoPolyChain(){

	if ("COHOMO" == modus) { 

		//! Two different homo-polymers glued together
		//! chemically feasible with different length ratios
		setLinearCoHomoPolyChain();
	} else if ("STAT" == modus) {
		
		//! Statistical Copolymer
		setLinearStatCoPolyChain();
	} else if ("BLOCK_1A" == modus || "BLOCK_1B" == modus) {
		
		//! Block-Co-Polymer with block lengths lBlockA, lBlockB
		//! chemically hard to synthesize
		setLinearBlockCoPolyChain();
	} else {
		
		int typBlock;
		if (0 == lBlockA) {
			typBlock = typBlockB;
		} else {
			typBlock = typBlockA;
		}

		addSingleMonomer(typBlock);	// Start-Monomer der Chain
		int idStartMono = ingredients.getMolecules().size()-1; 

		for (int l=0; l<chainLength-1; l++) {
			addMonomerToParent(idStartMono+l,typBlock); 
		}
	}
	
	return true;
}

template < class IngredientsType >
bool UpdaterCreateCrossLink<IngredientsType>::setLinearCoHomoPolyChain(){

	std::cout << "Updater CCL: setLinearCoHomoPolyChain()" << std::endl;

	int typB1;
	int typB2;
	//! Random definition of typBlock-order
	if (randomNumbers.r250_drand() < 0.5) {
		typB1 = typBlockA;
		typB2 = typBlockB;
	} else {
		typB1 = typBlockB;
		typB2 = typBlockA;
	}
	
	// Block lengths not defined
	int lBlock;
	if (0 == lBlockA) {
		lBlock = chainLength / 2;
	} else {
		lBlock = lBlockA;
	}
	
	addSingleMonomer(typB1);	// Start-Monomer der Chain
	int idStartMono = ingredients.getMolecules().size()-1; 
	
	//! Two different homo-polymers glued together
	for (int l=0; l<lBlock-1; l++) {
		addMonomerToParent(idStartMono+l,typB1); 
	}
	for (int l=lBlock-1; l<chainLength-1; l++) {
		addMonomerToParent(idStartMono+l,typB2); 
	}	
	return true;
}

template < class IngredientsType >
bool UpdaterCreateCrossLink<IngredientsType>::setLinearStatCoPolyChain(){

	std::cout << "Updater CCL: setLinearStatCoPolyChain()" << std::endl;
	//! Random order to beads 
	int typBlock;

	//! Random definition of Start-Monomer-Typ of Chain
	if (randomNumbers.r250_drand() < 0.5) {
		typBlock = typBlockA;
	} else {
		typBlock = typBlockB;
	}
	
	addSingleMonomer(typBlock);	// Start-monomer der Chain
	int idStartMono = ingredients.getMolecules().size()-1; 
	
	for (int l=0; l<chainLength-1; l++) {
		if (randomNumbers.r250_drand() < 0.5) {
			typBlock = typBlockA;
		} else {
			typBlock = typBlockB;
		}
		addMonomerToParent(idStartMono+l,typBlock); 
	}
	return true;
}

template < class IngredientsType >
bool UpdaterCreateCrossLink<IngredientsType>::setLinearBlockCoPolyChain(){

	std::cout << "Updater CCL: setLinearBlockCoPolyChain() - " << modus << std::endl;
	//! 1) 	fixed length of lBlock for both types 
	//!  a) repeated order (chemical resonable?)
	//!  a) random order of blocks thus "random lengths"
	//!  Req.: lBlock needs to be factor/divisor of chainLength
	
	
	if ("BLOCK_1B" == modus) {
		//! Random order to bead-blocks 
		int typBlock;

		//check for Requirement
		if (0 != (chainLength % lBlockA)) {
			throw std::runtime_error("UpdaterCreateCrossLink::setLinearBlockCoPolyChain(): lBlock not divisor of chainLength!");	
		}
		int lBlock = lBlockA;	 // for 1) lblockA == blockB

		//! Random definition of Start-Monomer-Typ of Chain
		if (randomNumbers.r250_drand() < 0.5) {
			typBlock = typBlockA;
		} else {
			typBlock = typBlockB;
		}
		
		addSingleMonomer(typBlock);	// Start-Monomer der Chain
		int idStartMono = ingredients.getMolecules().size()-1; 
		
		for (int l=0; l<chainLength-1; l++) {
			if ( 0 == (l+1) % lBlock ) {
				std::cout << "UpdaterCreateCrossLink::setLinearBlockCoPolyChain(): Index=" << l+1 << ", lBlock=" << lBlock << std::endl;
				if (randomNumbers.r250_drand() < 0.5) {
					typBlock = typBlockA;
				} else {
					typBlock = typBlockB;
				}
			}
			addMonomerToParent(idStartMono+l,typBlock); 
		}
		
	} else {	// "BLOCK_1A" == modus // -> Default modus

		//! Repeated order to bead-blocks
		int typBlock, typB1, typB2;

		//check for Requirement
		if (0 != (chainLength % lBlockA)) {
			throw std::runtime_error("UpdaterCreateCrossLink::setLinearBlockCoPolyChain(): lBlock not divisor of chainLength!");	
		}
		int lBlock = lBlockA;	 // for 1) lblockA == blockB

		//! Random Block-Type-Order
		if (randomNumbers.r250_drand() < 0.5) {
			typB1 = typBlockA;
			typB2 = typBlockB;
		} else {
			typB1 = typBlockB;
			typB2 = typBlockA;
		}
		
		typBlock=typB1;
		addSingleMonomer(typBlock);	// Start-Monomer der Chain
		int idStartMono = ingredients.getMolecules().size()-1; 
		
		for (int l=0; l<chainLength-1; l++) {
			if ( (0 == (l+1) % lBlock) && (typBlock == typB1) ) {
				std::cout << "UpdaterCreateCrossLink::setLinearBlockCoPolyChain(): Index=" << l+1 << ", lBlock=" << lBlock <<", typBlock=" << typBlock << std::endl;
				typBlock = typB2;
			} else {
//				std::cout << "UpdaterCreateCrossLink::setLinearBlockCoPolyChain(): Index=" << l+1 << ", lBlock=" << lBlock <<", typBlock=" << typBlock << std::endl;
				typBlock = typB1;
			}
			addMonomerToParent(idStartMono+l,typBlock); 
		}
	}	

	//! OPEN TO DO IF WANTED/NEEDED
	//! 2)	two different fixed length lBlockA, lBlockB
	//!  a) repeated order
	//!  b) random order
	//!		Problem for both: How to fit in fixed chainLength
	return true;
}


#endif /* LEMONADE_UPDATER_CREATE_CROSSLINK_H */
