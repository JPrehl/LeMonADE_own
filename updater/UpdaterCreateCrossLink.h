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
  double nnInteraction;
  
  int typAlO;
  int typPA6;
  int typBlockA;
  int typBlockB;
  
  bool isExecuted;
  
  bool setWallMonomer();
  bool setLinearCoPolyChain();
  
public:
  UpdaterCreateCrossLink(IngredientsType& ingredients_, uint32_t boxX_, uint32_t boxY_, uint32_t boxZ_, bool pX_, bool pY_, bool pZ_, uint32_t chainLength_, double nnInteraction_);
  
  virtual void initialize();
  virtual bool execute();
  virtual void cleanup();

};

/*****************************************************************************
 *  Necessary (virtual) Methods + Constructor
 * ***************************************************************************
 */

template < class IngredientsType >
UpdaterCreateCrossLink<IngredientsType>::UpdaterCreateCrossLink(IngredientsType& ingredients_, uint32_t boxX_, uint32_t boxY_, uint32_t boxZ_, bool pX_, bool pY_, bool pZ_, uint32_t chainLength_, double nnInteraction_):
  BaseClass(ingredients_), boxX(boxX_), boxY(boxY_), boxZ(boxZ_), pX(pX_), pY(pY_), pZ(pZ_), chainLength(chainLength_), nnInteraction(nnInteraction_), isExecuted(false)
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

	std::cout << "Updater CL: cleanup" << std::endl;
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

	addSingleMonomer(typBlockA);	// Start-Monomer der Chain
	int idStartMono = ingredients.getMolecules().size()-1; 
	for (int l=0; l<chainLength-1; l++) {
		
		// Hier kommt die Komposition der Kette rein
		addMonomerToParent(idStartMono+l,typBlockA); 
	}
	
	return true;
}

#endif /* LEMONADE_UPDATER_CREATE_CROSSLINK_H */
