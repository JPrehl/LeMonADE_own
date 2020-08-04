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
class UpdaterCreateCoDendrimersInSolvent: public UpdaterAbstractCreate<IngredientsType>
{
  typedef UpdaterAbstractCreate<IngredientsType> BaseClass;
  
public:
  UpdaterCreateCoDendrimersInSolvent(IngredientsType& ingredients_, uint32_t chainLength_, uint32_t NumSolv_, double nnInteraction_, bool interactionType_, bool setWall_);
  
  virtual void initialize();
  virtual bool execute();
  virtual void cleanup();
  
  //! getter for initialised bool
  const bool& getIsInitialized() const {return isInitialized;}
  //! getter for wall bool
  const bool& getSetWall() const {return setWall;}
  //! get linear Chain Length
  const uint32_t& getChainLength() const {return Generation;}
  //! get number of solvent monomers
  const uint32_t& getNumOfsolventMonos() const {return NumOfsolventMonos;}
  //! get the box size
  const uint32_t& getBox() const {return box;}
  
private:
  using BaseClass::ingredients;
  
  using BaseClass::addMonomerToParent;
  using BaseClass::addSingleMonomer;
  using BaseClass::linearizeSystem;
  
  //! dendrimers generation
  uint32_t Generation;
  //! dendrimers spacer length
  uint32_t SpacerLength;
  //! coDendrimers grafted Chain Length
  uint32_t GraftedChainLength;
  //! dendrimers functionality
  uint32_t Functionality;
  
  //! number of solvent monomers
  uint32_t NumOfsolventMonos;
  
  //! simulation box sizes
  uint32_t box;
  //! number of Co-Dendrimers in the box
  uint32_t NumOfDendrimers;
  //! interactions between sensitive speciem and solvent
  double nnInteraction;
  //! bool to set micelle (true=1) or helmet structure (false=0)
  bool interactionType;
  //! bool to set nonperiodic boundaries in Z direction
  bool setWall;
  //! bool to check if groups are initilized
  bool isInitialized;

  //! attribute tag of dendritic monomers
  int32_t tagDendritic;
  //! attribute tag of grafted monomers
  int32_t tagGraftedChain;
};

/** 
* @brief Constructor handling the new systems paramters
*
* @param ingredients_ a reference to the IngredientsType - mainly the system
* @param G_ generation of the dendrimers
* @param S_ spacer length of the dendrimers
* @param NGC_ grafted chain length
* @param F_ functionality of the dendrimers
* @param box_ size of the cubic box
* @param NumDendr_ number of dendrimers to be added
* @param NumSolv_ number of solvent monomers to be added
* @param nnInteraction_ interaction strength passed to feature nn_interaction
* @param interactionType_ bool to change between micelle and helmet structures
* @param setWall_ bool to set a wall in z direction
*/
template < class IngredientsType >
UpdaterCreateCoDendrimersInSolvent<IngredientsType>::UpdaterCreateCoDendrimersInSolvent(IngredientsType& ingredients_, uint32_t G_, uint32_t S_, uint32_t NGC_, uint32_t F_, uint32_t box_, uint32_t NumDendr_, uint32_t NumSolv_, double nnInteraction_, bool interactionType_, bool setWall_):
  BaseClass(ingredients_), Generation(G_), SpacerLength(S_), GraftedChainLength(NGC_), Functionality(F_),
  NumOfsolventMonos(NumSolv_),
  box(box_), NumOfDendrimers(NumDendr_), 
  nnInteraction(nnInteraction_), interactionType(interactionType_),
  setWall(setWall_), isInitialized(false)
{}

/**
* The initialize function handles the new systems information.
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template < class IngredientsType >
void UpdaterCreateCoDendrimersInSolvent<IngredientsType>::initialize(){
  if(!isInitialized){
    std::cout << "initialize UpdaterCreateCoDendrimersInSolvent" << std::endl;
    
    // supress std::cout output of addBondset and synchronize
    std::streambuf *old = std::cout.rdbuf(); // <-- save
    std::stringstream ss;
    std::cout.rdbuf (ss.rdbuf());       // <-- redirect
    
    //set box size
    ingredients.setBoxX(box);
    ingredients.setBoxY(box);
    ingredients.setBoxZ(box);
    
    //set periodicity
    ingredients.setPeriodicX(true);
    ingredients.setPeriodicY(true);
    if(setWall){
      ingredients.setPeriodicZ(false);
    }else{
      ingredients.setPeriodicZ(true);
    }
    
    //set Bondset
    ingredients.modifyBondset().addBFMclassicBondset();
    //set nn interactions
    // interactions in the gpu code are between attributes 3 (solvent) and 4,5 (sensitve speciem)
    ingredients.setNNInteraction(3,4,nnInteraction);
    ingredients.setNNInteraction(3,5,nnInteraction);

    // set the attribtue tags for the different spezies
    // interaction type == true -> micelles; interaction type == false -> helmets
    if(interactionType){
      tagDendritic=1;
      tagGraftedChain=4;
    }else{
      tagDendritic=4;
      tagGraftedChain=1;
    }
    
    // call synchronize
    ingredients.synchronize(ingredients);
    // redirect output
    std::cout.rdbuf (old);
    
    std::cout << "create "<<NumOfDendrimers<<" dendrimers in the box"<<std::endl;
    
    isInitialized=true;
    execute();
  }
}

/**
* Execution of the system creation
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template < class IngredientsType >
bool UpdaterCreateCoDendrimersInSolvent<IngredientsType>::execute(){
  std::cout << "execute UpdaterCreateCoDendrimersInSolvent" << std::endl;
  
  //check if system was already created
  if(ingredients.getMolecules().size()!=0)
    return true;
  
  // #######  create the dendrimers  ########## //
  // create all central monomers
  for(uint32_t i=0; i<NumOfDendrimers;i++){
    addSingleMonomer(tagDendritic) ? : throw std::runtime_error("UpdaterCreateCoDendrimersInSolvent: addSingleMonomer is not able to place a monomer!");
  }
  
  // create the dendrimers stepwise generation by generation
  for(uint32_t generation=0;generation<Generation;generation++){
    if(generation==0){
      //add 3 spacers to every monomer
      //loop over all initial monomers
      for(uint32_t i=0; i<NumOfDendrimers; i++){
	//loop over all branches
	for(uint32_t functionality=0;functionality<Functionality;functionality++){
	  //add first monomer in new branch
	  if(ingredients.getMolecules()[i].getAttributeTag()==tagDendritic)
	    addMonomerToParent(i,tagDendritic+1) ? : (throw std::runtime_error("UpdaterCreateCoDendrimersInSolvent: addMonomerToParent not able to place monomer"));
	  else
	    addMonomerToParent(i,tagDendritic) ? : (throw std::runtime_error("UpdaterCreateCoDendrimersInSolvent: addMonomerToParent not able to place monomer"));
	  //loop over all remaining spacers
	  for(uint32_t spacer=0;spacer<(SpacerLength-1);spacer++){
	    if(ingredients.getMolecules()[ingredients.getMolecules().size()-1].getAttributeTag()==tagDendritic)
	      addMonomerToParent(ingredients.getMolecules().size()-1, tagDendritic+1) ? : (throw std::runtime_error("UpdaterCreateCoDendrimersInSolvent: addMonomerToParent not able to place monomer"));
	    else
	      addMonomerToParent(ingredients.getMolecules().size()-1, tagDendritic) ? : (throw std::runtime_error("UpdaterCreateCoDendrimersInSolvent: addMonomerToParent not able to place monomer"));
	  }
	}
      }
    }else if(generation>0){
      //add 2 spacers to every monomer with only one bond partner
      //loop over all monomers having only one neighbor
      int32_t N_old(ingredients.getMolecules().size());
      for(uint32_t i=0; i<N_old; i++){
	      if(ingredients.getMolecules().getNumLinks(i)==1){
	        //loop over all remaining branches
	        for(uint32_t functionality=0;functionality<(Functionality-1);functionality++){
	          //add first monomer in new branch
	          if(ingredients.getMolecules()[i].getAttributeTag()==tagDendritic)
	            addMonomerToParent(i,tagDendritic+1) ? : (throw std::runtime_error("UpdaterCreateCoDendrimersInSolvent: addMonomerToParent not able to place monomer"));
	          else
	            addMonomerToParent(i,tagDendritic) ? : (throw std::runtime_error("UpdaterCreateCoDendrimersInSolvent: addMonomerToParent not able to place monomer"));
	          for(uint32_t spacer=0;spacer<(SpacerLength-1);spacer++){
	            if(ingredients.getMolecules()[ingredients.getMolecules().size()-1].getAttributeTag()==tagDendritic)
	              addMonomerToParent(ingredients.getMolecules().size()-1, tagDendritic+1) ? : (throw std::runtime_error("UpdaterCreateCoDendrimersInSolvent: addMonomerToParent not able to place monomer"));
	            else
	              addMonomerToParent(ingredients.getMolecules().size()-1, tagDendritic) ? : (throw std::runtime_error("UpdaterCreateCoDendrimersInSolvent: addMonomerToParent not able to place monomer"));
	          }
	        }
	      }
      }
    }
  }

  std::cout << "adding grafted chains with "<< GraftedChainLength<<" monomers"<<std::endl;

  // #######  add the grafted chains  ########## //
  int32_t previousSystemSize(ingredients.getMolecules().size());
  // add the grafted chains to monomers with only one neighbor
  for(uint32_t i=0; i<previousSystemSize;i++){
    if(ingredients.getMolecules().getNumLinks(i)==1){
      for(uint32_t j=0;j<GraftedChainLength;j++){
        if(j==0){
          addMonomerToParent(i,tagGraftedChain) ? : (throw std::runtime_error("UpdaterCreateCoDendrimersInSolvent: addMonomerToParent not able to place monomer"));
        }else{
          if(ingredients.getMolecules()[ingredients.getMolecules().size()-1].getAttributeTag()==tagGraftedChain){
            addMonomerToParent(ingredients.getMolecules().size()-1,tagGraftedChain+1) ? : (throw std::runtime_error("UpdaterCreateCoDendrimersInSolvent: addMonomerToParent not able to place monomer"));
          }else{
            addMonomerToParent(ingredients.getMolecules().size()-1,tagGraftedChain) ? : (throw std::runtime_error("UpdaterCreateCoDendrimersInSolvent: addMonomerToParent not able to place monomer"));
          }
        }
      }
    }
  }

  std::cout << "adding "<<NumOfsolventMonos <<"solvent monomers"<<std::endl;

  // #######  add the solvent   ########## //
  previousSystemSize=ingredients.getMolecules().size();
  for(uint32_t nMono=0; nMono<NumOfsolventMonos; nMono++){
    addSingleMonomer(3) ? : throw std::runtime_error("UpdaterCreateCoDendrimersInSolvent: addSingleMonomer is not able to place a monomer!");
  }


  // #######  check the system   ########## //
  int32_t theoretical_system_size(NumOfDendrimers*((Functionality*pow(Functionality-1,Generation-1)*GraftedChainLength)+(1+Functionality*SpacerLength*( (pow((Functionality-1),Generation)-1 ) / (Functionality-2))))+NumOfsolventMonos);
  if(ingredients.getMolecules().size() != theoretical_system_size){
    std::cerr << "ingredients size = "<< ingredients.getMolecules().size() << ", calculated size = "<<theoretical_system_size<<std::endl;
    throw std::runtime_error("UpdaterCreateCoDendrimersInSolvent: number of monomers in molecules does not match the calculated number of monomers!");
  }else{
    std::cout << "real lattice occupation density =" << (8*ingredients.getMolecules().size())/(pow(box,3))<<std::endl;
    
    //linearize system
    linearizeSystem();
    // set compressed indizees
    std::cout << "set compressed output indices for chains of length 1 for idx = " << ingredients.getMolecules().size()-NumOfsolventMonos << "-" << ingredients.getMolecules().size()-1 << std::endl;
    ingredients.setCompressedOutputIndices(ingredients.getMolecules().size()-NumOfsolventMonos, ingredients.getMolecules().size()-1);
  
    return true;
  }
  
}

/**
* empty clean up.
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template < class IngredientsType >
void UpdaterCreateCoDendrimersInSolvent<IngredientsType>::cleanup(){

}


#endif /* LEMONADE_UPDATER_CREATE_CROSSLINK_H */
