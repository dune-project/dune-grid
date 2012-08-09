// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <config.h>

#include <set>

#include <dune/grid/uggrid.hh>

/** \todo Remove the following include once getAllSubfaces... is gone */
#include <dune/common/container/sllist.hh>
#include <dune/common/stdstreams.hh>


using namespace Dune;

//***********************************************************************
//
// --UGGrid
// --Grid
//
//***********************************************************************

template<> int Dune::UGGrid<2>::numOfUGGrids = 0;
template<> int Dune::UGGrid<3>::numOfUGGrids = 0;

template<> unsigned int Dune::UGGrid<2>::heapSize_ = 500;
template<> unsigned int Dune::UGGrid<3>::heapSize_ = 500;


template <int dim>
Dune::UGGrid < dim >::UGGrid()
  : multigrid_(NULL),
    leafIndexSet_(*this),
    globalIdSet_(*this),
    localIdSet_(*this),
    refinementType_(LOCAL),
    closureType_(GREEN),
    someElementHasBeenMarkedForRefinement_(false),
    someElementHasBeenMarkedForCoarsening_(false),
    numBoundarySegments_(0)
{
  // If no UGGrid object exists yet start up UG for 2d and 3d
  if ( (UGGrid<2>::numOfUGGrids+UGGrid<3>::numOfUGGrids)==0) {

    /* Why do we start UG for 2d and 3d?  Why not just the one for the grid dimension
     * that we actually need?  The reason is a certain error condition that you get
     * when using 2d and 3d UGGrids together. This was FlySpray task 887.
     *
     * UG has an internal directory-tree-like structure for storing different types of
     * miscellaneous data.  It is called the 'environment heap'. The method InitUg
     * creates a directory 'BVP' there, which holds, among other things, a 'problem name'
     * (a string) for each existing grid.  When you call InitUg twice (i.e., once for 2d
     * and once for 3d), the second call will not notice that a directory 'BVP' already exists.
     * This is possible because these directories have a name and a 'type', (an integer
     * of which I have not really understood what you need it for), and the new 'BVP' directory
     * has a different type. The new one will appear in the list before the old one.
     * However, some methods that look for entries in the environment heap only look for the name,
     * and hence only find the first 'BVP' entry, which only contains the information about
     * the second grid.  Surprisingly, this doesn't appear to be a problem, at least the dune-grid
     * unit test for UG successfully tests a 2d and a 3d grid side by side.  However when trying
     * to delete grids the method DisposeMultiGrid tries to remove a grid's information from the
     * BVP directory.  However the information of the first grid is in the second BVP directory,
     * which is effectively hidden by the first one.  DisposeMultiGrid doesn't find the
     * grid information and reports an error.
     *
     * If both calls to InitUg are directly in a row, then all information will be stored
     * in the first instance of 'BVP', and there is no error.
     */
    int argc = 1;
    char* arg = strdup("dune.exe");
    char** argv = &arg;

    if (UG_NS<2>::InitUg(&argc, &argv))
      DUNE_THROW(GridError, "UG" << dim << "d::InitUg() returned an error code!");

    if (UG_NS<3>::InitUg(&argc, &argv))
      DUNE_THROW(GridError, "UG" << dim << "d::InitUg() returned an error code!");

    free(arg);
  }

  // Create a dummy problem
  typename UG_NS<dim>::CoeffProcPtr coeffs[1] = {NULL};
  typename UG_NS<dim>::UserProcPtr upp[1] = {NULL};

  // Create unique problem name
  std::stringstream numberAsAscii;
  numberAsAscii << numOfUGGrids;
  name_ = "DuneUGGrid_" + std::string((dim==2) ? "2" : "3") + std::string("d_") + numberAsAscii.str();

  std::string problemName = name_ + "_Problem";

  if (UG_NS<dim>::CreateBoundaryValueProblem(problemName.c_str(), 1,coeffs,1,upp) == NULL)
    DUNE_THROW(GridError, "UG" << dim << "d::CreateBoundaryValueProblem() returned an error code!");

  if (numOfUGGrids==0) {

    if (dim==2)
    {
      char* nfarg = strdup("newformat DuneFormat2d");
      if (UG_NS<dim>::CreateFormatCmd(1, &nfarg))
        DUNE_THROW(GridError, "UG" << dim << "d::CreateFormat() returned an error code!");
      free(nfarg);
    }
    if (dim==3)
    {
      char* newArgs[2];
      for (int i=0; i<2; i++)
        newArgs[i] = (char*)::malloc(50*sizeof(char));

      sprintf(newArgs[0], "newformat DuneFormat3d" );
      sprintf(newArgs[1], "V s1 : vt 1" ); // generates side vectors in 3D

      if (UG_NS<dim>::CreateFormatCmd(2, newArgs))
        DUNE_THROW(GridError, "UG" << dim << "d::CreateFormat() returned an error code!");

      for (int i=0; i<2; i++)
        free(newArgs[i]);

    }
  }

  numOfUGGrids++;

  dverb << "UGGrid<" << dim << "> with name " << name_ << " created!" << std::endl;

}


template < int dim >
Dune::UGGrid < dim >::~UGGrid()
{
  // Delete the UG multigrid if there is one (== createEnd() has been called)
  if (multigrid_) {
    // Set UG's currBVP variable to the BVP corresponding to this
    // grid.  This is necessary if we have more than one UGGrid in use.
    // DisposeMultiGrid will crash if we don't do this
    UG_NS<dim>::Set_Current_BVP(multigrid_->theBVP);
    if (UG_NS<dim>::DisposeMultiGrid(multigrid_) != 0)
      DUNE_THROW(GridError, "UG" << dim << "d::DisposeMultiGrid returned error code!");
  }

  // DisposeMultiGrid cleans up the BVP as well.  But if there was no
  // multigrid we have to take care of the BVP ourselves.
  std::string problemName = name_ + "_Problem";
  void** BVP = UG_NS<dim>::BVP_GetByName(problemName.c_str());

  if (BVP)
    if (UG_NS<dim>::BVP_Dispose(BVP))
      DUNE_THROW(GridError, "Couldn't dispose of UG boundary value problem!");



  numOfUGGrids--;

  // Shut down UG if this was the last existing UGGrid object
  if (UGGrid<2>::numOfUGGrids + UGGrid<3>::numOfUGGrids == 0)
    UG_NS<dim>::ExitUg();

  // Delete levelIndexSets
  for (unsigned int i=0; i<levelIndexSets_.size(); i++)
    if (levelIndexSets_[i])
      delete levelIndexSets_[i];
}

template < int dim >
int Dune::UGGrid < dim >::maxLevel() const
{
  if (!multigrid_)
    DUNE_THROW(GridError, "The grid has not been properly initialized!");

  return multigrid_->topLevel;
}



template<int dim>
template<int codim>
typename Dune::UGGrid<dim>::Traits::template Codim<codim>::LevelIterator
Dune::UGGrid<dim>::lbegin (int level) const
{
  if (!multigrid_)
    DUNE_THROW(GridError, "The grid has not been properly initialized!");

  if (!multigrid_->grids[level])
    DUNE_THROW(GridError, "LevelIterator in nonexisting level " << level << " requested!");

  return UGGridLevelIterator<codim, All_Partition, const UGGrid<dim> >(*this,level);
}

template<int dim>
template<int codim, Dune::PartitionIteratorType PiType>
typename Dune::UGGrid<dim>::Traits::template Codim<codim>::template Partition<PiType>::LevelIterator
Dune::UGGrid<dim>::lbegin (int level) const
{
  if (!multigrid_)
    DUNE_THROW(GridError, "The grid has not been properly initialized!");

  if (!multigrid_->grids[level])
    DUNE_THROW(GridError, "LevelIterator in nonexisting level " << level << " requested!");

  return UGGridLevelIterator<codim, PiType, const UGGrid<dim> >(*this,level);
}

template < int dim >
template<int codim>
typename Dune::UGGrid<dim>::Traits::template Codim<codim>::LevelIterator
Dune::UGGrid < dim >::lend (int level) const
{
  return UGGridLevelIterator<codim,All_Partition, const UGGrid<dim> >();
}

template < int dim >
template<int codim, Dune::PartitionIteratorType PiType>
typename Dune::UGGrid<dim>::Traits::template Codim<codim>::template Partition<PiType>::LevelIterator
Dune::UGGrid < dim >::lend (int level) const
{
  return UGGridLevelIterator<codim,PiType, const UGGrid<dim> >();
}

template < int dim >
int Dune::UGGrid < dim >::size (int level, int codim) const
{
  if(codim == 0)
  {
    return levelIndexSet(level).size(GeometryType(GeometryType::simplex,dim))
           + levelIndexSet(level).size(GeometryType(GeometryType::cube,dim))
           + levelIndexSet(level).size(GeometryType(GeometryType::pyramid,dim))
           + levelIndexSet(level).size(GeometryType(GeometryType::prism,dim));
  }
  if(codim == dim)
  {
    return this->levelIndexSet(level).size(GeometryType(0));
  }
  if (codim == dim-1)
  {
    return this->levelIndexSet(level).size(GeometryType(1));
  }
  if (codim == 1)
  {
    return levelIndexSet(level).size(GeometryType(GeometryType::cube,dim-1))
           + levelIndexSet(level).size(GeometryType(GeometryType::simplex,dim-1));
  }
  DUNE_THROW(GridError, "UGGrid<" << dim << ">::size(int level, int codim) is only implemented"
                                  << " for codim==0 and codim==dim!");
}


template < int dim >
bool Dune::UGGrid < dim >::mark(int refCount,
                                const typename Traits::template Codim<0>::Entity& e )
{
  typename UG_NS<dim>::Element* target = this->getRealImplementation(e).target_;

  // No refinement requested
  if (refCount==0) {
    if (UG_NS<dim>::MarkForRefinement(target,
                                      UG_NS<dim>::NO_REFINEMENT,      // unset
                                      0)      // Irrelevant if refinement rule is not BLUE
        ) DUNE_THROW(GridError, "UG" << dim << "d::MarkForRefinement returned error code!");
    return true;
  }

  // Check whether element can be marked for refinement
  if (!EstimateHere(target))
    return false;

  if (refCount==1) {
    if (UG_NS<dim>::MarkForRefinement(target,
                                      UG_NS<dim>::RED,      // red refinement rule
                                      0)      // Irrelevant if refinement rule is not BLUE
        ) DUNE_THROW(GridError, "UG" << dim << "d::MarkForRefinement returned error code!");

    someElementHasBeenMarkedForRefinement_ = true;
    return true;
  } else if (refCount==-1) {

    if (UG_NS<dim>::MarkForRefinement(target,
                                      UG_NS<dim>::COARSE,      // coarsen the element
                                      0)      // Irrelevant if refinement rule is not BLUE
        ) DUNE_THROW(GridError, "UG" << dim << "d::MarkForRefinement returned error code!");

    someElementHasBeenMarkedForCoarsening_ = true;
    return true;
  } else
    DUNE_THROW(GridError, "UGGrid only supports refCount values -1, 0, and 1 for mark()!");

}

template < int dim >
bool Dune::UGGrid < dim >::mark(const typename Traits::template Codim<0>::Entity& e,
                                typename UG_NS<dim>::RefinementRule rule,
                                int side)
{
  typename UG_NS<dim>::Element* target = this->getRealImplementation(e).target_;

  if (!UG_NS<dim>::isLeaf(target))
    return false;

  someElementHasBeenMarkedForRefinement_ = true;

  return UG_NS<dim>::MarkForRefinement(target, rule, side);

}

template <int dim>
int Dune::UGGrid<dim>::getMark(const typename Traits::template Codim<0>::Entity& e) const
{
  typename UG_NS<dim>::Element* target = this->getRealImplementation(e).target_;

  // Return -1 if element is marked for coarsening
  if (UG_NS<dim>::ReadCW(target,UG_NS<dim>::COARSEN_CE))
    return -1;

  // If the element is not regular, it's mark is actually stored on its regular
  // ancestor.  That's what we look for here.
  if (dim==2)
    target = (typename UG_NS<dim>::Element*) UG::D2::ELEMENT_TO_MARK((UG::D2::element*)target);
  else
    target = (typename UG_NS<dim>::Element*) UG::D3::ELEMENT_TO_MARK((UG::D3::element*)target);

  // Return 0 if element is not marked at all
  if (UG_NS<dim>::ReadCW(target,UG_NS<dim>::MARK_CE)==UG_NS<dim>::NO_REFINEMENT)
    return 0;

  // Else return 1
  return 1;
}

template <int dim>
bool Dune::UGGrid <dim>::preAdapt()
{
  if( closureType_ == GREEN )
  {
    // when conforming refinement is enabled
    // green closure has to be removed although not
    // marked for coarsening
    return someElementHasBeenMarkedForCoarsening_ ||
           someElementHasBeenMarkedForRefinement_ ;
  }
  else
  {
    // in non-conforming meshes only elements marked for coarsening
    // will be coarsened (hopefully)
    return someElementHasBeenMarkedForCoarsening_;
  }
}

template < int dim >
bool Dune::UGGrid < dim >::adapt()
{
  assert(multigrid_);

  // Set UG's currBVP variable to the BVP corresponding to this
  // grid.  This is necessary if we have more than one UGGrid in use.
  UG_NS<dim>::Set_Current_BVP(multigrid_->theBVP);

  int mode = UG_NS<dim>::GM_REFINE_TRULY_LOCAL;

  if (refinementType_==COPY)
    mode = mode | UG_NS<dim>::GM_COPY_ALL;

  if (closureType_==NONE)
    mode = mode | UG_NS<dim>::GM_REFINE_NOT_CLOSED;

  // I don't really know what this means
  int seq = UG_NS<dim>::GM_REFINE_PARALLEL;

  // I don't really know what this means either
  int mgtest = UG_NS<dim>::GM_REFINE_NOHEAPTEST;

  int rv = AdaptMultiGrid(multigrid_,mode,seq,mgtest);

  if (rv!=0)
    DUNE_THROW(GridError, "UG::adapt() returned with error code " << rv);

  // Renumber everything
  setIndices(false, NULL);

  // Return true iff the grid hierarchy changed
  //return !(bool)multigrid_->status;

  // grid has changed if at least one element was marked for refinement
  return someElementHasBeenMarkedForRefinement_;
}

template <int dim>
void Dune::UGGrid <dim>::postAdapt()
{
  for (int i=0; i<=maxLevel(); i++) {

    typename Traits::template Codim<0>::LevelIterator eIt    = lbegin<0>(i);
    typename Traits::template Codim<0>::LevelIterator eEndIt = lend<0>(i);

    for (; eIt!=eEndIt; ++eIt)
      UG_NS<dim>::WriteCW(this->getRealImplementation(*eIt).target_, UG_NS<dim>::NEWEL_CE, 0);

  }

  // reset marker flags
  someElementHasBeenMarkedForRefinement_ = false;
  someElementHasBeenMarkedForCoarsening_ = false;
}

template < int dim >
void Dune::UGGrid < dim >::globalRefine(int n)
{
  for (int i=0; i<n; i++) {

    // mark all entities for grid refinement
    typename Traits::template Codim<0>::LeafIterator iIt    = leafbegin<0>();
    typename Traits::template Codim<0>::LeafIterator iEndIt = leafend<0>();

    for (; iIt!=iEndIt; ++iIt)
      mark(1, *iIt);

    this->preAdapt();
    adapt();

  }

  this->postAdapt();

}

template <int dim>
void Dune::UGGrid<dim>::getChildrenOfSubface(const typename Traits::template Codim<0>::EntityPointer & e,
                                             int elementSide,
                                             int maxl,
                                             std::vector<typename Traits::template Codim<0>::EntityPointer>& childElements,
                                             std::vector<unsigned char>& childElementSides) const
{

  typedef std::pair<typename UG_NS<dim>::Element*,int> ListEntryType;

  SLList<ListEntryType> list;

  // //////////////////////////////////////////////////////////////////////
  //   Change the input face number from Dune numbering to UG numbering
  // //////////////////////////////////////////////////////////////////////

  elementSide = UGGridRenumberer<dim>::facesDUNEtoUG(elementSide, e->type());

  // ///////////////
  //   init list
  // ///////////////
  if (!e->isLeaf()   // Get_Sons_of_ElementSide returns GM_FATAL when called for a leaf !?!
      && e->level() < maxl) {

    typename UG_NS<dim>::Element* theElement = this->getRealImplementation(*e).target_;

    int Sons_of_Side = 0;
    typename UG_NS<dim>::Element* SonList[UG_NS<dim>::MAX_SONS];
    int SonSides[UG_NS<dim>::MAX_SONS];

    int rv = Get_Sons_of_ElementSide(theElement,
                                     elementSide,
                                     &Sons_of_Side,
                                     SonList,          // the output elements
                                     SonSides,         // Output element side numbers
                                     true,            // Element sons are not precomputed
                                     true);            // ioflag: I have no idea what this is supposed to do
    if (rv!=0)
      DUNE_THROW(GridError, "Get_Sons_of_ElementSide returned with error value " << rv);

    for (int i=0; i<Sons_of_Side; i++)
      list.push_back(ListEntryType(SonList[i],SonSides[i]));

  }

  // //////////////////////////////////////////////////
  //   Traverse and collect all children of the side
  // //////////////////////////////////////////////////

  typename SLList<ListEntryType>::iterator f = list.begin();
  for (; f!=list.end(); ++f) {

    typename UG_NS<dim>::Element* theElement = f->first;
    int side                                  = f->second;

    int Sons_of_Side = 0;
    typename UG_NS<dim>::Element* SonList[UG_NS<dim>::MAX_SONS];
    int SonSides[UG_NS<dim>::MAX_SONS];

    if (UG_NS<dim>::myLevel(theElement) < maxl) {

      Get_Sons_of_ElementSide(theElement,
                              side,         // Input element side number
                              &Sons_of_Side,       // Number of topological sons of the element side
                              SonList,            // Output elements
                              SonSides,           // Output element side numbers
                              true,
                              true);

      for (int i=0; i<Sons_of_Side; i++)
        list.push_back(ListEntryType(SonList[i],SonSides[i]));

    }

  }

  // //////////////////////////////
  //   Extract result from stack
  // //////////////////////////////

  // Use reserve / push_back since EntityPointer is not default constructable
  childElements.clear();
  childElements.reserve( list.size() );
  childElementSides.resize(list.size());

  int i=0;
  for (f = list.begin(); f!=list.end(); ++f, ++i)
  {

    // Set element
    typedef typename Traits::template Codim< 0 >::EntityPointer EntityPointer;
    childElements.push_back( EntityPointer( UGGridEntityPointer< 0, const UGGrid< dim > >( f->first, this ) ) );

    int side = f->second;

    // Dune numbers the faces of several elements differently than UG.
    // The following switch does the transformation
    childElementSides[i] = UGGridRenumberer<dim>::facesUGtoDUNE(side, childElements[i]->type());

  }

}

template < int dim >
bool Dune::UGGrid < dim >::loadBalance(int strategy, int minlevel, int depth, int maxLevel, int minelement)
{
  /** \todo Test for valid arguments */
  std::string argStrings[4];
  std::stringstream numberAsAscii[4];

  numberAsAscii[0] << strategy;
  argStrings[0] = "lb " + numberAsAscii[0].str();

  numberAsAscii[1] << minlevel;
  argStrings[1] = "c " + numberAsAscii[1].str();

  numberAsAscii[2] << depth;
  argStrings[2] = "d " + numberAsAscii[2].str();

  numberAsAscii[3] << minelement;
  argStrings[3] = "e " + numberAsAscii[3].str();

  const char* argv[4] = {argStrings[0].c_str(),
                         argStrings[1].c_str(),
                         argStrings[2].c_str(),
                         argStrings[3].c_str()};

  int errCode = UG_NS<dim>::LBCommand(4, argv);

  if (errCode)
    DUNE_THROW(GridError, "UG" << dim << "d::LBCommand returned error code " << errCode);

  // Renumber everything.
  setIndices(true, NULL);

  return true;
}


template < int dim >
void Dune::UGGrid < dim >::setPosition(const typename Traits::template Codim<dim>::EntityPointer& e,
                                       const FieldVector<double, dim>& pos)
{
  typename UG_NS<dim>::Node* target = this->getRealImplementation(*e).target_;

  for (int i=0; i<dim; i++)
    target->myvertex->iv.x[i] = pos[i];
}

template <int dim>
void Dune::UGGrid<dim>::saveState(const std::string& filename) const
{
  const char* type = "asc";
  const char* comment = "written by DUNE";

  if (dim==2)
    UG::D2::SaveMultiGrid((UG::D2::multigrid*)multigrid_,
                          filename.c_str(),
                          type,
                          comment,
                          0,      // autosave
                          0       // rename
                          );
  else
    UG::D3::SaveMultiGrid((UG::D3::multigrid*)multigrid_,
                          filename.c_str(),
                          type,
                          comment,
                          0,      // autosave
                          0       // rename
                          );
}


template <int dim>
void Dune::UGGrid<dim>::loadState(const std::string& filename)
{
  const char* type = "asc";
  std::string problemName = name_ + "_Problem";
  std::string formatName = "DuneFormat2d";

  if (dim==2) {
    std::string formatName = "DuneFormat2d";
    multigrid_ = (typename UG_NS<dim>::MultiGrid*) UG::D2::LoadMultiGrid(name_.c_str(),
                                                                         filename.c_str(),
                                                                         type,
                                                                         problemName.c_str(),
                                                                         formatName.c_str(),
                                                                         heapSize_,
                                                                         true, //force,
                                                                         true, //optimizedIO,
                                                                         false //autosave
                                                                         );
  } else {
    std::string formatName = "DuneFormat3d";
    multigrid_ = (typename UG_NS<dim>::MultiGrid*) UG::D3::LoadMultiGrid(name_.c_str(),
                                                                         filename.c_str(),
                                                                         type,
                                                                         problemName.c_str(),
                                                                         formatName.c_str(),
                                                                         heapSize_,
                                                                         true, //force,
                                                                         true, //optimizedIO,
                                                                         false //autosave
                                                                         );
  }

  if (multigrid_==NULL)
    DUNE_THROW(GridError, "In loadState()");
}

template < int dim >
void Dune::UGGrid < dim >::setIndices(bool setLevelZero,
                                      std::vector<unsigned int>* nodePermutation)
{
  /** \todo The code in the following if-clause contains two workarounds to fix FlySpray issues 170
      (inconsistent codim 1 subIndices) and 810 (Destructor of parallel UGGrid fails).

      UGGrid uses UG's SideVector data structure to store the codim 1 subIndices of elements.
      However the UG SideVector code is buggy and SideVectors are not correctly created during
      refinement and load balancing.  In particular
        - after refinement there may be _two_ SideVector objects for a single face,
          instead of only one,
        - After load balancing, and also after refinement, the VCOUNT field of the SideVector
          (which stores how many elements reference the SideVector) is sometimes not correct.

      I provide a workaround because it is outside of my capabilities to actually fix these
      issues in UG itself.  I'd probably screw up other things if I started to mess around
      in the UG refinement routines.

      If ever UGGrid stops relying on SideVectors to store indices the following
      if-clause can be deleted.
   */
  if (dim==3) {

    for (int i=0; i<=maxLevel(); i++) {

      typename Traits::template Codim<0>::LevelIterator eIt = lbegin<0>(i);
      typename Traits::template Codim<0>::LevelIterator eEndIt = lend<0>(i);

      for (; eIt!=eEndIt; ++eIt) {

        typename UG_NS<dim>::Element* elem0 = this->getRealImplementation(*eIt).target_;

        typename Traits::template Codim<0>::Entity::LevelIntersectionIterator nIt = eIt->ilevelbegin();
        typename Traits::template Codim<0>::Entity::LevelIntersectionIterator nEndIt = eIt->ilevelend();

        for (; nIt!=nEndIt; ++nIt) {

          int side0 = UGGridRenumberer<dim>::facesDUNEtoUG(nIt->indexInInside(), eIt->type());

          if (nIt->neighbor()) {

            typename UG_NS<dim>::Element* elem1 = this->getRealImplementation(*nIt->outside()).target_;

            int side1 = UGGridRenumberer<dim>::facesDUNEtoUG(nIt->indexInOutside(), nIt->outside()->type());

            // If there are two SideVector objects instead of only one (as there should be),
            // delete one of them.
            // Isn't it great that UG even provides a dedicated method for this?
            UG::D3::DisposeDoubledSideVector((typename UG_NS<3>::Grid*)multigrid_->grids[i],
                                             (typename UG_NS<3>::Element*)elem0,
                                             side0,
                                             (typename UG_NS<3>::Element*)elem1,
                                             side1);

          }

          // Set the correct value for the VCOUNT field:
          // the number of elements adjacent to this face.
          // This method may not be called before DisposeDoubledSideVector
          UG_NS<dim>::setVCount(elem0,side0, (nIt->neighbor() ? 2 : 1));

        }

      }

    }

  }

  // Create new level index sets if necessary
  for (int i=levelIndexSets_.size(); i<=maxLevel(); i++)
    levelIndexSets_.push_back(new UGGridLevelIndexSet<const UGGrid<dim> >());

  // Update the zero level LevelIndexSet.  It is updated only once, at the time
  // of creation of the coarse grid.  After that it is not touched anymore.
  if (setLevelZero)
    levelIndexSets_[0]->update(*this, 0, nodePermutation);

  // Update the remaining level index sets
  for (int i=1; i<=maxLevel(); i++)
    if (levelIndexSets_[i])
      levelIndexSets_[i]->update(*this, i);

  leafIndexSet_.update(nodePermutation);

  // id sets don't need updating
}

// /////////////////////////////////////////////////////////////////////////////////
//   Explicit instantiation of the dimensions that are actually supported by UG.
//   g++-4.0 wants them to be _after_ the method implementations.
// /////////////////////////////////////////////////////////////////////////////////

template class Dune::UGGrid<2>;
template class Dune::UGGrid<3>;

// ////////////////////////////////////////////////////////////////////////////////////
//   Explicitly instantiate the necessary member templates contained in UGGrid<2>
// ////////////////////////////////////////////////////////////////////////////////////
template Dune::UGGrid<2>::Traits::Codim<0>::LevelIterator Dune::UGGrid<2>::lbegin<0>(int level) const;
template Dune::UGGrid<2>::Traits::Codim<2>::LevelIterator Dune::UGGrid<2>::lbegin<2>(int level) const;

template Dune::UGGrid<2>::Traits::Codim<0>::LevelIterator Dune::UGGrid<2>::lend<0>(int level) const;
template Dune::UGGrid<2>::Traits::Codim<2>::LevelIterator Dune::UGGrid<2>::lend<2>(int level) const;

template Dune::UGGrid<2>::Traits::Codim<0>::LeafIterator Dune::UGGrid<2>::leafbegin<0>() const;
template Dune::UGGrid<2>::Traits::Codim<2>::LeafIterator Dune::UGGrid<2>::leafbegin<2>() const;

template Dune::UGGrid<2>::Traits::Codim<0>::LeafIterator Dune::UGGrid<2>::leafend<0>() const;
template Dune::UGGrid<2>::Traits::Codim<2>::LeafIterator Dune::UGGrid<2>::leafend<2>() const;


// Element level iterators
template Dune::UGGrid<2>::Traits::Codim<0>::Partition<Dune::Interior_Partition>::LevelIterator
Dune::UGGrid<2>::lbegin<0,Dune::Interior_Partition>(int level) const;
template Dune::UGGrid<2>::Traits::Codim<0>::Partition<Dune::InteriorBorder_Partition>::LevelIterator
Dune::UGGrid<2>::lbegin<0,Dune::InteriorBorder_Partition>(int level) const;
template Dune::UGGrid<2>::Traits::Codim<0>::Partition<Dune::Overlap_Partition>::LevelIterator
Dune::UGGrid<2>::lbegin<0,Dune::Overlap_Partition>(int level) const;
template Dune::UGGrid<2>::Traits::Codim<0>::Partition<Dune::OverlapFront_Partition>::LevelIterator
Dune::UGGrid<2>::lbegin<0,Dune::OverlapFront_Partition>(int level) const;
template Dune::UGGrid<2>::Traits::Codim<0>::Partition<Dune::All_Partition>::LevelIterator
Dune::UGGrid<2>::lbegin<0,Dune::All_Partition>(int level) const;
template Dune::UGGrid<2>::Traits::Codim<0>::Partition<Dune::Ghost_Partition>::LevelIterator
Dune::UGGrid<2>::lbegin<0,Dune::Ghost_Partition>(int level) const;

template Dune::UGGrid<2>::Traits::Codim<0>::Partition<Dune::Interior_Partition>::LevelIterator
Dune::UGGrid<2>::lend<0,Dune::Interior_Partition>(int level) const;
template Dune::UGGrid<2>::Traits::Codim<0>::Partition<Dune::InteriorBorder_Partition>::LevelIterator
Dune::UGGrid<2>::lend<0,Dune::InteriorBorder_Partition>(int level) const;
template Dune::UGGrid<2>::Traits::Codim<0>::Partition<Dune::Overlap_Partition>::LevelIterator
Dune::UGGrid<2>::lend<0,Dune::Overlap_Partition>(int level) const;
template Dune::UGGrid<2>::Traits::Codim<0>::Partition<Dune::OverlapFront_Partition>::LevelIterator
Dune::UGGrid<2>::lend<0,Dune::OverlapFront_Partition>(int level) const;
template Dune::UGGrid<2>::Traits::Codim<0>::Partition<Dune::All_Partition>::LevelIterator
Dune::UGGrid<2>::lend<0,Dune::All_Partition>(int level) const;
template Dune::UGGrid<2>::Traits::Codim<0>::Partition<Dune::Ghost_Partition>::LevelIterator
Dune::UGGrid<2>::lend<0,Dune::Ghost_Partition>(int level) const;

// Vertex level iterators
template Dune::UGGrid<2>::Traits::Codim<2>::Partition<Dune::Interior_Partition>::LevelIterator
Dune::UGGrid<2>::lbegin<2,Dune::Interior_Partition>(int level) const;
template Dune::UGGrid<2>::Traits::Codim<2>::Partition<Dune::InteriorBorder_Partition>::LevelIterator
Dune::UGGrid<2>::lbegin<2,Dune::InteriorBorder_Partition>(int level) const;
template Dune::UGGrid<2>::Traits::Codim<2>::Partition<Dune::Overlap_Partition>::LevelIterator
Dune::UGGrid<2>::lbegin<2,Dune::Overlap_Partition>(int level) const;
template Dune::UGGrid<2>::Traits::Codim<2>::Partition<Dune::OverlapFront_Partition>::LevelIterator
Dune::UGGrid<2>::lbegin<2,Dune::OverlapFront_Partition>(int level) const;
template Dune::UGGrid<2>::Traits::Codim<2>::Partition<Dune::All_Partition>::LevelIterator
Dune::UGGrid<2>::lbegin<2,Dune::All_Partition>(int level) const;
template Dune::UGGrid<2>::Traits::Codim<2>::Partition<Dune::Ghost_Partition>::LevelIterator
Dune::UGGrid<2>::lbegin<2,Dune::Ghost_Partition>(int level) const;

template Dune::UGGrid<2>::Traits::Codim<2>::Partition<Dune::Interior_Partition>::LevelIterator
Dune::UGGrid<2>::lend<2,Dune::Interior_Partition>(int level) const;
template Dune::UGGrid<2>::Traits::Codim<2>::Partition<Dune::InteriorBorder_Partition>::LevelIterator
Dune::UGGrid<2>::lend<2,Dune::InteriorBorder_Partition>(int level) const;
template Dune::UGGrid<2>::Traits::Codim<2>::Partition<Dune::Overlap_Partition>::LevelIterator
Dune::UGGrid<2>::lend<2,Dune::Overlap_Partition>(int level) const;
template Dune::UGGrid<2>::Traits::Codim<2>::Partition<Dune::OverlapFront_Partition>::LevelIterator
Dune::UGGrid<2>::lend<2,Dune::OverlapFront_Partition>(int level) const;
template Dune::UGGrid<2>::Traits::Codim<2>::Partition<Dune::All_Partition>::LevelIterator
Dune::UGGrid<2>::lend<2,Dune::All_Partition>(int level) const;
template Dune::UGGrid<2>::Traits::Codim<2>::Partition<Dune::Ghost_Partition>::LevelIterator
Dune::UGGrid<2>::lend<2,Dune::Ghost_Partition>(int level) const;

// Element leaf iterators
template Dune::UGGrid<2>::Traits::Codim<0>::Partition<Dune::Interior_Partition>::LeafIterator
Dune::UGGrid<2>::leafbegin<0,Dune::Interior_Partition>() const;
template Dune::UGGrid<2>::Traits::Codim<0>::Partition<Dune::InteriorBorder_Partition>::LeafIterator
Dune::UGGrid<2>::leafbegin<0,Dune::InteriorBorder_Partition>() const;
template Dune::UGGrid<2>::Traits::Codim<0>::Partition<Dune::Overlap_Partition>::LeafIterator
Dune::UGGrid<2>::leafbegin<0,Dune::Overlap_Partition>() const;
template Dune::UGGrid<2>::Traits::Codim<0>::Partition<Dune::OverlapFront_Partition>::LeafIterator
Dune::UGGrid<2>::leafbegin<0,Dune::OverlapFront_Partition>() const;
template Dune::UGGrid<2>::Traits::Codim<0>::Partition<Dune::All_Partition>::LeafIterator
Dune::UGGrid<2>::leafbegin<0,Dune::All_Partition>() const;
template Dune::UGGrid<2>::Traits::Codim<0>::Partition<Dune::Ghost_Partition>::LeafIterator
Dune::UGGrid<2>::leafbegin<0,Dune::Ghost_Partition>() const;

template Dune::UGGrid<2>::Traits::Codim<0>::Partition<Dune::Interior_Partition>::LeafIterator
Dune::UGGrid<2>::leafend<0,Dune::Interior_Partition>() const;
template Dune::UGGrid<2>::Traits::Codim<0>::Partition<Dune::InteriorBorder_Partition>::LeafIterator
Dune::UGGrid<2>::leafend<0,Dune::InteriorBorder_Partition>() const;
template Dune::UGGrid<2>::Traits::Codim<0>::Partition<Dune::Overlap_Partition>::LeafIterator
Dune::UGGrid<2>::leafend<0,Dune::Overlap_Partition>() const;
template Dune::UGGrid<2>::Traits::Codim<0>::Partition<Dune::OverlapFront_Partition>::LeafIterator
Dune::UGGrid<2>::leafend<0,Dune::OverlapFront_Partition>() const;
template Dune::UGGrid<2>::Traits::Codim<0>::Partition<Dune::All_Partition>::LeafIterator
Dune::UGGrid<2>::leafend<0,Dune::All_Partition>() const;
template Dune::UGGrid<2>::Traits::Codim<0>::Partition<Dune::Ghost_Partition>::LeafIterator
Dune::UGGrid<2>::leafend<0,Dune::Ghost_Partition>() const;

// Vertex leaf iterators
template Dune::UGGrid<2>::Traits::Codim<2>::Partition<Dune::Interior_Partition>::LeafIterator
Dune::UGGrid<2>::leafbegin<2,Dune::Interior_Partition>() const;
template Dune::UGGrid<2>::Traits::Codim<2>::Partition<Dune::InteriorBorder_Partition>::LeafIterator
Dune::UGGrid<2>::leafbegin<2,Dune::InteriorBorder_Partition>() const;
template Dune::UGGrid<2>::Traits::Codim<2>::Partition<Dune::Overlap_Partition>::LeafIterator
Dune::UGGrid<2>::leafbegin<2,Dune::Overlap_Partition>() const;
template Dune::UGGrid<2>::Traits::Codim<2>::Partition<Dune::OverlapFront_Partition>::LeafIterator
Dune::UGGrid<2>::leafbegin<2,Dune::OverlapFront_Partition>() const;
template Dune::UGGrid<2>::Traits::Codim<2>::Partition<Dune::All_Partition>::LeafIterator
Dune::UGGrid<2>::leafbegin<2,Dune::All_Partition>() const;
template Dune::UGGrid<2>::Traits::Codim<2>::Partition<Dune::Ghost_Partition>::LeafIterator
Dune::UGGrid<2>::leafbegin<2,Dune::Ghost_Partition>() const;

template Dune::UGGrid<2>::Traits::Codim<2>::Partition<Dune::Interior_Partition>::LeafIterator
Dune::UGGrid<2>::leafend<2,Dune::Interior_Partition>() const;
template Dune::UGGrid<2>::Traits::Codim<2>::Partition<Dune::InteriorBorder_Partition>::LeafIterator
Dune::UGGrid<2>::leafend<2,Dune::InteriorBorder_Partition>() const;
template Dune::UGGrid<2>::Traits::Codim<2>::Partition<Dune::Overlap_Partition>::LeafIterator
Dune::UGGrid<2>::leafend<2,Dune::Overlap_Partition>() const;
template Dune::UGGrid<2>::Traits::Codim<2>::Partition<Dune::OverlapFront_Partition>::LeafIterator
Dune::UGGrid<2>::leafend<2,Dune::OverlapFront_Partition>() const;
template Dune::UGGrid<2>::Traits::Codim<2>::Partition<Dune::All_Partition>::LeafIterator
Dune::UGGrid<2>::leafend<2,Dune::All_Partition>() const;
template Dune::UGGrid<2>::Traits::Codim<2>::Partition<Dune::Ghost_Partition>::LeafIterator
Dune::UGGrid<2>::leafend<2,Dune::Ghost_Partition>() const;


// ////////////////////////////////////////////////////////////////////////////////////
//   Explicitly instantiate the necessary member templates contained in UGGrid<3>
// ////////////////////////////////////////////////////////////////////////////////////
template Dune::UGGrid<3>::Traits::Codim<0>::LevelIterator Dune::UGGrid<3>::lbegin<0>(int level) const;
template Dune::UGGrid<3>::Traits::Codim<3>::LevelIterator Dune::UGGrid<3>::lbegin<3>(int level) const;

template Dune::UGGrid<3>::Traits::Codim<0>::LevelIterator Dune::UGGrid<3>::lend<0>(int level) const;
template Dune::UGGrid<3>::Traits::Codim<3>::LevelIterator Dune::UGGrid<3>::lend<3>(int level) const;

template Dune::UGGrid<3>::Traits::Codim<0>::LeafIterator Dune::UGGrid<3>::leafbegin<0>() const;
template Dune::UGGrid<3>::Traits::Codim<3>::LeafIterator Dune::UGGrid<3>::leafbegin<3>() const;

template Dune::UGGrid<3>::Traits::Codim<0>::LeafIterator Dune::UGGrid<3>::leafend<0>() const;
template Dune::UGGrid<3>::Traits::Codim<3>::LeafIterator Dune::UGGrid<3>::leafend<3>() const;

// element level iterators
template Dune::UGGrid<3>::Traits::Codim<0>::Partition<Dune::Interior_Partition>::LevelIterator
Dune::UGGrid<3>::lbegin<0,Dune::Interior_Partition>(int level) const;
template Dune::UGGrid<3>::Traits::Codim<0>::Partition<Dune::InteriorBorder_Partition>::LevelIterator
Dune::UGGrid<3>::lbegin<0,Dune::InteriorBorder_Partition>(int level) const;
template Dune::UGGrid<3>::Traits::Codim<0>::Partition<Dune::Overlap_Partition>::LevelIterator
Dune::UGGrid<3>::lbegin<0,Dune::Overlap_Partition>(int level) const;
template Dune::UGGrid<3>::Traits::Codim<0>::Partition<Dune::OverlapFront_Partition>::LevelIterator
Dune::UGGrid<3>::lbegin<0,Dune::OverlapFront_Partition>(int level) const;
template Dune::UGGrid<3>::Traits::Codim<0>::Partition<Dune::All_Partition>::LevelIterator
Dune::UGGrid<3>::lbegin<0,Dune::All_Partition>(int level) const;
template Dune::UGGrid<3>::Traits::Codim<0>::Partition<Dune::Ghost_Partition>::LevelIterator
Dune::UGGrid<3>::lbegin<0,Dune::Ghost_Partition>(int level) const;

template Dune::UGGrid<3>::Traits::Codim<0>::Partition<Dune::Interior_Partition>::LevelIterator
Dune::UGGrid<3>::lend<0,Dune::Interior_Partition>(int level) const;
template Dune::UGGrid<3>::Traits::Codim<0>::Partition<Dune::InteriorBorder_Partition>::LevelIterator
Dune::UGGrid<3>::lend<0,Dune::InteriorBorder_Partition>(int level) const;
template Dune::UGGrid<3>::Traits::Codim<0>::Partition<Dune::Overlap_Partition>::LevelIterator
Dune::UGGrid<3>::lend<0,Dune::Overlap_Partition>(int level) const;
template Dune::UGGrid<3>::Traits::Codim<0>::Partition<Dune::OverlapFront_Partition>::LevelIterator
Dune::UGGrid<3>::lend<0,Dune::OverlapFront_Partition>(int level) const;
template Dune::UGGrid<3>::Traits::Codim<0>::Partition<Dune::All_Partition>::LevelIterator
Dune::UGGrid<3>::lend<0,Dune::All_Partition>(int level) const;
template Dune::UGGrid<3>::Traits::Codim<0>::Partition<Dune::Ghost_Partition>::LevelIterator
Dune::UGGrid<3>::lend<0,Dune::Ghost_Partition>(int level) const;

// Vertex level iterators
template Dune::UGGrid<3>::Traits::Codim<3>::Partition<Dune::Interior_Partition>::LevelIterator
Dune::UGGrid<3>::lbegin<3,Dune::Interior_Partition>(int level) const;
template Dune::UGGrid<3>::Traits::Codim<3>::Partition<Dune::InteriorBorder_Partition>::LevelIterator
Dune::UGGrid<3>::lbegin<3,Dune::InteriorBorder_Partition>(int level) const;
template Dune::UGGrid<3>::Traits::Codim<3>::Partition<Dune::Overlap_Partition>::LevelIterator
Dune::UGGrid<3>::lbegin<3,Dune::Overlap_Partition>(int level) const;
template Dune::UGGrid<3>::Traits::Codim<3>::Partition<Dune::OverlapFront_Partition>::LevelIterator
Dune::UGGrid<3>::lbegin<3,Dune::OverlapFront_Partition>(int level) const;
template Dune::UGGrid<3>::Traits::Codim<3>::Partition<Dune::All_Partition>::LevelIterator
Dune::UGGrid<3>::lbegin<3,Dune::All_Partition>(int level) const;
template Dune::UGGrid<3>::Traits::Codim<3>::Partition<Dune::Ghost_Partition>::LevelIterator
Dune::UGGrid<3>::lbegin<3,Dune::Ghost_Partition>(int level) const;

template Dune::UGGrid<3>::Traits::Codim<3>::Partition<Dune::Interior_Partition>::LevelIterator
Dune::UGGrid<3>::lend<3,Dune::Interior_Partition>(int level) const;
template Dune::UGGrid<3>::Traits::Codim<3>::Partition<Dune::InteriorBorder_Partition>::LevelIterator
Dune::UGGrid<3>::lend<3,Dune::InteriorBorder_Partition>(int level) const;
template Dune::UGGrid<3>::Traits::Codim<3>::Partition<Dune::Overlap_Partition>::LevelIterator
Dune::UGGrid<3>::lend<3,Dune::Overlap_Partition>(int level) const;
template Dune::UGGrid<3>::Traits::Codim<3>::Partition<Dune::OverlapFront_Partition>::LevelIterator
Dune::UGGrid<3>::lend<3,Dune::OverlapFront_Partition>(int level) const;
template Dune::UGGrid<3>::Traits::Codim<3>::Partition<Dune::All_Partition>::LevelIterator
Dune::UGGrid<3>::lend<3,Dune::All_Partition>(int level) const;
template Dune::UGGrid<3>::Traits::Codim<3>::Partition<Dune::Ghost_Partition>::LevelIterator
Dune::UGGrid<3>::lend<3,Dune::Ghost_Partition>(int level) const;

// Element leaf iterators
template Dune::UGGrid<3>::Traits::Codim<0>::Partition<Dune::Interior_Partition>::LeafIterator
Dune::UGGrid<3>::leafbegin<0,Dune::Interior_Partition>() const;
template Dune::UGGrid<3>::Traits::Codim<0>::Partition<Dune::InteriorBorder_Partition>::LeafIterator
Dune::UGGrid<3>::leafbegin<0,Dune::InteriorBorder_Partition>() const;
template Dune::UGGrid<3>::Traits::Codim<0>::Partition<Dune::Overlap_Partition>::LeafIterator
Dune::UGGrid<3>::leafbegin<0,Dune::Overlap_Partition>() const;
template Dune::UGGrid<3>::Traits::Codim<0>::Partition<Dune::OverlapFront_Partition>::LeafIterator
Dune::UGGrid<3>::leafbegin<0,Dune::OverlapFront_Partition>() const;
template Dune::UGGrid<3>::Traits::Codim<0>::Partition<Dune::All_Partition>::LeafIterator
Dune::UGGrid<3>::leafbegin<0,Dune::All_Partition>() const;
template Dune::UGGrid<3>::Traits::Codim<0>::Partition<Dune::Ghost_Partition>::LeafIterator
Dune::UGGrid<3>::leafbegin<0,Dune::Ghost_Partition>() const;

template Dune::UGGrid<3>::Traits::Codim<0>::Partition<Dune::Interior_Partition>::LeafIterator
Dune::UGGrid<3>::leafend<0,Dune::Interior_Partition>() const;
template Dune::UGGrid<3>::Traits::Codim<0>::Partition<Dune::InteriorBorder_Partition>::LeafIterator
Dune::UGGrid<3>::leafend<0,Dune::InteriorBorder_Partition>() const;
template Dune::UGGrid<3>::Traits::Codim<0>::Partition<Dune::Overlap_Partition>::LeafIterator
Dune::UGGrid<3>::leafend<0,Dune::Overlap_Partition>() const;
template Dune::UGGrid<3>::Traits::Codim<0>::Partition<Dune::OverlapFront_Partition>::LeafIterator
Dune::UGGrid<3>::leafend<0,Dune::OverlapFront_Partition>() const;
template Dune::UGGrid<3>::Traits::Codim<0>::Partition<Dune::All_Partition>::LeafIterator
Dune::UGGrid<3>::leafend<0,Dune::All_Partition>() const;
template Dune::UGGrid<3>::Traits::Codim<0>::Partition<Dune::Ghost_Partition>::LeafIterator
Dune::UGGrid<3>::leafend<0,Dune::Ghost_Partition>() const;

// Vertex leaf iterators
template Dune::UGGrid<3>::Traits::Codim<3>::Partition<Dune::Interior_Partition>::LeafIterator
Dune::UGGrid<3>::leafbegin<3,Dune::Interior_Partition>() const;
template Dune::UGGrid<3>::Traits::Codim<3>::Partition<Dune::InteriorBorder_Partition>::LeafIterator
Dune::UGGrid<3>::leafbegin<3,Dune::InteriorBorder_Partition>() const;
template Dune::UGGrid<3>::Traits::Codim<3>::Partition<Dune::Overlap_Partition>::LeafIterator
Dune::UGGrid<3>::leafbegin<3,Dune::Overlap_Partition>() const;
template Dune::UGGrid<3>::Traits::Codim<3>::Partition<Dune::OverlapFront_Partition>::LeafIterator
Dune::UGGrid<3>::leafbegin<3,Dune::OverlapFront_Partition>() const;
template Dune::UGGrid<3>::Traits::Codim<3>::Partition<Dune::All_Partition>::LeafIterator
Dune::UGGrid<3>::leafbegin<3,Dune::All_Partition>() const;
template Dune::UGGrid<3>::Traits::Codim<3>::Partition<Dune::Ghost_Partition>::LeafIterator
Dune::UGGrid<3>::leafbegin<3,Dune::Ghost_Partition>() const;

template Dune::UGGrid<3>::Traits::Codim<3>::Partition<Dune::Interior_Partition>::LeafIterator
Dune::UGGrid<3>::leafend<3,Dune::Interior_Partition>() const;
template Dune::UGGrid<3>::Traits::Codim<3>::Partition<Dune::InteriorBorder_Partition>::LeafIterator
Dune::UGGrid<3>::leafend<3,Dune::InteriorBorder_Partition>() const;
template Dune::UGGrid<3>::Traits::Codim<3>::Partition<Dune::Overlap_Partition>::LeafIterator
Dune::UGGrid<3>::leafend<3,Dune::Overlap_Partition>() const;
template Dune::UGGrid<3>::Traits::Codim<3>::Partition<Dune::OverlapFront_Partition>::LeafIterator
Dune::UGGrid<3>::leafend<3,Dune::OverlapFront_Partition>() const;
template Dune::UGGrid<3>::Traits::Codim<3>::Partition<Dune::All_Partition>::LeafIterator
Dune::UGGrid<3>::leafend<3,Dune::All_Partition>() const;
template Dune::UGGrid<3>::Traits::Codim<3>::Partition<Dune::Ghost_Partition>::LeafIterator
Dune::UGGrid<3>::leafend<3,Dune::Ghost_Partition>() const;
