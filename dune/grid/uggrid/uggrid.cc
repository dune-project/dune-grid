// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <config.h>

#include <set>
#include <map>
#include <memory>

#include <dune/grid/uggrid.hh>

#if ModelP and DUNE_UGGRID_HAVE_PPIFCONTEXT
#  include <dune/uggrid/parallel/ppif/ppifcontext.hh>
#endif

/** \todo Remove the following two includes once getAllSubfaces... is gone */
#include <list>
#include <iterator>
#include <dune/common/stdstreams.hh>
#include <dune/grid/common/mcmgmapper.hh>

namespace Dune {

//***********************************************************************
//
// --UGGrid
// --Grid
//
//***********************************************************************

template<> int UGGrid<2>::numOfUGGrids = 0;
template<> int UGGrid<3>::numOfUGGrids = 0;


template <int dim>
UGGrid < dim >::UGGrid(UGCommunication comm)
  : multigrid_(nullptr),
    ccobj_(comm),
    leafIndexSet_(*this),
    idSet_(*this),
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
  typename UG_NS<dim>::CoeffProcPtr coeffs[1] = {nullptr};
  typename UG_NS<dim>::UserProcPtr upp[1] = {nullptr};

  // Create unique problem name
  std::stringstream numberAsAscii;
  numberAsAscii << numOfUGGrids;
  name_ = "DuneUGGrid_" + std::string((dim==2) ? "2" : "3") + std::string("d_") + numberAsAscii.str();

  std::string problemName = name_ + "_Problem";

  if (UG_NS<dim>::CreateBoundaryValueProblem(problemName.c_str(), 1,coeffs,1,upp) == nullptr)
    DUNE_THROW(GridError, "UG" << dim << "d::CreateBoundaryValueProblem() returned an error code!");

  numOfUGGrids++;

  dverb << "UGGrid<" << dim << "> with name " << name_ << " created!" << std::endl;

}


template < int dim >
UGGrid < dim >::~UGGrid() noexcept(false)
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

  // Shut down UG if this was the last existing UGGrid object.
  // Since we started both 2d and 3d versions of UG even if we don't need both,
  // we have to shut them down both, too.
  if (UGGrid<2>::numOfUGGrids + UGGrid<3>::numOfUGGrids == 0) {
    UG_NS<2>::ExitUg();
    UG_NS<3>::ExitUg();
  }
}

template < int dim >
int UGGrid < dim >::maxLevel() const
{
  if (!multigrid_)
    DUNE_THROW(GridError, "The grid has not been properly initialized!");

  return multigrid_->topLevel;
}


template < int dim >
int UGGrid < dim >::size (int level, int codim) const
{
  return levelIndexSet(level).size(codim);
}


template < int dim >
bool UGGrid < dim >::mark(int refCount,
                          const typename Traits::template Codim<0>::Entity& e )
{
  typename UG_NS<dim>::Element* target = e.impl().target_;

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
bool UGGrid < dim >::mark(const typename Traits::template Codim<0>::Entity& e,
                          typename UG_NS<dim>::RefinementRule rule,
                          int side)
{
  typename UG_NS<dim>::Element* target = e.impl().target_;

  if (!UG_NS<dim>::isLeaf(target))
    return false;

  someElementHasBeenMarkedForRefinement_ = true;

  return UG_NS<dim>::MarkForRefinement(target, rule, side);

}

template <int dim>
int UGGrid<dim>::getMark(const typename Traits::template Codim<0>::Entity& e) const
{
  typename UG_NS<dim>::Element* target = e.impl().target_;

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
bool UGGrid <dim>::preAdapt()
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
bool UGGrid < dim >::adapt()
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

  // Skip test whether we have enough memory available
  int mgtest = UG_NS<dim>::GM_REFINE_NOHEAPTEST;

  int rv = AdaptMultiGrid(multigrid_,mode,seq,mgtest);

  if (rv!=0)
    DUNE_THROW(GridError, "UG::adapt() returned with error code " << rv);

  // Renumber everything
  setIndices(false, nullptr);

  // Return true iff the grid hierarchy changed
  //return !(bool)multigrid_->status;

  // grid has changed if at least one element was marked for refinement
  return someElementHasBeenMarkedForRefinement_;
}

template <int dim>
void UGGrid <dim>::postAdapt()
{
  for (int i=0; i<=maxLevel(); i++)
    for (const auto& element : elements(this->levelGridView(i)))
      UG_NS<dim>::WriteCW(element.impl().target_, UG_NS<dim>::NEWEL_CE, 0);

  // reset marker flags
  someElementHasBeenMarkedForRefinement_ = false;
  someElementHasBeenMarkedForCoarsening_ = false;
}

template < int dim >
void UGGrid < dim >::globalRefine(int n)
{
  for (int i=0; i<n; i++) {

    // mark all entities for grid refinement
    for (const auto& element : elements(this->leafGridView()))
      mark(1, element);

    this->preAdapt();
    adapt();

  }

  this->postAdapt();

}

template <int dim>
void UGGrid<dim>::getChildrenOfSubface(const typename Traits::template Codim<0>::Entity & e,
                                       int elementSide,
                                       int maxl,
                                       std::vector<typename Traits::template Codim<0>::Entity>& childElements,
                                       std::vector<unsigned char>& childElementSides) const
{

  typedef std::pair<typename UG_NS<dim>::Element*,int> ListEntryType;

  std::list<ListEntryType> list;

  // //////////////////////////////////////////////////////////////////////
  //   Change the input face number from Dune numbering to UG numbering
  // //////////////////////////////////////////////////////////////////////

  elementSide = UGGridRenumberer<dim>::facesDUNEtoUG(elementSide, e.type());

  // ///////////////
  //   init list
  // ///////////////
  if (!e.isLeaf()   // Get_Sons_of_ElementSide returns GM_FATAL when called for a leaf !?!
      && e.level() < maxl) {

    typename UG_NS<dim>::Element* theElement = e.impl().target_;
    list.emplace_back(theElement, elementSide);

  }

  // //////////////////////////////////////////////////
  //   Traverse and collect all children of the side
  // //////////////////////////////////////////////////

  for (const auto& f : list)
  {
    typename UG_NS<dim>::Element* theElement = f.first;
    int side                                  = f.second;

    UG::INT Sons_of_Side = 0;
    typename UG_NS<dim>::Element* SonList[UG_NS<dim>::MAX_SONS];
    UG::INT SonSides[UG_NS<dim>::MAX_SONS];

    if (UG_NS<dim>::myLevel(theElement) < maxl) {

      int rv = Get_Sons_of_ElementSide(theElement,
                                       side,          // Input element side number
                                       &Sons_of_Side, // Number of topological sons of the element side
                                       SonList,       // Output elements
                                       SonSides,      // Output element side numbers
                                       true,
                                       true);

      if (rv != 0)
        DUNE_THROW(GridError, "Get_Sons_of_ElementSide returned with error value " << rv);

      for (int i=0; i<Sons_of_Side; i++)
        list.emplace_back(SonList[i], SonSides[i]);

    }
  }

  // Remove the element itself. We are only interested in the children.
  list.pop_front();

  // //////////////////////////////
  //   Extract result from stack
  // //////////////////////////////

  childElements.clear();
  childElements.reserve(std::distance(list.begin(), list.end()));
  childElementSides.resize(childElements.size());

  int i=0;
  for (const auto &f : list)
  {
    // Set element
    typedef typename Traits::template Codim< 0 >::Entity Entity;
    childElements.push_back( Entity( UGGridEntity< 0, dim, const UGGrid< dim > >( f.first, this ) ) );

    int side = f.second;

    // Dune numbers the faces of several elements differently than UG.
    // The following switch does the transformation
    childElementSides[i] = UGGridRenumberer<dim>::facesUGtoDUNE(side, childElements[i].type());
    ++i;
  }

}

template < int dim >
bool UGGrid < dim >::loadBalance(int minlevel)
{
  // Do nothing if we are on a single process
  if (comm().size()==1)
    return true;

  std::stringstream levelarg;
  levelarg << minlevel;
  UG_NS<dim>::lbs(levelarg.str().c_str(), multigrid_);

  // Renumber everything.
  // Note: this must not be called when on a single process, because it renumbers the zero-level
  // elements and vertices.
  setIndices(true, nullptr);

  return true;
}

template < int dim >
bool UGGrid < dim >::loadBalance(const std::vector<Rank>& targetProcessors, unsigned int fromLevel)
{
  // Do nothing if we are on a single process
  if (comm().size()==1)
    return true;

#ifdef ModelP  // The Partition field only exists if ModelP is set
  if (int(targetProcessors.size()) != this->leafGridView().size(0))
    DUNE_THROW(Exception, "targetProcessors argument does not have the correct size");

  // Get unique consecutive index across different element types
  typedef MultipleCodimMultipleGeomTypeMapper<typename Base::LeafGridView> ElementMapper;
  ElementMapper elementMapper(this->leafGridView(), mcmgElementLayout());

  // Loop over all elements of all level, in decreasing level number.
  // If the element is a leaf, take its target rank from the input targetProcessors array.
  // If it is not, assign it to the processor most of its children are assigned to.
  for (int i=maxLevel(); i>=0; i--) {
    for (const auto& element : elements(this->levelGridView(i), Partitions::interior)) {

      if (element.isLeaf()) {

        int targetRank = targetProcessors[elementMapper.index(element)];

        // sanity check
        if (targetRank >= comm().size())
          DUNE_THROW(GridError, "Requesting target processor " << targetRank <<
                     ", but only " << comm().size() << " processors are available.");

        UG_NS<dim>::Partition(element.impl().target_) = targetRank;
      } else {

        std::map<Rank,unsigned int> rank;
        Rank mostFrequentRank = 0;    // which rank occurred most often?
        unsigned int mostFrequentCount = 0;   // how often did it occur?

        // Loop over all children and collect the ranks they are assigned to
        for (const auto& child : descendantElements(element, element.level() + 1)) {

          auto childRank = UG_NS<dim>::Partition(child.impl().target_);

          if (rank.find(childRank) == rank.end())
            rank[childRank] = 1;
          else
            rank[childRank]++;

          if (rank[childRank] > mostFrequentCount) {
            mostFrequentRank = childRank;
            mostFrequentCount = rank[childRank];
          }

        }

        // Assign rank that occurred most often
        UG_NS<dim>::Partition(element.impl().target_) = mostFrequentRank;

      }
    }
  }
#endif

  int errCode = UG_NS<dim>::TransferGridFromLevel(multigrid_, fromLevel);

  if (errCode)
    DUNE_THROW(GridError, "UG" << dim << "d::TransferGridFromLevel returned error code " << errCode);

  // Renumber everything.
  // Note: this must not be called when on a single process, because it renumbers the zero-level
  // elements and vertices.
  setIndices(true, nullptr);

  return true;
}

#ifdef ModelP
template <int dim>
std::vector<typename UG_NS<dim>::DDD_IF> UGGrid<dim>::findDDDInterfaces(InterfaceType iftype,
                                                                        int codim) const
{
  std::vector<typename UG_NS<dim>::DDD_IF> dddIfaces;

  // UGGrid does not have overlap or front entities
  if (iftype == Overlap_OverlapFront_Interface || iftype == Overlap_All_Interface)
    return dddIfaces;

  if (codim == 0)
  {
    switch (iftype) {
    case InteriorBorder_InteriorBorder_Interface :
      // do not communicate anything: Elements cannot be in
      // the interior of two processes at the same time
      break;
    case InteriorBorder_All_Interface :
      dddIfaces.push_back(UG_NS<dim>::ElementVHIF(multigrid_->dddContext()));
      break;
    case All_All_Interface :
      dddIfaces.push_back(UG_NS<dim>::ElementSymmVHIF(multigrid_->dddContext()));
      break;
    default :
      DUNE_THROW(GridError,
                 "Element communication not supported for "
                 "interfaces of type  "
                 << iftype);
    }
  }
  else if (codim == dim)
  {
    switch (iftype)
    {
    case InteriorBorder_InteriorBorder_Interface :
      dddIfaces.push_back(UG_NS<dim>::BorderNodeSymmIF(multigrid_->dddContext()));
      break;
    case InteriorBorder_All_Interface :
      dddIfaces.push_back(UG_NS<dim>::NodeInteriorBorderAllIF(multigrid_->dddContext()));
      break;
    case All_All_Interface :
      dddIfaces.push_back(UG_NS<dim>::NodeAllIF(multigrid_->dddContext()));
      break;
    default :
      DUNE_THROW(GridError,
                 "Node communication not supported for "
                 "interfaces of type  "
                 << iftype);
    }
  }
  else if (codim == dim-1)
  {
    switch (iftype)
    {
    case InteriorBorder_InteriorBorder_Interface :
      dddIfaces.push_back(UG_NS<dim>::BorderEdgeSymmIF(multigrid_->dddContext()));
      break;
    case InteriorBorder_All_Interface :
      dddIfaces.push_back(UG_NS<dim>::EdgeVHIF(multigrid_->dddContext()));
      break;
    case All_All_Interface :
      dddIfaces.push_back(UG_NS<dim>::EdgeSymmVHIF(multigrid_->dddContext()));
      break;
    default :
      DUNE_THROW(GridError,
                 "Edge communication not supported for "
                 "interfaces of type  "
                 << iftype);
    }
  }
  else if (codim == 1)
  {
    switch (iftype)
    {
    case InteriorBorder_InteriorBorder_Interface :
      dddIfaces.push_back(UG_NS<dim>::BorderVectorSymmIF(multigrid_->dddContext()));
      break;
    case InteriorBorder_All_Interface :
      dddIfaces.push_back(UG_NS<dim>::FacetInteriorBorderAllIF(multigrid_->dddContext()));
      break;
    case All_All_Interface :
      dddIfaces.push_back(UG_NS<dim>::FacetAllAllIF(multigrid_->dddContext()));
      break;
    default :
      DUNE_THROW(GridError,
                 "Face communication not supported for "
                 "interfaces of type  "
                 << iftype);
    }
  }
  else
  {
    DUNE_THROW(GridError,
               "Communication codimension must be between 0 and " << dim << "!");
  }

  return dddIfaces;
};
#endif // ModelP


template < int dim >
void UGGrid < dim >::setPosition(const typename Traits::template Codim<dim>::Entity& e,
                                 const FieldVector<double, dim>& pos)
{
  typename UG_NS<dim>::Node* target = e.impl().target_;

  for (int i=0; i<dim; i++)
    target->myvertex->iv.x[i] = pos[i];
}

template <int dim>
void UGGrid<dim>::saveState(const std::string& filename) const
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
void UGGrid<dim>::loadState(const std::string& filename)
{
  const char* type = "asc";
  std::string problemName = name_ + "_Problem";
  std::string formatName = "DuneFormat2d";

  if (dim==2) {
    std::string formatName = "DuneFormat2d";
    multigrid_ = (typename UG_NS<dim>::MultiGrid*) UG::D2::LoadMultiGrid(
      name_.c_str(),
      filename.c_str(),
      type,
      problemName.c_str(),
      formatName.c_str(),
      0,    // dummy heap size
      true, //force,
      true, //optimizedIO,
      false //autosave
#if ModelP and DUNE_UGGRID_HAVE_PPIFCONTEXT
      , std::make_shared<PPIF::PPIFContext>(comm())
#endif
      );
  } else {
    std::string formatName = "DuneFormat3d";
    multigrid_ = (typename UG_NS<dim>::MultiGrid*) UG::D3::LoadMultiGrid(
      name_.c_str(),
      filename.c_str(),
      type,
      problemName.c_str(),
      formatName.c_str(),
      0,    // dummy heap size
      true, //force,
      true, //optimizedIO,
      false //autosave
#if ModelP and DUNE_UGGRID_HAVE_PPIFCONTEXT
      , std::make_shared<PPIF::PPIFContext>(comm())
#endif
      );
  }

  if (multigrid_==nullptr)
    DUNE_THROW(GridError, "In loadState()");
}

template < int dim >
void UGGrid < dim >::setIndices(bool setLevelZero,
                                      std::vector<unsigned int>* nodePermutation)
{
  // Create new level index sets if necessary
  for (int i=levelIndexSets_.size(); i<=maxLevel(); i++)
    levelIndexSets_.push_back(std::make_shared<UGGridLevelIndexSet<const UGGrid<dim> > >());

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

template class UGGrid<2>;
template class UGGrid<3>;

} /* namespace Dune */
