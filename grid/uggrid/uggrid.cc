// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <set>

#include <dune/grid/uggrid.hh>
#include "boundaryextractor.hh"

/** \todo Remove the following include once getAllSubfaces... is gone */
#include <dune/common/sllist.hh>
#include <dune/common/stdstreams.hh>


using namespace Dune;

static int boundarySegmentWrapper2d(void *data, double *param, double *result)
{
  const BoundarySegment<2>* boundarySegment = static_cast<const BoundarySegment<2>*>(data);

  FieldVector<double, 2> global = (*boundarySegment)(*((FieldVector<double,1>*)param));

  result[0] = global[0];
  result[1] = global[1];

  return 0;
}

static int boundarySegmentWrapper3d(void *data, double *param, double *result)
{
  const BoundarySegment<3>* boundarySegment = static_cast<const BoundarySegment<3>*>(data);

  FieldVector<double, 3> global = (*boundarySegment)(*((FieldVector<double,2>*)param));

  result[0] = global[0];
  result[1] = global[1];
  result[2] = global[2];

  return 0;
}

//***********************************************************************
//
// --UGGrid
// --Grid
//
//***********************************************************************

template<> int Dune::UGGrid<2>::numOfUGGrids = 0;
template<> int Dune::UGGrid<3>::numOfUGGrids = 0;


template <int dim>
Dune::UGGrid < dim >::UGGrid(unsigned int heapSize)
  : multigrid_(NULL),
    leafIndexSet_(*this),
    globalIdSet_(*this),
    localIdSet_(*this),
    refinementType_(LOCAL),
    closureType_(GREEN),
    someElementHasBeenMarkedForRefinement_(false)
{
  heapsize = heapSize;

  if (numOfUGGrids==0) {

    // Init the UG system
    int argc = 1;
    char* arg = strdup("dune.exe");
    char** argv = &arg;

    if (UG_NS<dim>::InitUg(&argc, &argv))
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

#ifndef UG_LGMDOMAIN
  if (UG_NS<dim>::CreateBoundaryValueProblem(problemName.c_str(), 1,coeffs,1,upp) == NULL)
    DUNE_THROW(GridError, "UG" << dim << "d::CreateBoundaryValueProblem() returned and error code!");
#endif

  if (numOfUGGrids==0) {

    if (dim==2)
    {
      char* nfarg = strdup("newformat DuneFormat2d");
      if (UG_NS<dim>::CreateFormatCmd(1, &nfarg))
        DUNE_THROW(GridError, "UG" << dim << "d::CreateFormat() returned and error code!");
      free(nfarg);
    }
    if (dim==3)
    {
      char* newArgs[2];
      for (int i=0; i<2; i++)
        newArgs[i] = (char*)::malloc(50*sizeof(char));

      sprintf(newArgs[0], "newformat DuneFormat3d" );
      sprintf(newArgs[1], "V s1 : vt 1" );             // generates side vectors in 3D

      if (UG_NS<dim>::CreateFormatCmd(2, newArgs))
        DUNE_THROW(GridError, "UG" << dim << "d::CreateFormat() returned and error code!");

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
  for (unsigned int i=0; i<boundarySegments_.size(); i++)
    delete boundarySegments_[i];

  // Delete the UG multigrid if there is one (== createEnd() has been called)
  if (multigrid_) {
    // Set UG's currBVP variable to the BVP corresponding to this
    // grid.  This is necessary if we have more than one UGGrid in use.
    // DisposeMultiGrid will crash if we don't do this
    UG_NS<dim>::Set_Current_BVP(multigrid_->theBVP);
    UG_NS<dim>::DisposeMultiGrid(multigrid_);
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

  const typename UG_NS<dim>::Grid* theGrid = multigrid_->grids[level];

  if (!theGrid)
    DUNE_THROW(GridError, "LevelIterator in nonexisting level " << level << " requested!");

  return UGGridLevelIterator<codim, All_Partition, const UGGrid<dim> >(theGrid);
}

template<int dim>
template<int codim, Dune::PartitionIteratorType PiType>
typename Dune::UGGrid<dim>::Traits::template Codim<codim>::template Partition<PiType>::LevelIterator
Dune::UGGrid<dim>::lbegin (int level) const
{
  if (!multigrid_)
    DUNE_THROW(GridError, "The grid has not been properly initialized!");

  typename UG_NS<dim>::Grid* theGrid = multigrid_->grids[level];

  if (!theGrid)
    DUNE_THROW(GridError, "LevelIterator in nonexisting level " << level << " requested!");

  return UGGridLevelIterator<codim, PiType, const UGGrid<dim> >(theGrid);
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
  typename UG_NS<dim>::Element* target = getRealImplementation(e).target_;

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

    someElementHasBeenMarkedForRefinement_ = true;
    return true;
  } else
    DUNE_THROW(GridError, "UGGrid only supports refCount values -1, 0, and 1 for mark()!");

}

template < int dim >
bool Dune::UGGrid < dim >::mark(const typename Traits::template Codim<0>::Entity& e,
                                typename UG_NS<dim>::RefinementRule rule,
                                int side)
{
  typename UG_NS<dim>::Element* target = getRealImplementation(e).target_;

  if (!UG_NS<dim>::isLeaf(target))
    return false;

  someElementHasBeenMarkedForRefinement_ = true;

  return UG_NS<dim>::MarkForRefinement(target, rule, side);

}

template <int dim>
int Dune::UGGrid<dim>::getMark(const typename Traits::template Codim<0>::Entity& e) const
{
  typename UG_NS<dim>::Element* target = getRealImplementation(e).target_;

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
  return someElementHasBeenMarkedForRefinement_;
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

  someElementHasBeenMarkedForRefinement_ = false;

  // Return true iff the grid hierarchy changed
  return !(bool)multigrid_->status;
}

template <int dim>
void Dune::UGGrid <dim>::postAdapt()
{
  for (int i=0; i<=maxLevel(); i++) {

    typename Traits::template Codim<0>::LevelIterator eIt    = lbegin<0>(i);
    typename Traits::template Codim<0>::LevelIterator eEndIt = lend<0>(i);

    for (; eIt!=eEndIt; ++eIt)
      UG_NS<dim>::WriteCW(getRealImplementation(*eIt).target_, UG_NS<dim>::NEWEL_CE, 0);

  }
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
void Dune::UGGrid<dim>::getChildrenOfSubface(typename Traits::template Codim<0>::EntityPointer & e,
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

    typename UG_NS<dim>::Element* theElement = getRealImplementation(*e).target_;

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
                              side,             // Input element side number
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

  /** \todo Initialized with a dummy value because EntityPointer isn't default constructible */
  childElements.resize(list.size(), lbegin<0>(0));
  childElementSides.resize(list.size());

  int i=0;
  for (f = list.begin(); f!=list.end(); ++f, ++i) {

    // Set element
    this->getRealImplementation(childElements[i]).setToTarget(f->first);

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
void Dune::UGGrid < dim >::createEnd()
{
#ifdef UG_LGMDOMAIN
  DUNE_THROW(GridError, "You cannot call createEnd() when your UGGrid has been configured for LGM!");
#else

  // ///////////////////////////////////////////
  //   Extract grid boundary segments
  // ///////////////////////////////////////////
  std::set<UGGridBoundarySegment<dim> > boundarySegments;
  typedef typename std::set<UGGridBoundarySegment<dim> >::iterator SetIterator;

  BoundaryExtractor::detectBoundarySegments(elementTypes_, elementVertices_, boundarySegments);
  if (boundarySegments.size() == 0)
    DUNE_THROW(GridError, "Couldn't extract grid boundary.");

  std::vector<int> isBoundaryNode;
  BoundaryExtractor::detectBoundaryNodes(boundarySegments, vertexPositions_.size(), isBoundaryNode);

  dverb << boundarySegments.size() << " boundary segments were found!" << std::endl;

  if (boundarySegments_.size() > boundarySegments.size())
    DUNE_THROW(GridError, "You have supplied " << boundarySegments_.size()
                                               << " parametrized boundary segments,  but the coarse grid has only "
                                               << boundarySegments.size() << " boundary faces!");

  // Count number of nodes on the boundary
  int noOfBNodes = 0;
  for (unsigned int i=0; i<isBoundaryNode.size(); i++) {
    if (isBoundaryNode[i] != -1)
      noOfBNodes++;
  }

  // ///////////////////////////////////////////////////////
  //   UG needs all boundary vertices first.  We set up
  //   up an array to keep track of the reordering.
  // ///////////////////////////////////////////////////////
  std::vector<unsigned int> nodePermutation;

  for (unsigned int i=0; i<isBoundaryNode.size(); ++i) {
    if (isBoundaryNode[i]!=-1)
      nodePermutation.push_back(i);
  }

  for (unsigned int i=0; i<isBoundaryNode.size(); ++i) {
    if (isBoundaryNode[i]==-1)
      nodePermutation.push_back(i);
  }

  // ///////////////////////////////////////////
  //   Create the domain data structure
  // ///////////////////////////////////////////
  unsigned int noOfBSegments = boundarySegments.size();
  std::string domainName = name_ + "_Domain";
  const double midPoint[3] = {0, 0, 0};

  if (UG_NS<dim>::CreateDomain(domainName.c_str(),     // The domain name
                               midPoint,               // Midpoint of a circle enclosing the grid, only needed for the UG graphics
                               1,                      // Radius of the enclosing circle
                               noOfBSegments,
                               noOfBNodes,
                               false) == NULL)                 // The domain is not convex
    DUNE_THROW(GridError, "Calling UG" << dim << "d::CreateDomain failed!");

  // ///////////////////////////////////////////
  //   Insert the boundary segments
  // ///////////////////////////////////////////
  unsigned int i;
  for (i=0; i<boundarySegments_.size(); i++) {

    /** \todo Due to some UG weirdness, in 3d, CreateBoundarySegment always expects
        this array to have four entries, even if only a triangular segment is
        inserted.  If not, undefined values are will be introduced. */
    int vertices_c_style[dim*2-2];
    for (int j=0; j<dim*2-2; j++)
      vertices_c_style[j] = boundarySegmentVertices_[i][j];

    // Create dummy parameter ranges
    const double alpha[2] = {0, 0};
    const double beta[2]  = {1, 1};

    // Create some boundary segment name
    char segmentName[20];
    if(sprintf(segmentName, "BS %d", i) < 0)
      DUNE_THROW(GridError, "sprintf returned error code!");

    // Actually create the segment
    if (UG_NS<dim>::CreateBoundarySegment(segmentName,                     // internal name of the boundary segment
                                          1,                               //  id of left subdomain
                                          2,                               //  id of right subdomain
                                          i,         // Index of the segment
                                          1,                        // Resolution, only for the UG graphics
                                          vertices_c_style,
                                          alpha,
                                          beta,
                                          (dim==2)
                                          ? boundarySegmentWrapper2d
                                          : boundarySegmentWrapper3d,
                                          const_cast<BoundarySegment<dim>*>(boundarySegments_[i]))==NULL) {
      DUNE_THROW(GridError, "Calling UG" << dim << "d::CreateBoundarySegment failed!");
    }

    // /////////////////////////////////////////////////////////////////////
    //   Remove this segment from the set of computed boundary segments.
    // /////////////////////////////////////////////////////////////////////

    UGGridBoundarySegment<dim> thisSegment;
    /** \todo Not nice: we need to copy because the array types are different */
    for (int j=0; j<2*dim-2; j++)
      thisSegment[j] = boundarySegmentVertices_[i][j];

    if (boundarySegments.erase(thisSegment)==0)
      DUNE_THROW(GridError, "You have provided a boundary parametrization for"
                 << " a segment which is not boundary segment in the grid!");
  }


  // ///////////////////////////////////////////////////////////////////////
  //   The boundary segments remaining in the std::set boundarySegments
  //   have not been provided with an explicit parametrization.  They are
  //   inserted into the domain as straight boundary segments.
  // ///////////////////////////////////////////////////////////////////////
  SetIterator it = boundarySegments.begin();

  for (; it != boundarySegments.end(); ++it, ++i) {

    const UGGridBoundarySegment<dim>& thisSegment = *it;

    // Copy the vertices into a C-style array
    int vertices_c_style[4];

    for (int j=0; j<thisSegment.numVertices(); j++)
      vertices_c_style[j] = isBoundaryNode[thisSegment[j]];

#ifndef UG_LGMDOMAIN
    // Create some boundary segment name
    char segmentName[20];
    if(sprintf(segmentName, "BS %d", i) < 0)
      DUNE_THROW(GridError, "sprintf returned error code!");

    if (dim==2) {

      double segmentCoordinates[2][2];
      for (int j=0; j<thisSegment.numVertices(); j++)
        for (int k=0; k<dim; k++)
          segmentCoordinates[j][k] = vertexPositions_[thisSegment[j]][k];

      if (UG::D2::CreateLinearSegment(segmentName,
                                      1,               /*id of left subdomain */
                                      2,              /*id of right subdomain*/
                                      i,                  /*id of segment*/
                                      thisSegment.numVertices(),          // Number of corners
                                      vertices_c_style,
                                      segmentCoordinates
                                      )==NULL)
        DUNE_THROW(IOError, "Error calling CreateLinearSegment");

    } else {

      double segmentCoordinates[4][3];
      for (int j=0; j<thisSegment.numVertices(); j++)
        for (int k=0; k<dim; k++)
          segmentCoordinates[j][k] = vertexPositions_[thisSegment[j]][k];

      if (UG::D3::CreateLinearSegment(segmentName,
                                      1,               /*id of left subdomain */
                                      2,              /*id of right subdomain*/
                                      i,                  /*id of segment*/
                                      thisSegment.numVertices(),                  // Number of corners
                                      vertices_c_style,
                                      segmentCoordinates
                                      )==NULL)
        DUNE_THROW(IOError, "Error calling CreateLinearSegment");
    }
#endif

  }


  // ///////////////////////////////////////////
  //   Call configureCommand and newCommand
  // ///////////////////////////////////////////

  //configure @PROBLEM $d @DOMAIN;
  std::string configureArgs[2] = {"configure " + name_ + "_Problem", "d " + name_ + "_Domain"};
  const char* configureArgs_c[2] = {configureArgs[0].c_str(), configureArgs[1].c_str()};

  if (UG_NS<dim>::ConfigureCommand(2, configureArgs_c))
    DUNE_THROW(GridError, "Calling UG" << dim << "d::ConfigureCommand failed!");

  //new @PROBLEM $b @PROBLEM $f @FORMAT $h @HEAP;
  char* newArgs[4];
  for (int i=0; i<4; i++)
    newArgs[i] = (char*)::malloc(50*sizeof(char));

  sprintf(newArgs[0], "new %s", name_.c_str());

  sprintf(newArgs[1], "b %s_Problem", name_.c_str());
  sprintf(newArgs[2], "f DuneFormat%dd", dim);
  sprintf(newArgs[3], "h %dM", heapsize);

  if (UG_NS<dim>::NewCommand(4, newArgs))
    DUNE_THROW(GridError, "UGGrid<" << dim << ">::makeNewMultigrid failed!");

  for (int i=0; i<4; i++)
    free(newArgs[i]);

  // Get a direct pointer to the newly created multigrid
  multigrid_ = UG_NS<dim>::GetMultigrid(name_.c_str());
  if (!multigrid_)
    DUNE_THROW(GridError, "UGGrid<" << dim << ">::GetMultigrid failed!");

  // ///////////////////////////////////////////////////////////////
  // If we are in a parallel setting and we are _not_ the master
  // process we can stop here.
  // ///////////////////////////////////////////////////////////////
  if (PPIF::me!=0) {
    // Complete the UG-internal grid data structure even if we are
    // not the master process. (CreateAlgebra communicates via MPI
    // so we would be out of sync if we don't do this here...)
    if (CreateAlgebra(multigrid_) != UG_NS<dim>::GM_OK)
      DUNE_THROW(IOError, "Call of 'UG::D" << dim << "::CreateAlgebra' failed!");

    /* here all temp memory since CreateMultiGrid is released */
    Release(multigrid_->theHeap, UG::FROM_TOP, multigrid_->MarkKey);
    multigrid_->MarkKey = 0;

    return;
  }

  // ////////////////////////////////////////////////
  //   Actually insert the interior vertices
  // ////////////////////////////////////////////////
  int nodeCounter = noOfBNodes;
  for (size_t i=0; i<vertexPositions_.size(); i++) {
    if (isBoundaryNode[i] != -1)
      continue;
    if (UG_NS<dim>::InsertInnerNode(multigrid_->grids[0], &((vertexPositions_[i])[0])) == NULL)
      DUNE_THROW(GridError, "Inserting a vertex into UGGrid failed!");

    isBoundaryNode[i] = nodeCounter++;
  }

  vertexPositions_.resize(0);

  // ////////////////////////////////////////////////
  //   Actually insert all the elements
  // ////////////////////////////////////////////////

  std::vector<const typename UG_NS<dim>::Node*> nodePointers(isBoundaryNode.size());
  for (typename UG_NS<dim>::Node* theNode=UG_NS<dim>::FirstNode(multigrid_->grids[0]); theNode!=NULL; theNode=theNode->succ)
    nodePointers[theNode->id] = theNode;

  int idx = 0;
  for (size_t i=0; i<elementTypes_.size(); i++) {

    const typename UG_NS<dim>::Node* vertices[elementTypes_[i]];
    for (size_t j=0; j<elementTypes_[i]; j++)
      vertices[j] = nodePointers[isBoundaryNode[elementVertices_[idx++]]];

    if (InsertElement(multigrid_->grids[0], elementTypes_[i],const_cast<typename UG_NS<dim>::Node**>(vertices),NULL,NULL,NULL)==NULL)
      DUNE_THROW(GridError, "Inserting element into UGGrid failed!");
  }

  // Not needed any more
  elementTypes_.resize(0);
  elementVertices_.resize(0);

  // Complete the UG-internal grid data structure
  if (CreateAlgebra(multigrid_) != UG_NS<dim>::GM_OK)
    DUNE_THROW(IOError, "Call of 'UG::D" << dim << "::CreateAlgebra' failed!");

  /* here all temp memory since CreateMultiGrid is released */
  Release(multigrid_->theHeap, UG::FROM_TOP, multigrid_->MarkKey);
  multigrid_->MarkKey = 0;

  // Set the local indices
  setIndices(true, &nodePermutation);

  // Clear refinement flags
  typename Traits::template Codim<0>::LevelIterator eIt    = lbegin<0>(0);
  typename Traits::template Codim<0>::LevelIterator eEndIt = lend<0>(0);

  for (; eIt!=eEndIt; ++eIt)
    UG_NS<dim>::WriteCW(getRealImplementation(*eIt).target_, UG_NS<dim>::NEWEL_CE, 0);

#endif
}

template < int dim >
void Dune::UGGrid < dim >::createLGMGrid(const std::string& name)
{
#ifndef UG_LGMDOMAIN
  DUNE_THROW(GridError, "You cannot call createLGMGrid() when your UGGrid hasn't been configured for LGM!");
#else
  // ////////////////////////////////////////////////
  //   Create a boundary value problem (BVP)
  // ////////////////////////////////////////////////
  if (dim==2) {
    if (!UG::D2::CreateProblem(name.c_str(),NULL,NULL,NULL,0,NULL,0,NULL))
      DUNE_THROW(GridError, "LGM::CreateProblem() returned error code!");
  } else {
    if (!UG::D3::CreateProblem(name.c_str(),NULL,NULL,NULL,0,NULL,0,NULL))
      DUNE_THROW(GridError, "LGM::CreateProblem() returned error code!");
  }

  // ///////////////////////////////////////////
  //   Call configureCommand and newCommand
  // ///////////////////////////////////////////
#if 0 // Do we need this?
      //configure @PROBLEM $d @DOMAIN;
  std::string configureArgs[2] = {"configure " + name_ + "_Problem", "d " + name_ + "_Domain"};
  const char* configureArgs_c[2] = {configureArgs[0].c_str(), configureArgs[1].c_str()};

  if (UG_NS<dim>::ConfigureCommand(2, configureArgs_c))
    DUNE_THROW(GridError, "Calling UG" << dim << "d::ConfigureCommand failed!");
#endif

  //new @PROBLEM $b @PROBLEM $f @FORMAT $h @HEAP;
  char* newArgs[4];
  for (int i=0; i<4; i++)
    newArgs[i] = (char*)::malloc(50*sizeof(char));

  sprintf(newArgs[0], "new %s", name_.c_str());

  sprintf(newArgs[1], "b %s", name.c_str());
  sprintf(newArgs[2], "f DuneFormat%dd", dim);
  sprintf(newArgs[3], "h %dM", heapsize);

  if (UG_NS<dim>::NewCommand(4, newArgs))
    DUNE_THROW(GridError, "UGGrid::makeNewMultigrid failed!");

  for (int i=0; i<4; i++)
    free(newArgs[i]);

  // Get a direct pointer to the newly created multigrid
  multigrid_ = UG_NS<dim>::GetMultigrid(name_.c_str());
  if (!multigrid_)
    DUNE_THROW(GridError, "UGGrid::makeNewMultigrid failed!");

  // Complete the UG-internal grid data structure
  if (CreateAlgebra(multigrid_) != UG_NS<dim>::GM_OK)
    DUNE_THROW(IOError, "Call of 'UG::CreateAlgebra' failed!");

  /* here all temp memory since CreateMultiGrid is released */
  Release(multigrid_->theHeap, UG::FROM_TOP, multigrid_->MarkKey);
  multigrid_->MarkKey = 0;

  // Set the local indices
  setIndices();

#endif
}

template <int dim>
void Dune::UGGrid<dim>::insertBoundarySegment(const std::vector<unsigned int> vertices,
                                              const BoundarySegment<dim>* boundarySegment)
{
  array<unsigned int, dim*2-2> segmentVertices;

  for (size_t i=0; i<vertices.size(); i++)
    segmentVertices[i] = vertices[i];

  for (size_t i=vertices.size(); i<dim*2-2; i++)
    segmentVertices[i] = -1;

  // DUNE --> UG vertex renumbering for quadrilateral boundary segments
  if (vertices.size()==4) {
    segmentVertices[2] = vertices[3];
    segmentVertices[3] = vertices[2];
  }

  boundarySegmentVertices_.push_back(segmentVertices);

  // Append boundary segment class to the boundary segment class list, so we can
  // delete them all in the destructor
  boundarySegments_.push_back(boundarySegment);

}

template <int dim>
void Dune::UGGrid<dim>::insertVertex(const FieldVector<double,dim>& pos)
{
  vertexPositions_.push_back(pos);
}

template <int dim>
void Dune::UGGrid<dim>::insertElement(GeometryType type,
                                      const std::vector<unsigned int>& vertices)
{
  if (dim!=type.dim())
    DUNE_THROW(GridError, "You cannot insert a " << type
                                                 << " into a UGGrid<" << dim << ">!");

  int newIdx = elementVertices_.size();

  elementTypes_.push_back(vertices.size());
  for (unsigned int i=0; i<vertices.size(); i++)
    elementVertices_.push_back(vertices[i]);

  if (type.isTriangle()) {
    // Everything alright
    if (vertices.size() != 3)
      DUNE_THROW(GridError, "You have requested to enter a triangle, but you"
                 << " have provided " << vertices.size() << " vertices!");

  } else if (type.isQuadrilateral()) {

    if (vertices.size() != 4)
      DUNE_THROW(GridError, "You have requested to enter a quadrilateral, but you"
                 << " have provided " << vertices.size() << " vertices!");

    // DUNE and UG numberings differ --> reorder the vertices
    elementVertices_[newIdx+2] = vertices[3];
    elementVertices_[newIdx+3] = vertices[2];

  } else if (type.isTetrahedron()) {

    if (vertices.size() != 4)
      DUNE_THROW(GridError, "You have requested to enter a tetrahedron, but you"
                 << " have provided " << vertices.size() << " vertices!");

  } else if (type.isPyramid()) {

    if (vertices.size() != 5)
      DUNE_THROW(GridError, "You have requested to enter a pyramid, but you"
                 << " have provided " << vertices.size() << " vertices!");

  } else if (type.isPrism()) {

    if (vertices.size() != 6)
      DUNE_THROW(GridError, "You have requested to enter a prism, but you"
                 << " have provided " << vertices.size() << " vertices!");

  } else if (type.isHexahedron()) {

    if (vertices.size() != 8)
      DUNE_THROW(GridError, "You have requested to enter a hexahedron, but you"
                 << " have provided " << vertices.size() << " vertices!");

    // DUNE and UG numberings differ --> reorder the vertices
    elementVertices_[newIdx+2] = vertices[3];
    elementVertices_[newIdx+3] = vertices[2];
    elementVertices_[newIdx+6] = vertices[7];
    elementVertices_[newIdx+7] = vertices[6];

  } else {
    DUNE_THROW(GridError, "You cannot insert a " << type
                                                 << " into a UGGrid<" << dim << ">!");
  }

}

template < int dim >
void Dune::UGGrid < dim >::setPosition(typename Traits::template Codim<dim>::EntityPointer& e,
                                       const FieldVector<double, dim>& pos)
{
  typename UG_NS<dim>::Node* target = getRealImplementation(*e).target_;

  for (int i=0; i<dim; i++)
    target->myvertex->iv.x[i] = pos[i];
}

template <int dim>
Dune::FieldVector<typename UGGrid<dim>::ctype,dim> Dune::UGGrid<dim>::getBoundaryPosition(const typename Traits::LevelIntersectionIterator& iIt,
                                                                                          const FieldVector<ctype,dim-1>& localPos) const
{
  if (!iIt->boundary())
    DUNE_THROW(GridError, "call getBoundaryPosition() only for boundary intersections!");

  FieldVector<ctype,dim> result;

  if (dim==2) {

    int ugEdgeNumber = UGGridRenumberer<dim>::edgesDUNEtoUG(iIt->numberInSelf(), iIt->inside()->type());

    // Get UG element from which this intersection iterator is taken
    const typename UG_NS<dim>::Element* target = getRealImplementation(*iIt->inside()).target_;

    const typename UG_NS<dim>::Vertex* v0 = UG_NS<dim>::Corner(target,UG_NS<dim>::Corner_Of_Edge(target, ugEdgeNumber,0))->myvertex;
    const typename UG_NS<dim>::Vertex* v1 = UG_NS<dim>::Corner(target,UG_NS<dim>::Corner_Of_Edge(target, ugEdgeNumber,1))->myvertex;

    typename UG_NS<dim>::BNDP* bndp = UG_NS<dim>::BNDP_CreateBndP(multigrid_->theHeap,v0->bv.bndp, v1->bv.bndp,localPos[0]);
    if (bndp == NULL)
      DUNE_THROW(GridError, "UG::D" << dim << "::BNDP_CreateBndP() returned NULL!");

    if (UG_NS<dim>::BNDP_Global(bndp,&result[0]))
      DUNE_THROW(GridError, "UG::D" << dim << "::BNDP_Global() returned nonzero error code!");

    // Get rid of the boundary point again
    UG_NS<dim>::BNDP_Dispose(multigrid_->theHeap, bndp);

  } else
    DUNE_THROW(NotImplemented, "getBoundaryPosition for 3d not implemented!");

  return result;
}


template <int dim>
void Dune::UGGrid<dim>::saveState(const std::string& filename) const
{
  const char* type = "asc";
  const char* comment = "written by DUNE";

  if (dim==2)
    UG::D2::SaveMultiGrid((UG::D2::multigrid*)multigrid_,
                          (char*)filename.c_str(),
                          (char*)type,
                          (char*)comment,
                          0,      // autosave
                          0       // rename
                          );
  else
    UG::D3::SaveMultiGrid((UG::D3::multigrid*)multigrid_,
                          (char*)filename.c_str(),
                          (char*)type,
                          (char*)comment,
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
    multigrid_ = (typename UG_NS<dim>::MultiGrid*) UG::D2::LoadMultiGrid((char*)name_.c_str(),
                                                                         (char*)filename.c_str(),
                                                                         (char*)type,
                                                                         (char*)problemName.c_str(),
                                                                         (char*)formatName.c_str(),
                                                                         heapsize,
                                                                         true, //force,
                                                                         true, //optimizedIO,
                                                                         false //autosave
                                                                         );
  } else {
    std::string formatName = "DuneFormat3d";
    multigrid_ = (typename UG_NS<dim>::MultiGrid*) UG::D3::LoadMultiGrid((char*)name_.c_str(),
                                                                         (char*)filename.c_str(),
                                                                         (char*)type,
                                                                         (char*)problemName.c_str(),
                                                                         (char*)formatName.c_str(),
                                                                         heapsize,
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
  /** \todo The code in the following if-clause is a workaround to fix FlySpray issue 170
      (inconsistent codim 1 subIndices).  UGGrid uses UG's SideVector data structure to
      store the codim 1 subIndices of elements.  However the UG SideVector code is buggy
      and SideVectors are not correctly created.  There are two reasons why I provide
      a workaround instead of fixing UG directly:
      - it is not trivial to fix in UG, and I'd probably screw up other things if I
      started to mess around in the UG refinement routines
      - AFAIK Stefan Lang is rewriting that part of UG anyways.
      As soon as UGGrid does not rely on SideVectors anymore to store indices the following
      if-clause can be deleted.
   */
  if (dim==3) {

    for (int i=1; i<=maxLevel(); i++) {

      typename Traits::template Codim<0>::LevelIterator eIt = lbegin<0>(i);
      typename Traits::template Codim<0>::LevelIterator eEndIt = lend<0>(i);

      for (; eIt!=eEndIt; ++eIt) {

        typename UG_NS<dim>::Element* elem0 = getRealImplementation(*eIt).target_;

        typename Traits::template Codim<0>::Entity::LevelIntersectionIterator nIt = eIt->ilevelbegin();
        typename Traits::template Codim<0>::Entity::LevelIntersectionIterator nEndIt = eIt->ilevelend();

        for (; nIt!=nEndIt; ++nIt) {

          if (!nIt->neighbor())
            continue;

          typename UG_NS<dim>::Element* elem1 = getRealImplementation(*nIt->outside()).target_;

          int side0 = UGGridRenumberer<dim>::facesDUNEtoUG(nIt->numberInSelf(), eIt->type());
          int side1 = UGGridRenumberer<dim>::facesDUNEtoUG(nIt->numberInNeighbor(), nIt->outside()->type());

          UG::D3::DisposeDoubledSideVector((typename UG_NS<3>::Grid*)multigrid_->grids[i],
                                           (typename UG_NS<3>::Element*)elem0,
                                           side0,
                                           (typename UG_NS<3>::Element*)elem1,
                                           side1);

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
