// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_YASPGRIDREVERSELEVELITERATOR_HH
#define DUNE_GRID_YASPGRIDREVERSELEVELITERATOR_HH

/** \file
 * \brief The YaspLevelIterator class
 */

namespace Dune {


  /** \brief Iterates over entities of one grid level
   */
  template<int codim, PartitionIteratorType pitype, class GridImp>
  class YaspReverseLevelIterator :
    public YaspEntityPointer<codim,GridImp>
  {
    //! know your own dimension
    enum { dim=GridImp::dimension };
    //! know your own dimension of world
    enum { dimworld=GridImp::dimensionworld };
    typedef typename GridImp::ctype ctype;
  public:
    typedef typename GridImp::template Codim<codim>::Entity Entity;
    typedef typename GridImp::YGridLevelIterator YGLI;
    typedef typename GridImp::YGrid::Iterator I;

    //! default constructor
    YaspReverseLevelIterator ()
    {}

    //! constructor
    YaspReverseLevelIterator (const YGLI & g, const I& it) :
      YaspEntityPointer<codim,GridImp>(g,it) {}

    //! copy constructor
    YaspReverseLevelIterator (const YaspReverseLevelIterator& i) :
      YaspEntityPointer<codim,GridImp>(i) {}

    //! increment
    void increment()
    {
      --(GridImp::getRealImplementation(this->_entity)._it);
    }
  };

}

#endif   // DUNE_GRID_YASPGRIDLEVELITERATOR_HH
