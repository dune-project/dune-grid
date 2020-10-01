// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_YASPGRIDLEVELITERATOR_HH
#define DUNE_GRID_YASPGRIDLEVELITERATOR_HH

/** \file
 * \brief The YaspLevelIterator class
 */

namespace Dune {


  /** \brief Iterates over entities of one grid level
   */
  template<int codim, PartitionIteratorType pitype, class GridImp>
  class YaspLevelIterator
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
    YaspLevelIterator ()
    {}

    //! constructor
    YaspLevelIterator (const YGLI & g, const I& it)
      : _entity(YaspEntity<codim, dim, GridImp>(g,it))
    {}

    //! increment
    void increment()
    {
      ++(_entity.impl()._it);
    }

    //! equality
    bool equals (const YaspLevelIterator& rhs) const
    {
      return (_entity == rhs._entity);
    }

    //! dereferencing
    const Entity& dereference() const
    {
      return _entity;
    }

  protected:
    Entity _entity; //!< entity
  };

}

#endif   // DUNE_GRID_YASPGRIDLEVELITERATOR_HH
