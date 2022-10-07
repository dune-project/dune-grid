// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_YASPGRIDPERSISTENTCONTAINER_HH
#define DUNE_GRID_YASPGRIDPERSISTENTCONTAINER_HH

/** \file
 * \brief Specialization of the PersistentContainer for YaspGrid
 */

#include <cassert>
#include <vector>

#include <dune/grid/utility/persistentcontainer.hh>
#include <dune/grid/utility/persistentcontainervector.hh>
#include "../yaspgrid.hh"

namespace Dune
{

  /** \internal
   *  \brief implement a consecutive index for all entities of given codim of a YaspGrid
   *
   * this class behaves similar to an IndexSet
   *
   * to make this work for a persistent vector we have to implement...
   *
   *      indexSet.size( codim )
   *      indexSet.index( entity );
   *      indexSet.subIndex( entity, subEntity, codimension() );
   */
  template<typename Grid>
  class YaspPersistentContainerIndex
  {

  public:
    YaspPersistentContainerIndex(const Grid & grid, int codim)
      : _grid(grid), _codim(codim)
    {
      recomputeOffsets();
    }

    /** \copydoc IndexSet::IndexType */
    typedef std::size_t IndexType;

    /** \copydoc IndexSet::index */
    template<class Entity>
    IndexType index (const Entity& e) const
    {
      static const int cc = Entity::codimension;
      std::size_t level = e.level();
      return _grid.indexsets[level]->template index<cc>(e) + _offsets[level];
    }

    /** \copydoc IndexSet::subIndex */
    template< class Entity >
    IndexType subIndex ( const Entity &e, int i, unsigned int codim ) const
    {

      static const int cc = Entity::codimension;
      std::size_t level = e.level();
      return _grid.indexsets[level]->template subIndex<cc>(e,i,codim) + _offsets[level];
    }

    /** \copydoc IndexSet::size */
    std::size_t size (int /* codim */) const
    {
      if (_grid.indexsets.size()+1 != _offsets.size())
        recomputeOffsets();
      return _offsets.back();
    }

  private:
    void recomputeOffsets() const
    {
      _offsets.resize(_grid.indexsets.size()+1,0);
      _offsets[0] = 0;
      for (std::size_t i=0; i<_grid.indexsets.size(); i++)
        _offsets[i+1] = _offsets[i] + _grid.indexsets[i]->size(_codim);
    }

    const Grid& _grid;
    int _codim;
    mutable std::vector<std::size_t> _offsets;
  };

  /** \brief Specialization of the PersistentContainer for YaspGrid */
  /**
     \note questions regarding the interface of PersistentContainer:
     - is it allowed to access the container after grid modification, without calling resize?
     - which assumptions on the indexset doe the PersistentContainerVector have?
   */
  template<int dim, class CoordCont, class T>
  class PersistentContainer< YaspGrid<dim, CoordCont>, T >
  /* We have to pass the reference to the IndexSet to the constructor
     of the PersistentContainerVector.  In order to have a valid
     indexset available, we inherit from a private indexset
   */
    : private YaspPersistentContainerIndex< const YaspGrid<dim, CoordCont> >,
      public PersistentContainerVector< YaspGrid<dim, CoordCont>,
                                        YaspPersistentContainerIndex< const YaspGrid<dim, CoordCont> >,
                                        std::vector<T> >
  {
    typedef YaspPersistentContainerIndex< const YaspGrid<dim, CoordCont> > IndexSet;
    typedef PersistentContainerVector< YaspGrid<dim, CoordCont>, IndexSet, std::vector<T> > Base;

  public:
    typedef typename Base::Grid Grid;
    typedef typename Base::Value Value;

    using Base::size;

    PersistentContainer ( const Grid &grid, int codim, const Value &value = Value() )
      : IndexSet(grid, codim),
        Base(*this, codim, value)
    {}
  };

} // end namespace Dune

#endif // end DUNE_GRID_YASPGRIDPERSISTENTCONTAINER_HH
