// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_INTERSECTIONITERATOR_HH
#define DUNE_ALBERTA_INTERSECTIONITERATOR_HH

#include <dune/grid/common/intersectioniterator.hh>

#include <dune/grid/albertagrid/intersection.hh>

/** \file
 *  \author Martin Nolte
 *  \brief  Implementation of the IntersectionIterator for AlbertaGrid
 */

#if HAVE_ALBERTA

namespace Dune
{

  // AlbertaGridLeafIntersectionIterator
  // -----------------------------------

  template< class GridImp >
  class AlbertaGridLeafIntersectionIterator
  {
    typedef AlbertaGridLeafIntersectionIterator< GridImp > This;

  public:
    typedef Dune::Intersection< GridImp, AlbertaGridLeafIntersection< GridImp > > Intersection;

    static const int dimension = Intersection::Entity::dimension;

    struct Begin {};
    struct End {};

  private:
    typedef AlbertaGridLeafIntersection< GridImp > IntersectionImp;

  public:
    AlbertaGridLeafIntersectionIterator ()
    {}

    template< class EntityImp >
    AlbertaGridLeafIntersectionIterator ( const EntityImp &entity, Begin )
      : intersection_( IntersectionImp( entity, 0 ) )
    {}

    template< class EntityImp >
    AlbertaGridLeafIntersectionIterator ( const EntityImp &entity, End )
      : intersection_( IntersectionImp( entity, dimension+1 ) )
    {}

    AlbertaGridLeafIntersectionIterator ( const This &other )
      : intersection_( other.intersection_.impl() )
    {}

    This &operator= ( const This &other )
    {
      intersection_.impl() = other.intersection_.impl();
      return *this;
    }

    const Intersection &dereference () const
    {
      return intersection_;
    }

    bool equals ( const This &other ) const
    {
      return (intersection_.impl() == other.intersection_.impl());
    }

    void increment ()
    {
      intersection_.impl().next();
    }

  private:
    Intersection intersection_;
  };

}

#endif // #if HAVE_ALBERTA

#endif // #ifndef DUNE_ALBERTA_INTERSECTIONITERATOR_HH
