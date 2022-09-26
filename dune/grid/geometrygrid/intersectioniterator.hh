// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_INTERSECTIONITERATOR_HH
#define DUNE_GEOGRID_INTERSECTIONITERATOR_HH

#include <dune/grid/geometrygrid/intersection.hh>

namespace Dune
{

  namespace GeoGrid
  {

    // IntersectionIterator
    // --------------------

    template< class Grid, class HostIntersectionIterator >
    class IntersectionIterator
    {
      typedef typename std::remove_const< Grid >::type::Traits Traits;

      typedef GeoGrid::Intersection< Grid, typename HostIntersectionIterator::Intersection > IntersectionImpl;

      typedef typename Traits::template Codim< 0 >::Geometry ElementGeometry;
      typedef typename Traits::template Codim< 0 >::GeometryImpl ElementGeometryImpl;

    public:
      typedef Dune::Intersection< Grid, IntersectionImpl > Intersection;

      IntersectionIterator()
      {}

      template< class Entity >
      IntersectionIterator ( const Entity &inside,
                             const HostIntersectionIterator &hostIterator )
        : hostIterator_( hostIterator )
        , insideGeo_( inside.geometry().impl() )
      {}

      IntersectionIterator ( const IntersectionIterator &other )
        : hostIterator_( other.hostIterator_ )
        , insideGeo_( other.insideGeo_ )
      {}

      IntersectionIterator ( IntersectionIterator&& other )
        : hostIterator_( std::move( other.hostIterator_ ) )
        , insideGeo_( std::move( other.insideGeo_ ) )
      {}

      IntersectionIterator &operator= ( const IntersectionIterator &other )
      {
        hostIterator_ = other.hostIterator_;
        insideGeo_ = other.insideGeo_;
        return *this;
      }

      IntersectionIterator &operator= ( IntersectionIterator&& other )
      {
        hostIterator_ = std::move( other.hostIterator_ );
        insideGeo_ = std::move( other.insideGeo_ );
        return *this;
      }

      bool equals ( const IntersectionIterator &other ) const
      {
        return (hostIterator_ == other.hostIterator_);
      }

      void increment ()
      {
        ++hostIterator_;
      }

      Intersection dereference () const
      {
        return IntersectionImpl( *hostIterator_, insideGeo_ );
      }

    private:

      HostIntersectionIterator hostIterator_;
      ElementGeometryImpl insideGeo_;

    };

  } // namespace GeoGrid

} // namespace Dune

#endif // #ifndef DUNE_GEOGRID_INTERSECTIONITERATOR_HH
