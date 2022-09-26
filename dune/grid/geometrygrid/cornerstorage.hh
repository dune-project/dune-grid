// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_CORNERSTORAGE_HH
#define DUNE_GEOGRID_CORNERSTORAGE_HH

#include <array>

#include <dune/grid/geometrygrid/coordfunctioncaller.hh>

namespace Dune
{

  namespace GeoGrid
  {

    // CoordVector
    // -----------

    template< int mydim, class Grid, bool fake >
    class CoordVector;


    template< int mydim, class Grid >
    class CoordVector< mydim, Grid, false >
    {
      typedef typename std::remove_const< Grid >::type::Traits Traits;

      typedef typename Traits::ctype ctype;

      static const int dimension = Traits::dimension;
      static const int mydimension = mydim;
      static const int codimension = dimension - mydimension;
      static const int dimensionworld = Traits::dimensionworld;

      typedef FieldVector< ctype, dimensionworld > Coordinate;

      typedef typename Traits::HostGrid HostGrid;
      typedef typename Traits::CoordFunction CoordFunction;

      typedef typename HostGrid::template Codim< codimension >::Entity HostEntity;

      typedef GeoGrid :: CoordFunctionCaller< HostEntity, typename CoordFunction::Interface >
      CoordFunctionCaller;

    public:
      CoordVector ( const HostEntity &hostEntity,
                    const CoordFunction &coordFunction )
        : coordFunctionCaller_( hostEntity, coordFunction )
      {}

      template< std::size_t size >
      void calculate ( std::array< Coordinate, size > (&corners) ) const
      {
        const std::size_t numCorners = coordFunctionCaller_.size();
        assert( size >= numCorners );
        for( std::size_t i = 0; i < numCorners; ++i )
          coordFunctionCaller_.evaluate( i, corners[ i ] );
      }

    private:
      const CoordFunctionCaller coordFunctionCaller_;
    };


    template< int mydim, class Grid >
    class CoordVector< mydim, Grid, true >
    {
      typedef typename std::remove_const< Grid > :: type :: Traits Traits;

      typedef typename Traits::ctype ctype;

      static const int dimension = Traits::dimension;
      static const int mydimension = mydim;
      static const int codimension = dimension - mydimension;
      static const int dimensionworld = Traits::dimensionworld;

      typedef FieldVector< ctype, dimensionworld > Coordinate;

      typedef typename Traits::HostGrid HostGrid;
      typedef typename Traits::CoordFunction CoordFunction;

      typedef typename HostGrid::template Codim< 0 >::Entity HostElement;

      typedef GeoGrid::CoordFunctionCaller< HostElement, typename CoordFunction::Interface >
      CoordFunctionCaller;

    public:
      CoordVector ( const HostElement &hostElement,
                    const unsigned int subEntity,
                    const CoordFunction &coordFunction )
        : coordFunctionCaller_( hostElement, coordFunction ),
          subEntity_( subEntity )
      {}

      template< std::size_t size >
      void calculate ( std::array< Coordinate, size > (&corners) ) const
      {
        const GeometryType type = coordFunctionCaller_.type();
        auto refElement = referenceElement< ctype, dimension >( type );
        const std::size_t numCorners = refElement.size( subEntity_, codimension, dimension );
        assert( size >= numCorners );
        for( std::size_t i = 0; i < numCorners; ++i )
        {
          const std::size_t j = refElement.subEntity( subEntity_, codimension, i, dimension );
          coordFunctionCaller_.evaluate( j, corners[ i ] );
        }
      }

    private:
      const CoordFunctionCaller coordFunctionCaller_;
      const unsigned int subEntity_;
    };



    // IntersectionCoordVector
    // -----------------------

    template< class Grid >
    class IntersectionCoordVector
    {
      typedef typename std::remove_const< Grid >::type::Traits Traits;

      typedef typename Traits::ctype ctype;

      static const int dimension = Traits::dimension;
      static const int codimension = 1;
      static const int mydimension = dimension-codimension;
      static const int dimensionworld = Traits::dimensionworld;

      typedef FieldVector< ctype, dimensionworld > Coordinate;

      typedef typename Traits::HostGrid HostGrid;

      typedef typename Traits::template Codim< 0 >::GeometryImpl ElementGeometryImpl;
      typedef typename Traits::template Codim< codimension >::LocalGeometry HostLocalGeometry;

    public:
      IntersectionCoordVector ( const ElementGeometryImpl &elementGeometry,
                                const HostLocalGeometry &hostLocalGeometry )
        : elementGeometry_( elementGeometry ),
          hostLocalGeometry_( hostLocalGeometry )
      {}

      template< std::size_t size >
      void calculate ( std::array< Coordinate, size > (&corners) ) const
      {
        const std::size_t numCorners = hostLocalGeometry_.corners();
        assert( size >= numCorners );
        for( std::size_t i = 0; i < numCorners; ++i )
          corners[ i ] = elementGeometry_.global( hostLocalGeometry_.corner( i ) );
      }

      template< unsigned int numCorners >
      void calculate ( Coordinate (&corners)[ numCorners ] ) const
      {
        assert( numCorners == hostLocalGeometry_.corners() );
      }

    private:
      const ElementGeometryImpl &elementGeometry_;
      HostLocalGeometry hostLocalGeometry_;
    };



    // CornerStorage
    // -------------

    template< int mydim, int cdim, class Grid >
    class CornerStorage
    {
      typedef typename std::remove_const< Grid >::type::Traits Traits;

      typedef typename Traits::ctype ctype;
      typedef FieldVector< ctype, cdim > Coordinate;

      typedef std::array< Coordinate, (1 << mydim) > Coords;

    public:
      typedef typename Coords::const_iterator const_iterator;

      template< bool fake >
      explicit CornerStorage ( const CoordVector< mydim, Grid, fake > &coords )
      {
        coords.calculate( coords_ );
      }

      explicit CornerStorage ( const IntersectionCoordVector< Grid > &coords )
      {
        coords.calculate( coords_ );
      }

      const Coordinate &operator[] ( unsigned int i ) const
      {
        return coords_[ i ];
      }

      const_iterator begin () const { return coords_.begin(); }
      const_iterator end () const { return coords_.end(); }

    private:
      Coords coords_;
    };

  } // namespace GeoGrid

} // namespace Dune

#endif // #ifndef DUNE_GEOGRID_CORNERSTORAGE_HH
