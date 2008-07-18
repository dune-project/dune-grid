// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_GEOMETRY_HH
#define DUNE_GENERICGEOMETRY_GEOMETRY_HH

#include <dune/grid/genericgeometry/mappings.hh>
#include <dune/grid/genericgeometry/subgeometry.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    template< class ctype, int cdim >
    struct CoordTraitsBase
    {
      typedef typename ctype FieldType;
      enum { dimCoord = cdim };

      template< int dim >
      struct Vector
      {
        typedef FieldVector< FieldType, dim > Type;
      };

      template< int dimR, int dim C >
      struct Matrix
      {
        typedef FieldMatrix< FieldType, dimR, dimC > Type;
      };

      typedef typename Vector< dimCoord > CoordinateType;
    };



    template< class GridImp >
    struct GeometryTraits;

    template< class GridImp >
    struct GeometryTraits< const GridImp >
      : public GeometryTraits< GridImp >
    {};



    template< int mydim, int cdim, class GridImp >
    class Geometry
    {
      typedef GeometryTraits< GridImp > Traits;

      enum { dimGrid = Traits :: dimGrid };

    public:
      enum { mydimension = mydim };
      enum { coorddimension = cdim };

      typedef typename GridImp :: ctype ctype;

      typedef FieldVector< ctype, mydimension > LocalCoordinate;
      typedef FieldVector< ctype, coorddimension > GlobalCoordinate;

    private:
      dune_static_assert( (0 <= mydimension) && (mydimension <= dimGrid),
                          "Invalid geometry dimension." );
      enum { codimension = dimGrid - mydimension };

      typedef typename Traits :: CoordTraits CoordTraits;

      template< bool >
      struct Hybrid
      {
        typedef HybridGeometry< dimGrid, CoordTraits, Traits :: template Caching >
        GeometryType;
      };

      template< bool >
      struct NonHybrid
      {
        typedef typename Convert< Traits :: dunetype, dimGrid > :: type Topology;
        typedef Geometry< Topology, CoordTraits, Traits :: template Caching >
        GeometryType;
      };

      typedef ProtectedIf< Traits :: hybrid, Hybrid, NonHybrid > :: GeometryType
      ElementGeometryType;
      typedef GenericGeometry :: GeometryProvider< ElementGeometryType, codimension >
      GeometryProvider;
      typedef typename GeometryProvider :: Geometry GeometryType;

      mutable GeometryType *geometry_;

    public:
      template< class CoordVector >
      explicit Geometry ( const CoordVector &coords,
                          const CachingType &cache = CachingType() )
        : geometry_( GeometryProvider :: geometry( coords, cache ) )
      {}

      Geometry ( const Geometry &other )
        : geometry_( other.geometry_ )
      {
        other.geometry_ = 0;
      }

      ~Geometry ()
      {
        if( geometry_ != 0 )
          delete geometry_;
      }

      Dune :: GeometryType type () const
      {
        return geometry().type();
      }

      int corners () const
      {
        return geometry().corners();
      }

      const GlobalCoordinate &operator[] ( int i ) const
      {
        return geometry()[ i ];
      }

      GlobalCoordinate global ( const LocalCoordinate &local ) const
      {
        return geometry().global( local );
      }

      LocalCoordinate local ( const GlobalCoordinate &global ) const
      {
        return geometry().local( global );
      }

      bool checkInside ( const LocalCoordinate &local ) const
      {
        return geometry().checkInside( local );
      }

      bool affine () const
      {
        return geometry().affine();
      }

      ctype integrationElement ( const LocalCoordinate &local ) const
      {
        return geometry().integrationElement( local );
      }

      ctype volume () const
      {
        return geometry().volume();
      }

      const Jacobian &jacobianInverseTransposed ( const LocalCoordinate &local ) const
      {
        return geometry().jacobianInverseTransposed( local );
      }

    private:
      const GeometryType &geometry () const
      {
        return *geometry_;
      }
    };
  }

}

#endif
