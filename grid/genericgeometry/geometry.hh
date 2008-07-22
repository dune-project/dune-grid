// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_GEOMETRY_HH
#define DUNE_GENERICGEOMETRY_GEOMETRY_HH

#include <dune/grid/genericgeometry/mappings.hh>
#include <dune/grid/genericgeometry/submapping.hh>
#include <dune/grid/genericgeometry/hybridmapping.hh>
#include <dune/grid/common/geometry.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    template< class ctype, int cdim,
        bool alwaysAffine = false,
        GeometryType::BasicType oneD = GeometryType:: simplex>
    struct DefaultCoordTraits
    {
      typedef ctype FieldType;
      enum { dimCoord = cdim };

      template< int dim >
      struct Vector
      {
        typedef FieldVector< FieldType, dim > Type;
      };

      template< int dimR, int dimC >
      struct Matrix
      {
        typedef FieldMatrix< FieldType, dimR, dimC > Type;
      };

      typedef typename Vector< dimCoord >::Type CoordinateType;

      enum {affine = alwaysAffine};
      enum {oneDType = oneD};
    };

    template< class GridImp >
    struct GeometryTraits;

    template< class GridImp >
    struct GeometryTraits< const GridImp >
      : public GeometryTraits< GridImp >
    {};

    template <class ctype,int gdim,int cdim>
    struct DefaultGeometryTraits {
      typedef DefaultCoordTraits<ctype,cdim> CoordTraits;
      template <class Traits>
      struct Caching {
        enum {jTCompute = geoCompute,
              jTInvCompute = geoCompute,
              intElCompute = geoCompute,
              normalCompute = geoCompute};
        void jacobianT(typename Traits::JacobianTransposedType& d) const {}
        void integrationElement(typename Traits::FieldType& intEl) const {}
        void jacobianInverseTransposed(typename Traits::JacobianType& dInv) const {}
        void normal(int face, typename Traits::GlobalCoordType& n) const {}
      };
      enum {dimGrid = gdim};
      //   hybrid   [ true if Codim 0 is hybrid ]
      enum {hybrid  = true};
      //   dunetype [ for Codim 0, needed for (hybrid=false) ]
      // enum {dunetype = GeometryType::simplex};
    };

    /*
       template <int dimGrid>
       struct GeometryTraits<MyGrid<dimGrid> > :
       public DefaultGeometryTraits<MyGrid<dimGrid>::ctype,
                                   dimGrid, dimGrid>
       {};
     */

    template< int mydim, int cdim, class GridImp >
    class Geometry :
      public GeometryDefaultImplementation <mydim, cdim, GridImp, Geometry>
    {
      typedef GeometryTraits< GridImp > Traits;

      enum { dimGrid = Traits :: dimGrid };

    public:
      enum { mydimension = mydim };
      enum { coorddimension = cdim };

      typedef typename GridImp :: ctype ctype;

      typedef FieldVector< ctype, mydimension > LocalCoordinate;
      typedef FieldVector< ctype, coorddimension > GlobalCoordinate;
      typedef FieldMatrix< ctype, coorddimension,mydimension > Jacobian;

    private:
      dune_static_assert( (0 <= mydimension) && (mydimension <= dimGrid),
                          "Invalid geometry dimension." );
      enum { codimension = dimGrid - mydimension };

      typedef typename Traits :: CoordTraits CoordTraits;

      template< bool >
      struct Hybrid
      {
        typedef HybridMapping< dimGrid, CoordTraits, Traits :: template Caching >
        GeometryType;
      };

      template< bool >
      struct NonHybrid
      {
        typedef typename Convert< Traits :: dunetype, dimGrid > :: type Topology;
        typedef CachedMapping< Topology, CoordTraits, Traits :: template Caching >
        GeometryType;
      };

      typedef typename ProtectedIf< Traits :: hybrid, Hybrid, NonHybrid > :: GeometryType
      ElementGeometryType;
      typedef typename ElementGeometryType :: template Codim <codimension> :: SubMapping
      GeometryType;

      mutable GeometryType *geometry_;

    public:
      template< class CoordVector >
      explicit Geometry ( const CoordVector& coord)
        : geometry_( new GeometryType(coord) )
      {}
      explicit Geometry ( GeometryType& geometry)
        : geometry_( &geometry )
      {}

      /*
         template< class GeoClass >
         explicit Geometry ( const GeoClass &geo,
                          const CachingType &cache = CachingType() )
         : geometry_( GeometryProvider :: geometry( geo,geo.type(),cache ) ) // ???
         {}
       */

      /*
         typedef typename Traits::Caching Caching;
         template< int fatherdim >
         explicit Geometry ( const Geometry<fatherdim,cdim,GridImp> &father,
                          int i,
                          const Caching &cache = Caching() )
         : geometry_( father.geometry().subMapping<fatherdim-mydim>(i,cache) )
         {}
       */

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

      template <int,int,class>
      friend class Geometry;

    private:
      const GeometryType &geometry () const
      {
        return *geometry_;
      }
    };
  }

}

#endif
