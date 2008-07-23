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
      typedef DefaultCoordTraits<ctype,gdim> LocalCoordTraits;
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



    // MappingProvider
    // ---------------

    template< class ElementMapping, unsigned int codim >
    struct MappingProvider
    {
      typedef GenericGeometry :: SubMappingTraits< ElementMapping, codim > SubMappingTraits;
      typedef typename SubMappingTraits :: SubMapping Mapping;
      typedef typename SubMappingTraits :: CachingType CachingType;

      template< GeometryType :: BasicType type, class CoordVector >
      static Mapping *virtualMapping ( const CoordVector &coords,
                                       const CachingType &cache )
      {
        typedef typename SubMappingTraits :: template VirtualMapping< type > :: type
        VirtualMapping;
        return new VirtualMapping( coords, cache );
      }

      template< bool > struct Virtual;
      template< bool > struct NonVirtual;

      template< class CoordVector >
      static Mapping *mapping ( const GeometryType &type,
                                const CoordVector &coords,
                                const CachingType &cache )
      {
        typedef ProtectedIf< SubMappingTraits :: isVirtual, Virtual, NonVirtual > Switch;
        return Switch :: mapping( type, coords, cache );
      }
    };

    template< class ElementMapping, unsigned int codim >
    template< bool >
    struct MappingProvider< ElementMapping, codim > :: Virtual
    {
      template< class CoordVector >
      static Mapping *
      mapping ( const GeometryType &type, const CoordVector &coords, const CachingType &cache )
      {
        assert( type.dim() == Mapping :: dimG );

        switch( type.basicType() )
        {
        case GeometryType :: simplex :
          return virtualMapping< GeometryType :: simplex, CoordVector >( coords, cache );

        case GeometryType :: cube :
          return virtualMapping< GeometryType :: cube, CoordVector >( coords, cache );
        /*
           case GeometryType :: prism:
           return virtualMapping< GeometryType :: prism, CoordVector >( coords, cache );

           case GeometryType :: pyramid:
           return virtualMapping< GeometryType :: pyramid, CoordVector >( coords, cache );
         */
        default :
          DUNE_THROW( RangeError, "Unknown basic geometry type: " << type.basicType() );
        }
      }
    };

    template< class ElementMapping, unsigned int codim >
    template< bool >
    struct MappingProvider< ElementMapping, codim > :: NonVirtual
    {
      template< class CoordVector >
      static Mapping *
      mapping ( const GeometryType &type, const CoordVector &coords, const CachingType &cache )
      {
        assert( type.dim() == Mapping :: dimG );
        return new Mapping( coords, cache );
      }
    };



    // Geometry
    // --------

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
      dune_static_assert( (0 <= mydimension) && (mydimension <= int(dimGrid)),
                          "Invalid geometry dimension." );
      enum { codimension = dimGrid - mydimension };

      typedef typename Traits :: CoordTraits CoordTraits;

      template< bool >
      struct Hybrid
      {
        typedef HybridMapping< dimGrid, CoordTraits, Traits :: template Caching >
        Mapping;
      };

      template< bool >
      struct NonHybrid
      {
        typedef typename Convert< GeometryType::BasicType(Traits :: dunetype), dimGrid > :: type Topology;
        typedef GenericGeometry :: CachedMapping< Topology, CoordTraits, Traits :: template Caching >
        Mapping;
      };

      typedef typename ProtectedIf< Traits :: hybrid, Hybrid, NonHybrid > :: Mapping
      ElementMapping;
      typedef GenericGeometry :: MappingProvider< ElementMapping, codimension > MappingProvider;
      typedef typename MappingProvider :: Mapping Mapping;
      typedef typename MappingProvider :: CachingType CachingType;

      mutable Mapping *mapping_;

    public:
      explicit Geometry ( Mapping &mapping )
        : mapping_( &mapping )
      {}

      template< class GeoType >
      explicit Geometry ( const GeoType &coords,
                          const CachingType &cache = CachingType() )
        : mapping_( MappingProvider :: mapping( coords.type(), coords, cache ) )
      {}
      template< class CoordVector >
      Geometry ( const GeometryType &type,
                 const CoordVector &coords,
                 const CachingType &cache = CachingType() )
        : mapping_( MappingProvider :: mapping( type, coords, cache ) )
      {}

      template< int fatherdim >
      Geometry ( const Geometry< fatherdim, cdim, GridImp > &father, int i,
                 const CachingType &cache = CachingType() )
        : mapping_( father.mapping().template subMapping< fatherdim-mydim >( i, cache ) )
      {}

      Geometry ( const Geometry &other )
        : mapping_( other.mapping_ )
      {
        other.mapping_ = 0;
      }

      ~Geometry ()
      {
        if( mapping_ != 0 )
          delete mapping_;
      }

      GeometryType type () const
      {
        return mapping().type();
      }

      int corners () const
      {
        return mapping().corners();
      }

      const GlobalCoordinate &operator[] ( int i ) const
      {
        return mapping()[ i ];
      }

      GlobalCoordinate global ( const LocalCoordinate &local ) const
      {
        return mapping().global( local );
      }

      LocalCoordinate local ( const GlobalCoordinate &global ) const
      {
        return mapping().local( global );
      }

      bool checkInside ( const LocalCoordinate &local ) const
      {
        return mapping().checkInside( local );
      }

      bool affine () const
      {
        return mapping().affine();
      }

      ctype integrationElement ( const LocalCoordinate &local ) const
      {
        return mapping().integrationElement( local );
      }

      ctype volume () const
      {
        return mapping().volume();
      }

      const Jacobian &jacobianInverseTransposed ( const LocalCoordinate &local ) const
      {
        return mapping().jacobianInverseTransposed( local );
      }

      GlobalCoordinate normal ( int face, const LocalCoordinate &local ) const
      {
        return mapping().normal( face , local );
      }

      template <int,int,class>
      friend class Geometry;

    private:
      const Mapping &mapping () const
      {
        return *mapping_;
      }
    };
    template< int mydim, int cdim, class GridImp >
    class LocalGeometry :
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
      dune_static_assert( (0 <= mydimension) && (mydimension <= int(dimGrid)),
                          "Invalid geometry dimension." );
      enum { codimension = dimGrid - mydimension };

      typedef typename Traits :: LocalCoordTraits CoordTraits;

      template< bool >
      struct Hybrid
      {
        typedef HybridMapping< dimGrid, CoordTraits, Traits :: template Caching >
        Mapping;
      };

      template< bool >
      struct NonHybrid
      {
        typedef typename Convert< GeometryType::BasicType(Traits :: dunetype), dimGrid > :: type Topology;
        typedef GenericGeometry :: CachedMapping< Topology, CoordTraits, Traits :: template Caching >
        Mapping;
      };

      typedef typename ProtectedIf< Traits :: hybrid, Hybrid, NonHybrid > :: Mapping
      ElementMapping;
      typedef GenericGeometry :: MappingProvider< ElementMapping, codimension > MappingProvider;
      typedef typename MappingProvider :: Mapping Mapping;
      typedef typename MappingProvider :: CachingType CachingType;

      mutable Mapping *mapping_;

    public:
      explicit LocalGeometry ( Mapping &mapping )
        : mapping_( &mapping )
      {}

      template< class GeoType >
      explicit LocalGeometry ( const GeoType &coords,
                               const CachingType &cache = CachingType() )
        : mapping_( MappingProvider :: mapping( coords.type(), coords, cache ) )
      {}
      template< class CoordVector >
      LocalGeometry ( const GeometryType &type,
                      const CoordVector &coords,
                      const CachingType &cache = CachingType() )
        : mapping_( MappingProvider :: mapping( type, coords, cache ) )
      {}

      template< int fatherdim >
      LocalGeometry ( const Geometry< fatherdim, cdim, GridImp > &father, int i,
                      const CachingType &cache = CachingType() )
        : mapping_( father.mapping().template subMapping< fatherdim-mydim >( i, cache ) )
      {}

      LocalGeometry ( const LocalGeometry &other )
        : mapping_( other.mapping_ )
      {
        other.mapping_ = 0;
      }

      ~LocalGeometry ()
      {
        if( mapping_ != 0 )
          delete mapping_;
      }

      GeometryType type () const
      {
        return mapping().type();
      }

      int corners () const
      {
        return mapping().corners();
      }

      const GlobalCoordinate &operator[] ( int i ) const
      {
        return mapping()[ i ];
      }

      GlobalCoordinate global ( const LocalCoordinate &local ) const
      {
        return mapping().global( local );
      }

      LocalCoordinate local ( const GlobalCoordinate &global ) const
      {
        return mapping().local( global );
      }

      bool checkInside ( const LocalCoordinate &local ) const
      {
        return mapping().checkInside( local );
      }

      bool affine () const
      {
        return mapping().affine();
      }

      ctype integrationElement ( const LocalCoordinate &local ) const
      {
        return mapping().integrationElement( local );
      }

      ctype volume () const
      {
        return mapping().volume();
      }

      const Jacobian &jacobianInverseTransposed ( const LocalCoordinate &local ) const
      {
        return mapping().jacobianInverseTransposed( local );
      }

      GlobalCoordinate normal ( int face, const LocalCoordinate &local ) const
      {
        return mapping().normal( face , local );
      }

      template <int,int,class>
      friend class Geometry;

    private:
      const Mapping &mapping () const
      {
        return *mapping_;
      }
    };
  }

}

#endif
