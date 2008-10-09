// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_GEOMETRY_HH
#define DUNE_GENERICGEOMETRY_GEOMETRY_HH

#include <dune/grid/genericgeometry/mappingprovider.hh>
#include <dune/grid/common/geometry.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    template< class ct, int dimW >
    struct DuneCoordTraits
    {
      typedef ct ctype;

      static const int dimWorld = dimW;

      template< int dim >
      struct Vector
      {
        typedef FieldVector< ctype, dim > type;
      };

      template< int rows, int cols >
      struct Matrix
      {
        typedef FieldMatrix< ctype, rows, cols > type;
      };

      typedef typename Vector< dimWorld > :: type GlobalCoordinate;
    };



    template< class ct, int dimW, bool alwaysAffine = false >
    struct DefaultCoordTraits
      : public DuneCoordTraits< ct, dimW >
    {
      static const bool affine = alwaysAffine;

      typedef typename DuneCoordTraits< ct, dimW > :: GlobalCoordinate
      GlobalCoordinate;

      template< class Topology >
      struct CornerStorage
      {
        typedef CoordPointerStorage< Topology, GlobalCoordinate > Type;
      };
    };



    template< class ctype, int dimG, int dimW >
    struct DefaultGeometryTraits
    {
      typedef DefaultCoordTraits< ctype, dimW > GlobalCoordTraits;
      typedef DefaultCoordTraits< ctype, dimG > LocalCoordTraits;

      // are all global geometries affine?
      static const bool globallyAffine = false;
      // are all local geometries affine?
      static const bool locallyAffine = false;

      static const int dimGrid = dimG;

      //   hybrid   [ true if Codim 0 is hybrid ]
      static const bool hybrid = true;
      //   dunetype [ for Codim 0, needed for (hybrid=false) ]
      // static const GeometryType :: BasicType dunetype = GeometryType :: simplex;

      // what basic geometry type shall the line be considered?
      static const GeometryType :: BasicType linetype = GeometryType :: simplex;

      template< class Traits >
      struct Caching
        : public ComputeAll< Traits >
      {};
    };



    /** \struct  GeometryTraits
     *  \ingroup GenericGeometry
     *  \brief   grid specific information required by the generic geometries
     *
     *  Every implementation of a DUNE Geometry is required to have the same
     *  template parameter list:
     *  \code
     *  template< int mydim, int cdim, class Grid >
     *  \endcode
     *  Consequently, there is no direct way to pass compile time static
     *  information to a unified implementation such as the generic geometries.
     *  The structure GeometryTraits realizes an indirect way to do this.
     *
     *  For every grid implementation using the generic geometries, this
     *  structure must be specialized. The following default implementation
     *  can be used (via derivation) to provide the necessary information. It
     *  contains exactly the fields that are necessary:
     *  \code
     *  template< class ctype, int dimG, int dimW >
     *  struct DefaultGeometryTraits
     *  {
     *    typedef DefaultCoordTraits< ctype, dimW > GlobalCoordTraits;
     *    typedef DefaultCoordTraits< ctype, dimG > LocalCoordTraits;
     *
     *    // are all global geometries affine?
     *    static const bool globallyAffine = false;
     *    // are all local geometries affine?
     *    static const bool locallyAffine = false;
     *
     *    static const int dimGrid = dimG;
     *
     *    // hybrid   [ true if Codim 0 is hybrid ]
     *    static const bool hybrid = true;
     *    // dunetype [ for Codim 0, needed for (hybrid=false) ]
     *    // static const GeometryType :: BasicType dunetype = GeometryType :: simplex;
     *
     *    // what basic geometry type shall the line be considered?
     *    static const GeometryType :: BasicType linetype = GeometryType :: simplex;
     *
     *    template< class Traits >
     *    struct Caching
     *    : public ComputeAll< Traits >
     *    {};
     *  };
     *  \endcode
     */
    template< class Grid >
    struct GeometryTraits;

    template< class Grid >
    struct GeometryTraits< const Grid >
      : public GeometryTraits< Grid >
    {};



    template< class CoordTraits, template< class > class Caching >
    struct CornerMappingTraits
    {
      typedef CoordTraits CoordinateTraits;

      template< unsigned int dimension >
      struct Traits
        : public MappingTraits< dimension, CoordTraits >
      {
        typedef Caching< MappingTraits< dimension, CoordTraits > > CachingType;
      };

      template< class Topology >
      struct Mapping
      {
        typedef typename CoordTraits :: template CornerStorage< Topology > :: Type
        CornerStorage;
        typedef CornerMapping< Topology, CoordTraits, CornerStorage > Type;
      };
    };



    // BasicGeometry
    // -------------

    template< int mydim, int cdim, class Grid, class CoordTraits,
        GeometryType :: BasicType linetype >
    class BasicGeometry
    {
      typedef GeometryTraits< Grid > Traits;

      static const int dimGrid = Traits :: dimGrid;

      template< int, int, class, class, GeometryType :: BasicType >
      friend class BasicGeometry;

    public:
      static const int mydimension = mydim;
      static const int coorddimension = cdim;

      typedef typename Grid :: ctype ctype;

      typedef FieldVector< ctype, mydimension > LocalCoordinate;
      typedef FieldVector< ctype, coorddimension > GlobalCoordinate;
      typedef FieldMatrix< ctype, coorddimension, mydimension > Jacobian;

    private:
      dune_static_assert( (0 <= mydimension) && (mydimension <= dimGrid),
                          "Invalid geometry dimension." );

      static const int codimension = dimGrid - mydimension;

      typedef CornerMappingTraits< CoordTraits, Traits :: template Caching >
      GeometricMappingTraits;

      template< bool >
      struct Hybrid
      {
        typedef HybridMapping< dimGrid, GeometricMappingTraits > Mapping;
      };

      template< bool >
      struct NonHybrid
      {
        typedef typename Convert< Traits :: dunetype, dimGrid > :: type Topology;
        typedef GenericGeometry :: CachedMapping< Topology, GeometricMappingTraits >
        Mapping;
      };

      typedef GenericGeometry :: DuneGeometryTypeProvider< mydimension, linetype >
      DuneGeometryTypeProvider;

      typedef typename ProtectedIf< Traits :: hybrid, Hybrid, NonHybrid > :: Mapping
      ElementMapping;
      typedef GenericGeometry :: MappingProvider< ElementMapping, codimension > MappingProvider;

    protected:
      typedef typename MappingProvider :: Mapping Mapping;
      typedef typename MappingProvider :: CachingType CachingType;

    private:
      Mapping *mapping_;

    public:
      BasicGeometry ()
        : mapping_( 0 )
      {}

#if 0
      explicit BasicGeometry ( Mapping &mapping )
        : mapping_( &mapping )
      {
        ++mapping_->referenceCount;
      }
#endif

      template< class CoordVector >
      BasicGeometry ( const GeometryType &type,
                      const CoordVector &coords,
                      const CachingType &cache )
        : mapping_( MappingProvider :: mapping( type, coords, cache ) )
      {
        mapping_->referenceCount = 1;
      }

      template< int fatherdim >
      BasicGeometry ( const BasicGeometry< fatherdim, cdim, Grid, CoordTraits, linetype > &father,
                      int i,
                      const CachingType &cache )
        : mapping_( subMapping( father, i, cache ) )
      {
        mapping_->referenceCount = 1;
      }

      BasicGeometry ( const BasicGeometry &other )
        : mapping_( other.mapping_ )
      {
        if( mapping_ != 0 )
          ++(mapping_->referenceCount);
      }

      ~BasicGeometry ()
      {
        if( (mapping_ != 0) && ((--mapping_->referenceCount) == 0) )
          delete mapping_;
      }

      BasicGeometry &operator= ( const BasicGeometry &other )
      {
        if( other.mapping_ != 0 )
          ++(other.mapping_->referenceCount);
        if( (mapping_ != 0) && (--(mapping_->referenceCount) == 0) )
          delete mapping_;
        mapping_ = other.mapping_;
        return *this;
      }

    public:
      bool operator! () const
      {
        return (mapping_ == 0);
      }

      GeometryType type () const
      {
        return DuneGeometryTypeProvider :: type( mapping().topologyId() );
      }

      int corners () const
      {
        return mapping().numCorners();
      }

      const GlobalCoordinate &operator[] ( int i ) const
      {
        return mapping().corner( i );
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
        const unsigned int tid = mapping().topologyId();
        const unsigned int i = MapNumberingProvider< mydimension >
                               :: template dune2generic< 1 >( tid, face );
        return mapping().normal( i, local );
      }

    private:
      const Mapping &mapping () const
      {
        assert( mapping_ != 0 );
        return *mapping_;
      }

      template< int fatherdim >
      Mapping *
      subMapping ( const BasicGeometry< fatherdim, cdim, Grid, CoordTraits, linetype > &father,
                   int i, const CachingType &cache )
      {
        const unsigned int codim = fatherdim - mydim;
        const unsigned int ftid = father.mapping().topologyId();
        const unsigned int j = MapNumberingProvider< fatherdim >
                               :: template dune2generic< codim >( ftid, i );
        return father.mapping().template subMapping< codim >( j, cache );
      }
    };



    // Geometry
    // --------

    template< int mydim, int cdim, class Grid >
    class Geometry
      : public BasicGeometry
        < mydim, cdim, Grid,
            typename GeometryTraits< Grid > :: GlobalCoordTraits,
            //GeometryTraits< Grid > :: globallyAffine,
            GeometryTraits< Grid > :: linetype >
    {
      typedef typename GeometryTraits< Grid > :: GlobalCoordTraits CoordTraits;
      typedef BasicGeometry
      < mydim, cdim, Grid, CoordTraits, GeometryTraits< Grid > :: linetype >
      Base;

    protected:
      typedef typename Base :: CachingType CachingType;
      typedef typename Base :: Mapping Mapping;

    public:
      Geometry ()
        : Base()
      {}

      explicit Geometry ( Mapping &mapping )
        : Base( mapping )
      {}

      template< class Geo >
      explicit Geometry ( const Geo &geo,
                          const CachingType &cache = CachingType() )
        : Base( geo.type(), geo, cache )
      {}

      template< class CoordVector >
      Geometry ( const GeometryType &type,
                 const CoordVector &coords,
                 const CachingType &cache = CachingType() )
        : Base( type, coords, cache )
      {}

      template< int fatherdim >
      Geometry ( const Geometry< fatherdim, cdim, Grid > &father, int i,
                 const CachingType &cache = CachingType() )
        : Base( father, i, cache )
      {}
    };



    // LocalGeometry
    // -------------

    template< int mydim, int cdim, class Grid >
    class LocalGeometry
      : public BasicGeometry
        < mydim, cdim, Grid, typename GeometryTraits< Grid > :: LocalCoordTraits,
            // GeometryTraits< Grid > :: locallyAffine,
            GeometryTraits< Grid > :: linetype >
    {
      typedef typename GeometryTraits< Grid > :: LocalCoordTraits CoordTraits;
      typedef BasicGeometry
      < mydim, cdim, Grid, CoordTraits, GeometryTraits< Grid > :: linetype >
      Base;

    protected:
      typedef typename Base :: CachingType CachingType;
      typedef typename Base :: Mapping Mapping;

    public:
      LocalGeometry ()
        : Base()
      {}

      explicit LocalGeometry ( Mapping &mapping )
        : Base( mapping )
      {}

      template< class Geo >
      explicit LocalGeometry ( const Geo &geo,
                               const CachingType &cache = CachingType() )
        : Base( geo.type(), geo, cache )
      {}

      template< class CoordVector >
      LocalGeometry ( const GeometryType &type,
                      const CoordVector &coords,
                      const CachingType &cache = CachingType() )
        : Base( type, coords, cache )
      {}

      template< int fatherdim >
      LocalGeometry ( const Geometry< fatherdim, cdim, Grid > &father, int i,
                      const CachingType &cache = CachingType() )
        : Base( father, i, cache )
      {}
    };

  }

}

#endif
