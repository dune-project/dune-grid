// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_GEOMETRY_HH
#define DUNE_GEOGRID_GEOMETRY_HH

#include <dune/common/nullptr.hh>
#include <dune/common/typetraits.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/genericgeometry/mappingprovider.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/geometrygrid/cornerstorage.hh>

namespace Dune
{

  namespace GeoGrid
  {

    // GeometryTraits
    // --------------

    template< int cdim, class Grid >
    struct GeometryTraits
    {
      typedef typename remove_const< Grid >::type::Traits Traits;

      typedef GenericGeometry::DuneCoordTraits< typename Traits::ctype > CoordTraits;

      static const int dimWorld = cdim;

      template< class Topology >
      struct Mapping
      {
        typedef GeoGrid::CornerStorage< Topology, const Grid > CornerStorage;
        typedef GenericGeometry::CornerMapping< CoordTraits, Topology, dimWorld, CornerStorage > type;
      };

      struct Caching
      {
        static const GenericGeometry::EvaluationType evaluateJacobianTransposed = GenericGeometry::ComputeOnDemand;
        static const GenericGeometry::EvaluationType evaluateJacobianInverseTransposed = GenericGeometry::ComputeOnDemand;
        static const GenericGeometry::EvaluationType evaluateIntegrationElement = GenericGeometry::ComputeOnDemand;
        static const GenericGeometry::EvaluationType evaluateNormal = GenericGeometry::ComputeOnDemand;
      };
    };



    // Geometry
    // --------

    template< int mydim, int cdim, class Grid >
    class Geometry
    {
      typedef Geometry< mydim, cdim, Grid > This;

      typedef typename remove_const< Grid >::type::Traits Traits;
      typedef GeoGrid::GeometryTraits< cdim, Grid > GeometryTraits;

      template< int, int, class > friend class Geometry;

    public:
      static const int mydimension = mydim;
      static const int coorddimension = cdim;
      static const int dimension = Traits::dimension;
      static const int codimension = dimension - mydimension;

      typedef typename GeometryTraits::CoordTraits::ctype ctype;

      typedef FieldVector< ctype, mydimension > LocalCoordinate;
      typedef FieldVector< ctype, coorddimension > GlobalCoordinate;

    private:
      typedef typename Traits::HostGrid HostGrid;

      template< bool >
      struct Hybrid
      {
        typedef GenericGeometry::HybridMapping< dimension, GeometryTraits > Mapping;
      };

      template< bool >
      struct NonHybrid
      {
        static const unsigned int topologyId = Capabilities::hasSingleGeometryType< HostGrid >::topologyId;
        typedef typename GenericGeometry::Topology< topologyId, dimension >::type Topology;
        typedef GenericGeometry::CachedMapping< Topology, GeometryTraits > Mapping;
      };

      typedef typename SelectType< Capabilities::hasSingleGeometryType< HostGrid >::v, NonHybrid< true >, Hybrid< false > >::Type::Mapping ElementMapping;
      typedef GenericGeometry::MappingProvider< ElementMapping, codimension > MappingProvider;

    protected:
      typedef typename MappingProvider::Mapping Mapping;

    public:
      typedef typename Mapping::JacobianTransposed JacobianTransposed;
      typedef typename Mapping::JacobianInverseTransposed JacobianInverseTransposed;

      // for cenvencience: Jacobian is the name of the type in the geometry interface
      typedef JacobianInverseTransposed Jacobian;

      Geometry ()
        : mapping_( nullptr )
      {}

      template< class CoordVector >
      Geometry ( const GeometryType &type, const CoordVector &coords )
      {
        mapping_ = MappingProvider::construct( type.id(), coords, mappingStorage_ );
      }

      template< int fatherdim >
      Geometry ( const Geometry< fatherdim, cdim, Grid > &father, int i )
      {
        const unsigned int codim = fatherdim - mydim;
        mapping_ = father.mapping_->template trace< codim >( i, mappingStorage_ );
      }

      Geometry ( const This &other )
        : mapping_( other.mapping_ ? other.mapping_->clone( mappingStorage_ ) : nullptr )
      {}

      ~Geometry ()
      {
        if( mapping_ )
          mapping_->~Mapping();
      }

      const This &operator= ( const This &other )
      {
        if( mapping_ )
          mapping_->~Mapping();
        mapping_ = (other.mapping_) ? other.mapping_->clone( mappingStorage_ ) : nullptr;
        return *this;
      }

      operator bool () const { return bool( mapping_ ); }

      bool affine () const { return mapping_->affine(); }
      GeometryType type () const { return mapping_->type(); }

      int corners () const { return mapping_->numCorners(); }
      GlobalCoordinate corner ( const int i ) const { return mapping_->corner( i ); }
      GlobalCoordinate center () const { return mapping_->center(); }

      GlobalCoordinate global ( const LocalCoordinate &local ) const { return mapping_->global( local ); }
      LocalCoordinate local ( const GlobalCoordinate &global ) const { return mapping_->local( global ); }

      ctype integrationElement ( const LocalCoordinate &local ) const { return mapping_->integrationElement( local ); }
      ctype volume () const { return mapping_->volume(); }

      const JacobianTransposed &jacobianTransposed ( const LocalCoordinate &local ) const { return mapping_->jacobianTransposed( local ); }
      const JacobianInverseTransposed &jacobianInverseTransposed ( const LocalCoordinate &local ) const { return mapping_->jacobianInverseTransposed( local ); }

    private:
      Mapping* mapping_;
      char mappingStorage_[ MappingProvider::maxMappingSize ];
    };

  } // namespace GeoGrid

} // namespace Dune

#endif // #ifndef DUNE_GEOGRID_GEOMETRY_HH
