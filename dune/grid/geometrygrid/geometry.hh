// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_GEOMETRY_HH
#define DUNE_GEOGRID_GEOMETRY_HH

#include <dune/common/nullptr.hh>
#include <dune/common/typetraits.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/mapping/cornermapping.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/geometrygrid/cornerstorage.hh>

namespace Dune
{

  namespace GeoGrid
  {

    // MappingTraits
    // -------------

    template< class Grid >
    struct MappingTraits
    {
      typedef typename remove_const< Grid >::type::Traits::ctype ctype;

      typedef GenericGeometry::MatrixHelper< GenericGeometry::DuneCoordTraits< ctype > > MatrixHelper;

      static ctype tolerance () { return 16 * std::numeric_limits< ctype >::epsilon(); }

      template< int mydim, int cdim >
      struct CornerStorage
      {
        typedef GeoGrid::CornerStorage< mydim, cdim, Grid > Type;
      };

      struct UserData
      {
        UserData () : refCount_( 0 ) {}

        void addReference () { ++refCount_; }
        bool removeReference () { return (--refCount_ == 0); }

      private:
        unsigned int refCount_;
      };
    };



    // Geometry
    // --------

    template< int mydim, int cdim, class Grid >
    class Geometry
    {
      typedef Geometry< mydim, cdim, Grid > This;

      typedef typename remove_const< Grid >::type::Traits Traits;

      template< int, int, class > friend class Geometry;

    public:
      typedef typename Traits::ctype ctype;

      static const int mydimension = mydim;
      static const int coorddimension = cdim;
      static const int dimension = Traits::dimension;
      static const int codimension = dimension - mydimension;

    protected:
      typedef CachedCornerMapping< ctype, mydimension, coorddimension, MappingTraits< Grid > > Mapping;

    public:
      typedef typename Mapping::LocalCoordinate LocalCoordinate;
      typedef typename Mapping::GlobalCoordinate GlobalCoordinate;

      typedef typename Mapping::JacobianTransposed JacobianTransposed;
      typedef typename Mapping::JacobianInverseTransposed JacobianInverseTransposed;

      Geometry ( const Grid &grid )
        : grid_( &grid ),
          mapping_( nullptr )
      {}

      template< class CoordVector >
      Geometry ( const Grid &grid, const GeometryType &type, const CoordVector &coords )
        : grid_( &grid )
      {
        assert( int( type.dim() ) == mydimension );
        void *mappingStorage = grid.allocateStorage( sizeof( Mapping ) );
        mapping_ = new( mappingStorage ) Mapping( type, coords );
        mapping_->userData().addReference();
      }

      Geometry ( const This &other )
        : grid_( other.grid_ ),
          mapping_( other.mapping_ )
      {
        if( mapping_ )
          mapping_->userData().addReference();
      }

      ~Geometry ()
      {
        if( mapping_ && mapping_->userData().removeReference() )
          destroyMapping();
      }

      const This &operator= ( const This &other )
      {
        if( other.mapping_ )
          other.mapping_->userData().addReference();
        if( mapping_ && mapping_->userData().removeReference() )
          destroyMapping();
        grid_ = other.grid_;
        mapping_ = other.mapping_;
        return *this;
      }

      operator bool () const { return bool( mapping_ ); }

      bool affine () const { return mapping_->affine(); }
      GeometryType type () const { return mapping_->type(); }

      int corners () const { return mapping_->corners(); }
      GlobalCoordinate corner ( const int i ) const { return mapping_->corner( i ); }
      GlobalCoordinate center () const { return mapping_->center(); }

      GlobalCoordinate global ( const LocalCoordinate &local ) const { return mapping_->global( local ); }
      LocalCoordinate local ( const GlobalCoordinate &global ) const { return mapping_->local( global ); }

      ctype integrationElement ( const LocalCoordinate &local ) const { return mapping_->integrationElement( local ); }
      ctype volume () const { return mapping_->volume(); }

      const JacobianTransposed &jacobianTransposed ( const LocalCoordinate &local ) const { return mapping_->jacobianTransposed( local ); }
      const JacobianInverseTransposed &jacobianInverseTransposed ( const LocalCoordinate &local ) const { return mapping_->jacobianInverseTransposed( local ); }

      const Grid &grid () const { return *grid_; }

    private:
      void destroyMapping ()
      {
        mapping_->~Mapping();
        grid().deallocateStorage( mapping_, sizeof( Mapping ) );
      }

      const Grid *grid_;
      Mapping* mapping_;
    };

  } // namespace GeoGrid



  // FacadeOptions
  // -------------

  namespace FacadeOptions
  {

    template< int mydim, int cdim, class Grid >
    struct StoreGeometryReference< mydim, cdim, Grid, GeoGrid::Geometry >
    {
      static const bool v = false;
    };

  } // namespace FacadeOptions

} // namespace Dune

#endif // #ifndef DUNE_GEOGRID_GEOMETRY_HH
