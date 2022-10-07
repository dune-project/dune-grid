// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_GEOMETRY_HH
#define DUNE_GEOGRID_GEOMETRY_HH

#include <utility>

#include <dune/common/typetraits.hh>

#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/geometrygrid/cornerstorage.hh>

namespace Dune
{

  namespace GeoGrid
  {

    // InferHasSingleGeometryType
    // --------------------------

    template< class hasSingleGeometryType, int dim, int mydim >
    struct InferHasSingleGeometryType
    {
    private:
      static const unsigned int id = hasSingleGeometryType::topologyId;
      static const unsigned int idMask = (1u << mydim) - 1u;

    public:
      static const bool v = hasSingleGeometryType::v && ((mydim == dim) || ((id | 1u) == 1u) || ((id | 1u) == idMask));
      static const unsigned int topologyId = (v ? id & idMask : ~0u);
    };

    template< class hasSingleGeometryType, int dim >
    struct InferHasSingleGeometryType< hasSingleGeometryType, dim, 1 >
    {
      static const bool v = true;
      static const unsigned int topologyId = GeometryTypes::cube(1).id();
    };

    template< class hasSingleGeometryType, int dim >
    struct InferHasSingleGeometryType< hasSingleGeometryType, dim, 0 >
    {
      static const bool v = true;
      static const unsigned int topologyId = GeometryTypes::cube(1).id();
    };



    // GeometryTraits
    // --------------

    template< class Grid >
    struct GeometryTraits
    {
      typedef typename std::remove_const< Grid >::type::Traits Traits;

      typedef typename Traits::ctype ctype;

      typedef Impl::FieldMatrixHelper< ctype > MatrixHelper;

      static ctype tolerance () { return 16 * std::numeric_limits< ctype >::epsilon(); }

      template< int mydim, int cdim >
      struct CornerStorage
      {
        typedef GeoGrid::CornerStorage< mydim, cdim, Grid > Type;
      };

      template< int mydim >
      struct hasSingleGeometryType
        : public InferHasSingleGeometryType< Capabilities::hasSingleGeometryType< Grid >, Traits::dimension, mydim >
      {};
    };



    // Geometry
    // --------

    template< int mydim, int cdim, class Grid >
    class Geometry
    {
      typedef Geometry< mydim, cdim, Grid > This;

      typedef typename std::remove_const< Grid >::type::Traits Traits;

      template< int, int, class > friend class Geometry;

    public:
      typedef typename Traits::ctype ctype;

      static const int mydimension = mydim;
      static const int coorddimension = cdim;
      static const int dimension = Traits::dimension;
      static const int codimension = dimension - mydimension;

    protected:
      typedef CachedMultiLinearGeometry< ctype, mydimension, coorddimension, GeometryTraits< Grid > > BasicMapping;

      struct Mapping
        : public BasicMapping
      {
        template< class CoordVector >
        Mapping ( const GeometryType &type, const CoordVector &coords )
          : BasicMapping( type, coords ),
            refCount_( 0 )
        {}

        void addReference () { ++refCount_; }
        bool removeReference () { return (--refCount_ == 0); }

      private:
        unsigned int refCount_;
      };

    public:
      typedef typename Mapping::LocalCoordinate LocalCoordinate;
      typedef typename Mapping::GlobalCoordinate GlobalCoordinate;

      typedef typename Mapping::JacobianTransposed JacobianTransposed;
      typedef typename Mapping::JacobianInverseTransposed JacobianInverseTransposed;
      typedef typename Mapping::Jacobian Jacobian;
      typedef typename Mapping::JacobianInverse JacobianInverse;


      Geometry () : grid_( nullptr ), mapping_( nullptr ) {}

      explicit Geometry ( const Grid &grid ) : grid_( &grid ), mapping_( nullptr ) {}

      template< class CoordVector >
      Geometry ( const Grid &grid, const GeometryType &type, const CoordVector &coords )
        : grid_( &grid )
      {
        assert( int( type.dim() ) == mydimension );
        void *mappingStorage = grid.allocateStorage( sizeof( Mapping ) );
        mapping_ = new( mappingStorage ) Mapping( type, coords );
        mapping_->addReference();
      }

      Geometry ( const This &other )
        : grid_( other.grid_ ),
          mapping_( other.mapping_ )
      {
        if( mapping_ )
          mapping_->addReference();
      }

      Geometry ( This&& other )
        : grid_( other.grid_ ),
          mapping_( other.mapping_ )
      {
        other.grid_ = nullptr;
        other.mapping_ = nullptr;
      }

      ~Geometry ()
      {
        if( mapping_ && mapping_->removeReference() )
          destroyMapping();
      }

      const This &operator= ( const This &other )
      {
        if( other.mapping_ )
          other.mapping_->addReference();
        if( mapping_ && mapping_->removeReference() )
          destroyMapping();
        grid_ = other.grid_;
        mapping_ = other.mapping_;
        return *this;
      }

      const This &operator= ( This&& other )
      {
        using std::swap;
        swap( grid_, other.grid_ );
        swap( mapping_, other.mapping_ );
        return *this;
      }

      explicit operator bool () const { return bool( mapping_ ); }

      bool affine () const { return mapping_->affine(); }
      GeometryType type () const { return mapping_->type(); }

      int corners () const { return mapping_->corners(); }
      GlobalCoordinate corner ( const int i ) const { return mapping_->corner( i ); }
      GlobalCoordinate center () const { return mapping_->center(); }

      GlobalCoordinate global ( const LocalCoordinate &local ) const { return mapping_->global( local ); }
      LocalCoordinate local ( const GlobalCoordinate &global ) const { return mapping_->local( global ); }

      ctype integrationElement ( const LocalCoordinate &local ) const { return mapping_->integrationElement( local ); }
      ctype volume () const { return mapping_->volume(); }

      JacobianTransposed jacobianTransposed ( const LocalCoordinate &local ) const { return mapping_->jacobianTransposed( local ); }
      JacobianInverseTransposed jacobianInverseTransposed ( const LocalCoordinate &local ) const { return mapping_->jacobianInverseTransposed( local ); }

      Jacobian jacobian ( const LocalCoordinate &local ) const { return mapping_->jacobian( local ); }
      JacobianInverse jacobianInverse ( const LocalCoordinate &local ) const { return mapping_->jacobianInverse( local ); }

      const Grid &grid () const { assert( grid_ ); return *grid_; }

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

} // namespace Dune

#endif // #ifndef DUNE_GEOGRID_GEOMETRY_HH
