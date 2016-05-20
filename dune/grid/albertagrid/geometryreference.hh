// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_ALBERTAGRID_GEOMETRYREFERENCE_HH
#define DUNE_GRID_ALBERTAGRID_GEOMETRYREFERENCE_HH

/** \file
    \brief Wrapper and interface classes for element geometries
 */

#include <dune/common/typetraits.hh>

#include <dune/geometry/type.hh>

#include <dune/grid/common/geometry.hh>

namespace Dune
{

  // GeometryReference
  // -----------------

  template< class Implementation >
  class GeometryReference
  {
    typedef GeometryReference< Implementation > This;

  public:
    static const int mydimension = Implementation::mydimension;
    static const int coorddimension = Implementation::coorddimension;

    typedef typename Implementation::ctype ctype;

    typedef typename Implementation::LocalCoordinate LocalCoordinate;
    typedef typename Implementation::GlobalCoordinate GlobalCoordinate;

    typedef typename Implementation::JacobianInverseTransposed JacobianInverseTransposed;
    typedef typename Implementation::JacobianTransposed JacobianTransposed;

    explicit GeometryReference ( const Implementation &impl )
      : impl_( &impl )
    {}

    GeometryType type () const { return impl().type(); }

    bool affine() const { return impl().affine(); }

    int corners () const { return impl().corners(); }
    GlobalCoordinate corner ( int i ) const { return impl().corner( i ); }
    GlobalCoordinate center () const { return impl().center(); }

    GlobalCoordinate global ( const LocalCoordinate &local ) const
    {
      return impl().global( local );
    }

    LocalCoordinate local ( const GlobalCoordinate &global ) const
    {
      return impl().local( global );
    }

    ctype integrationElement ( const LocalCoordinate &local ) const
    {
      return impl().integrationElement( local );
    }

    ctype volume () const { return impl().volume(); }

    JacobianTransposed jacobianTransposed ( const LocalCoordinate &local ) const
    {
      return impl().jacobianTransposed( local );
    }

    JacobianInverseTransposed jacobianInverseTransposed ( const LocalCoordinate &local ) const
    {
      return impl().jacobianInverseTransposed( local );
    }

    const Implementation &impl () const { return *impl_; }

  private:
    const Implementation *impl_;
  };


  // LocalGeometryReference
  // -----------------------

  template< int mydim, int cdim, class Grid >
  class LocalGeometryReference
    : public GeometryReference< typename std::remove_const< Grid >::type::Traits::template Codim< std::remove_const< Grid >::type::dimension - mydim >::LocalGeometryImpl >
  {
    typedef typename std::remove_const< Grid >::type::Traits::template Codim< std::remove_const< Grid >::type::dimension - mydim >::LocalGeometryImpl Implementation;

  public:
    LocalGeometryReference ( const Implementation &impl )
      : GeometryReference< Implementation >( impl )
    {}
  };



  // Definitions of GeometryReference
  // --------------------------------

  template< class Implementation >
  const int GeometryReference< Implementation >::mydimension;

  template< class Implementation >
  const int GeometryReference< Implementation >::coorddimension;

} // namespace Dune

#endif // #ifndef DUNE_GRID_ALBERTAGRID_GEOMETRYREFERENCE_HH
