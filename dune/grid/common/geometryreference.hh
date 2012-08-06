// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_GEOMETRYREFERENCE_HH
#define DUNE_GRID_GEOMETRYREFERENCE_HH

/** \file
    \brief Wrapper and interface classes for element geometries
 */

#include <dune/common/typetraits.hh>

#include <dune/geometry/type.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< int mydim, int cdim, class Grid >
  class GlobalGeometryReference;

  template< int mydim, int cdim, class Grid >
  class LocalGeometryReference;



  // FacadeOptions
  // -------------

  namespace FacadeOptions
  {

    template< int mydim, int cdim, class GridImp, template< int, int, class > class GeometryImp >
    struct StoreGeometryReference;

    template< int mydim, int cdim, class Grid >
    struct StoreGeometryReference< mydim, cdim, Grid, GlobalGeometryReference >
    {
      static const bool v = false;
    };

    template< int mydim, int cdim, class Grid >
    struct StoreGeometryReference< mydim, cdim, Grid, LocalGeometryReference >
    {
      static const bool v = false;
    };

  } // namespace FacadeOptions



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

    typedef typename Implementation::Jacobian Jacobian;
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

    const JacobianTransposed &jacobianTransposed ( const LocalCoordinate &local ) const
    {
      return impl().jacobianTransposed( local );
    }

    const Jacobian &jacobianInverseTransposed ( const LocalCoordinate &local ) const
    {
      return impl().jacobianInverseTransposed( local );
    }

    const Implementation &impl () const { return *impl_; }

  private:
    const Implementation *impl_;
  };



  // GlobalGeometryReference
  // -----------------------

  template< int mydim, int cdim, class Grid >
  class GlobalGeometryReference
    : public GeometryReference< typename remove_const< Grid >::type::Traits::template Codim< remove_const< Grid >::type::dimension - mydim >::GeometryImpl >
  {
    typedef typename remove_const< Grid >::type::Traits::template Codim< remove_const< Grid >::type::dimension - mydim >::GeometryImpl Implementation;

  public:
    GlobalGeometryReference ( const Implementation &impl )
      : GeometryReference< Implementation >( impl )
    {}
  };



  // LocalGeometryReference
  // -----------------------

  template< int mydim, int cdim, class Grid >
  class LocalGeometryReference
    : public GeometryReference< typename remove_const< Grid >::type::Traits::template Codim< remove_const< Grid >::type::dimension - mydim >::LocalGeometryImpl >
  {
    typedef typename remove_const< Grid >::type::Traits::template Codim< remove_const< Grid >::type::dimension - mydim >::LocalGeometryImpl Implementation;

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

#endif // #ifndef DUNE_GRID_GEOMETRYREFERENCE_HH
