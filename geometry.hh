// -*- tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PYTHON_GRID_GEOMETRY_HH
#define DUNE_PYTHON_GRID_GEOMETRY_HH

#include <string>
#include <type_traits>
#include <utility>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/visibility.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/python/common/vector.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/numpy.h>

namespace Dune
{

  namespace CorePy
  {

    // PyCornerRange
    // -------------

    template<class Geometry>
    struct PyCorners
    {
      PyCorners(const Geometry &geometry, pybind11::object ref)
        : geometry_(geometry), ref_(ref)
      {}

      const Geometry& geometry() { return geometry_; }

    private:
      const Geometry &geometry_;
      pybind11::object ref_;
    };


    // PyCornerIterator
    // ----------------

    template<class Geometry>
    struct PyCornerIterator
    {
      PyCornerIterator(const PyCorners<Geometry> corners) : corners_(corners) {}

      typename Geometry::GlobalCoordinate next()
      {
        if(index_ == corners_.geometry().corners())
          throw pybind11::stop_iteration();

        return corners_.geometry().corner(index_++);
      }

    private:
      PyCorners<Geometry> corners_;
      int index_ = 0;
    };



    namespace detail
    {

      // registerGridGeometry
      // --------------------

      template< class Geometry, class... options >
      void registerGridGeometry ( pybind11::handle scope,
           pybind11::class_<Geometry, options...> cls )
      {
        const int mydimension = Geometry::mydimension;
        const int coorddimension = Geometry::coorddimension;

        typedef typename Geometry::ctype ctype;
        typedef FieldVector< ctype, mydimension > LocalCoordinate;
        typedef FieldVector< ctype, coorddimension > GlobalCoordinate;
        typedef FieldMatrix< ctype, mydimension, coorddimension > JacobianTransposed;
        typedef FieldMatrix< ctype, coorddimension, mydimension > JacobianInverseTransposed;

        typedef pybind11::array_t< ctype > Array;

        auto itCls = insertClass< PyCornerIterator< Geometry > > ( cls, "CornerIterator",
            GenerateTypeName("PyCornerIterator", cls) ).first;
        itCls.def( "__iter__", [] ( PyCornerIterator< Geometry > &it ) -> PyCornerIterator< Geometry > & { return it; } );
        itCls.def( "__next__", &PyCornerIterator< Geometry >::next );

        auto cCls = insertClass< PyCorners< Geometry > > ( cls, "Corners",
            GenerateTypeName("PyCorners",cls) ).first;
        cCls.def( "__iter__", [] ( const PyCorners< Geometry > &c ) { return PyCornerIterator< Geometry >( c ); } );

        cls.def_property_readonly( "corners", [] ( pybind11::object geo ) {
            return PyCorners< Geometry >( geo.cast< const Geometry & >(), geo );
          } );

        cls.def_property_readonly( "center", &Geometry::center );
        cls.def_property_readonly( "volume", &Geometry::volume );

        cls.def_property_readonly( "affine", &Geometry::affine );

        cls.def( "position", [] ( const Geometry &self, const LocalCoordinate &x ) { return self.global( x ); } );
        cls.def( "position", [] ( const Geometry &self, Array x ) {
            return vectorize( [ &self ] ( const LocalCoordinate &x ) { return self.global( x ); }, x );
          } );

        cls.def( "localPosition", [] ( const Geometry &self, const GlobalCoordinate &y ) { return self.local( y ); } );
        cls.def( "localPosition", [] ( const Geometry &self, Array y ) {
            return vectorize( [ &self ] ( const GlobalCoordinate &y ) { return self.local( y ); }, y );
          } );

        cls.def( "integrationElement", [] ( const Geometry &self, const LocalCoordinate &x ) { return self.integrationElement( x ); } );
        cls.def( "integrationElement", [] ( const Geometry &self, Array x ) {
            return vectorize( [ &self ] ( const LocalCoordinate &x ) { return self.integrationElement( x ); }, x );
          } );

        cls.def( "jacobianTransposed", [] ( const Geometry &self, const LocalCoordinate &x ) {
            return static_cast< JacobianTransposed >( self.jacobianTransposed( x ) );
          } );
        cls.def( "jacobianTransposed", [] ( const Geometry &self, Array x ) {
            return vectorize( [ &self ] ( const LocalCoordinate &x ) { return static_cast< JacobianTransposed >( self.jacobianTransposed( x ) ); }, x );
          } );

        cls.def( "jacobianInverseTransposed", [] ( const Geometry &self, const LocalCoordinate &x ) {
            return static_cast< JacobianInverseTransposed >( self.jacobianInverseTransposed( x ) );
          } );
        cls.def( "jacobianInverseTransposed", [] ( const Geometry &self, Array x ) {
            return vectorize( [ &self ] ( const LocalCoordinate &x ) { return static_cast< JacobianInverseTransposed >( self.jacobianInverseTransposed( x ) ); }, x );
          } );
      }

    } // namespace detail

    template< class Base, class Geometry = typename Base::Geometry >
    pybind11::class_<Geometry> registerGridGeometry ( pybind11::handle scope )
    {
      auto entry = insertClass<Geometry>(scope, "Geometry",
         GenerateTypeName(scope, "Geometry"));
      if ( entry.second )
        detail::registerGridGeometry(scope,entry.first);
      return entry.first;
    }

  } // namespace CorePy

} // namespace Dune

#endif // #ifndef DUNE_PYTHON_GRID_GEOMETRY_HH
