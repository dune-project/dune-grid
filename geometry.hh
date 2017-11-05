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

  namespace Python
  {

    namespace detail
    {

      // registerGridGeometry
      // --------------------

      template <class Geometry, class Array>
      pybind11::array_t<double> pushForwardGradients( const Geometry &geo, Array xVec,
                                pybind11::array_t<double> gVec )
      {
        // x = (localCoord,nofQuad)
        auto x = xVec.unchecked();
        // g = (dimRange,localCoord,nofQuad)
        auto g = gVec.unchecked();
        // ret = (dimRange,globalCoord,nofQuad)
        pybind11::array_t<double> ret( std::vector<size_t>{g.shape(0),(size_t)Geometry::GlobalCoordinate::size(),g.shape(2)} );
        auto y = ret.template mutable_unchecked< 3 >();
        if (x.shape(1) != g.shape(2))
          std::cout << x.shape(1) << " " << g.shape(2) << std::endl;
        if (x.shape(0) != Geometry::LocalCoordinate::size())
          std::cout << x.shape(0) << " " << Geometry::LocalCoordinate::size() << std::endl;

        for (size_t p=0;p<g.shape(2);++p)
        {
          typename Geometry::LocalCoordinate loc;
          for (size_t l=0;l<loc.size();++l)
            loc[l] = x(l,p);
          auto jit = geo.jacobianInverseTransposed( loc );
          for (size_t range=0;range<g.shape(0);++range)
            for (size_t r=0;r<Geometry::GlobalCoordinate::size();++r)
            {
              y(range,r,p) = 0;
              for (size_t c=0;c<Geometry::LocalCoordinate::size();++c)
                y(range,r,p) += jit[r][c]*g(range,c,p);
            }
        }
        return ret;
      }

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

        using pybind11::operator""_a;

        cls.def( "corner", [] ( const Geometry &self, int i ) {
            const int size = self.corners();
            if( (i < 0) || (i >= size) )
              throw pybind11::value_error( "Invalid index: " + std::to_string( i ) + " (must be in [0, " + std::to_string( size ) + "))." );
            return self.corner( i );
          }, "index"_a );
        cls.def( "corner", [] ( const Geometry &self ) {
            const int size = self.corners();
            pybind11::array_t< ctype > cornersArray( { static_cast< ssize_t >( coorddimension ), static_cast< ssize_t >( size ) } );
            auto corners = cornersArray.template mutable_unchecked< 2 >();
            for( int i = 0; i < size; ++i )
            {
              const auto corner = self.corner( i );
              for( int j = 0; j < coorddimension; ++j )
                corners( j, i ) = corner[ j ];
            }
            return cornersArray;
          } );
        cls.def_property_readonly( "corners", [] ( const Geometry &self ) {
            const int size = self.corners();
            pybind11::tuple corners( size );
            for( int i = 0; i < size; ++i )
              corners[ i ] = pybind11::cast( self.corner( i ) );
            return corners;
          } );


        cls.def_property_readonly( "center", &Geometry::center );
        cls.def_property_readonly( "volume", &Geometry::volume );

        cls.def_property_readonly( "affine", &Geometry::affine );
        cls.def_property_readonly( "domain", []( const Geometry &self) { return referenceElement<double,mydimension>(self.type()); },
            pybind11::keep_alive<0,1>() );

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

        cls.def( "pushForwardGradients", [] ( const Geometry &self, Array x, pybind11::array_t<double> g ) {
            return pushForwardGradients(self,x,g);
          } );
      }

    } // namespace detail



    // registerGridGeometry
    // --------------------

    template< class Base, class Geometry = typename Base::Geometry >
    inline static pybind11::class_< Geometry > registerGridGeometry ( pybind11::handle scope )
    {
      auto entry = insertClass< Geometry >( scope, "Geometry", GenerateTypeName( scope, "Geometry" ) );
      if ( entry.second )
        detail::registerGridGeometry( scope, entry.first );
      return entry.first;
    }

  } // namespace Python

} // namespace Dune

#endif // #ifndef DUNE_PYTHON_GRID_GEOMETRY_HH
