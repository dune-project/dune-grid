// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PYTHON_GRID_GEOMETRY_HH
#define DUNE_PYTHON_GRID_GEOMETRY_HH

#include <cstddef>

#include <array>
#include <string>
#include <type_traits>
#include <utility>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/visibility.hh>

#include <dune/python/common/fvecmatregistry.hh>

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

      template< class Geometry, class Array >
      inline static pybind11::array_t< double >
      pushForwardGradients ( const Geometry &geo, Array xVec, pybind11::array_t< double > gVec )
      {
        typedef typename Geometry::LocalCoordinate LocalCoordinate;
        typedef typename Geometry::GlobalCoordinate GlobalCoordinate;

        // x = (localCoord,nofQuad)
        auto x = xVec.unchecked();
        // g = (dimRange,localCoord,nofQuad)
        auto g = gVec.unchecked();
        // ret = (dimRange,globalCoord,nofQuad)
        pybind11::array_t< double > ret( std::array< ssize_t, 3 >{{ g.shape( 0 ), static_cast< ssize_t >( Geometry::GlobalCoordinate::size() ), g.shape( 2 ) }} );
        auto y = ret.template mutable_unchecked< 3 >();
        if( x.shape( 1 ) != g.shape( 2 ) )
          std::cout << x.shape( 1 ) << " " << g.shape( 2 ) << std::endl;
        if( x.shape( 0 ) != ssize_t(LocalCoordinate::size()) )
          std::cout << x.shape( 0 ) << " " << Geometry::LocalCoordinate::size() << std::endl;

        for( ssize_t p = 0; p < g.shape( 2 ); ++p )
        {
          LocalCoordinate xLocal;
          for( std::size_t l = 0; l < LocalCoordinate::size(); ++l )
            xLocal[ l ] = x( l, p );
          const auto jit = geo.jacobianInverseTransposed( xLocal );

          for( ssize_t range = 0; range < g.shape( 0 ); ++range )
          {
            // Performance Issue:
            // The copies gradLocal and gradGlobal can be avoided by providing
            // a DenseVector implementation based on a single axis of the
            // pybind11::array accessor, because the `jit.mv` method is required
            // to take arbitrary implementations of the dense vector interface.
            LocalCoordinate gradLocal;
            for( std::size_t l = 0; l < LocalCoordinate::size(); ++l )
              gradLocal[ l ] = g(range, l, p );

            GlobalCoordinate gradGlobal;
            jit.mv( gradLocal, gradGlobal );

            for( std::size_t r = 0; r < GlobalCoordinate::size(); ++r )
              y( range, r, p ) = gradGlobal[ r ];
          }
        }
        return ret;
      }

      template< class Geometry, class... options >
      void registerGridGeometry ( pybind11::handle scope, pybind11::class_<Geometry, options...> cls )
      {
        const int mydimension = Geometry::mydimension;
        const int coorddimension = Geometry::coorddimension;

        typedef typename Geometry::ctype ctype;
        typedef FieldVector< ctype, mydimension > LocalCoordinate;
        typedef FieldVector< ctype, coorddimension > GlobalCoordinate;
        typedef FieldMatrix< ctype, mydimension, coorddimension > JacobianTransposed;
        typedef FieldMatrix< ctype, coorddimension, mydimension > JacobianInverseTransposed;
        registerFieldVecMat<LocalCoordinate>::apply();
        registerFieldVecMat<GlobalCoordinate>::apply();
        registerFieldVecMat<JacobianTransposed>::apply();
        registerFieldVecMat<JacobianInverseTransposed>::apply();

        typedef pybind11::array_t< ctype > Array;

        using pybind11::operator""_a;

        pybind11::options opts;
        opts.disable_function_signatures();

        cls.doc() = R"doc(
            A geometry describes a map from the reference domain into a
            Euclidian space, where the reference domain is given by the
            reference element.
            The mapping is required to be one-to-one.

            We refer to points within the reference domain as "local" points.
            The image of a local point is called its (global) position.

            Note: The image of the mapping may be a submanifold of the
                  Euclidian space.
          )doc";

        cls.def( "corner", [] ( const Geometry &self, int i ) {
            const int size = self.corners();
            if( (i < 0) || (i >= size) )
              throw pybind11::value_error( "Invalid index: " + std::to_string( i ) + " (must be in [0, " + std::to_string( size ) + "))." );
            return self.corner( i );
          }, "index"_a,
          R"doc(
            get global position a reference corner

            Args:
                index:    index of the reference corner

            Returns:
                global position of the reference corner

            Note: If the argument "index" is omitted, this method returns a
                  NumPy array containing the global position of all corners.
                  This version may be used to vectorize the code.
          )doc" );
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
          },
          R"doc(
            get global positions of all reference corners

            Note:
                This function differs from the vectorized version of 'corner'
                in the way the corners are returned. This method returns a
                tuple of global positions of type FieldVector.

            Returns:
                tuple of global positions, in the order given by the reference element
          )doc" );


        cls.def_property_readonly( "center", [] ( const Geometry &self ) { return self.center(); },
          R"doc(
            global position of the barycenter of the reference element
          )doc" );
        cls.def_property_readonly( "volume", [] ( const Geometry &self ) { return self.volume(); },
          R"doc(
            volume of the map's image

            The volume is measured using the Hausdorff measure of the corresponding dimension.
          )doc" );

        cls.def_property_readonly( "affine", [] ( const Geometry &self ) { return self.affine(); },
          R"doc(
            True, if the map is affine-linear, False otherwise
          )doc" );
        cls.def_property_readonly( "referenceElement", []( const Geometry &self ) {
            return referenceElement< double, mydimension >( self.type() );
          }, pybind11::keep_alive< 0, 1 >(),
          R"doc(
            corresponding reference element, describing the domain of the map
          )doc" );

        cls.def( "toGlobal", [] ( const Geometry &self, const LocalCoordinate &x ) { return self.global( x ); }, "x"_a,
          R"doc(
            obtain global position of a local point

            Args:
                x:    local point

            Returns:
                global position of x

            Note: This method may be used in vectorized form by passing in a
                  NumPy array for 'x'.
          )doc" );
        cls.def( "toGlobal", [] ( const Geometry &self, Array x ) {
            return vectorize( [ &self ] ( const LocalCoordinate &x ) { return self.global( x ); }, x );
          } );

        cls.def( "toLocal", [] ( const Geometry &self, const GlobalCoordinate &y ) { return self.local( y ); }, "y"_a,
          R"doc(
            obtain local point mapped to a global position

            Args:
                y:    global position

            Returns:
                local point mapped to y

            Note: This method may be used in vectorized form by passing in a
                  NumPy array for 'y'.
          )doc" );
        cls.def( "toLocal", [] ( const Geometry &self, Array y ) {
            return vectorize( [ &self ] ( const GlobalCoordinate &y ) { return self.local( y ); }, y );
          } );

        cls.def( "integrationElement", [] ( const Geometry &self, const LocalCoordinate &x ) { return self.integrationElement( x ); }, "x"_a,
          R"doc(
            obtain integration element in a local point

            The integration element is the factor appearing in the integral
            transformation formula.
            It describes the weight factor when transforming a quadrature rule
            on the reference element into a quadrature rule on the image of this
            map.

            Args:
                x:    local point

            Returns:
                integration element in x

            Note: This method may be used in vectorized form by passing in a
                  NumPy array for 'x'.
          )doc" );
        cls.def( "integrationElement", [] ( const Geometry &self, Array x ) {
            return vectorize( [ &self ] ( const LocalCoordinate &x ) { return self.integrationElement( x ); }, x );
          } );

        cls.def( "jacobianTransposed", [] ( const Geometry &self, const LocalCoordinate &x ) {
            return static_cast< JacobianTransposed >( self.jacobianTransposed( x ) );
          }, "x"_a,
          R"doc(
            obtain transposed of the Jacobian of this mapping in a local point

            The rows of the returned matrix describe the tangential vectors in
            the global position of the local point.

            The Jacobian itself describes the push-forward for tangential
            vectors from the reference domain to the image of this map.

            Args:
                x:    local point

            Returns:
                transposed of the Jacobian matrix in x

            Note: This method may be used in vectorized form by passing in a
                  NumPy array for 'x'.
          )doc" );
        cls.def( "jacobianTransposed", [] ( const Geometry &self, Array x ) {
            return vectorize( [ &self ] ( const LocalCoordinate &x ) { return static_cast< JacobianTransposed >( self.jacobianTransposed( x ) ); }, x );
          } );

        cls.def( "jacobianInverseTransposed", [] ( const Geometry &self, const LocalCoordinate &x ) {
            return static_cast< JacobianInverseTransposed >( self.jacobianInverseTransposed( x ) );
          }, "x"_a,
          R"doc(
            obtain transposed of the inverse Jacobian of this mapping in a local point

            This matrix describes the push-forward for local function gradients
            to the image of this map.

            The inverse Jacobian itself describes the push-forward of cotangential
            vectors from the reference domain to the image of this map.

            Args:
                x:    local point

            Returns:
                transposed of the Jacobian matrix in x

            Note: This method may be used in vectorized form by passing in a
                  NumPy array for 'x'.
          )doc" );
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
