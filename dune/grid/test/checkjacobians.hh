// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_TEST_CHECKJACOBIANS_HH
#define DUNE_GRID_TEST_CHECKJACOBIANS_HH

#include <dune/common/dynvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/typetraits.hh>

/*  Interface check for Jacobian related return types
 *  -------------------------------------------------
 *
 *  The Dune grid interface allows for implementation-defined return types
 *  of the Dune::Geometry member functions:
 *    JacobianInverseTransposed &jacobianInverseTransposed ( const LocalCoordinate &local );
 *    JacobianTransposed &jacobianTransposed ( const LocalCoordinate &local );
 *
 *  The return types 'JacobianInverseTransposed', 'JacobianTransposed'
 *  are expected to implement the following interface:
\code
  struct Jacobian
  {
    // field type
    typedef ImplementationDefined field_type;
    // size type
    typedef ImplementationDefined size_type;

    // number of rows
    static const int rows = implementationDefined;
    // number of cols
    static const int cols = implementationDefined;

    // linear operations
    template< class X, class Y >
    void mv ( const X &x, Y &y ) const;
    template< class X, class Y >
    void mtv ( const X &x, Y &y ) const;
    template< class X, class Y >
    void umv ( const X &x, Y &y ) const;
    template< class X, class Y >
    void umtv ( const X &x, Y &y ) const;
    template< class X, class Y >
    void umhv ( const X &x, Y &y ) const;
    template< class X, class Y >
    void mmv ( const X &x, Y &y ) const;
    template< class X, class Y >
    void mmtv ( const X &x, Y &y ) const;
    template< class X, class Y >
    void mmhv ( const X &x, Y &y ) const;
    template< class X, class Y >
    void usmv ( const field_type &alpha, const X &x, Y &y ) const;
    template< class X, class Y >
    void usmtv ( const field_type &alpha, const X &x, Y &y ) const;
    template< class X, class Y >
    void usmhv ( const field_type &alpha, const X &x, Y &y ) const;

    // norms
    typename FieldTraits< field_type >::real_type frobenius_norm () const;
    typename FieldTraits< field_type >::real_type frobenius_norm2 () const;
    typename FieldTraits< field_type >::real_type infinity_norm () const;
    typename FieldTraits< field_type >::real_type infinity_norm_real () const;

    // cast to FieldMatrix
    operator FieldMatrix< field_type, rows, cols > () const;
  };
\endcode
 *
 *  Note that a FieldMatrix itself is a valid return type.
 */

namespace Dune
{

  namespace
  {

    // CheckJacobianInterface
    // ----------------------

    template< class ctype, int dimRange, int dimDomain, class Jacobian,
              bool isFieldMatrix = std::is_same< Jacobian, FieldMatrix< ctype, dimRange, dimDomain > >::value
            >
    struct CheckJacobianInterface
    {
      // the interface check always holds if jacobian type is a FieldMatrix
      static void apply ( const Jacobian &jacobian ) {}
    };

    template< class ctype, int dimRange, int dimDomain, class Jacobian >
    struct CheckJacobianInterface< ctype, dimRange, dimDomain, Jacobian, false >
    {
      // field type
      typedef typename Jacobian::field_type field_type;
      typedef typename FieldTraits< field_type >::real_type real_type;
      static_assert(( std::is_convertible<ctype, field_type>::value &&
                      std::is_convertible<field_type, ctype>::value ),
                    "Field type not compatible with geometry's coordinate type");
      // size type
      typedef typename Jacobian::size_type size_type;

      // number of rows
      static const int rows = Jacobian::rows;
      static_assert((rows == dimRange), "Number of rows and range dimension do not coincide");
      // number of cols
      static const int cols = Jacobian::cols;
      static_assert((cols == dimDomain), "Number of columns and domain dimension do not coincide");

      static void apply ( const Jacobian &jacobian )
      {
        // call matrix-vector operations; methods must accept any dense
        // vector implementation
        {
          FieldVector< field_type, cols > x( field_type ( 0  ) );
          FieldVector< field_type, rows > y( field_type( 0 ) );
          callMappings( jacobian, x, y );
        }
        {
          DynamicVector< field_type > x( cols, field_type( 0 ) );
          DynamicVector< field_type > y( rows, field_type( 0 ) );
          callMappings( jacobian, x, y );
        }

        // call norms
        [[maybe_unused]] real_type frobenius_norm
          = jacobian.frobenius_norm();
        [[maybe_unused]] real_type frobenius_norm2
          = jacobian.frobenius_norm2();
        [[maybe_unused]] real_type infinity_norm
          = jacobian.infinity_norm();
        [[maybe_unused]] real_type infinity_norm_real
          = jacobian.infinity_norm_real();

        // cast to FieldMatrix
        FieldMatrix< field_type, rows, cols > A( jacobian );
        FieldMatrix< field_type, rows, cols > B;
        B = jacobian;

        // check consistency with matrix-vector multiplication
        assemble( jacobian, B );
        A -= B;
        if( A.frobenius_norm() > std::numeric_limits< real_type >::epsilon() )
          DUNE_THROW( MathError, "Cast to field matrix is inconsistent with matrix-vector multiplication" );
      }

    private:
      // call matrix-vector multiplication methods
      template< class Domain, class Range >
      static void callMappings ( const Jacobian &jacobian, Domain &x, Range &y )
      {
        field_type alpha( 1 );
        jacobian.mv( x, y );
        jacobian.mtv( y, x );
        jacobian.umv( x, y );
        jacobian.umtv( y, x );
        jacobian.umhv( y, x );
        jacobian.mmv( x, y );
        jacobian.mmtv( y, x );
        jacobian.mmhv( y, x );
        jacobian.usmv( alpha, x, y );
        jacobian.usmtv( alpha, y, x );
        jacobian.usmhv( alpha, y, x );
      }

      // use matrix-vector multiplication to assemble field matrix
      static void assemble ( const Jacobian &jacobian, FieldMatrix< field_type, rows, cols > &fieldMatrix )
      {
        for( int j = 0; j < cols; ++j )
        {
          // i-th standard basis vector
          FieldVector< field_type, cols > e( 0 );
          e[ j ] = 1;

          // compute i-th row
          FieldVector< field_type, rows > x;
          jacobian.mv( e, x );
          for( int i = 0; i < rows; ++i )
            fieldMatrix[ i ][ j ] = x[ i ];
        }
      }
    };

  } // namespace



  // checkJacobians
  // --------------

  template< class Geometry >
  void checkJacobians ( const Geometry &geometry )
  {
    typedef typename Geometry::ctype ctype;
    const int mydimension = Geometry::mydimension;
    const int coorddimension = Geometry::coorddimension;

    const typename Geometry::LocalCoordinate local = geometry.local( geometry.center() );

    // check jacobian inverse transposed
    typedef typename Geometry::JacobianInverseTransposed JacobianInverseTransposed;
    CheckJacobianInterface< ctype, coorddimension, mydimension, JacobianInverseTransposed >
      ::apply( geometry.jacobianInverseTransposed( local ) );

    // check jacobian transposed
    typedef typename Geometry::JacobianTransposed JacobianTransposed;
    CheckJacobianInterface< ctype, mydimension, coorddimension, JacobianTransposed >
      ::apply( geometry.jacobianTransposed( local ) );
  }

} // namespace Dune

#endif // #ifndef DUNE_GRID_TEST_CHECKJACOBIANS_HH
