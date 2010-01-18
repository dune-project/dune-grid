// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LABATTOQUADRATURE_HH
#define DUNE_LABATTOQUADRATURE_HH

#include <iostream>

#if HAVE_ALGLIB
#include <alglib/gqgenlobatto.h>
#endif

#include <dune/common/static_assert.hh>

#include <dune/grid/common/quadraturerules/genericquadrature.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // LobattoPoints
    // -----------

    /**
     * @brief Lobatto quadrature points and weights in 1d -  only available if ALGLib is used
     **/
    template< class F>
    struct LobattoPoints
    {
      dune_static_assert(sizeof(F)==0,"Lobatto Points only implemented for ALGLib ampf type");
      explicit LobattoPoints ( unsigned int n )
      {}
    };

#if HAVE_ALGLIB
    template< unsigned int precision >
    class LobattoPoints< amp::ampf< precision > >
      : public std::vector< QuadraturePoint<amp::ampf< precision>,1> >
    {
      typedef amp::ampf< precision > F;
      typedef std::vector< QuadraturePoint<F,1> > Base;
      typedef LobattoPoints< F > This;

    public:
      typedef F Field;

      explicit LobattoPoints ( unsigned int n )
      {
        Base::reserve( n );

        typedef ap::template_1d_array< Field > AlgLibVector;
        AlgLibVector p,w;
        p.setbounds( 0, n-1 );
        w.setbounds( 0, n-1 );
        AlgLibVector alpha,beta;
        alpha.setbounds( 0, n-1 );
        beta.setbounds( 0, n-1 );

        for( unsigned int i = 0; i < n; ++i )
        {
          alpha( i ) = 0;
          beta( i ) = Field( i*i ) / Field( 4*(i*i)-1 );
        }

        bool succ = gqgenlobatto::generategausslobattoquadrature< precision >( alpha, beta, Field( 2 ), Field(-1),Field(1), n, p, w );
        if (!succ)
        {
          std::cout << "Problem with computing Lobatto points!" << std::endl;
        }

        const Field half = Field( 1 ) / Field( 2 );
        for( unsigned int i = 0; i < n; ++i )
        {
          QuadraturePoint<Field,1> q( (p( i ) + Field( 1 )) * half, w( i )*half );
          Base::push_back( q );
        }
      }
      // order is the maximal order of polynomials which can be exactly integrated
      static unsigned int minPoints( unsigned int order )
      {
        return (order+4)/2;
      }
    };
#endif // #if HAVE_ALGLIB

    /**
     * @brief Singleton provider for Lobatto quadratures
     *
     * \tparam dim dimension of the reference elements contained in the factory
     * \tparam F field in which weight and point of the quadrature are stored
     * \tparam CF the compute field for the points and weights
     **/
    template< int dim, class F, class CF=F >
    struct LobattoQuadratureProvider
      : public TopologySingletonFactory< GenericQuadratureFactory< dim, F, LobattoPoints<CF> > >
    {};
  }

}

#endif // #ifndef DUNE_LOBATTOQUADRATURE_HH
