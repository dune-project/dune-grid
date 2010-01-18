// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GAUSSQUADRATURE_HH
#define DUNE_GAUSSQUADRATURE_HH

#if HAVE_ALGLIB
#include <alglib/gqgengauss.h>
#endif

#include <dune/grid/common/quadraturerules.hh>
#include <dune/grid/common/quadraturerules/genericquadrature.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // GaussPoints
    // -----------

    /**
     * @brief Gauss quadrature points and weights in 1d.
     *
     * if ALGLib is found higher precision Gauss points can be used by
     * prescribing the amp::ampf field type; otherwise
     * the dune-grid quadrature is used.
     **/
    template< class F>
    class GaussPoints
      : public std::vector< QuadraturePoint<F,1> >
    {
      typedef std::vector< QuadraturePoint<F,1> > Base;
      typedef GaussPoints< F > This;

    public:
      typedef F Field;
      // n is number of points required
      explicit GaussPoints ( unsigned int n )
      {
        Base::reserve( n );
        const QuadratureRule<Field,1>& points =
          QuadratureRules<Field,1>::rule(GeometryType(GeometryType::cube,1),
                                         2*n-2, QuadratureType::Gauss);
        for( unsigned int i = 0; i < n; ++i )
        {
          QuadraturePoint<Field,1> q( points[i].position()[0], points[i].weight() );
          Base::push_back( q );
        }
      }
      // order is the maximal order of polynomials which can be exactly integrated
      static unsigned int minPoints( unsigned int order )
      {
        return (order+2)/2;
      }
    };

#if HAVE_ALGLIB
    template< unsigned int precision >
    class GaussPoints< amp::ampf< precision > >
      : public std::vector< QuadraturePoint<amp::ampf< precision>,1> >
    {
      typedef amp::ampf< precision > F;
      typedef std::vector< QuadraturePoint<F,1> > Base;
      typedef GaussPoints< F > This;

    public:
      typedef F Field;

      explicit GaussPoints ( unsigned int n )
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

        gqgengauss::generategaussquadrature< precision >( alpha, beta, Field( 2 ), n, p, w );

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
        return (order+2)/2;
      }
    };
#endif

    /**
     * @brief Singleton provider for Gauss quadratures
     *
     * \tparam dim dimension of the reference elements contained in the factory
     * \tparam F field in which weight and point of the quadrature are stored
     * \tparam CF the compute field for the points and weights
     **/
    template< int dim, class F, class CF=F >
    struct GaussQuadratureProvider
      : public TopologySingletonFactory< GenericQuadratureFactory< dim, F, GaussPoints<CF> > >
    {};

  }

}

#endif // #ifndef DUNE_GAUSSQUADRATURE_HH
