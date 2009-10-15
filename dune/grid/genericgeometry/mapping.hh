// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_MAPPING_HH
#define DUNE_GENERICGEOMETRY_MAPPING_HH

#include <dune/grid/genericgeometry/topologytypes.hh>
#include <dune/grid/genericgeometry/referenceelements.hh>
#include <dune/grid/genericgeometry/matrix.hh>
#include <dune/grid/genericgeometry/geometrytraits.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // Mapping
    // -------

    /** \class Mapping
     *  \ingroup GenericGeometry
     *  \brief interface for a mapping
     *
     *  \tparam  CoordTraits  coordinate traits
     *  \tparam  Topology     topology of the reference domain
     *  \tparam  dimW         dimension of the world
     *  \tparam  Impl         implementation of the mapping
     */
    template< class CoordTraits, class Topo, int dimW, class Impl >
    class Mapping
    {
      typedef Mapping< CoordTraits, Topo, dimW, Impl > This;

      typedef Impl Implementation;

    public:
      typedef Topo Topology;
      typedef MappingTraits< CoordTraits, Topology :: dimension, dimW > Traits;

      static const unsigned int dimension = Traits :: dimension;
      static const unsigned int dimWorld = Traits :: dimWorld;

      typedef typename Traits :: FieldType FieldType;
      typedef typename Traits :: LocalCoordType LocalCoordType;
      typedef typename Traits :: GlobalCoordType GlobalCoordType;
      typedef typename Traits :: JacobianType JacobianType;
      typedef typename Traits :: JacobianTransposedType JacobianTransposedType;

      typedef typename Traits :: MatrixHelper MatrixHelper;

      typedef GenericGeometry :: ReferenceElement< Topology, FieldType > ReferenceElement;

      template< unsigned int codim, unsigned int i >
      struct SubTopology
      {
        typedef typename GenericGeometry :: SubTopology< Topo, codim, i > :: type Topology;
        typedef typename Implementation :: template SubTopology< codim, i > :: Trace TraceImpl;
        typedef Mapping< CoordTraits, Topology, dimWorld, TraceImpl > Trace;
      };

      static const bool alwaysAffine = Implementation :: alwaysAffine;

    protected:
      Implementation impl_;

    public:
      template< class CoordVector >
      explicit Mapping ( const CoordVector &coords )
        : impl_( coords )
      {}

      Mapping ( const Implementation &implementation )
        : impl_( implementation )
      {}

      const GlobalCoordType &corner ( int i ) const
      {
        return implementation().corner( i );
      }

      void global ( const LocalCoordType &x, GlobalCoordType &y ) const
      {
        implementation().global( x, y );
      }

      void local ( const GlobalCoordType &y, LocalCoordType &x ) const
      {
        x = ReferenceElement :: template baryCenter< 0 >( 0 );
        LocalCoordType dx;
        do
        { // DF^n dx^n = F^n, x^{n+1} -= dx^n
          JacobianTransposedType JT;
          jacobianTransposed( x, JT );
          GlobalCoordType z;
          global( x, z );
          z -= y;
          MatrixHelper :: template xTRightInvA< dimension, dimWorld >( JT, z, dx );
          x -= dx;
        } while( dx.two_norm2() > 1e-12 ); // use some numeric_limit eps?
      }

      bool jacobianTransposed ( const LocalCoordType &x,
                                JacobianTransposedType &JT ) const
      {
        return implementation().jacobianTransposed( x, JT );
      }

      FieldType
      jacobianInverseTransposed ( const LocalCoordType &x, JacobianType &JTInv ) const
      {
        JacobianTransposedType JT;
        jacobianTransposed( x, JT );
        return MatrixHelper :: template rightInvA< dimension, dimWorld >( JT, JTInv );
      }

      FieldType integrationElement ( const LocalCoordType &x ) const
      {
        JacobianTransposedType JT;
        jacobianTransposed( x, JT );
        return MatrixHelper :: template detAAT< dimension, dimWorld >( JT );
      }

      const Implementation &implementation () const
      {
        return impl_;
      }

      template< unsigned int codim, unsigned int i >
      typename SubTopology< codim, i > :: Trace trace () const
      {
        return impl_.template trace< codim, i >();
      }
    };

  }

}

#endif
