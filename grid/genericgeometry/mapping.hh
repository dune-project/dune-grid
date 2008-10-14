// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_MAPPING_HH
#define DUNE_GENERICGEOMETRY_MAPPING_HH

#include <dune/common/smallobject.hh>

#include <dune/grid/genericgeometry/misc.hh>
#include <dune/grid/genericgeometry/topologytypes.hh>
#include <dune/grid/genericgeometry/referenceelements.hh>
#include <dune/grid/genericgeometry/matrix.hh>
#include <dune/grid/genericgeometry/submapping.hh>
#include <dune/grid/genericgeometry/hybridmapping.hh>

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
    template< class CoordTraits, class Topology, int dimW, class Impl >
    class Mapping
    {
      typedef Mapping< CoordTraits, Topology, dimW, Impl > This;

      typedef Impl Implementation;

    public:
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

      static const bool alwaysAffine = Implementation :: alwaysAffine;

    protected:
      Implementation impl_;

    public:
      template< class CoordVector >
      explicit Mapping ( const CoordVector &coords )
        : impl_( coords )
      {}

      const GlobalCoordType &corner ( int i ) const
      {
        return implementation().corner( i );
      }

      void global ( const LocalCoordType &x, GlobalCoordType &y ) const
      {
        implementation().global( x, y );
      }

      bool jacobianTransposed ( const LocalCoordType &x,
                                JacobianTransposedType &JT ) const
      {
        return implementation().jacobianTransposed( x, JT );
      }

      void local ( const GlobalCoordType &y, LocalCoordType &x ) const
      {
        x = baryCenter();
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
        } while( dx.two_norm2() > 1e-12 );
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

    protected:
      static const LocalCoordType &baryCenter ()
      {
        return ReferenceElement :: template baryCenter< 0 >( 0 );
      }
    };

  }

}

#endif
