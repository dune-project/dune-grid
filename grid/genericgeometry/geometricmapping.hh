// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_GEOMETRICMAPPING_HH
#define DUNE_GENERICGEOMETRY_GEOMETRICMAPPING_HH

#include <dune/common/smallobject.hh>

#include <dune/grid/genericgeometry/misc.hh>
#include <dune/grid/genericgeometry/topologytypes.hh>
#include <dune/grid/genericgeometry/referenceelements.hh>
#include <dune/grid/genericgeometry/matrix.hh>
#include <dune/grid/genericgeometry/submapping.hh>
#include <dune/grid/genericgeometry/cornermapping.hh>
#include <dune/grid/genericgeometry/hybridmapping.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // GeometricMapping
    // ----------------

    /** \class GeometricMapping
     *  \ingroup GenericGeometry
     *  \brief add geometric functionality to a mapping
     *
     *  \tparam  Topology                topology of the reference domain
     *  \tparam  GeometricMappingTraits  structure containing required types
     *
     *  The GeometricMappingTraits for a mapping named <tt>MyMapping</tt> could
     *  look as follows:
     *  \code
     *  struct MyMappingTraits
     *  {
     *    typedef MyCoordTraits CoordTraits;
     *
     *    template< unsigned int dimension >
     *    struct Traits
     *    : public MappingTraits< dimension, CoordTraits >
     *    {
     *      typedef MyCaching< dimension, CoordTraits > CachingType;
     *    };
     *
     *    template< class Topology >
     *    struct Mapping
     *    {
     *      typedef MyMapping< Topology, CoordTraits > Type;
     *    };
     *  };
     *  \endcode
     *
     *  The mapping (called <tt>MyMapping</tt> here) must provide the following
     *  types and methods:
     *  \code
     *  template< unsigned int codim, unsigned int i >
     *  struct SubTopology
     *  {
     *    typedef MyTrace< MyMapping, codim, i > Trace;
     *  };
     *
     *  template< class CoordVector >
     *  explicit MyMapping ( const CoordVector &coords );
     *
     *  const GlobalCoordType &corner ( int i ) const;
     *
     *  void global ( const LocalCoordType &x, GlobalCoordType &ret ) const;
     *  bool jacobianTransposed ( const LocalCoordType &x, JacobianTransposedType &ret ) const;
     *
     *  template< unsigned int codim, unsigned int i >
     *  typename SubTopology< codim, i > :: Trace trace () const;
     *  \endcode
     */
    template< class Topology, class GeometricMappingTraits >
    class GeometricMapping
    {
      typedef GeometricMapping< Topology, GeometricMappingTraits > This;

    public:
      typedef typename GeometricMappingTraits :: template Mapping< Topology > :: type
      Mapping;
      typedef typename Mapping :: Traits Traits;

      static const unsigned int dimension = Traits :: dimension;
      static const unsigned int dimWorld = Traits :: dimWorld;

      typedef typename Traits :: FieldType FieldType;
      typedef typename Traits :: LocalCoordType LocalCoordType;
      typedef typename Traits :: GlobalCoordType GlobalCoordType;
      typedef typename Traits :: JacobianType JacobianType;
      typedef typename Traits :: JacobianTransposedType JacobianTransposedType;

      typedef typename Traits :: MatrixHelper MatrixHelper;

      typedef GenericGeometry :: ReferenceElement< Topology, FieldType > ReferenceElement;

      static const bool alwaysAffine = Mapping :: alwaysAffine;

    protected:
      Mapping mapping_;

    public:
      template< class CoordVector >
      explicit GeometricMapping ( const CoordVector &coords )
        : mapping_( coords )
      {}

      const GlobalCoordType &corner ( int i ) const
      {
        return mapping_.corner( i );
      }

      void global ( const LocalCoordType &x, GlobalCoordType &y ) const
      {
        mapping_.global( x, y );
      }

      bool jacobianTransposed ( const LocalCoordType &x,
                                JacobianTransposedType &JT ) const
      {
        return mapping_.jacobianTransposed( x, JT );
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

    protected:
      static const LocalCoordType &baryCenter ()
      {
        return ReferenceElement :: template baryCenter< 0 >( 0 );
      }
    };

  }

}

#endif
