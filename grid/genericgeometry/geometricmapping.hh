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

      static const unsigned int dimension = Topology :: dimension;

    public:
      typedef typename GeometricMappingTraits :: template Mapping< Topology > :: Type
      Mapping;
      typedef typename Mapping :: Traits Traits;

      static const unsigned int dimG = Traits :: dimG;
      static const unsigned int dimW = Traits :: dimW;

      typedef typename Traits :: FieldType FieldType;
      typedef typename Traits :: LocalCoordType LocalCoordType;
      typedef typename Traits :: GlobalCoordType GlobalCoordType;
      typedef typename Traits :: JacobianType JacobianType;
      typedef typename Traits :: JacobianTransposedType JacobianTransposedType;

      typedef typename Traits :: MatrixHelper MatrixHelper;

      typedef GenericGeometry :: ReferenceElement< Topology, FieldType > ReferenceElement;

    protected:
      static const unsigned int numNormals = ReferenceElement :: numNormals;

      Mapping mapping_;

      mutable JacobianTransposedType jacobianTransposed_;
      mutable JacobianType jTInv_;
      mutable FieldType intEl_;
      mutable GlobalCoordType faceNormal_[ numNormals ];

      mutable bool affine_;
      mutable bool jacobianTransposedComputed_;
      mutable bool jTInvComputed, intElComputed;
      mutable FieldVector< bool, numNormals > normalComputed;

    public:
      template< class CoordVector >
      explicit GeometricMapping ( const CoordVector &coords )
        : mapping_( coords ),
          jacobianTransposedComputed_( false ),
          jTInvComputed( false ),
          intElComputed( false ),
          normalComputed( false )
      {
        if( !Mapping :: alwaysAffine )
        {
          jacobianT( baryCenter() );
          affine_ = jacobianTransposedComputed_;
        }
        else
          affine_ = true;
      }

      bool affine () const
      {
        return affine_;
      }

      unsigned int topologyId () const
      {
        return ReferenceElement :: topologyId;
      }

      const GlobalCoordType &corner ( int i ) const
      {
        return mapping_.corner( i );
      }

      int numCorners () const
      {
        return ReferenceElement :: numCorners;
      }

      GlobalCoordType global ( const LocalCoordType &x ) const
      {
        GlobalCoordType p;
        if( jacobianTransposedComputed_ )
        {
          MatrixHelper :: template ATx< dimG, dimW >( jacobianTransposed_, x, p );
          p += corner( 0 );
        }
        else
          mapping_.global( x, p );
        return p;
      }

      static bool checkInside ( const LocalCoordType &x )
      {
        return ReferenceElement :: checkInside( x );
      }

      const JacobianTransposedType &jacobianT ( const LocalCoordType &x ) const
      {
        if( !jacobianTransposedComputed_ )
        {
          jacobianTransposedComputed_
            = mapping_.jacobianTransposed( x, jacobianTransposed_ );
        }
        return jacobianTransposed_;
      }

      // additional methods
      LocalCoordType local ( const GlobalCoordType &p ) const
      {
        LocalCoordType x;
        GlobalCoordType y = p - corner( 0 );
        if( jTInvComputed )
          MatrixHelper :: template ATx< dimW, dimG >( jTInv_, y, x );
        else if( affine() )
          local_affine( baryCenter(), y, x );
        else
        {
          x = baryCenter();
          LocalCoordType dx;
          do
          { // DF^n dx^n = -F^n, x^{n+1} += dx^n
            y = p - global( x );
            local_affine( x, y, dx );
            x += dx;
          } while( dx.two_norm2() > 1e-12 );
        }
        return x;
      }

      FieldType
      jacobianInverseTransposed ( const LocalCoordType &x, JacobianType &jTInv ) const
      {
        const JacobianTransposedType &jT = jacobianT( x );
        return MatrixHelper :: template rightInvA< dimG, dimW >( jT, jTInv );
      }

      const JacobianType &jacobianInverseTransposed ( const LocalCoordType &x ) const
      {
        if( !jTInvComputed )
        {
          intEl_ = jacobianInverseTransposed( x, jTInv_ );
          jTInvComputed = affine();
          intElComputed = affine();
        }
        return jTInv_;
      }

      FieldType integrationElement ( const LocalCoordType &x ) const
      {
        if( !intElComputed )
        {
          const JacobianTransposedType &JT = jacobianT( x );
          intEl_ = MatrixHelper :: template detAAT< dimG, dimW >( JT );
          intElComputed = affine();
        }
        return intEl_;
      }

      const GlobalCoordType &normal ( int face, const LocalCoordType &x ) const
      {
        if( !normalComputed[ face ] )
        {
          const JacobianType &JT = jacobianInverseTransposed( x );
          const LocalCoordType &refNormal
            =  ReferenceElement :: integrationOuterNormal( face );
          MatrixHelper :: template Ax< dimW, dimG >( JT, refNormal, faceNormal_[ face ] );
          faceNormal_[ face ] *= intEl_;
          normalComputed[ face ] = affine();
        }
        return faceNormal_[ face ];
      }

      FieldType volume () const
      {
        const FieldType refVolume = ReferenceElement :: volume();
        return refVolume * integrationElement( baryCenter() );
      }

    protected:
      static const LocalCoordType &baryCenter ()
      {
        return ReferenceElement :: template baryCenter< 0 >( 0 );
      }

      // affine local to global mapping
      void local_affine ( const LocalCoordType &x,
                          const GlobalCoordType &p,
                          LocalCoordType &y ) const
      {
        const JacobianTransposedType &JT = jacobianT( x );
        MatrixHelper :: template xTRightInvA< dimG, dimW >( JT, p, y );
      }
    };

  }

}

#endif
