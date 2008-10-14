// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_MAPPINGS_HH
#define DUNE_GENERICGEOMETRY_MAPPINGS_HH

#include <dune/common/smallobject.hh>

#include <dune/grid/genericgeometry/misc.hh>
#include <dune/grid/genericgeometry/topologytypes.hh>
#include <dune/grid/genericgeometry/referenceelements.hh>
#include <dune/grid/genericgeometry/matrix.hh>
#include <dune/grid/genericgeometry/submapping.hh>
#include <dune/grid/genericgeometry/cornermapping.hh>
#include <dune/grid/genericgeometry/geometricmapping.hh>
#include <dune/grid/genericgeometry/hybridmapping.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    template< class Topology, class GeometricMappingTraits >
    class CachedMapping;



    // CachedMapping
    // -------------

    template< class Topology, class GeometricMappingTraits >
    class CachedMapping
      : public GeometricMapping< Topology, GeometricMappingTraits >,
        public SmallObject
    {
      typedef GeometricMapping< Topology, GeometricMappingTraits > Base;
      typedef CachedMapping< Topology, GeometricMappingTraits > This;

      typedef MappingTraits
      < typename GeometricMappingTraits :: CoordTraits,
          Topology :: dimension, GeometricMappingTraits :: dimWorld >
      Traits;

      typedef typename Traits :: MatrixHelper MatrixHelper;

      static const bool alwaysAffine = Base :: alwaysAffine;

    public:
      static const unsigned int dimension = Traits :: dimension;
      static const unsigned int dimWorld = Traits :: dimWorld;

      typedef typename Traits :: FieldType FieldType;
      typedef typename Traits :: LocalCoordType LocalCoordType;
      typedef typename Traits :: GlobalCoordType GlobalCoordType;
      typedef typename Traits :: JacobianType JacobianType;
      typedef typename Traits :: JacobianTransposedType JacobianTransposedType;

      typedef typename GeometricMappingTraits :: Caching CachingType;
      typedef typename Base :: ReferenceElement ReferenceElement;

      template< unsigned int codim >
      struct Codim
      {
        typedef typename SubMappingTraits< This, codim > :: SubMapping SubMapping;
      };

      template< unsigned int codim, unsigned int i >
      struct SubTopology
      {
        typedef typename GenericGeometry :: SubTopology< Topology, codim, i > :: type type;
        typedef CachedMapping< type, GeometricMappingTraits > Trace;
      };

    public:
      unsigned int referenceCount;

    protected:
      using Base :: mapping_;

    private:
      static const unsigned int numNormals = ReferenceElement :: numNormals;

      mutable JacobianTransposedType jacobianTransposed_;
      mutable JacobianType jacobianInverseTransposed_;
      mutable FieldType integrationElement_;
      mutable array< GlobalCoordType, numNormals > normal_;

      mutable bool affine_;

      mutable bool jacobianTransposedComputed_;
      mutable bool jacobianInverseTransposedComputed_;
      mutable bool integrationElementComputed_;
      mutable array< bool, numNormals > normalComputed_;

    public:
      template< class CoordVector >
      explicit CachedMapping ( const CoordVector &coords )
        : Base( coords ),
          jacobianTransposedComputed_( false ),
          jacobianInverseTransposedComputed_( false ),
          integrationElementComputed_( false )
      {
        for( unsigned int i = 0; i < numNormals; ++i )
          normalComputed_[ i ] = false;

        if( alwaysAffine )
          affine_ = true;
        else
          computeJacobianTransposed( baryCenter() );

        if( affine_ )
        {
          if( (CachingType :: evaluateJacobianTransposed == PreCompute) && !jacobianTransposedComputed_ )
            computeJacobianTransposed( baryCenter() );

          if( CachingType :: evaluateJacobianInverseTransposed == PreCompute )
            computeJacobianInverseTransposed( baryCenter() );
          else if( CachingType :: evaluateIntegrationElement == PreCompute )
            computeIntegrationElement( baryCenter() );
        }
      }

      unsigned int topologyId () const
      {
        return ReferenceElement :: topologyId;
      }

      using Base :: corner;

      int numCorners () const
      {
        return ReferenceElement :: numCorners;
      }

      static bool checkInside ( const LocalCoordType &x )
      {
        return ReferenceElement :: checkInside( x );
      }

      bool affine () const
      {
        return affine_;
      }

      GlobalCoordType global ( const LocalCoordType &x ) const
      {
        GlobalCoordType y;
        if( jacobianTransposedComputed_ )
        {
          MatrixHelper :: template ATx< dimension, dimWorld >( jacobianTransposed_, x, y );
          y += corner( 0 );
        }
        else
          Base :: global( x, y );
        return y;
      }

      LocalCoordType local ( const GlobalCoordType &y ) const
      {
        LocalCoordType x;
        if( jacobianInverseTransposedComputed_ )
        {
          GlobalCoordType z = y - corner( 0 );
          MatrixHelper :: template ATx< dimWorld, dimension >( jacobianInverseTransposed_, z, x );
        }
        else if( affine_ )
        {
          const JacobianTransposedType &JT = jacobianTransposed( baryCenter() );
          GlobalCoordType z = y - corner( 0 );
          MatrixHelper :: template xTRightInvA< dimension, dimWorld >( JT, z, x );
        }
        else
          Base :: local( y, x );
        return x;
      }

      const JacobianTransposedType &jacobianTransposed ( const LocalCoordType &x ) const
      {
        const EvaluationType evaluate = CachingType :: evaluateJacobianTransposed;
        if( (evaluate == PreCompute) && alwaysAffine )
          return jacobianTransposed_;

        if( !jacobianTransposedComputed_ )
          computeJacobianTransposed( x );
        return jacobianTransposed_;
      }

      // additional methods
      FieldType integrationElement ( const LocalCoordType &x ) const
      {
        const EvaluationType evaluateI = CachingType :: evaluateIntegrationElement;
        const EvaluationType evaluateJ = CachingType :: evaluateJacobianInverseTransposed;
        if( ((evaluateI == PreCompute) || (evaluateJ == PreCompute)) && alwaysAffine )
          return integrationElement_;

        if( !integrationElementComputed_ )
          computeIntegrationElement( x );
        return integrationElement_;
      }

      const JacobianType &jacobianInverseTransposed ( const LocalCoordType &x ) const
      {
        const EvaluationType evaluate = CachingType :: evaluateJacobianInverseTransposed;
        if( (evaluate == PreCompute) && alwaysAffine )
          return jacobianInverseTransposed_;

        if( !jacobianInverseTransposedComputed_ )
          computeJacobianInverseTransposed( x );
        return jacobianInverseTransposed_;
      }

      FieldType volume () const
      {
        const FieldType refVolume = ReferenceElement :: volume();
        return refVolume * integrationElement( baryCenter() );
      }

      const GlobalCoordType &normal ( int face, const LocalCoordType &x ) const
      {
        assert( (unsigned int)face < numNormals );
        if( !normalComputed_[ face ] )
        {
          const JacobianType &JT = jacobianInverseTransposed( x );
          const LocalCoordType &refNormal =  ReferenceElement :: integrationOuterNormal( face );
          MatrixHelper :: template Ax< dimWorld, dimension >( JT, refNormal, normal_[ face ] );
          normal_[ face ] *= integrationElement_;
          normalComputed_[ face ] = affine_;
        }
        return normal_[ face ];
      }

      template< unsigned int codim >
      typename Codim< codim > :: SubMapping *subMapping ( unsigned int i ) const
      {
        return SubMappingProvider< This, codim > :: subMapping( *this, i );
      }

      template< unsigned int codim, unsigned int i >
      typename SubTopology< codim, i > :: Trace
      trace () const
      {
        typedef typename SubTopology< codim, i > :: Trace Trace;
        return Trace( mapping_.trace() );
      }

    private:
      static const LocalCoordType &baryCenter ()
      {
        return ReferenceElement :: template baryCenter< 0 >( 0 );
      }

      void computeJacobianTransposed ( const LocalCoordType &x ) const
      {
        affine_ = Base :: jacobianTransposed( x, jacobianTransposed_ );
        jacobianTransposedComputed_ = affine_;
      }

      void computeJacobianInverseTransposed ( const LocalCoordType &x ) const
      {
        integrationElement_ = Base :: jacobianInverseTransposed( x, jacobianInverseTransposed_ );
        integrationElementComputed_ = jacobianInverseTransposedComputed_ = affine_;
      }

      void computeIntegrationElement ( const LocalCoordType &x ) const
      {
        integrationElement_ = Base :: integrationElement( x );
        integrationElementComputed_ = affine_;
      }
    };

  }

}

#endif
