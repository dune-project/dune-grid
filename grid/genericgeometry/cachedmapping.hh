// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_MAPPINGS_HH
#define DUNE_GENERICGEOMETRY_MAPPINGS_HH
#include <dune/common/smallobject.hh>

#include <dune/grid/genericgeometry/misc.hh>
#include <dune/grid/genericgeometry/topologytypes.hh>
#include <dune/grid/genericgeometry/referenceelements.hh>
#include <dune/grid/genericgeometry/matrix.hh>
#include <dune/grid/genericgeometry/mapping.hh>
#include <dune/grid/genericgeometry/traceprovider.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // CachedMapping
    // -------------

    template< class Topology, class GeometryTraits >
    class CachedMapping
      : public SmallObject
    {
      typedef CachedMapping< Topology, GeometryTraits > This;

      typedef typename GeometryTraits :: template Mapping< Topology > :: type
      MappingImpl;

    public:
      typedef MappingTraits
      < typename GeometryTraits :: CoordTraits,
          Topology :: dimension, GeometryTraits :: dimWorld >
      Traits;

      typedef GenericGeometry :: Mapping
      < typename GeometryTraits :: CoordTraits,
          Topology, GeometryTraits :: dimWorld, MappingImpl >
      Mapping;

      static const unsigned int dimension = Traits :: dimension;
      static const unsigned int dimWorld = Traits :: dimWorld;

      typedef typename Traits :: FieldType FieldType;
      typedef typename Traits :: LocalCoordType LocalCoordType;
      typedef typename Traits :: GlobalCoordType GlobalCoordType;
      typedef typename Traits :: JacobianType JacobianType;
      typedef typename Traits :: JacobianTransposedType JacobianTransposedType;

      typedef GenericGeometry :: ReferenceElement< Topology, FieldType > ReferenceElement;

      static const bool alwaysAffine = Mapping :: alwaysAffine;

      template< unsigned int codim >
      struct Codim
      {
        typedef typename TraceProvider< Topology, GeometryTraits, codim, false > :: Trace
        Trace;
      };

      typedef typename GeometryTraits::Caching Caching;

    private:
      typedef typename Traits :: MatrixHelper MatrixHelper;

    public:
      unsigned int referenceCount;

    private:
      static const unsigned int numNormals = ReferenceElement :: numNormals;

      Mapping mapping_;

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
        : mapping_( coords ),
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

        if( affine() )
        {
          if( (Caching :: evaluateJacobianTransposed == PreCompute) && !jacobianTransposedComputed_ )
            computeJacobianTransposed( baryCenter() );

          if( Caching :: evaluateJacobianInverseTransposed == PreCompute )
            computeJacobianInverseTransposed( baryCenter() );
          else if( Caching :: evaluateIntegrationElement == PreCompute )
            computeIntegrationElement( baryCenter() );
        }
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

      static bool checkInside ( const LocalCoordType &x )
      {
        return ReferenceElement :: checkInside( x );
      }

      bool affine () const
      {
        return (alwaysAffine || affine_);
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
          mapping_.global( x, y );
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
        else if( affine() )
        {
          const JacobianTransposedType &JT = jacobianTransposed( baryCenter() );
          GlobalCoordType z = y - corner( 0 );
          MatrixHelper :: template xTRightInvA< dimension, dimWorld >( JT, z, x );
        }
        else
          mapping_.local( y, x );
        return x;
      }

      const JacobianTransposedType &jacobianTransposed ( const LocalCoordType &x ) const
      {
        const EvaluationType evaluate = Caching :: evaluateJacobianTransposed;
        if( (evaluate == PreCompute) && alwaysAffine )
          return jacobianTransposed_;

        if( !jacobianTransposedComputed_ )
          computeJacobianTransposed( x );
        return jacobianTransposed_;
      }

      FieldType integrationElement ( const LocalCoordType &x ) const
      {
        const EvaluationType evaluateI = Caching :: evaluateIntegrationElement;
        const EvaluationType evaluateJ = Caching :: evaluateJacobianInverseTransposed;
        if( ((evaluateI == PreCompute) || (evaluateJ == PreCompute)) && alwaysAffine )
          return integrationElement_;

        if( !integrationElementComputed_ )
          computeIntegrationElement( x );
        return integrationElement_;
      }

      const JacobianType &jacobianInverseTransposed ( const LocalCoordType &x ) const
      {
        const EvaluationType evaluate = Caching :: evaluateJacobianInverseTransposed;
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
          normalComputed_[ face ] = affine();
        }
        return normal_[ face ];
      }

      template< unsigned int codim, bool hybrid >
      typename TraceProvider< Topology, GeometryTraits, codim, hybrid > :: Trace *
      trace ( unsigned int i ) const
      {
        typedef TraceProvider< Topology, GeometryTraits, codim, hybrid > Provider;
        return Provider :: trace( mapping_, i );
      }

      template< unsigned int codim >
      typename Codim< codim > :: Trace *trace ( unsigned int i ) const
      {
        return trace< codim, false >( i );
      }

    private:
      static const LocalCoordType &baryCenter ()
      {
        return ReferenceElement :: template baryCenter< 0 >( 0 );
      }

      void computeJacobianTransposed ( const LocalCoordType &x ) const
      {
        affine_ = mapping_.jacobianTransposed( x, jacobianTransposed_ );
        jacobianTransposedComputed_ = affine_;
      }

      void computeJacobianInverseTransposed ( const LocalCoordType &x ) const
      {
        integrationElement_ = mapping_.jacobianInverseTransposed( x, jacobianInverseTransposed_ );
        integrationElementComputed_ = jacobianInverseTransposedComputed_ = affine();
      }

      void computeIntegrationElement ( const LocalCoordType &x ) const
      {
        integrationElement_ = mapping_.integrationElement( x );
        integrationElementComputed_ = affine();
      }
    };

  }

}

#endif
