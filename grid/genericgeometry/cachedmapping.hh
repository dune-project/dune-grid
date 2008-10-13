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
      using Base :: baryCenter;
      using Base :: jacobianTransposed_;
      using Base :: jTInv_;
      using Base :: intEl_;
      using Base :: jacobianTransposedComputed_;
      using Base :: jTInvComputed;
      using Base :: intElComputed;

    public:
      template< class CoordVector >
      explicit CachedMapping ( const CoordVector &coords )
        : Base( coords )
      {
        if( affine() )
        {
          if( CachingType :: evaluateJacobianTransposed == PreCompute )
            Base :: jacobianT( baryCenter() );

          if( CachingType :: evaluateJacobianInverseTransposed == PreCompute )
            Base :: jacobianInverseTransposed( baryCenter() );

          if( CachingType :: evaluateIntegrationElement == PreCompute )
            integrationElement( baryCenter() );
        }
      }

      using Base :: affine;
      using Base :: corner;
      using Base :: global;
      using Base :: local;
      using Base :: volume;
      using Base :: normal;
      using Base :: jacobianInverseTransposed;

      const JacobianTransposedType &jacobianT ( const LocalCoordType &x ) const
      {
        const EvaluationType evaluate = CachingType :: evaluateJacobianTransposed;
        if( (evaluate == ComputeOnDemand) || !affine() )
          Base :: jacobianT( x );
        return jacobianTransposed_;
      }

      // additional methods
      FieldType integrationElement ( const LocalCoordType &x ) const
      {
        const EvaluationType evaluate = CachingType :: evaluateIntegrationElement;
        if( (evaluate == ComputeOnDemand) || !affine() )
          Base :: integrationElement(x);
        return this->intEl_;
      }

      const JacobianType &jacobianInverseTransposed ( const LocalCoordType &x ) const
      {
        const EvaluationType evaluate = CachingType :: evaluateJacobianInverseTransposed;
        if( (evaluate == ComputeOnDemand) || !affine() )
          Base :: jacobianInverseTransposed(x);
        return this->jTInv_;
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
    };

  }

}

#endif
