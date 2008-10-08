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

    template< class CoordTraits >
    struct CornerMappingTraits
    {
      template< class Topology >
      struct Mapping
      {
        typedef MappingTraits< Topology :: dimension, CoordTraits > Traits;
        typedef CornerMapping< Topology, CoordTraits > Type;
      };
    };



    // CachedMapping
    // -------------

    template< class Topology, class CoordTraits,
        template< class > class Caching = ComputeAll >
    class CachedMapping
      : public GeometricMapping< Topology, CornerMappingTraits< CoordTraits > >,
        public SmallObject
    {
      typedef GeometricMapping< Topology, CornerMappingTraits< CoordTraits > >
      BaseType;
      typedef CachedMapping< Topology, CoordTraits, Caching > ThisType;

      typedef CachedMappingTraits< Topology :: dimension, CoordTraits, Caching > Traits;

    public:
      enum { dimG = Traits :: dimG };
      enum { dimW = Traits :: dimW };
      typedef typename Traits :: FieldType FieldType;
      typedef typename Traits :: LocalCoordType LocalCoordType;
      typedef typename Traits :: GlobalCoordType GlobalCoordType;
      typedef typename Traits :: JacobianType JacobianType;
      typedef typename Traits :: JacobianTransposedType JacobianTransposedType;

      typedef typename Traits :: CachingType CachingType;
      typedef typename BaseType :: ReferenceElement ReferenceElement;

      template< unsigned int codim >
      struct Codim
      {
        typedef typename SubMappingTraits< ThisType, codim > :: SubMapping SubMapping;
        typedef typename SubMappingTraits< ThisType, codim > :: CachingType CachingType;
      };

    public:
      unsigned int referenceCount;

    protected:
      using BaseType :: baryCenter;
      using BaseType :: jacobianTransposed_;
      using BaseType :: jTInv_;
      using BaseType :: intEl_;
      using BaseType :: jacobianTransposedComputed_;
      using BaseType :: jTInvComputed;
      using BaseType :: intElComputed;

    public:
      template< class CoordVector >
      inline explicit CachedMapping ( const CoordVector &coords,
                                      const CachingType &cache = CachingType() );

      using BaseType :: affine;
      using BaseType :: operator[];
      using BaseType :: global;
      using BaseType :: local;
      using BaseType :: volume;
      using BaseType :: normal;
      using BaseType :: jacobianInverseTransposed;

      const JacobianTransposedType &jacobianT ( const LocalCoordType &x ) const
      {
        const EvaluationType evaluate = CachingType :: evaluateJacobianTransposed;
        if( (evaluate == ComputeOnDemand) || !affine() )
          BaseType :: jacobianT( x );
        return jacobianTransposed_;
      }

      // additional methods
      FieldType integrationElement ( const LocalCoordType &x ) const
      {
        const EvaluationType evaluate = CachingType :: evaluateIntegrationElement;
        if( (evaluate == ComputeOnDemand) || !affine() )
          BaseType :: integrationElement(x);
        return this->intEl_;
      }

      const JacobianType &jacobianInverseTransposed ( const LocalCoordType &x ) const
      {
        const EvaluationType evaluate = CachingType :: evaluateJacobianInverseTransposed;
        if( (evaluate == ComputeOnDemand) || !affine() )
          BaseType :: jacobianInverseTransposed(x);
        return this->jTInv_;
      }

      template< unsigned int codim >
      typename Codim< codim > :: SubMapping *
      subMapping ( unsigned int i,
                   const typename Codim< codim > :: CachingType &cache ) const
      {
        return SubMappingProvider< ThisType, codim > :: subMapping( *this, i, cache );
      }
    };



    template< class Topology, class CoordTraits, template< class > class Caching >
    template< class CoordVector >
    inline CachedMapping< Topology, CoordTraits, Caching >
    :: CachedMapping ( const CoordVector &coords, const CachingType &cache )
      : BaseType( coords )
    {
      if( affine() )
      {
        switch( CachingType :: evaluateJacobianTransposed )
        {
        case IsComputed :
          cache.jacobianT( jacobianTransposed_ );
          jacobianTransposedComputed_ = true;
          break;

        case PreCompute :
          BaseType :: jacobianT( baryCenter() );
          break;

        case ComputeOnDemand :
          break;
        }

        switch( CachingType :: evaluateJacobianInverseTransposed )
        {
        case IsComputed :
          cache.jacobianInverseTransposed( jTInv_ );
          jTInvComputed = true;
          break;

        case PreCompute :
          BaseType :: jacobianInverseTransposed( baryCenter() );
          break;

        case ComputeOnDemand :
          break;
        }

        switch( CachingType :: evaluateIntegrationElement )
        {
        case IsComputed :
          cache.integrationElement( intEl_ );
          intElComputed = true;
          break;

        case PreCompute :
          integrationElement( baryCenter() );
          break;

        case ComputeOnDemand :
          break;
        }
      }
    }

  }

}

#endif
