// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_MAPPINGS_HH
#define DUNE_GENERICGEOMETRY_MAPPINGS_HH

#include <dune/grid/genericgeometry/misc.hh>
#include <dune/grid/genericgeometry/topologytypes.hh>
#include <dune/grid/genericgeometry/conversion.hh>
#include <dune/grid/genericgeometry/referenceelements.hh>
#include <dune/grid/genericgeometry/matrix.hh>
#include <dune/grid/genericgeometry/submapping.hh>
#include <dune/grid/genericgeometry/cornermapping.hh>
#include <dune/grid/genericgeometry/hybridmapping.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // Mapping
    // -------

    template< class BaseMapping >
    class Mapping;

    template< class Topology, class CoordTraits >
    class Mapping< CornerMapping< Topology, CoordTraits > >
      : public CornerMapping< Topology, CoordTraits >
    {
      typedef CornerMapping< Topology, CoordTraits > BaseType;
      typedef Mapping< BaseType > ThisType;

    public:
      typedef typename BaseType :: Traits Traits;

      enum { dimG = Traits :: dimG };
      enum { dimW = Traits :: dimW };
      typedef typename Traits :: FieldType FieldType;
      typedef typename Traits :: LocalCoordType LocalCoordType;
      typedef typename Traits :: GlobalCoordType GlobalCoordType;
      typedef typename Traits :: JacobianType JacobianType;
      typedef typename Traits :: JacobianTransposedType JacobianTransposedType;

      typedef typename BaseType :: ReferenceElement ReferenceElement;

    protected:
      mutable JacobianType jTInv_;
      mutable FieldType intEl_;
      mutable GlobalCoordType faceNormal_[Size<Topology,1>::value];
      mutable bool jTInvComputed, intElComputed;
      mutable FieldVector<bool,Size<Topology,1>::value> normalComputed;

      using BaseType :: jT_;
      using BaseType :: jTComputed;

    public:
      template< class CoordVector >
      explicit Mapping ( const CoordVector &coords )
        : BaseType( coords ),
          jTInvComputed( false ),
          intElComputed( false ),
          normalComputed( false )
      {}

      using BaseType :: operator[];
      using BaseType :: corners;
      using BaseType :: affine;
      using BaseType :: global;
      using BaseType :: jacobianT;

      GeometryType type () const
      {
        return DuneGeometryType< Topology,
            Dune::GeometryType::BasicType(CoordTraits::oneDType) > :: type();
      }

      // additional methods
      LocalCoordType local ( const GlobalCoordType &p ) const
      {
        LocalCoordType x;
        GlobalCoordType y = p - (*this)[ 0 ];
        if( jTInvComputed )
          MatrixHelper< CoordTraits > :: template ATx< dimW, dimG >( jTInv_, y, x );
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

      const JacobianType &jacobianInverseTransposed ( const LocalCoordType &x ) const
      {
        if( !jTInvComputed )
        {
          const JacobianTransposedType &JT = jacobianT( x );
          intEl_ = MatrixHelper< CoordTraits >
                   :: template rightInvA< dimG, dimW >( JT, jTInv_ );
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
          intEl_ = MatrixHelper< CoordTraits > :: template detAAT< dimG, dimW >( JT );
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
          MatrixHelper< CoordTraits >
          :: template Ax< dimW, dimG >( JT, refNormal, faceNormal_[ face ] );
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
        MatrixHelper< CoordTraits > :: template xTRightInvA< dimG, dimW >( JT, p, y );
      }
    };


    // CachedMapping
    // -------------

    template< class Topology, class CoordTraits,
        template< class > class Caching = ComputeAll >
    class CachedMapping
      : public Mapping< CornerMapping< Topology, CoordTraits > >
    {
      typedef Mapping< CornerMapping< Topology,CoordTraits > > BaseType;
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

    protected:
      using BaseType :: baryCenter;
      using BaseType :: jT_;
      using BaseType :: jTInv_;
      using BaseType :: intEl_;
      using BaseType :: jTComputed;
      using BaseType :: jTInvComputed;
      using BaseType :: intElComputed;

    public:
      template< class CoordVector >
      explicit CachedMapping ( const CoordVector &coords,
                               const CachingType &cache = CachingType() )
        : BaseType( coords )
      {
        if( affine() )
        {
          if( (int)CachingType :: jTCompute == (int)geoIsComputed )
          {
            cache.jacobianT( jT_ );
            jTComputed = true;
          }
          else if( (int)CachingType :: jTCompute == (int)geoPreCompute )
            BaseType :: jacobianT( baryCenter() );

          if( (int)CachingType :: jTInvCompute == (int)geoIsComputed )
          {
            cache.jacobianInverseTransposed( jTInv_ );
            jTInvComputed = true;
          }
          else if( (int)CachingType :: jTInvCompute == (int)geoPreCompute )
            BaseType :: jacobianInverseTransposed( baryCenter() );

          if( (int)CachingType :: intElCompute == (int)geoIsComputed )
          {
            cache.integrationElement( intEl_ );
            intElComputed = true;
          }
          else if( (int)CachingType :: intElCompute == (int)geoPreCompute )
            integrationElement( baryCenter() );
        }
      }

      using BaseType :: affine;
      using BaseType :: operator[];
      using BaseType :: global;
      using BaseType :: local;
      using BaseType :: volume;
      using BaseType :: normal;

      const JacobianTransposedType &jacobianT ( const LocalCoordType &x ) const
      {
        if( ((int)CachingType :: jTCompute == (int)geoCompute) || !affine() )
          BaseType :: jacobianT( x );
        return jT_;
      }

      // additional methods
      FieldType integrationElement ( const LocalCoordType &x ) const
      {
        if( ((int)CachingType :: intElCompute == (int)geoCompute) || !affine() )
          BaseType :: integrationElement(x);
        return this->intEl_;
      }

      const JacobianType &jacobianInverseTransposed ( const LocalCoordType &x ) const
      {
        if( ((int)CachingType :: jTInvCompute == geoCompute) || !affine() )
        {
          BaseType::jacobianInverseTransposed(x);
        }
        return this->jTInv_;
      }

      template< unsigned int codim >
      typename Codim< codim > :: SubMapping *
      subMapping ( unsigned int i,
                   const typename Codim< codim > :: CachingType &cache
                     = typename Codim< codim > :: CachingType() ) const
      {
        return SubMappingProvider< ThisType, codim > :: subMapping( *this, i, cache );
      }
    };

  }

}

#endif
