// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMTRY_HYBRIDMAPPING_HH
#define DUNE_GENERICGEOMTRY_HYBRIDMAPPING_HH

#include <dune/common/smallobject.hh>

#include <dune/grid/genericgeometry/submapping.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // External Forward Declarations
    // -----------------------------

    template< class Topology, class GeometricMappingTraits >
    class CachedMapping;



    // HybridMapping
    // -------------

    template< unsigned int dim, class GeometricMappingTraits >
    class HybridMapping
      : public SmallObject
    {
      typedef HybridMapping< dim, GeometricMappingTraits > This;

    protected:
      typedef MappingTraits
      < typename GeometricMappingTraits :: CoordTraits,
          dim, GeometricMappingTraits :: dimWorld >
      Traits;

    public:
      static const unsigned int dimG = Traits :: dimG;
      static const unsigned int dimW = Traits :: dimW;

      typedef typename Traits :: FieldType FieldType;
      typedef typename Traits :: LocalCoordType LocalCoordType;
      typedef typename Traits :: GlobalCoordType GlobalCoordType;
      typedef typename Traits :: JacobianType JacobianType;

      template< unsigned int codim >
      struct Codim
      {
        typedef typename SubMappingTraits< This, codim > :: SubMapping SubMapping;
      };

      unsigned int referenceCount;

      virtual ~HybridMapping ()
      {}

      virtual unsigned int topologyId () const = 0;

      virtual const GlobalCoordType &corner ( int i ) const = 0;

      virtual int numCorners () const = 0;

      virtual GlobalCoordType global ( const LocalCoordType &local ) const = 0;

      virtual LocalCoordType local ( const GlobalCoordType &global ) const = 0;

      virtual bool checkInside ( const LocalCoordType &local ) const = 0;

      virtual bool affine () const = 0;

      virtual FieldType integrationElement ( const LocalCoordType &local ) const = 0;

      virtual FieldType volume () const = 0;

      virtual const JacobianType &
      jacobianInverseTransposed ( const LocalCoordType &local ) const = 0;

      virtual GlobalCoordType
      normal ( int face, const LocalCoordType &local ) const = 0;

      template< unsigned int codim >
      typename Codim< codim > :: SubMapping *
      subMapping ( unsigned int i ) const
      {
        typedef typename Codim< codim > :: SubMapping SubMapping;
        return reinterpret_cast< SubMapping * >( subMapping( codim, i ) );
      }

    protected:
      virtual void *
      subMapping ( unsigned int codim, unsigned int i ) const = 0;
    };



    template< class Topology, class GeometricMappingTraits >
    class VirtualMapping
      : public HybridMapping< Topology :: dimension, GeometricMappingTraits >
    {
      typedef HybridMapping< Topology :: dimension, GeometricMappingTraits > Base;
      typedef VirtualMapping< Topology, GeometricMappingTraits > This;

      typedef typename Base :: Traits Traits;

      typedef CachedMapping< Topology, GeometricMappingTraits > Mapping;

    public:
      enum { dimG = Traits :: dimG };
      enum { dimW = Traits :: dimW };

      typedef typename Traits :: FieldType FieldType;
      typedef typename Traits :: LocalCoordType LocalCoordType;
      typedef typename Traits :: GlobalCoordType GlobalCoordType;
      typedef typename Traits :: JacobianType JacobianType;

      typedef typename Mapping :: ReferenceElement ReferenceElement;

      template< unsigned int codim >
      struct Codim
      {
        typedef typename SubMappingTraits< This, codim > :: SubMapping SubMapping;
      };

    private:
      Mapping mapping_;

    public:
      template< class CoordVector >
      explicit VirtualMapping ( const CoordVector &coordVector )
        : mapping_( coordVector )
      {}

      virtual unsigned int topologyId () const
      {
        return mapping_.topologyId();
      }

      virtual const GlobalCoordType &corner ( int i ) const
      {
        return mapping_.corner( i );
      }

      virtual int numCorners () const
      {
        return mapping_.numCorners();
      }

      virtual GlobalCoordType global ( const LocalCoordType &local ) const
      {
        return mapping_.global( local );
      }

      virtual LocalCoordType local ( const GlobalCoordType &global ) const
      {
        return mapping_.local( global );
      }

      virtual bool checkInside ( const LocalCoordType &local ) const
      {
        return mapping_.checkInside( local );
      }

      virtual bool affine () const
      {
        return mapping_.affine();
      }

      virtual FieldType integrationElement ( const LocalCoordType &local ) const
      {
        return mapping_.integrationElement( local );
      }

      virtual FieldType volume () const
      {
        return mapping_.volume();
      }

      virtual const JacobianType &
      jacobianInverseTransposed ( const LocalCoordType &local ) const
      {
        return mapping_.jacobianInverseTransposed( local );
      }

      virtual GlobalCoordType
      normal ( int face, const LocalCoordType &local ) const
      {
        return mapping_.normal( face , local );
      }

      template< unsigned int codim >
      typename Codim< codim > :: SubMapping *
      subMapping ( unsigned int i ) const
      {
        return SubMappingProvider< This, codim > :: subMapping( *this, i );
      }

    private:
      class CodimCaller;

    protected:
      void *subMapping ( unsigned int codim, unsigned int i ) const
      {
        return CodimCaller :: subMapping ( *this, codim, i );
      }
    };


    template< class Topology, class GeometricMappingTraits >
    class VirtualMapping< Topology, GeometricMappingTraits > :: CodimCaller
    {
      typedef VirtualMapping< Topology, GeometricMappingTraits > Mapping;

      struct CallerInterface;
      template< int codim > struct CallerImplementation;

      CallerInterface *caller_[ dimG+1 ];

      CodimCaller ()
      {
        ForLoop< CallerImplementation, 0, dimG > :: apply( caller_ );
      }

      static const CodimCaller &instance ()
      {
        static CodimCaller inst;
        return inst;
      }

    public:
      static void *
      subMapping ( const Mapping &mapping, unsigned int codim, unsigned int i )
      {
        assert( codim <= dimG );
        return instance().caller_[ codim ]->subMapping( mapping, i );
      }
    };


    template< class Topology, class GeometricMappingTraits >
    struct VirtualMapping< Topology, GeometricMappingTraits > :: CodimCaller :: CallerInterface
    {
      typedef VirtualMapping< Topology, GeometricMappingTraits > Mapping;

      virtual ~CallerInterface ()
      {}

      virtual void *subMapping ( const Mapping &mapping, unsigned int i ) const = 0;
    };

    template< class Topology, class GeometricMappingTraits >
    template< int codim >
    struct VirtualMapping< Topology, GeometricMappingTraits > :: CodimCaller :: CallerImplementation
      : public CallerInterface
    {
      virtual void *
      subMapping ( const Mapping &mapping, unsigned int i ) const
      {
        return mapping.template subMapping< codim >( i );
      }

      static void apply ( CallerInterface *(&caller)[ dimG+1 ] )
      {
        caller[ codim ] = new CallerImplementation< codim >;
      }
    };

  }

}

#endif
