// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMTRY_HYBRIDMAPPING_HH
#define DUNE_GENERICGEOMTRY_HYBRIDMAPPING_HH

#include <dune/common/smallobject.hh>

#include <dune/grid/genericgeometry/geometrytraits.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // External Forward Declarations
    // -----------------------------

    template< class Topology, class GeometryTraits >
    class CachedMapping;



    // Internal Forward Declarations
    // -----------------------------

    template< unsigned int dim, class GeometryTraits >
    class HybridMapping;

    template< class Topology, class GeometryTraits >
    class VirtualMapping;



    // HybridMappingBase
    // -----------------

    template< unsigned int dim, class GeometryTraits, unsigned int codim = dim >
    class HybridMappingBase;

    template< unsigned int dim, class GeometryTraits, unsigned int codim >
    class HybridMappingBase
      : public virtual HybridMappingBase< dim, GeometryTraits, codim-1 >
    {
      typedef HybridMapping< dim, GeometryTraits > Mapping;

    public:
      virtual ~HybridMappingBase() {}

    protected:
      using HybridMappingBase< dim, GeometryTraits, codim-1 > :: trace;

      virtual HybridMapping< dim - codim, GeometryTraits > *
      trace ( Int2Type< codim >, unsigned int i ) const = 0;
    };

    template< unsigned int dim, class GeometryTraits >
    class HybridMappingBase< dim, GeometryTraits, 0 >
    {
      typedef HybridMapping< dim, GeometryTraits > Mapping;

    public:
      virtual ~HybridMappingBase() {}
    protected:
      virtual HybridMapping< dim, GeometryTraits > *
      trace ( Int2Type< 0 >, unsigned int i ) const = 0;
    };



    // HybridMapping
    // -------------

    template< unsigned int dim, class GeometryTraits >
    class HybridMapping
      : public virtual HybridMappingBase< dim, GeometryTraits >,
        public SmallObject
    {
      typedef HybridMapping< dim, GeometryTraits > This;

    protected:
      typedef MappingTraits
      < typename GeometryTraits :: CoordTraits, dim, GeometryTraits :: dimWorld >
      Traits;

    public:
      static const unsigned int dimension = Traits :: dimension;
      static const unsigned int dimWorld = Traits :: dimWorld;

      typedef typename Traits :: FieldType FieldType;
      typedef typename Traits :: LocalCoordType LocalCoordType;
      typedef typename Traits :: GlobalCoordType GlobalCoordType;
      typedef typename Traits :: JacobianType JacobianType;

      template< int codim >
      struct Codim
      {
        typedef HybridMapping< dimension - codim, GeometryTraits > Trace;
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

      template< int codim >
      typename Codim< codim > :: Trace *trace ( unsigned int i ) const
      {
        Int2Type< codim > codimVariable;
        return trace( codimVariable, i );
      }

    protected:
      using HybridMappingBase< dim, GeometryTraits > :: trace;
    };



    // VirtualMappingBase
    // ------------------

    template< class Topology, class GeometryTraits, unsigned int codim = Topology :: dimension >
    class VirtualMappingBase;

    template< class Topology, class GeometryTraits, unsigned int codim >
    class VirtualMappingBase
      : public VirtualMappingBase< Topology, GeometryTraits, codim-1 >,
        public virtual HybridMappingBase< Topology :: dimension, GeometryTraits, codim >
    {
      typedef GenericGeometry :: VirtualMapping< Topology, GeometryTraits >
      VirtualMapping;

    protected:
      using VirtualMappingBase< Topology, GeometryTraits, codim-1 > :: trace;

      virtual HybridMapping< Topology :: dimension - codim, GeometryTraits > *
      trace ( Int2Type< codim >, unsigned int i ) const
      {
        const VirtualMapping &impl = static_cast< const VirtualMapping & >( *this );
        return impl.template trace< codim >( i );
      }
    };

    template< class Topology, class GeometryTraits >
    class VirtualMappingBase< Topology, GeometryTraits, 0 >
      : public virtual HybridMappingBase< Topology :: dimension, GeometryTraits, 0 >
    {
      typedef GenericGeometry :: VirtualMapping< Topology, GeometryTraits >
      VirtualMapping;

    protected:
      virtual HybridMapping< Topology :: dimension, GeometryTraits > *
      trace ( Int2Type< 0 >, unsigned int i ) const
      {
        const VirtualMapping &impl = static_cast< const VirtualMapping & >( *this );
        return impl.template trace< 0 >( i );
      }
    };



    template< class Topology, class GeometryTraits >
    class VirtualMapping
      : public HybridMapping< Topology :: dimension, GeometryTraits >,
        public VirtualMappingBase< Topology, GeometryTraits >
    {
      typedef HybridMapping< Topology :: dimension, GeometryTraits > Base;
      typedef VirtualMapping< Topology, GeometryTraits > This;

      typedef typename Base :: Traits Traits;

      typedef CachedMapping< Topology, GeometryTraits > Mapping;

    public:
      static const unsigned int dimension = Traits :: dimension;
      static const unsigned int dimWorld = Traits :: dimWorld;

      typedef typename Traits :: FieldType FieldType;
      typedef typename Traits :: LocalCoordType LocalCoordType;
      typedef typename Traits :: GlobalCoordType GlobalCoordType;
      typedef typename Traits :: JacobianType JacobianType;

      typedef typename Mapping :: ReferenceElement ReferenceElement;

      template< int codim >
      struct Codim
      {
        typedef HybridMapping< dimension - codim, GeometryTraits > Trace;
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

      virtual GlobalCoordType normal ( int face, const LocalCoordType &local ) const
      {
        return mapping_.normal( face , local );
      }

      template< int codim >
      typename Codim< codim > :: Trace *trace ( unsigned int i ) const
      {
        return mapping_.template trace< codim, true >( i );
      }

    protected:
      using VirtualMappingBase< Topology, GeometryTraits > :: trace;
    };

  }

}

#endif
