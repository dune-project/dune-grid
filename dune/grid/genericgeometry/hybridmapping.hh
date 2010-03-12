// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMTRY_HYBRIDMAPPING_HH
#define DUNE_GENERICGEOMTRY_HYBRIDMAPPING_HH

#include <dune/common/smallobject.hh>

#include <dune/grid/genericgeometry/geometrytraits.hh>
#include <dune/grid/genericgeometry/cachedmapping.hh>

namespace Dune
{

  namespace GenericGeometry
  {


    // Internal Forward Declarations
    // -----------------------------

    template< unsigned int dim, class GeometryTraits >
    class HybridMapping;

    template< class Topology, class GeometryTraits >
    class VirtualMapping;



    // HybridMappingBase
    // -----------------

    /** \cond */
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
    /** \endcond */



    // HybridMapping
    // -------------

    /** \class   HybridMapping
     *  \ingroup GenericGeometry
     *  \brief   abstract base class for generic mapping
     *
     *  This is the user-visible class of the generic geometries if the
     *  topology type for each codimension is not unique. It is the abstract
     *  base class of VirtualMapping, which implements all methods by
     *  forwarding them to a CachedMapping for the corresponding topology.
     */
    template< unsigned int dim, class GeometryTraits >
    class HybridMapping
    /** \cond */
      : public virtual HybridMappingBase< dim, GeometryTraits >,
        public SmallObject
        /** \endcond */
    {
      typedef HybridMapping< dim, GeometryTraits > This;

    protected:
      typedef MappingTraits
      < typename GeometryTraits::CoordTraits, dim, GeometryTraits::dimWorld >
      Traits;

    public:
      static const unsigned int dimension = Traits::dimension;
      static const unsigned int dimWorld = Traits::dimWorld;

      typedef typename Traits::FieldType FieldType;
      typedef typename Traits::LocalCoordType LocalCoordType;
      typedef typename Traits::GlobalCoordType GlobalCoordType;
      typedef typename Traits::JacobianType JacobianType;
      typedef typename Traits::JacobianTransposedType JacobianTransposedType;

      template< int codim >
      struct Codim
      {
        typedef HybridMapping< dimension - codim, GeometryTraits > Trace;
      };

      typedef typename GeometryTraits::Caching Caching;

      unsigned int referenceCount;

      virtual ~HybridMapping ()
      {}

      /** \copydoc CachedMapping::topologyId */
      virtual unsigned int topologyId () const = 0;

      /** \copydoc CachedMapping::corner */
      virtual const GlobalCoordType &corner ( int i ) const = 0;

      /** \copydoc CachedMapping::numCorners */
      virtual int numCorners () const = 0;

      /** \copydoc CachedMapping::center */
      virtual GlobalCoordType center () const = 0;

      /** \copydoc CachedMapping::global */
      virtual GlobalCoordType global ( const LocalCoordType &x ) const = 0;

      /** \copydoc CachedMapping::local */
      virtual LocalCoordType local ( const GlobalCoordType &y ) const = 0;

      /** \copydoc CachedMapping::checkInside */
      virtual bool checkInside ( const LocalCoordType &x ) const = 0;

      /** \copydoc CachedMapping::affine */
      virtual bool affine () const = 0;

      /** \copydoc CachedMapping::integrationElement */
      virtual FieldType integrationElement ( const LocalCoordType &x ) const = 0;

      /** \copydoc CachedMapping::volume */
      virtual FieldType volume () const = 0;

      /** \copydoc CachedMapping::jacobianTransposed */
      virtual const JacobianTransposedType &
      jacobianTransposed ( const LocalCoordType &x ) const = 0;

      /** \copydoc CachedMapping::jacobianInverseTransposed */
      virtual const JacobianType &
      jacobianInverseTransposed ( const LocalCoordType &x ) const = 0;

      /** \copydoc CachedMapping::normal */
      virtual GlobalCoordType
      normal ( int face, const LocalCoordType &x ) const = 0;

    protected:
      using HybridMappingBase< dim, GeometryTraits >::trace;

    public:
      /** \copydoc CachedMapping::trace */
      template< int codim >
      typename Codim< codim > :: Trace *trace ( unsigned int i ) const
      {
        Int2Type< codim > codimVariable;
        return trace( codimVariable, i );
      }
    };



    // VirtualMappingBase
    // ------------------

    /** \cond */
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
    /** \endcond */



    template< class Topology, class GeometryTraits >
    class VirtualMapping
      : public HybridMapping< Topology :: dimension, GeometryTraits >,
        /** \cond */
        public VirtualMappingBase< Topology, GeometryTraits >
        /**\endcond*/
    {
      typedef HybridMapping< Topology :: dimension, GeometryTraits > Base;
      typedef VirtualMapping< Topology, GeometryTraits > This;

      typedef typename Base :: Traits Traits;

      typedef CachedMapping< Topology, GeometryTraits > Mapping;

    public:
      static const unsigned int dimension = Traits :: dimension;
      static const unsigned int dimWorld = Traits :: dimWorld;

      typedef typename Traits::FieldType FieldType;
      typedef typename Traits::LocalCoordType LocalCoordType;
      typedef typename Traits::GlobalCoordType GlobalCoordType;
      typedef typename Traits::JacobianType JacobianType;
      typedef typename Traits::JacobianTransposedType JacobianTransposedType;

      typedef typename Mapping :: ReferenceElement ReferenceElement;

      template< int codim >
      struct Codim
      {
        typedef HybridMapping< dimension - codim, GeometryTraits > Trace;
      };

      typedef typename GeometryTraits::Caching Caching;

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

      virtual GlobalCoordType center () const
      {
        return mapping_.center();
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

      virtual const JacobianTransposedType &
      jacobianTransposed ( const LocalCoordType &local ) const
      {
        return mapping_.jacobianTransposed( local );
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
