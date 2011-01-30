// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GENERICGEOMETRY_HYBRIDMAPPING_HH
#define DUNE_GENERICGEOMETRY_HYBRIDMAPPING_HH

#include <dune/common/typetraits.hh>

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

      virtual HybridMapping< dim - codim, GeometryTraits >*
      trace ( integral_constant< int, codim >, unsigned int i, char *mappingStorage ) const = 0;
    };

    template< unsigned int dim, class GeometryTraits >
    class HybridMappingBase< dim, GeometryTraits, 0 >
    {
      typedef HybridMapping< dim, GeometryTraits > Mapping;

    public:
      virtual ~HybridMappingBase() {}

    protected:
      virtual HybridMapping< dim, GeometryTraits >*
      trace ( integral_constant< int, 0 >, unsigned int i, char *mappingStorage ) const = 0;
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
      : public virtual HybridMappingBase< dim, GeometryTraits >
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
      typedef typename Traits::LocalCoordinate LocalCoordinate;
      typedef typename Traits::GlobalCoordinate GlobalCoordinate;

      typedef CachedJacobianTransposed< dimension, GeometryTraits > JacobianTransposed;
      typedef CachedJacobianInverseTransposed< dimension, GeometryTraits > JacobianInverseTransposed;

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
      virtual const GlobalCoordinate &corner ( int i ) const = 0;

      /** \copydoc CachedMapping::numCorners */
      virtual int numCorners () const = 0;

      /** \copydoc CachedMapping::center */
      virtual GlobalCoordinate center () const = 0;

      /** \copydoc CachedMapping::global */
      virtual GlobalCoordinate global ( const LocalCoordinate &x ) const = 0;

      /** \copydoc CachedMapping::local */
      virtual LocalCoordinate local ( const GlobalCoordinate &y ) const = 0;

      /** \copydoc CachedMapping::checkInside */
      virtual bool checkInside ( const LocalCoordinate &x ) const = 0;

      /** \copydoc CachedMapping::affine */
      virtual bool affine () const = 0;

      /** \copydoc CachedMapping::integrationElement */
      virtual FieldType integrationElement ( const LocalCoordinate &x ) const = 0;

      /** \copydoc CachedMapping::volume */
      virtual FieldType volume () const = 0;

      /** \copydoc CachedMapping::jacobianTransposed */
      virtual const JacobianTransposed &
      jacobianTransposed ( const LocalCoordinate &x ) const = 0;

      /** \copydoc CachedMapping::jacobianInverseTransposed */
      virtual const JacobianInverseTransposed &
      jacobianInverseTransposed ( const LocalCoordinate &x ) const = 0;

    protected:
      using HybridMappingBase< dim, GeometryTraits >::trace;

    public:
      virtual This *clone () const = 0;
      virtual This *clone ( char *mappingStorage ) const = 0;

      template< int codim >
      typename Codim< codim >::Trace* trace ( unsigned int i, char *mappingStorage ) const
      {
        integral_constant< int, codim > codimVariable;
        return trace( codimVariable, i, mappingStorage );
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

      virtual HybridMapping< Topology::dimension - codim, GeometryTraits >*
      trace ( integral_constant< int, codim >, unsigned int i, char *mappingStorage) const
      {
        return static_cast< const VirtualMapping & >( *this ).template trace< codim >( i, mappingStorage );
      }
    };

    template< class Topology, class GeometryTraits >
    class VirtualMappingBase< Topology, GeometryTraits, 0 >
      : public virtual HybridMappingBase< Topology :: dimension, GeometryTraits, 0 >
    {
      typedef GenericGeometry :: VirtualMapping< Topology, GeometryTraits >
      VirtualMapping;

    protected:
      virtual HybridMapping< Topology::dimension, GeometryTraits >*
      trace ( integral_constant< int, 0 >, unsigned int i,
              char *mappingStorage ) const
      {
        return static_cast< const VirtualMapping & >( *this ).template trace< 0 >( i, mappingStorage );
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
      typedef typename Traits::LocalCoordinate LocalCoordinate;
      typedef typename Traits::GlobalCoordinate GlobalCoordinate;

      typedef typename Base::JacobianTransposed JacobianTransposed;
      typedef typename Base::JacobianInverseTransposed JacobianInverseTransposed;

      typedef typename Mapping::ReferenceElement ReferenceElement;

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

      virtual const GlobalCoordinate &corner ( int i ) const
      {
        return mapping_.corner( i );
      }

      virtual int numCorners () const
      {
        return mapping_.numCorners();
      }

      virtual GlobalCoordinate center () const
      {
        return mapping_.center();
      }

      virtual GlobalCoordinate global ( const LocalCoordinate &local ) const
      {
        return mapping_.global( local );
      }

      virtual LocalCoordinate local ( const GlobalCoordinate &global ) const
      {
        return mapping_.local( global );
      }

      virtual bool checkInside ( const LocalCoordinate &local ) const
      {
        return mapping_.checkInside( local );
      }

      virtual bool affine () const
      {
        return mapping_.affine();
      }

      virtual FieldType integrationElement ( const LocalCoordinate &local ) const
      {
        return mapping_.integrationElement( local );
      }

      virtual FieldType volume () const
      {
        return mapping_.volume();
      }

      virtual const JacobianTransposed &
      jacobianTransposed ( const LocalCoordinate &local ) const
      {
        return mapping_.jacobianTransposed( local );
      }

      virtual const JacobianInverseTransposed &
      jacobianInverseTransposed ( const LocalCoordinate &local ) const
      {
        return mapping_.jacobianInverseTransposed( local );
      }

      virtual Base *clone () const
      {
        return new This( *this );
      }

      virtual Base* clone ( char *mappingStorage ) const
      {
        return new( mappingStorage ) This( *this );
      }

      template< int codim >
      typename Codim< codim >::Trace* trace ( unsigned int i, char *mappingStorage ) const
      {
        return mapping_.template trace< codim, true >( i, mappingStorage );
      }

    protected:
      using VirtualMappingBase< Topology, GeometryTraits > :: trace;
    };

  }

}

#endif // #ifndef DUNE_GENERICGEOMETRY_HYBRIDMAPPING_HH
