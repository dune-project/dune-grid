// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_CACHED_MAPPING_HH
#define DUNE_GENERICGEOMETRY_CACHED_MAPPING_HH

#include <dune/grid/genericgeometry/topologytypes.hh>
#include <dune/grid/genericgeometry/referenceelements.hh>
#include <dune/grid/genericgeometry/matrixhelper.hh>
#include <dune/grid/genericgeometry/mapping.hh>
#include <dune/grid/genericgeometry/traceprovider.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // CachedMapping
    // -------------

    /** \class   CachedMapping
     *  \ingroup GenericGeometry
     *  \brief   caching implementation of a geometric mapping
     *
     *  This is the first user-visible class of the generic geometry
     *  implementation and the last class that explicitly depends on the
     *  topology.
     *
     *  All functions required for a mapping (that should be used as a Geometry)
     *  are implemented.
     *  Moreover, a caching mechanism is added to speed up affine mappings.
     */
    template< class Topology, class GeometryTraits >
    class CachedMapping
    {
      typedef CachedMapping< Topology, GeometryTraits > This;

      typedef typename GeometryTraits::template Mapping< Topology >::type
      MappingImpl;

    public:
      typedef MappingTraits
      < typename GeometryTraits::CoordTraits, Topology::dimension, GeometryTraits::dimWorld >
      Traits;

      typedef GenericGeometry::Mapping
      < typename GeometryTraits::CoordTraits, Topology, GeometryTraits::dimWorld, MappingImpl >
      Mapping;

      static const unsigned int dimension = Traits::dimension;
      static const unsigned int dimWorld = Traits::dimWorld;

      typedef typename Traits::FieldType FieldType;
      typedef typename Traits::LocalCoordinate LocalCoordinate;
      typedef typename Traits::GlobalCoordinate GlobalCoordinate;
      typedef typename Traits::JacobianType JacobianType;
      typedef typename Traits::JacobianTransposedType JacobianTransposedType;

      typedef GenericGeometry::ReferenceElement< Topology, FieldType > ReferenceElement;

      //! can we safely assume that this mapping is always affine?
      static const bool alwaysAffine = Mapping::alwaysAffine;

      template< unsigned int codim >
      struct Codim
      {
        typedef typename TraceProvider< Topology, GeometryTraits, codim, false > :: Trace
        Trace;
      };

      typedef typename GeometryTraits::Caching Caching;

    private:
      typedef typename Traits::MatrixHelper MatrixHelper;

    public:
      template< class CoordVector >
      explicit CachedMapping ( const CoordVector &coords )
        : mapping_( coords ),
          jacobianTransposedComputed_( false ),
          jacobianInverseTransposedComputed_( false ),
          integrationElementComputed_( false )
      {
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

      /** \brief obtain topology id of the corresponding reference element */
      unsigned int topologyId () const
      {
        return ReferenceElement::topologyId;
      }

      // do we still require this mehtod?
      const GlobalCoordinate &corner ( int i ) const
      {
        return mapping_.corner( i );
      }

      /** \brief obtain number of corners of the corresponding reference element */
      int numCorners () const
      {
        return ReferenceElement::numCorners;
      }

      /** \brief obtain the centroid of the mapping's image
       *
       *  \note Currently, this method is defined to return the image
       *        of the reference element's barycenter.
       */
      GlobalCoordinate center () const
      {
        return global( ReferenceElement::template baryCenter< 0 >( 0 ) );
      }

      /** \brief check whether a point lies within the reference element
       *
       *  \param[in]  x  local coorinate of point to check
       *
       *  \note Historically, this method was part of the geometry interface.
       *        It is still required for the GenericReferenceElement.
       */
      static bool checkInside ( const LocalCoordinate &x )
      {
        return ReferenceElement::checkInside( x );
      }

      /** \brief is this mapping affine? */
      bool affine () const
      {
        return (alwaysAffine || affine_);
      }

      /** \brief evaluate the mapping
       *
       *  \param[in]  x  local coordinate to map
       *
       *  \returns corresponding global coordinate
       */
      GlobalCoordinate global ( const LocalCoordinate &x ) const
      {
        GlobalCoordinate y;
        if( jacobianTransposedComputed_ )
        {
          MatrixHelper::template ATx< dimension, dimWorld >( jacobianTransposed_, x, y );
          y += corner( 0 );
        }
        else
          mapping_.global( x, y );
        return y;
      }

      /** \brief evaluate the inverse mapping
       *
       *  \param[in]  y  global coorindate to map
       *
       *  \return corresponding local coordinate
       *
       *  \note The returned local coordinate y minimizes
       *  \code
       *  (global( x ) - y).two_norm()
       *  \endcode
       */
      LocalCoordinate local ( const GlobalCoordinate &y ) const
      {
        LocalCoordinate x;
        if( jacobianInverseTransposedComputed_ )
        {
          GlobalCoordinate z = y - corner( 0 );
          MatrixHelper :: template ATx< dimWorld, dimension >( jacobianInverseTransposed_, z, x );
        }
        else if( affine() )
        {
          const JacobianTransposedType &JT = jacobianTransposed( baryCenter() );
          GlobalCoordinate z = y - corner( 0 );
          MatrixHelper :: template xTRightInvA< dimension, dimWorld >( JT, z, x );
        }
        else
          mapping_.local( y, x );
        return x;
      }

      /** \brief obtain the transposed of the Jacobian
       *
       *  \param[in]  x  local coordinate to evaluate Jacobian in
       *
       *  \returns a reference to the transposed of the Jacobian
       *
       *  \note The returned reference is reused on the next call to
       *        JacobianTransposed, destroying the previous value.
       */
      const JacobianTransposedType &jacobianTransposed ( const LocalCoordinate &x ) const
      {
        const EvaluationType evaluate = Caching::evaluateJacobianTransposed;
        if( (evaluate == PreCompute) && alwaysAffine )
          return jacobianTransposed_;

        if( !jacobianTransposedComputed_ )
          computeJacobianTransposed( x );
        return jacobianTransposed_;
      }

      /** \brief obtain the integration element
       *
       *  If the Jacobian of the mapping is denoted by $J(x)$, the integration
       *  integration element \f$\mu(x)\f$ is given by
       *  \f[ \mu(x) = \sqrt{|\det (J^T(x) J(x))|}.\f]
       *
       *  \param[in]  x  local coordinate to evaluate the integration element in
       *
       *  \returns the integration element \f$\mu(x)\f$.
       *
       *  \note For affine mappings, it is more efficient to call
       *        jacobianInverseTransposed before integrationElement, if both
       *        are required.
       */
      FieldType integrationElement ( const LocalCoordinate &x ) const
      {
        const EvaluationType evaluateI = Caching::evaluateIntegrationElement;
        const EvaluationType evaluateJ = Caching::evaluateJacobianInverseTransposed;
        if( ((evaluateI == PreCompute) || (evaluateJ == PreCompute)) && alwaysAffine )
          return integrationElement_;

        if( !integrationElementComputed_ )
          computeIntegrationElement( x );
        return integrationElement_;
      }

      /** \brief obtain the transposed of the Jacobian's inverse
       *
       *  The Jacobian's inverse is defined as a pseudo-inverse. If we denote
       *  the Jacobian by \f$J(x)\f$, the following condition holds:
       *  \f[J^{-1}(x) J(x) = I.\f]
       */
      const JacobianType &jacobianInverseTransposed ( const LocalCoordinate &x ) const
      {
        const EvaluationType evaluate = Caching::evaluateJacobianInverseTransposed;
        if( (evaluate == PreCompute) && alwaysAffine )
          return jacobianInverseTransposed_;

        if( !jacobianInverseTransposedComputed_ )
          computeJacobianInverseTransposed( x );
        return jacobianInverseTransposed_;
      }

      /** \brief obtain the volume of the mapping's image
       *
       *  \note The current implementation just returns
       *  \code
       *  integrationElement( baryCenter() ) * ReferenceElement::volume()
       *  \endcode
       *  which is wrong for n-linear surface maps and other nonlinear maps.
       */
      FieldType volume () const
      {
        // do we need a quadrature of higher order, here?
        const FieldType refVolume = ReferenceElement::volume();
        return refVolume * integrationElement( baryCenter() );
      }

      template< unsigned int codim, bool hybrid >
      typename TraceProvider< Topology, GeometryTraits, codim, hybrid >::Trace *
      trace ( unsigned int i, typename GeometryTraits::Allocator &allocator ) const
      {
        typedef TraceProvider< Topology, GeometryTraits, codim, hybrid > Provider;
        return Provider::trace( mapping_, i, allocator );
      }

      /** \brief obtain a trace of this mapping to some subentity
       *
       *  \tparam  codim   codimension of the subentity
       *
       *  \param[in]  i  number of the subentity
       *
       *  \returns a pointer to a mapping representing the trace
       *
       *  \note The caller is responsible for deleting the returned mapping.
       */
      template< unsigned int codim >
      typename Codim< codim >::Trace *trace ( unsigned int i ) const
      {
        return trace< codim, false >( i );
      }

    private:
      static const LocalCoordinate &baryCenter ()
      {
        return ReferenceElement::template baryCenter< 0 >( 0 );
      }

      void computeJacobianTransposed ( const LocalCoordinate &x ) const
      {
        affine_ = mapping_.jacobianTransposed( x, jacobianTransposed_ );
        jacobianTransposedComputed_ = affine_;
      }

      void computeJacobianInverseTransposed ( const LocalCoordinate &x ) const
      {
        // jacobian is computed here instead of using the cached value or returning the jacobian
        integrationElement_ = mapping_.jacobianInverseTransposed( x, jacobianInverseTransposed_ );
        integrationElementComputed_ = jacobianInverseTransposedComputed_ = affine();
      }

      void computeIntegrationElement ( const LocalCoordinate &x ) const
      {
        integrationElement_ = mapping_.integrationElement( x );
        integrationElementComputed_ = affine();
      }

    public:
      unsigned int referenceCount;

    private:
      Mapping mapping_;

      mutable JacobianTransposedType jacobianTransposed_;
      mutable JacobianType jacobianInverseTransposed_;
      mutable FieldType integrationElement_;

      mutable bool affine_;

      mutable bool jacobianTransposedComputed_;
      mutable bool jacobianInverseTransposedComputed_;
      mutable bool integrationElementComputed_;
    };

  }

}

#endif // #ifndef DUNE_GENERICGEOMETRY_CACHED_MAPPING_HH
