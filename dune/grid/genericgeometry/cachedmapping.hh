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

    // Internal Forward Declarations
    // -----------------------------

    template< unsigned int, class >
    class CachedJacobianTransposed;

    template< unsigned int, class >
    class CachedJacobianInverseTransposed;



    // CachedStorage
    // -------------

    template< unsigned int dim, class GeometryTraits >
    class CachedStorage
    {
      friend class CachedJacobianTransposed< dim, GeometryTraits >;

    public:
      static const unsigned int dimension = dim;
      static const unsigned int dimWorld = GeometryTraits::dimWorld;

      typedef MappingTraits< typename GeometryTraits::CoordTraits, dimension, dimWorld > Traits;

      typedef typename GeometryTraits::Caching Caching;

      typename Traits::JacobianTransposedType jacobianTransposed;
      typename Traits::JacobianType jacobianInverseTransposed;
      typename Traits::FieldType integrationElement;

      CachedStorage ()
        : affine( false ),
          jacobianTransposedComputed( false ),
          jacobianInverseTransposedComputed( false ),
          integrationElementComputed( false )
      {}

      bool affine;

      bool jacobianTransposedComputed;        // = affine, if jacobian transposed was computed
      bool jacobianInverseTransposedComputed; // = affine, if jacobian inverse transposed was computed
      bool integrationElementComputed;        // = affine, if integration element was computed
    };



    // CachedJacobianTranposed
    // -----------------------

    template< unsigned int dim, class GeometryTraits >
    class CachedJacobianTransposed
    {
      friend class CachedJacobianInverseTransposed< dim, GeometryTraits >;

      typedef CachedStorage< dim, GeometryTraits > Storage;
      typedef typename Storage::Traits Traits;

      typedef typename Traits::MatrixHelper MatrixHelper;

    public:
      typedef typename Traits::FieldType ctype;

      static const int rows = Traits::dimension;
      static const int cols = Traits::dimWorld;

      typedef typename Traits::JacobianTransposedType FieldMatrix;

      operator bool () const
      {
        return storage().jacobianTransposedComputed;
      }

      operator const FieldMatrix & () const
      {
        return storage().jacobianTransposed;
      }

      template< class X, class Y >
      void mv ( const X &x, Y &y ) const
      {
        static_cast< const FieldMatrix & >( *this ).mv( x, y );
      }

      template< class X, class Y >
      void mtv ( const X &x, Y &y ) const
      {
        static_cast< const FieldMatrix & >( *this ).mtv( x, y );
      }

      template< class X, class Y >
      void umv ( const X &x, Y &y ) const
      {
        static_cast< const FieldMatrix & >( *this ).umv( x, y );
      }

      template< class X, class Y >
      void umtv ( const X &x, Y &y ) const
      {
        static_cast< const FieldMatrix & >( *this ).umtv( x, y );
      }

      template< class X, class Y >
      void mmv ( const X &x, Y &y ) const
      {
        static_cast< const FieldMatrix & >( *this ).mmv( x, y );
      }

      template< class X, class Y >
      void mmtv ( const X &x, Y &y ) const
      {
        static_cast< const FieldMatrix & >( *this ).mmtv( x, y );
      }

      ctype det () const
      {
        if( !storage().integrationElementComputed )
        {
          storage().integrationElement = MatrixHelper::template sqrtDetAAT< rows, cols >( storage().jacobianTransposed );
          storage().integrationElementComputed = storage().affine;
        }
        return storage().integrationElement;
      }

    private:
      Storage &storage () const { return storage_; }

      mutable Storage storage_;
    };



    // CachedJacobianInverseTransposed
    // -------------------------------

    template< unsigned int dim, class GeometryTraits >
    class CachedJacobianInverseTransposed
    {
      template< class, class > friend class CachedMapping;

      typedef CachedJacobianTransposed< dim, GeometryTraits > JacobianTransposed;
      typedef typename JacobianTransposed::Storage Storage;
      typedef typename JacobianTransposed::Traits Traits;

      typedef typename Traits::MatrixHelper MatrixHelper;

    public:
      typedef typename Traits::FieldType ctype;

      static const int rows = Traits::dimWorld;
      static const int cols = Traits::dimension;

      typedef typename Traits::JacobianType FieldMatrix;

      operator bool () const
      {
        return storage().jacobianInverseTransposedComputed;
      }

      operator const FieldMatrix & () const
      {
        return storage().jacobianInverseTransposed;
      }

      template< class X, class Y >
      void mv ( const X &x, Y &y ) const
      {
        static_cast< const FieldMatrix & >( *this ).mv( x, y );
      }

      template< class X, class Y >
      void mtv ( const X &x, Y &y ) const
      {
        static_cast< const FieldMatrix & >( *this ).mtv( x, y );
      }

      template< class X, class Y >
      void umv ( const X &x, Y &y ) const
      {
        static_cast< const FieldMatrix & >( *this ).umv( x, y );
      }

      template< class X, class Y >
      void umtv ( const X &x, Y &y ) const
      {
        static_cast< const FieldMatrix & >( *this ).umtv( x, y );
      }

      template< class X, class Y >
      void mmv ( const X &x, Y &y ) const
      {
        static_cast< const FieldMatrix & >( *this ).mmv( x, y );
      }

      template< class X, class Y >
      void mmtv ( const X &x, Y &y ) const
      {
        static_cast< const FieldMatrix & >( *this ).mmtv( x, y );
      }

      ctype det () const
      {
        // integrationElement is always computed with jacobianInverseTransposed
        return ctype( 1 ) / storage().integrationElement;
      }

    private:
      JacobianTransposed &jacobianTransposed () { return jacobianTransposed_; }
      const JacobianTransposed &jacobianTransposed () const { return jacobianTransposed_; }

      Storage &storage () const { return jacobianTransposed().storage(); }

      JacobianTransposed jacobianTransposed_;
    };



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

      typedef CachedStorage< dimension, GeometryTraits > Storage;
      typedef CachedJacobianTransposed< dimension, GeometryTraits > JacobianTransposed;
      typedef CachedJacobianInverseTransposed< dimension, GeometryTraits > JacobianInverseTransposed;

      typedef GenericGeometry::ReferenceElement< Topology, FieldType > ReferenceElement;

      //! can we safely assume that this mapping is always affine?
      static const bool alwaysAffine = Mapping::alwaysAffine;

      template< unsigned int codim >
      struct Codim
      {
        typedef typename TraceProvider< Topology, GeometryTraits, codim, false >::Trace Trace;
      };

      typedef typename GeometryTraits::Caching Caching;

    private:
      typedef typename Traits::MatrixHelper MatrixHelper;

    public:
      template< class CoordVector >
      explicit CachedMapping ( const CoordVector &coords )
        : mapping_( coords )
      {
        if( alwaysAffine )
          storage().affine = true;
        else
          computeJacobianTransposed( baryCenter() );
        preCompute();
      }

      template< class CoordVector >
      explicit CachedMapping ( const std::pair< const CoordVector &, bool > &coords )
        : mapping_( coords.first )
      {
        storage().affine = coords.second;
        preCompute();
      }

      /** \brief obtain topology id of the corresponding reference element */
      unsigned int topologyId () const
      {
        return ReferenceElement::topologyId;
      }

      /** \brief obtain coordinates of the i-th corner */
      GlobalCoordinate corner ( int i ) const
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
        return (alwaysAffine || storage().affine);
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
        if( jacobianTransposed() )
        {
          y = corner( 0 );
          jacobianTransposed().umtv( x, y );
          //MatrixHelper::template ATx< dimension, dimWorld >( jacobianTransposed_, x, y );
          //y += corner( 0 );
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
        if( jacobianInverseTransposed() )
        {
          GlobalCoordinate z = y - corner( 0 );
          jacobianInverseTransposed().mtv( z, x );
          // MatrixHelper::template ATx< dimWorld, dimension >( jacobianInverseTransposed(), z, x );
        }
        else if( affine() )
        {
          const JacobianTransposed &JT = jacobianTransposed( baryCenter() );
          GlobalCoordinate z = y - corner( 0 );
          MatrixHelper::template xTRightInvA< dimension, dimWorld >( JT, z, x );
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
      const JacobianTransposed &jacobianTransposed ( const LocalCoordinate &x ) const
      {
        const EvaluationType evaluate = Caching::evaluateJacobianTransposed;
        if( (evaluate == PreCompute) && alwaysAffine )
          return jacobianTransposed();

        if( !jacobianTransposed() )
          computeJacobianTransposed( x );
        return jacobianTransposed();
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
          return storage().integrationElement;
        else
          return jacobianTransposed( x ).det();
      }

      /** \brief obtain the transposed of the Jacobian's inverse
       *
       *  The Jacobian's inverse is defined as a pseudo-inverse. If we denote
       *  the Jacobian by \f$J(x)\f$, the following condition holds:
       *  \f[J^{-1}(x) J(x) = I.\f]
       */
      const JacobianInverseTransposed &
      jacobianInverseTransposed ( const LocalCoordinate &x ) const
      {
        const EvaluationType evaluate = Caching::evaluateJacobianInverseTransposed;
        if( (evaluate == PreCompute) && alwaysAffine )
          return jacobianInverseTransposed();

        if( !jacobianInverseTransposed() )
          computeJacobianInverseTransposed( x );
        return jacobianInverseTransposed();
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

      This *clone () const
      {
        return new This( *this );
      }

      This* clone ( char *mappingStorage ) const
      {
        return new( mappingStorage ) This( *this );
      }

      template< unsigned int codim, bool hybrid >
      typename TraceProvider< Topology, GeometryTraits, codim, hybrid >::Trace*
      trace ( unsigned int i, char *mappingStorage ) const
      {
        return TraceProvider< Topology, GeometryTraits, codim, hybrid >::construct( mapping_, i, mappingStorage );
      }

    private:
      static const LocalCoordinate &baryCenter ()
      {
        return ReferenceElement::template baryCenter< 0 >( 0 );
      }

      Storage &storage () const
      {
        return jacobianInverseTransposed().storage();
      }

      const JacobianTransposed &jacobianTransposed () const
      {
        return jacobianInverseTransposed().jacobianTransposed();
      }

      const JacobianInverseTransposed &jacobianInverseTransposed () const
      {
        return jacobianInverseTransposed_;
      }

      void preCompute ()
      {
        assert( affine() == mapping_.jacobianTransposed( baryCenter(), storage().jacobianTransposed ) );
        if( !affine() )
          return;

        if( (Caching::evaluateJacobianTransposed == PreCompute) && !jacobianTransposed() )
          computeJacobianTransposed( baryCenter() );

        if( Caching::evaluateJacobianInverseTransposed == PreCompute )
          computeJacobianInverseTransposed( baryCenter() );
        else if( Caching::evaluateIntegrationElement == PreCompute )
          jacobianTransposed().det();
      }

      void computeJacobianTransposed ( const LocalCoordinate &x ) const
      {
        storage().affine = mapping_.jacobianTransposed( x, storage().jacobianTransposed );
        storage().jacobianTransposedComputed = affine();
      }

      void computeJacobianInverseTransposed ( const LocalCoordinate &x ) const
      {
        storage().integrationElement
          = MatrixHelper::template rightInvA< dimension, dimWorld >( jacobianTransposed( x ), storage().jacobianInverseTransposed );
        storage().integrationElementComputed = affine();
        storage().jacobianInverseTransposedComputed = affine();
      }

    private:
      Mapping mapping_;
      JacobianInverseTransposed jacobianInverseTransposed_;
    };

  } // namespace GenericGeometry

} // namespace Dune

#endif // #ifndef DUNE_GENERICGEOMETRY_CACHED_MAPPING_HH
