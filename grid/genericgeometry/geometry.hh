// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_GEOMETRY_HH
#define DUNE_GENERICGEOMETRY_GEOMETRY_HH

#include <dune/grid/common/geometry.hh>

#include <dune/grid/genericgeometry/mappingprovider.hh>
#include <dune/grid/genericgeometry/geometrytraits.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    /** \addtogroup GenericGeometry
     *
     *  \section General
     *
     *  Based on a recursive definition of the reference elements, a generic
     *  implementation of Dune::Geometry is provided. The class used for the
     *  implementation of the Dune::Geometry engine, is
     *  GenericGeometry::BasicGeometry.
     *
     *  The BasicGeometry class takes a template argument Traits specifying
     *  details of the reference mapping implementation and some performance
     *  settings. A default implementation for this class is
     *  GenericGeometry::DefaultGeometryTraits. The traits class must contain
     *  the same types as this default implementation.
     *
     *  To conform with the Dune::Geometry engine, two further classes are
     *  provided: GenericGeometry::Geometry and GenericGeometry::LocalGeometry.
     *  To use these classes instead of GenericGeometry::BasicGeometry, the
     *  traits classes
     *  \code
     *  template< class Grid> GenericGeometry::GlobalGeometryTraits<Grid>
     *  template< class Grid> GenericGeometry::LocalGeometryTraits<Grid>
     *  \endcode
     *  have to be specialized. These classes are simply passed as Traits
     *  argument to GenericGeometry::BasicGeometry.
     *
     *  The reference mapping for a given topology type is given by
     *  Mapping<Topology>::type in the traits class. Here, Topology is one of
     *  the generic topology classes GenericGeometry::Point,
     *  GenericGeometry::Prism, GenericGeometry::Pyramid.
     *  An interface for the mapping is provided by GenericGeometry::Mapping.
     *  The implementation of this interface must have constructors taking a
     *  single argument. The constructor of GenericGeometry::BasicGeometry
     *  looks as follows:
     *  \code
     *  template< class CoordVector >
     *  BasicGeometry ( const GeometryType &type, const CoordVector &coords );
     *  \endcode
     *  Its first argument, <em>type</em>, specifies the type of the reference
     *  element (as a Dune::GeometryType). The second argument, <em>coords</em>
     *  is passed directly to the constructor of the mapping implementation.
     *  The most prominent implementation of GenericGeometry::Mapping is
     *  GenericGeometry::CornerMapping. It provides a polynomial interpolation
     *  of the entity's corners with minimal degree. In this case,
     *  <em>coords</em> represents the entity's corners.
     *
     *  \section Simple Usage
     *  To add first order Lagrange type geometries to a grid implementation
     *  the following steps suffice:
     *  - Overload the traits classes
     *    \code
            template<>
            struct GenericGeometry::GlobalGeometryTraits< MyGrid >
            : public GenericGeometry::DefaultGeometryTraits
                     <MyGrid::cytpe,MyGrid::dimension,MyGrid::dimworld>
            {};
            template<>
            struct GenericGeometry::LocalGeometryTraits< MyGrid >
            : public GenericGeometry::DefaultGeometryTraits
                     <MyGrid::cytpe,MyGrid::dimension,MyGrid::dimworld>
            {};
     *    \endcode
     *    Note that these classes are default implementation which should cover all
     *    cases but
     *    are in a specific situation far from the optimal choice.
     *    For example, an increase
     *    in efficiency can be achieved for grids with a fixed element type
     *    (set hyprid to false and set the dunetype variable)
     *    or for grids with only affine transformations - which in the case of
     *    the local geometries is often true the variable; the last template
     *    argument (default false) can be used to switch to mappings which are
     *    assumed to always be affined (no checking done).
     *  - Add to the GridFamily::Traits::Codim<codim> structure:
     *    \code
            typedef Dune :: Geometry
              < dimension-codim, dimensionworld, const MyGrid,
                Dune :: GenericGeometry :: Geometry > Geometry;
            typedef Dune :: Geometry
              < dimension-codim, dimension, const MyGrid,
                Dune :: GenericGeometry :: LocalGeometry > LocalGeometry;
     *    \endcode
     *  - Both geometries can be build by calling the constructor taking
     *    a DUNE grid type and an instance of an arbitrary class with a method
            \code
              const FieldVector<  ctype, dimensionworld >& operator[](unsigned int i);
            \endcode
     *    The references returned must remain valid during the whole life span
     *    of the geometry.
     *  - In MyGrid::Entity<0> the following methods can then be easily implemented:
     *    - geometry(): this requires the knowledge of the dune geometry type of the entity
     *      and the coordinates of the corner points.
     *    - geometryInFather(): The corner points for each child in the
     *      reference element of the father can be used to construct the local geometry -
     *      note that this geometry is mostly affine and these geometries can be
     *      precomputed and stored.
     *    .
     *  - For the Dune::Intersection class the geometries the following implementations for the geometries can be used:
     *    - intersectionGlobal(): can be implemented in the same way as the geometry of the
     *      entity using the coordinates of the corners of the intersection.
     *      Alternatively, in the case of a conform intersection,
     *      the class GenericGeometry::Geometry provides a possibility
     *      to construct traces of a given geometry, e.g., a reference mapping
     *      restricted to a codimension one subentities of the reference
     *      element. This is achieved by calling the constructor on the
     *      GenericGeometry::Geometry class (with the codim template equal to
     *      one) passing a codimension zero geometry implementation and the number of the
     *      codimension one subentity (in DUNE numbering).
     *      \code
            GenericGeometry::Geometry<myGridDim-1,myWorldDim,MyGrid>
                             (inside->geometry(),numberInSelf());
     *      \endcode
     *    - intersectionInSelf()/intersectionInNeighbor():
     *      A similar strategy as described above for the intersectionGlobal
     *      can also be used for the geometry mapping to the codimension zero
     *      reference element. Either the corners of the intersection in
     *      local coordinates can be used in the construction of the local
     *      geometries, or (for conform intersections) the traces can be used,
     *      passing an identity mapping as codimension zero geometry.
     *      The GenericGeometry::GenericReferenceElement provides these
     *      mappings directly via the template method
     *      GenericGeometry::GenericReferenceElement::mapping.
     *      The return value of this method can be directly used to construct
     *      a GenericGeometry::Geometry instance:
     *      \code
            typedef GenericReferenceElementContainer<ctype,myGridDim> RefElementContType;
            RefElementContType refElemCont;
            const RefElementContType::value_type& refElem=refElemCont(insideGeometryType);
            GenericGeometry::Geometry<myGridDim-1,myGridDim,MyGrid>(refElem.mapping(numberInSelf()));
     *      \endcode
     *    - integrationOuterNormal(): the generic geometry implementation provides a method
     *      to compute the integration outer normals, so that the following code
     *      fragment can be used:
            \code
               typedef typename Grid :: template Codim< 0 > :: Geometry Geometry;
               const Geometry &geo = inside()->geometry();
               FieldVector< ctype, dimension > x( intersectionSelfLocal().global( local ) );
               return Grid :: getRealImplementation( geo ).normal( numberInSelf(), x );
            \endcode
     *    .
     *  - To add geometries for subentitiies of codim>0
     *    given a entity en of codimension zero and the subentity number subNr
     *    in DUNE numbering:
     *    - geometry: the geometry can be constructed by the following line of code
            \code
            GenericGeometry::Geometry<myGridDim-codim,myWorldDim,MyGrid>
                             (en.geometry(),subNr);
            \endcode
     *    .
     *  .
     *
     */

    // BasicGeometry
    // -------------

    /** \class   BasicGeometry
     *  \ingroup GenericGeometry
     *  \brief   generic implementation of DUNE geometries
     *
     *  This class is provides a generic implementation of a DUNE geometry.
     *
     *  Parameters shared by all codimensions are summarized into one class
     *  parameter called Traits. The following default implementation can be
     *  used (via subclassing) to provide the necessary information. It contains
     *  exactly the required fields:
     *  \code
     *  template< class ctype, int dimG, int dimW >
     *  struct DefaultGeometryTraits
     *  {
     *    typedef DuneCoordTraits< ctype > CoordTraits;
     *
     *    static const int dimGrid = dimG;
     *    static const int dimWorld = dimW;
     *
     *    //   hybrid   [ true if Codim 0 is hybrid ]
     *    static const bool hybrid = true;
     *    //   dunetype [ for Codim 0, needed for (hybrid=false) ]
     *    // static const GeometryType :: BasicType dunetype = GeometryType :: simplex;
     *
     *    // what basic geometry type shall the line be considered?
     *    static const GeometryType :: BasicType linetype = GeometryType :: simplex;
     *
     *    // explained below
     *    template< class Topology >
     *    struct Mapping
     *    {
     *      typedef MappingTraits< CoordTraits, Topology :: dimension, dimWorld > Traits;
     *      typedef CoordPointerStorage< Topology, typename Traits :: GlobalCoordType >
     *        CornerStorage;
     *      typedef CornerMapping< Topology, Traits, CornerStorage > type;
     *    };
     *
     *    // explained below
     *    struct Caching
     *    {
     *      static const EvaluationType evaluateJacobianTransposed = ComputeOnDemand;
     *      static const EvaluationType evaluateJacobianInverseTransposed = ComputeOnDemand;
     *      static const EvaluationType evaluateIntegrationElement = ComputeOnDemand;
     *      static const EvaluationType evaluateNormal = ComputeOnDemand;
     *    };
     *  };
     *  \endcode
     *
     *  The structure specifying the reference mapping is
     *  Traits::Mapping::type. An example implementation
     *  is the GenericGeometry::CornerMapping which defines
     *  the simple mapping taking corners of the reference
     *  elements to corner of the entity in space.
     *  The class given by Traits::Mapping::CornerStorage
     *  (an example is given by GenericGeometry::CoordPointerStorage).
     *  is a container for the coordinates of the corners
     *  returned by the Dune::Geometry (also required for
     *  non-linear reference mappings).
     *  The third type in Traits::Mapping is a traits class
     *  for the reference mapping, following the structure of
     *  GenericGeometry::MappingTraits.
     *
     *  The central reference mapping specified by Traits::Mapping::type
     *  requires
     *  a constructure taking a single argument. The
     *  GenericGeometry::BasicGeometry has a constructure
     *  with one template argument which is passed on to
     *  the constructure of the used provided reference mapping.
     *  The interface for the this class is
     *  GenericGeometry::Mapping.
     *
     *  To increase the efficiency of the geometry
     *  implementation, different strategies for
     *  the caching of parts of the geometry data
     *  is provided. The specifics are given
     *  by the structure Traits::Caching. Possible
     *  values are:
     *  - ComputeOnDemand:    use caching if method called using barycenter
     *  - PreCompute:         use caching in constructor using barycenter
     *  .
     *
     *  \note This class cannot be used directly as an implementation of
     *        Dune::Geometry. Its template parameter list differs from what
     *        is expected there from the engine.
     *        One of the following derived classes
     *        can be used instead:
     *        - Dune::GenericGeometry::Geometry
     *        - Dune::GenericGeometry::LocalGeometry
     *        .
     */
    template< int mydim, class Traits >
    class BasicGeometry
    {
      typedef typename Traits :: CoordTraits CoordTraits;

      static const int dimGrid = Traits :: dimGrid;

      template< int, class > friend class BasicGeometry;

    public:

      /** \brief The dimension of the parameter space of this geometry */
      static const int mydimension = mydim;

      /** \brief The dimension of the world space of this geometry */
      static const int coorddimension = Traits :: dimWorld;

      /** \brief Type used for coordinate components */
      typedef typename CoordTraits :: ctype ctype;

      /** \brief Type used for parameter coordinates */
      typedef FieldVector< ctype, mydimension > LocalCoordinate;

      /** \brief Type used for world coordinates */
      typedef FieldVector< ctype, coorddimension > GlobalCoordinate;

      /** \brief Type used for Jacobian matrices */
      typedef FieldMatrix< ctype, coorddimension, mydimension > Jacobian;
      typedef FieldMatrix< ctype, mydimension, coorddimension > JacobianTransposed;

    private:
      dune_static_assert( (0 <= mydimension) && (mydimension <= dimGrid),
                          "Invalid geometry dimension." );

      static const int codimension = dimGrid - mydimension;

      template< bool >
      struct Hybrid
      {
        typedef HybridMapping< dimGrid, Traits > Mapping;
      };

      template< bool >
      struct NonHybrid
      {
        typedef typename Convert< Traits :: dunetype, dimGrid > :: type Topology;
        typedef GenericGeometry :: CachedMapping< Topology, Traits > Mapping;
      };

      typedef GenericGeometry :: DuneGeometryTypeProvider< mydimension, Traits :: linetype >
      DuneGeometryTypeProvider;

      typedef typename ProtectedIf< Traits :: hybrid, Hybrid, NonHybrid > :: Mapping
      ElementMapping;
      typedef GenericGeometry :: MappingProvider< ElementMapping, codimension >
      MappingProvider;

    protected:
      typedef typename MappingProvider :: Mapping Mapping;

    private:
      Mapping *mapping_;

    public:

      /** \brief Default constructor */
      BasicGeometry ()
        : mapping_( 0 )
      {}

#if 0
      /** \brief Constructor taking a mapping */
      explicit BasicGeometry ( Mapping &mapping )
        : mapping_( &mapping )
      {
        ++mapping_->referenceCount;
      }
#endif

      /** \brief Constructor using a GeometryType and a list of corner coordinates */
      template< class CoordVector >
      BasicGeometry ( const GeometryType &type, const CoordVector &coords )
        : mapping_( MappingProvider :: mapping( type, coords ) )
      {
        mapping_->referenceCount = 1;
      }

      /** \brief obtain a geometry for a subentity
       *
       *  Assume that we have a geometry for some entity d-dimensional E.
       *  This method can provide a geometry for the i-th subentity of E
       *  (of codimension d - mydimension).
       *
       *  \note This method can be more efficient than just building up the
       *        geometry for the subentity. For example, the subgeometry
       *        automatically inherits affinity.
       *
       *  \param[in]  father  geometry of entity \em E
       *  \param[in]  i       number of the subentity (in generic numbering)
       */
      template< int fatherdim >
      BasicGeometry ( const BasicGeometry< fatherdim, Traits > &father, int i )
        : mapping_( subMapping( father, i ) )
      {
        mapping_->referenceCount = 1;
      }

      /** \brief Copy constructor */
      BasicGeometry ( const BasicGeometry &other )
        : mapping_( other.mapping_ )
      {
        if( mapping_ != 0 )
          ++(mapping_->referenceCount);
      }

      /** \brief Destructor */
      ~BasicGeometry ()
      {
        if( (mapping_ != 0) && ((--mapping_->referenceCount) == 0) )
          delete mapping_;
      }

      /** \brief Assignment from other BasicGeometry */
      BasicGeometry &operator= ( const BasicGeometry &other )
      {
        if( other.mapping_ != 0 )
          ++(other.mapping_->referenceCount);
        if( (mapping_ != 0) && (--(mapping_->referenceCount) == 0) )
          delete mapping_;
        mapping_ = other.mapping_;
        return *this;
      }

      /** \brief Test whether this BasicGeometry is properly set up
          \todo Please doc me better!
       */
      bool operator! () const
      {
        return (mapping_ == 0);
      }

      /** \brief Return the topological type of this geometry */
      GeometryType type () const
      {
        return DuneGeometryTypeProvider :: type( mapping().topologyId() );
      }

      /** \brief Return the number of corners */
      int corners () const
      {
        return mapping().numCorners();
      }

      /** \brief Return the world coordinates of the i-th corner */
      const GlobalCoordinate &operator[] ( int i ) const
      {
        return mapping().corner( i );
      }

      /** \brief Return the world coordinates of the i-th corner */
      GlobalCoordinate corner ( const int i ) const
      {
        return mapping().corner( i );
      }

      /** \brief Map local to global coordinates */
      GlobalCoordinate global ( const LocalCoordinate &local ) const
      {
        return mapping().global( local );
      }

      /** \brief Map global to local coordinates */
      LocalCoordinate local ( const GlobalCoordinate &global ) const
      {
        return mapping().local( global );
      }

#if 0
      /** \brief Return true if a given point is within the parameter domain */
      bool checkInside ( const LocalCoordinate &local ) const
      {
        return mapping().checkInside( local );
      }
#endif

      /** \brief Return true if this is an affine geometry */
      bool affine () const
      {
        return mapping().affine();
      }

      /** \brief Return the factor \$|det F|\$ that appears in the integral transformation formula */
      ctype integrationElement ( const LocalCoordinate &local ) const
      {
        return mapping().integrationElement( local );
      }

      /** \brief Return the volume of the element */
      ctype volume () const
      {
        return mapping().volume();
      }

      const JacobianTransposed &jacobianTransposed ( const LocalCoordinate &local ) const
      {
        return mapping().jacobianTransposed( local );
      }

      /** \brief Compute the transpose of the inverse Jacobian matrix of the transformation
          from the reference element into the world space */
      const Jacobian &jacobianInverseTransposed ( const LocalCoordinate &local ) const
      {
        return mapping().jacobianInverseTransposed( local );
      }

      /** \brief compute an outer normal
       *
       *  \param[in]  face   number of the face (in generic numbering)
       *  \param[in]  local  point to compute the normal in (in local coordinates)
       *
       *  \returns an outer normal to the given face at the given point
       *
       *  \note Thouogh the local coordinates are given with respect to geometry's
       *        reference domain, the point is required to be on the given face.
       */
      GlobalCoordinate normal ( int face, const LocalCoordinate &local ) const
      {
        return mapping().normal( face, local );
      }

    private:
      const Mapping &mapping () const
      {
        assert( mapping_ != 0 );
        return *mapping_;
      }

      template< int fatherdim >
      Mapping *
      subMapping ( const BasicGeometry< fatherdim, Traits > &father, int i )
      {
        const unsigned int codim = fatherdim - mydim;
        return father.mapping().template trace< codim >( i );
      }
    };



    // Geometry
    // --------

    /** \class   Geometry
     *  \ingroup GenericGeometry
     *  \brief   generic implementation of a DUNE (global) geometry
     *
     *  Geometry inherits all its features from Geometry. It only add
     *  GlobalGeometryTraits< Grid > as Traits parameter to the template
     *  parameter list.
     *
     * \tparam mydim Dimension of the entity
     * \tparam cdom Dimension of the coordinate space
     * \tparam Grid The grid this geometry will be used in
     */
    template< int mydim, int cdim, class Grid >
    class Geometry
      : public BasicGeometry< mydim, GlobalGeometryTraits< Grid > >
    {
      typedef BasicGeometry< mydim, GlobalGeometryTraits< Grid > > Base;

    protected:
      typedef typename Base :: Mapping Mapping;

    public:
      /** \brief Default constructor */
      Geometry ()
        : Base()
      {}

      /** \brief Constructor accepting a mapping */
      explicit Geometry ( Mapping &mapping )
        : Base( mapping )
      {}

      /** \brief Copy constructor from another geometry */
      template< class Geo >
      explicit Geometry ( const Geo &geo )
        : Base( geo.type(), geo )
      {}

      /** \brief Constructor with a GeometryType and a set of coordinates */
      template< class CoordVector >
      Geometry ( const GeometryType &type,
                 const CoordVector &coords )
        : Base( type, coords )
      {}

      /** \todo Please doc me! */
      template< int fatherdim >
      Geometry ( const Geometry< fatherdim, cdim, Grid > &father, int i )
        : Base( father, i )
      {}
    };



    // LocalGeometry
    // -------------

    /** \class   LocalGeometry
     *  \ingroup GenericGeometry
     *  \brief   generic implementation of a DUNE (local) geometry
     *
     *  Geometry inherits all its features from Geometry. It only adds
     *  LocalGeometryTraits< Grid > as Traits parameter to the template
     *  parameter list.
     *
     * \tparam mydim Dimension of the entity
     * \tparam cdom Dimension of the coordinate space
     * \tparam Grid The grid this geometry will be used in
     */
    template< int mydim, int cdim, class Grid >
    class LocalGeometry
      : public BasicGeometry< mydim, LocalGeometryTraits< Grid > >
    {
      typedef BasicGeometry< mydim, LocalGeometryTraits< Grid > > Base;

    protected:
      typedef typename Base :: Mapping Mapping;

    public:
      /** \brief Default constructor */
      LocalGeometry ()
        : Base()
      {}

      /** \brief Constructor accepting a mapping */
      explicit LocalGeometry ( Mapping &mapping )
        : Base( mapping )
      {}

      /** \brief Copy constructor from another geometry */
      template< class Geo >
      explicit LocalGeometry ( const Geo &geo )
        : Base( geo.type(), geo )
      {}

      /** \brief Constructor with a GeometryType and a set of coordinates */
      template< class CoordVector >
      LocalGeometry ( const GeometryType &type, const CoordVector &coords )
        : Base( type, coords )
      {}

      /** \todo Please doc me! */
      template< int fatherdim >
      LocalGeometry ( const Geometry< fatherdim, cdim, Grid > &father, int i )
        : Base( father, i )
      {}
    };

  }

}

#endif
