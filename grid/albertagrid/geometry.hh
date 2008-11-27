// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_GEOMETRY
#define DUNE_ALBERTA_GEOMETRY

#include <dune/grid/common/geometry.hh>

namespace Dune
{

  // Forward Declarations
  // --------------------

  template< int dim, int dimworld >
  class AlbertaGrid;



  // AlbertaGridGeometry
  // -------------------

  /** \class AlbertaGridGeometry
   *  \brief geometry implementation for AlbertaGrid
   *
   *  Defines the geometry part of a mesh entity. Works for all dimensions,
   *  element types and dim of world. Provides reference element and mapping
   *  between local and global coordinates.
   *
   *  \tparam  mydim    dimension of the element (0 <= dim <= 3)
   *  \tparam  cdim     dimension of global coordinates
   *  \tparam  GridImp  grid implementation
   *                    (always const AlbertaGrid< dim, dimworld >)
   */
  template< int mydim, int cdim, class GridImp >
  class AlbertaGridGeometry
  //: public GeometryDefaultImplementation<mydim,cdim,GridImp,AlbertaGridGeometry>
  {
    typedef AlbertaGridGeometry< mydim, cdim, GridImp > This;

    // remember type of the grid
    typedef GridImp Grid;

    // dimension of barycentric coordinates
    static const int dimbary = mydim + 1;

  public:
    //! type of coordinates
    typedef albertCtype ctype;

    static const int dimension = Grid :: dimension;
    static const int mydimension = mydim;
    static const int codimension = dimension - mydimension;
    static const int coorddimension = cdim;

    //! Default constructor
    AlbertaGridGeometry();

    //! constructor building geometry in father
    AlbertaGridGeometry( const int child, const int orientation );

    //! return the element type identifier
    //! line , triangle or tetrahedron, depends on dim
    GeometryType type () const;

    /** \brief obtain the number of corners of this element */
    int corners () const;

    //! access to coordinates of corners. Index is the number of the corner
    const FieldVector< ctype, cdim > &operator[] (int i) const;

    //! maps a local coordinate within reference element to
    //! global coordinate in element
    FieldVector< ctype, cdim >
    global ( const FieldVector< ctype, mydim> &local ) const;

    //! maps a global coordinate within the element to a
    //! local coordinate in its reference element
    FieldVector< ctype, mydim >
    local ( const FieldVector< ctype, cdim  > &global ) const;

    //! returns true if the point in local coordinates is inside reference element
    bool checkInside( const FieldVector< ctype, mydim > &local ) const;

    /*!
       Copy from sgrid.hh:

       Integration over a general element is done by integrating over the reference element
       and using the transformation from the reference element to the global element as follows:
       \f[\int\limits_{\Omega_e} f(x) dx = \int\limits_{\Omega_{ref}} f(g(l)) A(l) dl \f] where
       \f$g\f$ is the local to global mapping and \f$A(l)\f$ is the integration element.

       For a general map \f$g(l)\f$ involves partial derivatives of the map (surface element of
       the first kind if \f$d=2,w=3\f$, determinant of the Jacobian of the transformation for
       \f$d=w\f$, \f$\|dg/dl\|\f$ for \f$d=1\f$).

       For linear elements, the derivatives of the map with respect to local coordinates
       do not depend on the local coordinates and are the same over the whole element.

       For a structured mesh where all edges are parallel to the coordinate axes, the
       computation is the length, area or volume of the element is very simple to compute.

       Each grid module implements the integration element with optimal efficieny. This
       will directly translate in substantial savings in the computation of finite element
       stiffness matrices.
     */

    // A(l)
    ctype integrationElement ( const FieldVector< ctype, mydim > &local ) const;

    // volume if geometry
    ctype volume () const;

    //! can only be called for dim=dimworld!
    //! Note that if both methods are called on the same element, then
    //! call jacobianInverseTransposed first because integration element is calculated
    //! during calculation of the transposed of the jacobianInverse
    const FieldMatrix< ctype, cdim, mydim > &
    jacobianInverseTransposed ( const FieldVector< ctype, mydim > &local ) const;

    //***********************************************************************
    //  Methods that not belong to the Interface, but have to be public
    //***********************************************************************
    //! generate the geometry for the ALBERTA EL_INFO
    //! no interface method
    bool builtGeom( const Grid &grid, ALBERTA EL_INFO *elInfo, int subEntity );

    //! build geometry for intersectionSelfLocal and
    //! intersectionNeighborLocal
    template <class GeometryType, class LocalGeomType >
    bool builtLocalGeom(const GeometryType & geo , const LocalGeomType & lg,
                        ALBERTA EL_INFO *elInfo, int face);

    // init geometry with zeros
    //! no interface method
    void initGeom();

    //! print internal data
    //! no interface method
    void print (std::ostream& ss) const;

  private:
    // build geometry with local coords of child in reference element
    void buildGeomInFather(const int child, const int orientation);

    // calculate Matrix for Mapping from reference element to actual element
    void calcElMatrix () const;

    //! build the transposed of the jacobian inverse and store the volume
    void buildJacobianInverseTransposed () const;

    // template method for map the vertices of EL_INFO to the actual
    // coords needed for operator []
    static int mapVertices ( int i, int face, int edge, int vertex );
    static int mapVertices ( int subEntity, int i );

    // calculates the volume of the element
    ctype elDeterminant () const;

    //! the vertex coordinates
    FieldMatrix< ctype, mydim+1, cdim > coord_;

    mutable FieldMatrix< ctype, cdim, mydim > Jinv_; //!< storage for inverse of jacobian

    mutable FieldMatrix< ctype, cdim, mydim > elMat_; //!< storage for mapping matrix

    //! is true if elMat_ was calced
    mutable bool builtElMat_;
    //! is true if Jinv_ and volume_ is calced
    mutable bool builtinverse_;


    mutable bool calcedDet_; //! true if determinant was calculated
    mutable ctype elDet_; //!< storage of element determinant
  };



  // AlbertaGridLocalGeometryProvider
  // --------------------------------

  template< int dim, int dimworld >
  class AlbertaGridLocalGeometryProvider
  {
    typedef AlbertaGridLocalGeometryProvider< dim, dimworld > This;

    typedef AlbertaGrid< dim, dimworld > Grid;

  public:
    template< int codim >
    struct Codim
    {
      typedef Geometry< dim-codim, dim, const Grid, AlbertaGridGeometry >
      LocalGeometry;
    };

    typedef typename Codim< 0 > :: LocalGeometry LocalElementGeometry;

    static const int numChildren = 2;

  private:
    const LocalElementGeometry *geometryInFather_[ numChildren ][ 2 ];

    AlbertaGridLocalGeometryProvider ()
    {
      typedef MakeableInterfaceObject< LocalElementGeometry > LocalGeoObject;
      typedef typename LocalGeoObject :: ImplementationType LocalGeoImp;

      for( int child = 0; child < numChildren; ++child )
      {
        geometryInFather_[ child ][ 0 ]
          = new LocalGeoObject( LocalGeoImp( child, -1 ) );
        geometryInFather_[ child ][ 1 ]
          = new LocalGeoObject( LocalGeoImp( child, 1 ) );
      }
    }

    ~AlbertaGridLocalGeometryProvider ()
    {
      for( int child = 0; child < numChildren; ++child )
      {
        delete geometryInFather_[ child ][ 0 ];
        delete geometryInFather_[ child ][ 1 ];
      }
    }

  public:
    const LocalElementGeometry &
    geometryInFather ( int child, const int orientation = 1 ) const
    {
      assert( (child >= 0) && (child < numChildren) );
      assert( (orientation == 1) || (orientation == -1) );
      return *geometryInFather_[ child ][ (orientation + 1) / 2 ];
    }

    static const This &instance ()
    {
      static This theInstance;
      return theInstance;
    }
  };

}

#endif
