// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_GEOMETRY_HH
#define DUNE_ALBERTA_GEOMETRY_HH

#include <dune/grid/common/geometry.hh>
#include <dune/grid/albertagrid/misc.hh>
#include <dune/grid/albertagrid/elementinfo.hh>

#if HAVE_ALBERTA

namespace Dune
{

  // Forward Declarations
  // --------------------

  template< int dim, int dimworld >
  class AlbertaGrid;



  // AlbertaGridCoordinateReader
  // ---------------------------

  template< int codim, class GridImp >
  struct AlbertaGridCoordinateReader
  {
    typedef typename std::remove_const< GridImp >::type Grid;

    static constexpr int dimension = Grid::dimension;
    static constexpr int codimension = codim;
    static constexpr int mydimension = dimension - codimension;
    static constexpr int coorddimension = Grid::dimensionworld;

    typedef Alberta::Real ctype;

    typedef Alberta::ElementInfo< dimension > ElementInfo;
    typedef FieldVector< ctype, coorddimension > Coordinate;

    AlbertaGridCoordinateReader ( const GridImp &grid,
                                  const ElementInfo &elementInfo,
                                  int subEntity )
      : grid_( grid ),
        elementInfo_( elementInfo ),
        subEntity_( subEntity )
    {}

    const ElementInfo &elementInfo () const
    {
      return elementInfo_;
    }

    void coordinate ( int i, Coordinate &x ) const
    {
      assert( !elementInfo_ == false );
      assert( (i >= 0) && (i <= mydimension) );

      const int k = mapVertices( subEntity_, i );
      const Alberta::GlobalVector &coord = grid_.getCoord( elementInfo_, k );
      for( int j = 0; j < coorddimension; ++j )
        x[ j ] = coord[ j ];
    }

    bool hasDeterminant () const
    {
      return false;
    }

    ctype determinant () const
    {
      assert( hasDeterminant() );
      return ctype( 0 );
    }

  private:
    static int mapVertices ( int subEntity, int i )
    {
      return Alberta::MapVertices< dimension, codimension >::apply( subEntity, i );
    }

    const Grid &grid_;
    const ElementInfo &elementInfo_;
    const int subEntity_;
  };



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
  {
    typedef AlbertaGridGeometry< mydim, cdim, GridImp > This;

    // remember type of the grid
    typedef GridImp Grid;

    // dimension of barycentric coordinates
    static constexpr int dimbary = mydim + 1;

  public:
    //! type of coordinates
    typedef Alberta::Real ctype;

    static constexpr int dimension = Grid :: dimension;
    static constexpr int mydimension = mydim;
    static constexpr int codimension = dimension - mydimension;
    static constexpr int coorddimension = cdim;

    typedef FieldVector< ctype, mydimension > LocalCoordinate;
    typedef FieldVector< ctype, coorddimension > GlobalCoordinate;

    typedef FieldMatrix< ctype, mydimension, coorddimension > JacobianTransposed;
    typedef FieldMatrix< ctype, coorddimension, mydimension > JacobianInverseTransposed;

    typedef FieldMatrix< ctype, coorddimension, mydimension > Jacobian;
    typedef FieldMatrix< ctype, mydimension, coorddimension > JacobianInverse;

  private:
    static constexpr int numCorners = mydimension + 1;

    typedef FieldMatrix< ctype, numCorners, coorddimension > CoordMatrix;

  public:
    AlbertaGridGeometry ()
    {
      invalidate();
    }

    template< class CoordReader >
    AlbertaGridGeometry ( const CoordReader &coordReader )
    {
      build( coordReader );
    }

    /** \brief obtain the type of reference element */
    GeometryType type () const
    {
      return GeometryTypes::simplex( mydimension );
    }

    //! returns always true since we only have simplices
    bool affine () const { return true; }

    /** \brief number of corner the geometry */
    int corners () const
    {
      return numCorners;
    }

    /** \brief obtain the i-th corner of this geometry */
    GlobalCoordinate corner ( const int i ) const
    {
      assert( (i >= 0) && (i < corners()) );
      return coord_[ i ];
    }

    /** \brief return center of geometry */
    GlobalCoordinate center () const
    {
      return centroid_;
    }

    /** \brief map a point from the reference element to the geometry */
    GlobalCoordinate global ( const LocalCoordinate &local ) const;

    /** \brief map a point from the geometry to the reference element */
    LocalCoordinate local ( const GlobalCoordinate &global ) const;

    /** \brief integration element of the geometry mapping
     *
     *  \note This method is not part of the geometry interface. It is used
     *        internally only.
     */
    ctype integrationElement () const
    {
      assert( calcedDet_ );
      return elDet_;
    }

    /** \brief integration element of the geometry mapping */
    ctype integrationElement ( [[maybe_unused]] const LocalCoordinate &local ) const
    {
      return integrationElement();
    }

    /** \brief volume of geometry */
    ctype volume () const
    {
      return integrationElement() / ctype( factorial(mydimension) );
    }

    /** \brief transposed of the geometry mapping's Jacobian
     *
     *  \note This method is not part of the geometry interface. It is used
     *        internally only.
     */
    const JacobianTransposed &jacobianTransposed () const;

    /** \brief transposed of the geometry mapping's Jacobian */
    const JacobianTransposed &
    jacobianTransposed ( [[maybe_unused]] const LocalCoordinate &local ) const
    {
      return jacobianTransposed();
    }

    /** \brief transposed inverse of the geometry mapping's Jacobian
     *
     *  \note This method is not part of the geometry interface. It is used
     *        internally only.
     */
    const JacobianInverseTransposed &jacobianInverseTransposed () const;

    /** \brief transposed inverse of the geometry mapping's Jacobian */
    const JacobianInverseTransposed &
    jacobianInverseTransposed ( [[maybe_unused]] const LocalCoordinate &local ) const
    {
      return jacobianInverseTransposed();
    }

    /** \brief geometry mapping's Jacobian */
    Jacobian jacobian ( const LocalCoordinate &local ) const
    {
      return jacobianTransposed(local).transposed();
    }

    /** \brief inverse of the geometry mapping's Jacobian */
    JacobianInverse jacobianInverse ( const LocalCoordinate &local ) const
    {
      return jacobianInverseTransposed(local).transposed();
    }

    //***********************************************************************
    //  Methods that not belong to the Interface, but have to be public
    //***********************************************************************

    /** \brief invalidate the geometry */
    void invalidate ()
    {
      builtJT_ = false;
      builtJTInv_ = false;
      calcedDet_ = false;
    }

    /** \brief build the geometry from a coordinate reader */
    template< class CoordReader >
    void build ( const CoordReader &coordReader );

    void print ( std::ostream &out ) const;

  private:
    // calculates the volume of the element
    ctype elDeterminant () const
    {
      return std::abs( Alberta::determinant( jacobianTransposed() ) );
    }

    //! the vertex coordinates
    CoordMatrix coord_;

    //! the center/centroid
    GlobalCoordinate centroid_;

    // storage for the transposed of the jacobian
    mutable JacobianTransposed jT_;

    // storage for the transposed inverse of the jacboian
    mutable JacobianInverseTransposed jTInv_;

    // has jT_ been computed, yet?
    mutable bool builtJT_;
    // has jTInv_ been computed, yet?
    mutable bool builtJTInv_;

    mutable bool calcedDet_; //!< true if determinant was calculated
    mutable ctype elDet_; //!< storage of element determinant
  };



  // AlbertaGridGlobalGeometry
  // -------------------------

  template< int mydim, int cdim, class GridImp >
  class AlbertaGridGlobalGeometry
    : public AlbertaGridGeometry< mydim, cdim, GridImp >
  {
    typedef AlbertaGridGlobalGeometry< mydim, cdim, GridImp > This;
    typedef AlbertaGridGeometry< mydim, cdim, GridImp > Base;

  public:
    AlbertaGridGlobalGeometry ()
      : Base()
    {}

    template< class CoordReader >
    AlbertaGridGlobalGeometry ( const CoordReader &coordReader )
      : Base( coordReader )
    {}
  };


#if !DUNE_ALBERTA_CACHE_COORDINATES
  template< int dim, int cdim >
  class AlbertaGridGlobalGeometry< dim, cdim, const AlbertaGrid< dim, cdim > >
  {
    typedef AlbertaGridGlobalGeometry< dim, cdim, const AlbertaGrid< dim, cdim > > This;

    // remember type of the grid
    typedef AlbertaGrid< dim, cdim > Grid;

    // dimension of barycentric coordinates
    static constexpr int dimbary = dim + 1;

    typedef Alberta::ElementInfo< dim > ElementInfo;

  public:
    //! type of coordinates
    typedef Alberta::Real ctype;

    static constexpr int dimension = Grid::dimension;
    static constexpr int mydimension = dimension;
    static constexpr int codimension = dimension - mydimension;
    static constexpr int coorddimension = cdim;

    typedef FieldVector< ctype, mydimension > LocalCoordinate;
    typedef FieldVector< ctype, coorddimension > GlobalCoordinate;

    typedef FieldMatrix< ctype, mydimension, coorddimension > JacobianTransposed;
    typedef FieldMatrix< ctype, coorddimension, mydimension > JacobianInverseTransposed;

    typedef FieldMatrix< ctype, coorddimension, mydimension > Jacobian;
    typedef FieldMatrix< ctype, mydimension, coorddimension > JacobianInverse;

  private:
    static constexpr int numCorners = mydimension + 1;

    typedef FieldMatrix< ctype, numCorners, coorddimension > CoordMatrix;

  public:
    AlbertaGridGlobalGeometry ()
    {
      invalidate();
    }

    template< class CoordReader >
    AlbertaGridGlobalGeometry ( const CoordReader &coordReader )
    {
      build( coordReader );
    }

    /** \brief obtain the type of reference element */
    GeometryType type () const
    {
      return GeometryTypes::simplex( mydimension );
    }

    /** \brief number of corner the geometry */
    int corners () const
    {
      return numCorners;
    }

    /** \brief obtain the i-th corner of this geometry */
    GlobalCoordinate corner ( const int i ) const
    {
      assert( (i >= 0) && (i < corners()) );
      const Alberta::GlobalCoordinate &x = elementInfo_.coordinate( i );
      GlobalCoordinate y;
      for( int j = 0; j < coorddimension; ++j )
        y[ j ] = x[ j ];
      return y;
    }

    /** \brief return center of geometry */
    GlobalCoordinate center () const
    {
      GlobalCoordinate centroid_ = corner( 0 );
      for( int i = 1; i < numCorners; ++i )
        centroid_ += corner( i );
      centroid_ *= ctype( 1 ) / ctype( numCorners );
      return centroid_;
    }

    /** \brief map a point from the reference element to the geometry */
    GlobalCoordinate global ( const LocalCoordinate &local ) const;

    /** \brief map a point from the geometry to the reference element */
    LocalCoordinate local ( const GlobalCoordinate &global ) const;

    /** \brief integration element of the geometry mapping
     *
     *  \note This method is not part of the geometry interface. It is used
     *        internally only.
     */
    ctype integrationElement () const
    {
      return elementInfo_.geometryCache().integrationElement();
    }

    /** \brief integration element of the geometry mapping */
    ctype integrationElement ( const LocalCoordinate &local ) const
    {
      return integrationElement();
    }

    /** \brief volume of geometry */
    ctype volume () const
    {
      return integrationElement() / ctype( factorial(mydimension) );
    }

    /** \brief transposed of the geometry mapping's Jacobian
     *
     *  \note This method is not part of the geometry interface. It is used
     *        internally only.
     */
    const JacobianTransposed &jacobianTransposed () const
    {
      return elementInfo_.geometryCache().jacobianTransposed();
    }

    /** \brief transposed of the geometry mapping's Jacobian */
    const JacobianTransposed &
    jacobianTransposed ( const LocalCoordinate &local ) const
    {
      return jacobianTransposed();
    }

    /** \brief transposed inverse of the geometry mapping's Jacobian
     *
     *  \note This method is not part of the geometry interface. It is used
     *        internally only.
     */
    const JacobianInverseTransposed &jacobianInverseTransposed () const
    {
      return elementInfo_.geometryCache().jacobianInverseTransposed();
    }

    /** \brief transposed inverse of the geometry mapping's Jacobian */
    const JacobianInverseTransposed &
    jacobianInverseTransposed ( const LocalCoordinate &local ) const
    {
      return jacobianInverseTransposed();
    }

    /** \brief geometry mapping's Jacobian */
    Jacobian jacobian ( const LocalCoordinate &local ) const
    {
      return jacobianTransposed(local).transposed();
    }

    /** \brief inverse of the geometry mapping's Jacobian */
    JacobianInverse jacobianInverse ( const LocalCoordinate &local ) const
    {
      return jacobianInverseTransposed(local).transposed();
    }

    //***********************************************************************
    //  Methods that not belong to the Interface, but have to be public
    //***********************************************************************

    /** \brief invalidate the geometry */
    void invalidate ()
    {
      elementInfo_ = ElementInfo();
    }

    /** \brief build the geometry from a coordinate reader */
    template< class CoordReader >
    void build ( const CoordReader &coordReader )
    {
      elementInfo_ = coordReader.elementInfo();
    }

  private:
    ElementInfo elementInfo_;
  };
#endif // #if !DUNE_ALBERTA_CACHE_COORDINATES



  // AlbertaGridLocalGeometryProvider
  // --------------------------------

  template< class Grid >
  class AlbertaGridLocalGeometryProvider
  {
    typedef AlbertaGridLocalGeometryProvider< Grid > This;

  public:
    typedef typename Grid::ctype ctype;

    static constexpr int dimension = Grid::dimension;

    template< int codim >
    struct Codim
    {
      typedef AlbertaGridGeometry< dimension-codim, dimension, Grid > LocalGeometry;
    };

    typedef typename Codim< 0 >::LocalGeometry LocalElementGeometry;
    typedef typename Codim< 1 >::LocalGeometry LocalFaceGeometry;

    static constexpr int numChildren = 2;
    static constexpr int numFaces = dimension + 1;

    static constexpr int minFaceTwist = Alberta::Twist< dimension, dimension-1 >::minTwist;
    static constexpr int maxFaceTwist = Alberta::Twist< dimension, dimension-1 >::maxTwist;
    static constexpr int numFaceTwists = maxFaceTwist - minFaceTwist + 1;

  private:
    struct GeoInFatherCoordReader;
    struct FaceCoordReader;

    AlbertaGridLocalGeometryProvider ()
    {
      buildGeometryInFather();
      buildFaceGeometry();
    }

    ~AlbertaGridLocalGeometryProvider ()
    {
      for( int child = 0; child < numChildren; ++child )
      {
        delete geometryInFather_[ child ][ 0 ];
        delete geometryInFather_[ child ][ 1 ];
      }

      for( int i = 0; i < numFaces; ++i )
      {
        for( int j = 0; j < numFaceTwists; ++j )
          delete faceGeometry_[ i ][ j ];
      }
    }

    void buildGeometryInFather();
    void buildFaceGeometry();

  public:
    const LocalElementGeometry &
    geometryInFather ( int child, const int orientation = 1 ) const
    {
      assert( (child >= 0) && (child < numChildren) );
      assert( (orientation == 1) || (orientation == -1) );
      return *geometryInFather_[ child ][ (orientation + 1) / 2 ];
    }

    const LocalFaceGeometry &
    faceGeometry ( int face, int twist = 0 ) const
    {
      assert( (face >= 0) && (face < numFaces) );
      assert( (twist >= minFaceTwist) && (twist <= maxFaceTwist) );
      return *faceGeometry_[ face ][ twist - minFaceTwist ];
    }

    static const This &instance ()
    {
      static This theInstance;
      return theInstance;
    }

  private:
    template< int codim >
    static int mapVertices ( int subEntity, int i )
    {
      return Alberta::MapVertices< dimension, codim >::apply( subEntity, i );
    }

    const LocalElementGeometry *geometryInFather_[ numChildren ][ 2 ];
    const LocalFaceGeometry *faceGeometry_[ numFaces ][ numFaceTwists ];
  };

} // namespace Dune

#endif // #if HAVE_ALBERTA

#endif // #ifndef DUNE_ALBERTA_GEOMETRY_HH
