// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_ALBERTA_GRIDFACTORY_HH
#define DUNE_ALBERTA_GRIDFACTORY_HH

/** \file
 *  \author Martin Nolte
 *  \brief  specialization of the generic GridFactory for AlbertaGrid
 */

#include <algorithm>
#include <array>
#include <limits>
#include <map>
#include <memory>

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/gridfactory.hh>

#include <dune/grid/albertagrid/agrid.hh>

#if HAVE_ALBERTA

namespace Dune
{

  /** \brief specialization of the generic GridFactory for AlbertaGrid
   *
   *  \ingroup GridFactory
   *
   *  The GridFactory for AlbertaGrid adds some extensions to the standard
   *  GridFactoryInterface. It provides the following additional features:
   *  - It allows to set boundary ids via insertBoundary. For ALBERTA 1.2,
   *    these boundary ids are ignored, though.
   *  - For ALBERTA 3.0 and above, you can add face transformation to identify
   *    faces. This allows the construction of periodic grids.
   *  - The grid can be written in ALBERTA's native format for macro
   *    triangulations via write (both ASCII and XDR format are supported).
   *  - On grid creation, a name can be assigned to the grid.
   *  - On grid creation, the grid can be reordered such that ALBERTA uses
   *    the longest edge as refinement edge during recursive bisection.
   *  .
   */
  template< int dim, int dimworld >
  class GridFactory< AlbertaGrid< dim, dimworld > >
    : public GridFactoryInterface< AlbertaGrid< dim, dimworld > >
  {
    typedef GridFactory< AlbertaGrid< dim, dimworld > > This;

  public:
    //! type of grid this factory is for
    typedef AlbertaGrid< dim, dimworld > Grid;

    //! type of (scalar) coordinates
    typedef typename Grid::ctype ctype;

    //! dimension of the grid
    static const int dimension = Grid::dimension;
    //! dimension of the world
    static const int dimensionworld = Grid::dimensionworld;

    //! type of vector for world coordinates
    typedef FieldVector< ctype, dimensionworld > WorldVector;
    //! type of matrix from world coordinates to world coordinates
    typedef FieldMatrix< ctype, dimensionworld, dimensionworld > WorldMatrix;

    typedef DuneBoundaryProjection< dimensionworld > DuneProjection;
    typedef std::shared_ptr< const DuneProjection > DuneProjectionPtr;
    typedef Dune::BoundarySegment< dimension, dimensionworld > BoundarySegment;

    template< int codim >
    struct Codim
    {
      typedef typename Grid::template Codim< codim >::Entity Entity;
    };

  private:
    typedef Dune::BoundarySegmentWrapper< dimension, dimensionworld > BoundarySegmentWrapper;

    static const int numVertices
      = Alberta::NumSubEntities< dimension, dimension >::value;

    typedef Alberta::MacroElement< dimension > MacroElement;
    typedef Alberta::ElementInfo< dimension > ElementInfo;
    typedef Alberta::MacroData< dimension > MacroData;
    typedef Alberta::NumberingMap< dimension, Alberta::Dune2AlbertaNumbering > NumberingMap;
    typedef Alberta::DuneBoundaryProjection< dimension > Projection;

    typedef std::array< unsigned int, dimension > FaceId;
    typedef std::map< FaceId, size_t > BoundaryMap;

    class ProjectionFactory;

  public:
    //! are boundary ids supported by this factory?
    static const bool supportsBoundaryIds = true;
    //! is the factory able to create periodic meshes?
    static const bool supportPeriodicity = MacroData::supportPeriodicity;

    /** default constructor */
    GridFactory ()
      : globalProjection_( (const DuneProjection *) 0 )
    {
      macroData_.create();
    }

    virtual ~GridFactory ();

    /** \brief insert a vertex into the macro grid
     *
     *  \param[in]  pos  position of the vertex (in world coordinates)
     */
    virtual void insertVertex ( const WorldVector &pos )
    {
      macroData_.insertVertex( pos );
    }

    /** \brief insert an element into the macro grid
     *
     *  \param[in]  type      GeometryType of the new element
     *  \param[in]  vertices  indices of the element vertices (in DUNE numbering)
     */
    virtual void insertElement ( const GeometryType &type,
                                 const std::vector< unsigned int > &vertices )
    {
      if( (int)type.dim() != dimension )
        DUNE_THROW( AlbertaError, "Inserting element of wrong dimension: " << type.dim() );
      if( !type.isSimplex() )
        DUNE_THROW( AlbertaError, "Alberta supports only simplices." );

      if( vertices.size() != (size_t)numVertices )
        DUNE_THROW( AlbertaError, "Wrong number of vertices passed: " << vertices.size() << "." );

      int array[ numVertices ];
      for( int i = 0; i < numVertices; ++i )
        array[ i ] = vertices[ numberingMap_.alberta2dune( dimension, i ) ];
      macroData_.insertElement( array );
    }

    /** \brief mark a face as boundary (and assign a boundary id)
     *
     *  \internal
     *
     *  \param[in]  element  index of the element, the face belongs to
     *  \param[in]  face     local number of the face within the element
     *  \param[in]  id       boundary id to assign to the face
     *
     *  \note ALBERTA supports only boundary id in the range 1,...,127.
     */
    virtual void insertBoundary ( int element, int face, int id )
    {
      if( (id <= 0) || (id > 127) )
        DUNE_THROW( AlbertaError, "Invalid boundary id: " << id << "." );
      macroData_.boundaryId( element, numberingMap_.dune2alberta( 1, face ) ) = id;
    }

    /** \brief insert a boundary projection into the macro grid
     *
     *  \internal
     *
     *  \param[in]  type        geometry type of boundary face
     *  \param[in]  vertices    vertices of the boundary face
     *  \param[in]  projection  boundary projection
     *
     *  \note The grid takes control of the projection object.
     */
    virtual void
    insertBoundaryProjection ( const GeometryType &type,
                               const std::vector< unsigned int > &vertices,
                               const DuneProjection *projection )
    {
      if( (int)type.dim() != dimension-1 )
        DUNE_THROW( AlbertaError, "Inserting boundary face of wrong dimension: " << type.dim() );
      if( !type.isSimplex() )
        DUNE_THROW( AlbertaError, "Alberta supports only simplices." );

      FaceId faceId;
      if( vertices.size() != faceId.size() )
        DUNE_THROW( AlbertaError, "Wrong number of face vertices passed: " << vertices.size() << "." );
      for( size_t i = 0; i < faceId.size(); ++i )
        faceId[ i ] = vertices[ i ];
      std::sort( faceId.begin(), faceId.end() );

      typedef std::pair< typename BoundaryMap::iterator, bool > InsertResult;
      const InsertResult result = boundaryMap_.insert( std::make_pair( faceId, boundaryProjections_.size() ) );
      if( !result.second )
        DUNE_THROW( GridError, "Only one boundary projection can be attached to a face." );
      boundaryProjections_.push_back( DuneProjectionPtr( projection ) );
    }


    /** \brief insert a global (boundary) projection into the macro grid
     *
     *  \internal
     *
     *  \param[in]  projection  global (boundary) projection
     *
     *  \note The grid takes control of the projection object.
     */
    virtual void insertBoundaryProjection ( const DuneProjection *projection )
    {
      if( globalProjection_ )
        DUNE_THROW( GridError, "Only one global boundary projection can be attached to a grid." );
      globalProjection_ = DuneProjectionPtr( projection );
    }

    /** \brief insert a boundary segment into the macro grid
     *
     *  Only influences the ordering of the boundary segments
     *  \param[in]  vertices         vertex indices of boundary face
     */
    virtual void
    insertBoundarySegment ( const std::vector< unsigned int >& vertices )
    {
      insertBoundaryProjection( GeometryTypes::simplex( dimension-1 ), vertices, 0 );
    }

    /** \brief insert a shaped boundary segment into the macro grid
     *
     *  \param[in]  vertices         vertex indices of boundary face
     *  \param[in]  boundarySegment  geometric realization of shaped boundary
     */
    virtual void
    insertBoundarySegment ( const std::vector< unsigned int > &vertices,
                            const std::shared_ptr< BoundarySegment > &boundarySegment )
    {
      auto refSimplex = ReferenceElements< ctype, dimension-1 >::simplex();

      if( !boundarySegment )
        DUNE_THROW( GridError, "Trying to insert null as a boundary segment." );
      if( (int)vertices.size() != refSimplex.size( dimension-1 ) )
        DUNE_THROW( GridError, "Wrong number of face vertices passed: " << vertices.size() << "." );

      std::vector< WorldVector > coords( refSimplex.size( dimension-1 ) );
      for( int i = 0; i < dimension; ++i )
      {
        Alberta::GlobalVector &x = macroData_.vertex( vertices[ i ] );
        for( int j = 0; j < dimensionworld; ++j )
          coords[ i ][ j ] = x[ j ];
        if( ((*boundarySegment)( refSimplex.position( i, dimension-1 ) ) - coords[ i ]).two_norm() > 1e-6 )
          DUNE_THROW( GridError, "Boundary segment does not interpolate the corners." );
      }

      const GeometryType gt = refSimplex.type( 0, 0 );
      const DuneProjection *prj = new BoundarySegmentWrapper( gt, coords, boundarySegment );
      insertBoundaryProjection( gt, vertices, prj );
    }

    /** \brief add a face transformation (for periodic identification)
     *
     *  A face transformation is an affine mapping T from world coordinates
     *  to world coordinates. ALBERTA periodically identifies two faces f and g
     *  if T( f ) = g or T( g ) = f.
     *
     *  \param[in]  matrix  matrix describing the linear part of T
     *  \param[in]  shift   vector describing T( 0 )
     *
     *  \note ALBERTA requires the matrix to be orthogonal.
     *
     *  \note ALBERTA automatically adds the inverse transformation.
     */
    void insertFaceTransformation ( const WorldMatrix &matrix, const WorldVector &shift );

    /** \brief mark the longest edge as refinemet edge
     *
     *  Marking the longest edge avoids cycles in the recursive bisection
     *  algorithm, if the longest edge of each element is unique. It also
     *  makes sure the angles degenerate least. It can, hoowever, produce
     *  more nonlocal refinements than necessary. Therefore this feature is
     *  disabled by default.
     */
    void markLongestEdge ()
    {
      macroData_.markLongestEdge();
    }

    /** \brief finalize grid creation and hand over the grid
     *
     *  This version of createGrid is original to the AlbertaGrid grid factroy,
     *  allowing to specity a grid name.
     *
     *  \returns a pointer to the newly created grid
     *
     *  \note The caller takes responsibility of freeing the memory allocated
     *        for the grid.
     *  \note ALBERTA's grid factory provides a static method for freeing the
     *        grid (destroyGrid).
     */
    std::unique_ptr<Grid> createGrid ()
    {
      macroData_.finalize();
      if( macroData_.elementCount() == 0 )
        DUNE_THROW( GridError, "Cannot create empty AlbertaGrid." );
      if( dimension < 3 )
        macroData_.setOrientation( Alberta::Real( 1 ) );
      assert( macroData_.checkNeighbors() );
      macroData_.checkCycles();
      ProjectionFactory projectionFactory( *this );
      return std::make_unique<Grid>( macroData_, projectionFactory );
    }

    /** \brief destroy a grid previously obtain from this factory
     *
     *  \param[in]  grid  pointer to the grid to destroy
     */
    static void destroyGrid ( Grid *grid )
    {
      delete grid;
    }

    /** \brief write out the macro triangulation in native grid file format
     *
     *  \param[in]  filename  name of the file to write to
     *
     *  \returns \c true on success
     */
    bool write ( const std::string &filename )
    {
      macroData_.finalize();
      if( dimension < 3 )
        macroData_.setOrientation( Alberta::Real( 1 ) );
      assert( macroData_.checkNeighbors() );
      return macroData_.write( filename, false );
    }

    virtual unsigned int
    insertionIndex ( const typename Codim< 0 >::Entity &entity ) const
    {
      return insertionIndex( entity.impl().elementInfo() );
    }

    virtual unsigned int
    insertionIndex ( const typename Codim< dimension >::Entity &entity ) const
    {
      const int elIndex = insertionIndex( entity.impl().elementInfo() );
      const typename MacroData::ElementId &elementId = macroData_.element( elIndex );
      return elementId[ entity.impl().subEntity() ];
    }

    virtual unsigned int
    insertionIndex ( const typename Grid::LeafIntersection &intersection ) const
    {
      const Grid &grid = intersection.impl().grid();
      const ElementInfo &elementInfo = intersection.impl().elementInfo();
      const int face = grid.generic2alberta( 1, intersection.indexInInside() );
      return insertionIndex( elementInfo, face );
    }

    virtual bool
    wasInserted ( const typename Grid::LeafIntersection &intersection ) const
    {
      return (insertionIndex( intersection ) < std::numeric_limits< unsigned int >::max());
    }

  private:
    unsigned int insertionIndex ( const ElementInfo &elementInfo ) const;
    unsigned int insertionIndex ( const ElementInfo &elementInfo, const int face ) const;

    FaceId faceId ( const ElementInfo &elementInfo, const int face ) const;

    MacroData macroData_;
    NumberingMap numberingMap_;
    DuneProjectionPtr globalProjection_;
    BoundaryMap boundaryMap_;
    std::vector< DuneProjectionPtr > boundaryProjections_;
  };


  template< int dim, int dimworld >
  GridFactory< AlbertaGrid< dim, dimworld > >::~GridFactory ()
  {
    macroData_.release();
  }


  template< int dim, int dimworld >
  inline void
  GridFactory< AlbertaGrid< dim, dimworld > >
  ::insertFaceTransformation ( const WorldMatrix &matrix, const WorldVector &shift )
  {
    // make sure the matrix is orthogonal
    for( int i = 0; i < dimworld; ++i )
      for( int j = 0; j < dimworld; ++j )
      {
        const ctype delta = (i == j ? ctype( 1 ) : ctype( 0 ));
        const ctype epsilon = (8*dimworld)*std::numeric_limits< ctype >::epsilon();

        if( std::abs( matrix[ i ] * matrix[ j ] - delta ) > epsilon )
        {
          DUNE_THROW( AlbertaError,
                      "Matrix of face transformation is not orthogonal." );
        }
      }

    // copy matrix
    Alberta::GlobalMatrix M;
    for( int i = 0; i < dimworld; ++i )
      for( int j = 0; j < dimworld; ++j )
        M[ i ][ j ] = matrix[ i ][ j ];

    // copy shift
    Alberta::GlobalVector t;
    for( int i = 0; i < dimworld; ++i )
      t[ i ] = shift[ i ];

    // insert into ALBERTA macro data
    macroData_.insertWallTrafo( M, t );
  }


  template< int dim, int dimworld >
  inline unsigned int
  GridFactory< AlbertaGrid< dim, dimworld > >
  ::insertionIndex ( const ElementInfo &elementInfo ) const
  {
    const MacroElement &macroElement = elementInfo.macroElement();
    const unsigned int index = macroElement.index;

#ifndef NDEBUG
    const typename MacroData::ElementId &elementId = macroData_.element( index );
    for( int i = 0; i <= dimension; ++i )
    {
      const Alberta::GlobalVector &x = macroData_.vertex( elementId[ i ] );
      const Alberta::GlobalVector &y = macroElement.coordinate( i );
      for( int j = 0; j < dimensionworld; ++j )
      {
        if( x[ j ] != y[ j ] )
          DUNE_THROW( GridError, "Vertex in macro element does not coincide with same vertex in macro data structure." );
      }
    }
#endif // #ifndef NDEBUG

    return index;
  }


  template< int dim, int dimworld >
  inline unsigned int
  GridFactory< AlbertaGrid< dim, dimworld > >
  ::insertionIndex ( const ElementInfo &elementInfo, const int face ) const
  {
    typedef typename BoundaryMap::const_iterator Iterator;
    const Iterator it = boundaryMap_.find( faceId( elementInfo, face ) );
    if( it != boundaryMap_.end() )
      return it->second;
    else
      return std::numeric_limits< unsigned int >::max();
  }


  template< int dim, int dimworld >
  inline typename GridFactory< AlbertaGrid< dim, dimworld > >::FaceId
  GridFactory< AlbertaGrid< dim, dimworld > >
  ::faceId ( const ElementInfo &elementInfo, const int face ) const
  {
    const unsigned int index = insertionIndex( elementInfo );
    const typename MacroData::ElementId &elementId = macroData_.element( index );

    FaceId faceId;
    for( size_t i = 0; i < faceId.size(); ++i )
    {
      const int k = Alberta::MapVertices< dimension, 1 >::apply( face, i );
      faceId[ i ] = elementId[ k ];
    }
    std::sort( faceId.begin(), faceId.end() );
    return faceId;
  }



  // GridFactory::ProjectionFactory
  // ------------------------------

  template< int dim, int dimworld >
  class GridFactory< AlbertaGrid< dim, dimworld > >::ProjectionFactory
    : public Alberta::ProjectionFactory< Alberta::DuneBoundaryProjection< dim >, ProjectionFactory >
  {
    typedef ProjectionFactory This;
    typedef Alberta::ProjectionFactory< Alberta::DuneBoundaryProjection< dim >, ProjectionFactory > Base;

    typedef typename Dune::GridFactory< AlbertaGrid< dim, dimworld > > Factory;

  public:
    typedef typename Base::Projection Projection;
    typedef typename Base::ElementInfo ElementInfo;

    typedef typename Projection::Projection DuneProjection;

    ProjectionFactory( const GridFactory &gridFactory )
      : gridFactory_( gridFactory )
    {}

    bool hasProjection ( const ElementInfo &elementInfo, const int face ) const
    {
      if( gridFactory().globalProjection_ )
        return true;

      const unsigned int index = gridFactory().insertionIndex( elementInfo, face );
      if( index < std::numeric_limits< unsigned int >::max() )
        return bool( gridFactory().boundaryProjections_[ index ] );
      else
        return false;
    }

    bool hasProjection ( const ElementInfo & /* elementInfo */ ) const
    {
      return bool( gridFactory().globalProjection_ );
    }

    Projection projection ( const ElementInfo &elementInfo, const int face ) const
    {
      const unsigned int index = gridFactory().insertionIndex( elementInfo, face );
      if( index < std::numeric_limits< unsigned int >::max() )
      {
        const DuneProjectionPtr &projection = gridFactory().boundaryProjections_[ index ];
        if( projection )
          return Projection( projection );
      }

      assert( gridFactory().globalProjection_ );
      return Projection( gridFactory().globalProjection_ );
    };

    Projection projection ( const ElementInfo & /* elementInfo */ ) const
    {
      assert( gridFactory().globalProjection_ );
      return Projection( gridFactory().globalProjection_ );
    };

    const GridFactory &gridFactory () const
    {
      return gridFactory_;
    }

  private:
    const GridFactory &gridFactory_;
  };

}

#endif // #if HAVE_ALBERTA

#endif // #ifndef DUNE_ALBERTA_GRIDFACTORY_HH
