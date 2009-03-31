// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_INTERSECTION_CC
#define DUNE_ALBERTA_INTERSECTION_CC

#include <dune/grid/albertagrid/intersection.hh>

// set to 1 to use ElementInfo::leafNeighbor instead of faking leaf neighbor
#define TRAVERSE_LEAFNEIGHBOR 1

namespace Dune
{

  // AlbertaGridIntersection
  // -----------------------

  template< class GridImp >
  inline AlbertaGridIntersection< GridImp >
  ::AlbertaGridIntersection ( const GridImp &grid, int level )
    : grid_( grid ),
      neighborCount_( dimension+1 ),
      elementInfo_(),
      fakeNeighObj_( LocalGeometryImp() ),
      fakeSelfObj_ ( LocalGeometryImp() ),
      neighGlobObj_( GeometryImp() ),
      neighborInfo_()
  {}


  template< class GridImp >
  inline AlbertaGridIntersection< GridImp >
  ::AlbertaGridIntersection ( const This &other )
    : grid_( other.grid_ ),
      neighborCount_( other.neighborCount_ ),
      elementInfo_( other.elementInfo_ ),
      fakeNeighObj_( LocalGeometryImp() ),
      fakeSelfObj_ ( LocalGeometryImp() ),
      neighGlobObj_( GeometryImp() ),
      neighborInfo_()
  {}


  template< class GridImp >
  inline void
  AlbertaGridIntersection< GridImp >
  ::first ( const EntityImp &entity, int level )
  {
    neighborCount_ = 0;
    neighborInfo_ = ElementInfo();
    elementInfo_ = entity.elementInfo_;

    assert( !elementInfo_ == false );
    assert( elementInfo_.level() == level );
  }


  template< class GridImp >
  inline void AlbertaGridIntersection< GridImp >::done ()
  {
    neighborCount_ = dimension+1;
    neighborInfo_ = ElementInfo();
    elementInfo_ = ElementInfo();
  }


  // assignment operator
  template< class GridImp >
  inline void
  AlbertaGridIntersection< GridImp >::assign ( const This &other )
  {
    // only assign iterators from the same grid
    assert( &(this->grid_) == &(other.grid_) );

    neighborCount_ = other.neighborCount_;
    elementInfo_ = other.elementInfo_;
    neighborInfo_ = ElementInfo();
  }


  template< class GridImp >
  inline bool
  AlbertaGridIntersection< GridImp >::equals ( const This &other ) const
  {
    const ALBERTA EL *e1 = elementInfo_.el();
    const ALBERTA EL *e2 = other.elementInfo_.el();
    return (e1 == e2) && (neighborCount_ == other.neighborCount_);
  }

  template< class GridImp >
  inline void AlbertaGridIntersection<GridImp>::increment()
  {
    neighborInfo_ = ElementInfo();
    ++neighborCount_;
    if( neighborCount_ > dimension )
      this->done();
  }

  template< class GridImp >
  inline typename AlbertaGridIntersection< GridImp >::EntityPointer
  AlbertaGridIntersection< GridImp >::outside () const
  {
    typedef AlbertaGridEntityPointer< 0, GridImp > EntityPointerImp;

    if( !neighborInfo_ )
    {
      assert( !elementInfo_ == false );

#if TRAVERSE_LEAFNEIGHBOR
      neighborInfo_ = elementInfo_.leafNeighbor( neighborCount_ );
#else
      neighborInfo_ = ElementInfo::createFake( grid_.meshPointer(), NULL, 0 );
      std::memcpy( &(neighborInfo_.elInfo()), &(elementInfo_.elInfo()), sizeof( ALBERTA EL_INFO ) );

      setupVirtEn();
#endif
    }

    assert( !neighborInfo_ == false );
    assert( neighborInfo_.el() != NULL );
    return EntityPointerImp( grid_, neighborInfo_, 0 );
  }

  template< class GridImp >
  inline typename AlbertaGridIntersection< GridImp >::EntityPointer
  AlbertaGridIntersection< GridImp >::inside () const
  {
    typedef AlbertaGridEntityPointer< 0, GridImp > EntityPointerImp;
    assert( !elementInfo_ == false );
    return EntityPointerImp( grid_, elementInfo_, 0 );
  }

  template< class GridImp >
  inline int
  AlbertaGridIntersection<GridImp>::boundaryId () const
  {
    assert( !elementInfo_ == false );

    // id of interior intersections is 0
    if( !boundary() )
      return 0;

    const int id = elementInfo_.boundaryId( neighborCount_ );
    assert( id != 0 );
    return id;
  }

  template< class GridImp >
  inline bool AlbertaGridIntersection<GridImp>::conforming() const
  {
    return true;
  }

  template< class GridImp >
  inline bool AlbertaGridIntersection< GridImp >::boundary() const
  {
    assert( !!elementInfo_ );
    return elementInfo_.isBoundary( neighborCount_ );
  }

  template< class GridImp >
  inline bool AlbertaGridIntersection< GridImp >::neighbor () const
  {
    assert( !!elementInfo_ );
    const ALBERTA EL_INFO &elInfo = elementInfo_.elInfo();
    return (elInfo.neigh[ neighborCount_ ] != NULL);
  }

  template<class GridImp>
  inline const typename AlbertaGridIntersection< GridImp >::NormalVector
  AlbertaGridIntersection< GridImp >
  ::unitOuterNormal ( const LocalCoordType &local ) const
  {
    NormalVector normal;
    calcOuterNormal( normal );
    normal *= (1.0 / normal.two_norm());
    return normal;
  }

  template< class GridImp >
  inline const typename AlbertaGridIntersection< GridImp >::NormalVector
  AlbertaGridIntersection< GridImp >
  ::integrationOuterNormal ( const LocalCoordType &local ) const
  {
    NormalVector normal;
    calcOuterNormal( normal );
    return normal;
  }

  template< class GridImp >
  inline const typename AlbertaGridIntersection< GridImp >::NormalVector
  AlbertaGridIntersection< GridImp >
  ::outerNormal( const LocalCoordType &local ) const
  {
    return integrationOuterNormal( local );
  }

  template< class GridImp >
  inline void
  AlbertaGridIntersection<GridImp>::calcOuterNormal( NormalVector &n ) const
  {
    DUNE_THROW( NotImplemented, "AlbertaGrid: outer normal for dim != dimworld "
                "has not been implemented, yet." );
  }

  template<>
  inline void
  AlbertaGridIntersection< const AlbertaGrid< 1, 1 > >
  ::calcOuterNormal ( NormalVector &n ) const
  {
    assert( !!elementInfo_ );
    const Alberta::GlobalVector &oppCoord = grid_.getCoord( elementInfo_, neighborCount_ );
    const Alberta::GlobalVector &myCoord = grid_.getCoord( elementInfo_, 1-neighborCount_ );
    n[ 0 ] = (myCoord[ 0 ] > oppCoord[ 0 ] ? ctype( 1 ) : -ctype( 1 ));
  }

  template<>
  inline void
  AlbertaGridIntersection< const AlbertaGrid< 2, 2 > >
  ::calcOuterNormal ( NormalVector &n ) const
  {
    assert( !!elementInfo_ );
    const Alberta::GlobalVector &coordOne = grid_.getCoord( elementInfo_, (neighborCount_+1)%3 );
    const Alberta::GlobalVector &coordTwo = grid_.getCoord( elementInfo_, (neighborCount_+2)%3 );

    n[ 0 ] = -(coordOne[ 1 ] - coordTwo[ 1 ]);
    n[ 1 ] =   coordOne[ 0 ] - coordTwo[ 0 ];
  }

  template<>
  inline void
  AlbertaGridIntersection< const AlbertaGrid< 3, 3 > >
  ::calcOuterNormal ( NormalVector &n ) const
  {
    assert( !!elementInfo_ );

    // in this case the orientation is negative, multiply by -1
#if (DUNE_ALBERTA_VERSION >= 0x200) || (DIM == 3)
    const ALBERTA EL_INFO &elInfo = elementInfo_.elInfo();
    const ctype val = (elInfo.orientation > 0) ? 1.0 : -1.0;
#else
    const ctype val = 1.0;
#endif

    // neighborCount_ is the local face number
    static const int faceVertices[ 4 ][ 3 ]
      = { {1,3,2}, {0,2,3}, {0,3,1}, {0,1,2} };
    const int *localFaces = faceVertices[ neighborCount_ ];

    const Alberta::GlobalVector &coord0 = grid_.getCoord( elementInfo_, localFaces[ 0 ] );
    const Alberta::GlobalVector &coord1 = grid_.getCoord( elementInfo_, localFaces[ 1 ] );
    const Alberta::GlobalVector &coord2 = grid_.getCoord( elementInfo_, localFaces[ 2 ] );

    FieldVector< ctype, dimensionworld > u;
    FieldVector< ctype, dimensionworld > v;
    for( int i = 0; i < dimension; ++i )
    {
      v[ i ] = coord1[ i ] - coord0[ i ];
      u[ i ] = coord2[ i ] - coord1[ i ];
    }

    // outNormal_ has length 3
    for( int i = 0; i < dimension; ++i )
    {
      const int j = (i+1)%dimension;
      const int k = (i+2)%dimension;
      n[ i ] = val * (u[ j ] * v[ k ] - u[ k ] * v[ j ]);
    }
  }


  template< class GridImp >
  inline AlbertaTransformation
  AlbertaGridIntersection< GridImp >::transformation () const
  {
    return AlbertaTransformation( elementInfo_.transformation( neighborCount_ ) );
  }


  template< class GridImp >
  inline const typename AlbertaGridIntersection< GridImp >::LocalGeometry &
  AlbertaGridIntersection< GridImp >::geometryInInside () const
  {
    assert( !!elementInfo_ );

    LocalGeometryImp &geo = GridImp::getRealImplementation( fakeSelfObj_ );
    const LocalCoordReader coordReader( inside()->geometry(), geometry() );
    geo.build( coordReader );
    return fakeSelfObj_;
  }


  template< class GridImp >
  inline const typename AlbertaGridIntersection< GridImp >::LocalGeometry &
  AlbertaGridIntersection< GridImp >::geometryInOutside () const
  {
    assert( neighbor() );

    LocalGeometryImp &geo = GridImp::getRealImplementation( fakeNeighObj_ );
    const LocalCoordReader coordReader( outside()->geometry(), geometry() );
    geo.build( coordReader );
    return fakeNeighObj_;
  }


  template< class GridImp >
  inline const typename AlbertaGridIntersection< GridImp >::Geometry &
  AlbertaGridIntersection< GridImp >::geometry () const
  {
    assert( !elementInfo_ == false );

    GeometryImp &geo = GridImp::getRealImplementation( neighGlobObj_ );
    const GlobalCoordReader coordReader( grid_, elementInfo_, neighborCount_ );
    geo.build( coordReader );
    return neighGlobObj_;
  }


  template< class GridImp >
  inline GeometryType AlbertaGridIntersection< GridImp >::type () const
  {
    return GeometryType( GeometryType::simplex, dimension-1 );
  }


  template< class GridImp >
  inline int AlbertaGridIntersection< GridImp >::numberInInside () const
  {
    const int oppVertex = neighborCount_;
    const int number = (dimension > 1 ? grid_.alberta2dune( 1, oppVertex ) : 1-oppVertex);

    typedef GenericGeometry::MapNumberingProvider< dimension > Numbering;
    const unsigned int tid = GenericGeometry::topologyId( inside()->type() );
    return Numbering::template dune2generic< dimension >( tid, number );
  }

  template< class GridImp >
  inline int AlbertaGridIntersection< GridImp >::numberInOutside () const
  {
    assert( !!elementInfo_ );
    const ALBERTA EL_INFO &elInfo = elementInfo_.elInfo();
    const int oppVertex = elInfo.opp_vertex[ neighborCount_ ];
    const int number = (dimension > 1 ? grid_.alberta2dune( 1, oppVertex ) : 1-oppVertex);

    typedef GenericGeometry::MapNumberingProvider< dimension > Numbering;
    const unsigned int tid = GenericGeometry::topologyId( inside()->type() );
    return Numbering::template dune2generic< dimension >( tid, number );
  }


  template< class GridImp >
  inline int
  AlbertaGridIntersection< GridImp >::twistInSelf () const
  {
    return elementInfo_.template twist< 1 >( neighborCount_ );
  }


  template< class GridImp >
  inline int
  AlbertaGridIntersection< GridImp >::twistInNeighbor () const
  {
    return elementInfo_.twistInNeighbor( neighborCount_ );
  }


  // setup neighbor element with the information of elInfo_
  template< class GridImp >
  inline bool AlbertaGridIntersection<GridImp>::neighborHasSameLevel () const
  {
    assert( neighbor() );
    assert( !elementInfo_ == false );
    const ALBERTA EL_INFO &elInfo = elementInfo_.elInfo();

    const ALBERTA EL * myEl    = elementInfo_.el();
    const ALBERTA EL * neighEl = elInfo.neigh[ neighborCount_ ];

    return (grid_.levelProvider() ( myEl ) == grid_.levelProvider() ( neighEl ));
  }


#if !TRAVERSE_LEAFNEIGHBOR
  //*****************************************
  //  setup for 2d
  //*****************************************
  template <class GridImp, int dimworld , int dim >
  struct SetupVirtualNeighbour;


  #if DIM == 2
  template <class GridImp>
  struct SetupVirtualNeighbour<GridImp,2,2>
  {
    enum { dim      = 2 };
    enum { dimworld = 2 };
    static int setupNeighInfo(GridImp & grid, const ALBERTA EL_INFO * elInfo,
                              const int vx, const int nb, ALBERTA EL_INFO * neighInfo)
    {

  #if CALC_COORD
      // vx is the face number in the neighbour
      const int (& neighmap)[dim] = ALBERTA AlbertHelp :: localTriangleFaceNumber[vx];
      // neighborCount is the face number in the actual element
      const int (& selfmap) [dim] = ALBERTA AlbertHelp :: localTriangleFaceNumber[nb];

      // copy the two edge vertices
      for(int i=0; i<dim; i++)
      {
        const ALBERTA REAL_D & coord = elInfo->coord[ selfmap[i] ];
        // here the twist is simple, just swap the vertices, and do this by hand
        ALBERTA REAL_D & newcoord    = neighInfo->coord[ neighmap[(dim-1) - i] ];
        for(int j=0; j<dimworld; j++) newcoord[j] = coord[j];
      }
  #endif

      //****************************************
      // twist is always 1
      return 1;
    }
  };
  #endif

  //******************************************
  //  setup for 3d
  //******************************************
  #if DIM == 3
  template <class GridImp>
  struct SetupVirtualNeighbour<GridImp,3,3>
  {
    enum { dim      = 3 };
    enum { dimworld = 3 };
    static int setupNeighInfo(GridImp & grid, const ALBERTA EL_INFO * elInfo,
                              const int vx, const int nb, ALBERTA EL_INFO * neighInfo)
    {
      // the face might be twisted when look from different elements
      // default is no, and then rthe orientation is -1
      int facemap[dim]   = {0,1,2};
      bool rightOriented = false;
      {
        const typename GridImp::HierarchicIndexSet &hIndexSet = grid.hierarchicIndexSet();

        int myvx[dim];
        int neighvx[dim];

        const int * vxmap = ALBERTA AlbertHelp :: localAlbertaFaceNumber[ vx ];
        const int * nbmap = ALBERTA AlbertHelp :: localAlbertaFaceNumber[ nb ];

        bool allRight = true;
        for(int i=0; i<dim; i++)
        {
          myvx[ i ] = hIndexSet.template subIndex( elInfo->el, dim, nbmap[ i ] );
          neighvx[ i ] = hIndexSet.template subIndex( neighInfo->el, dim, vxmap[ i ] );
          if( myvx[ i ] != neighvx[ i ] )
            allRight = false;
        }

        // note: if the vertices are equal then the face in the neighbor
        // is not oriented right, because all face are oriented math. pos when
        // one looks from the outside of the element.
        // if the vertices are the same, the face in the neighbor is therefore
        // wrong oriented
        if( !allRight )
        {
          for(int i=0; i<dim; i++)
          {
            if(myvx[i] != neighvx[i])
            {
              for(int k=1; k<dim; k++)
              {
                int newvx = (i+k) % dim;
                if( myvx[i] == neighvx[newvx] ) facemap[i] = newvx;
              }
            }
          }
          rightOriented = true;
        }
      }

      // TODO check infulence of orientation
      // is used when calculation the outerNormal
      neighInfo->orientation = ( rightOriented ) ? elInfo->orientation : -elInfo->orientation;

  #if CALC_COORD
      // vx is the face number in the neighbour
      const int * neighmap = ALBERTA AlbertHelp :: localAlbertaFaceNumber[vx];
      // neighborCount is the face number in the actual element
      const int * selfmap  = ALBERTA AlbertHelp :: localAlbertaFaceNumber[nb];

      // copy the three face vertices
      for(int i=0; i<dim; i++)
      {
        const ALBERTA REAL_D & coord = elInfo->coord[selfmap[i]];
        // here consider possible twist
        ALBERTA REAL_D & newcoord    = neighInfo->coord[ neighmap[ facemap[i] ] ];
        for(int j=0; j<dimworld; j++) newcoord[j] = coord[j];
      }
  #endif

      //****************************************
      if (facemap[1] == (facemap[0]+1)%3) {
        return facemap[0];
      }
      // twist
      return facemap[1]-3;
    }
  };
  #endif

  // setup neighbor element with the information of elInfo_
  template< class GridImp >
  inline void AlbertaGridIntersection< GridImp >::setupVirtEn () const
  {
    // if this assertion fails then outside was called without checking
    // neighbor first
    assert( neighbor() );
    assert( neighborCount_ < dimension+1 );

    assert( !elementInfo_ == false );
    const ALBERTA EL_INFO &elInfo = elementInfo_.elInfo();

    // set the neighbor element as element
    // use ALBERTA macro to get neighbour
    assert( !neighborInfo_ == false );
    ALBERTA EL_INFO &nbInfo = neighborInfo_.elInfo();

    nbInfo.el = elInfo.neigh[ neighborCount_ ];
    assert( nbInfo.el != NULL );
    nbInfo.level = grid_.levelProvider() ( nbInfo.el );

    const int vx = elInfo.opp_vertex[ neighborCount_ ];
    assert( (vx >= 0) && (vx < ElementInfo::maxNeighbors) );

    // reset neighbor information
    for( int i = 0; i <= dimension; ++i )
    {
      nbInfo.neigh[ i ] = NULL;
      nbInfo.opp_vertex[ i ] = 127;
    }

    // set origin
    nbInfo.neigh[ vx ] = elInfo.el;
    nbInfo.opp_vertex[ vx ] = neighborCount_;


  #if CALC_COORD
    // copy the one opposite vertex
    // the same for 2d and 3d
    {
      const ALBERTA REAL_D &coord = elInfo.opp_coord[ neighborCount_ ];
      ALBERTA REAL_D &newcoord = nbInfo.coord[ vx ];
      for( int j = 0; j < dimensionworld; ++j )
        newcoord[ j ] = coord[ j ];
    }
  #endif

    // setup coordinates of neighbour elInfo
    SetupVirtualNeighbour<GridImp,dimensionworld,dimension>::
    setupNeighInfo( this->grid_, &elInfo, vx, neighborCount_, &nbInfo );
  }
#endif // #if !TRAVERSE_LEAFNEIGHBOR



  // AlbertaGridIntersection::GlobalCoordReader
  // --------------------------------------------------

  template< class GridImp >
  struct AlbertaGridIntersection< GridImp >::GlobalCoordReader
  {
    typedef typename remove_const< GridImp >::type Grid;

    static const int dimension = Grid::dimension;
    static const int codimension = 1;
    static const int mydimension = dimension - codimension;
    static const int coorddimension = Grid::dimensionworld;

    typedef Alberta::Real ctype;

    typedef Alberta::ElementInfo< dimension > ElementInfo;
    typedef FieldVector< ctype, coorddimension > Coordinate;

    typedef Alberta::Twist< dimension, mydimension > Twist;

  private:
    const Grid &grid_;
    const ElementInfo &elementInfo_;
    const int subEntity_;
    const int twist_;

  public:
    GlobalCoordReader ( const GridImp &grid,
                        const ElementInfo &elementInfo,
                        int subEntity )
      : grid_( grid ),
        elementInfo_( elementInfo ),
        subEntity_( subEntity ),
        twist_( Twist::twist( elementInfo_.el(), subEntity ) )
    {}

    void coordinate ( int i, Coordinate &x ) const
    {
      assert( !elementInfo_ == false );
      assert( (i >= 0) && (i <= mydimension) );

      const int ti = Alberta::applyInverseTwist< mydimension >( twist_, i );
      const int k = mapVertices( subEntity_, ti );
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
      assert( false );
      return ctype( 0 );
    }

  private:
    static int mapVertices ( int subEntity, int i )
    {
      return Alberta::MapVertices< dimension, codimension >::apply( subEntity, i );
    }
  };




  // AlbertaGridIntersection::LocalCoordReader
  // -----------------------------------------

  template< class GridImp >
  struct AlbertaGridIntersection< GridImp >::LocalCoordReader
  {
    typedef typename remove_const< GridImp >::type Grid;

    static const int dimension = Grid::dimension;
    static const int codimension = 1;
    static const int mydimension = dimension - codimension;
    static const int coorddimension = dimension;

    typedef Alberta::Real ctype;

    typedef FieldVector< ctype, coorddimension > Coordinate;

    typedef typename Grid::template Codim< 0 >::Geometry ElementGeometry;
    typedef typename Grid::template Codim< 1 >::Geometry FaceGeometry;

  private:
    const ElementGeometry &elementGeometry_;
    const FaceGeometry &faceGeometry_;

  public:
    LocalCoordReader ( const ElementGeometry &elementGeometry,
                       const FaceGeometry &faceGeometry )
      : elementGeometry_( elementGeometry ),
        faceGeometry_( faceGeometry )
    {}

    void coordinate ( int i, Coordinate &x ) const
    {
      x = elementGeometry_.local( faceGeometry_.corner( i ) );
    }

    bool hasDeterminant () const
    {
      return false;
    }

    ctype determinant () const
    {
      return ctype( 0 );
    }
  };

}

#endif
