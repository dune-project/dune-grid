// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_TEST_CHECKINDEXSET_HH
#define DUNE_GRID_TEST_CHECKINDEXSET_HH

#include <algorithm>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/stdstreams.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/exceptions.hh>

/** @file
   @author Robert Kloefkorn
   @brief Provides a check of the grids index set.
 */

template< class Grid >
struct EnableLevelIntersectionIteratorCheck
{
  static const bool v = true;
};

template< class Grid >
struct EnableLevelIntersectionIteratorCheck< const Grid >
{
  static const bool v = EnableLevelIntersectionIteratorCheck< Grid >::v;
};


namespace Dune
{

  // Check whether two FieldVectors are equal in the max norm
  template <typename ctype, int dim>
  bool compareVec(const FieldVector<ctype,dim> & vx1 , const FieldVector<ctype,dim> & vx2 )
  {
    const ctype eps = 1e5 * std::numeric_limits<ctype>::epsilon();
    return (vx1-vx2).infinity_norm() < eps;
  }

  /** \brief Check various features of codim-0 entities (elements)
      \param grid The grid we are testing
      \param en The grid element we are testing
      \param lset The corresponding level set
      \param sout The output stream that status (non-error) messages will go to
      \param setOfVerticesPerSubEntity A map that associates a vector of vertex indices to a subEntity
      \param subEntityPerSetOfVertices The reverse: a map that associates a subEntity to a vector of vertex indices
      \param vertexCoordsMap Map that associates vertex positions to indices (integers)
   */
  template <int codim, class GridType,
      class IndexSetType, class OutputStreamImp,
      class MapType1 , class MapType2 , class MapType3 >
  void checkSubEntity ( const GridType &,
                        const typename GridType::template Codim<0>::Entity &en , const IndexSetType & lset,
                        OutputStreamImp & sout , MapType1 & setOfVerticesPerSubEntity ,
                        MapType2 & subEntityPerSetOfVertices,
                        const MapType3 & vertexCoordsMap )
  {
    enum { dim = GridType::dimension };
    const int dimworld = GridType::dimensionworld;
    typedef typename GridType::ctype coordType;

    const GeometryType type = en.type();
    assert( type == en.geometry().type() );

    if( type.isNone() )
      return;

    auto refElem = referenceElement( en.geometry() );

    // check whether the element has the number of codim-subentities mandated by the reference element
    if(int(en.subEntities(codim)) != refElem.size(0,0,codim))
    {
      std::cerr << "entity index = " << lset.index(en)
                << ", type = " << type
                << std::endl
                << "codim = " << codim
                << std::endl
                << "subEntities(codim) = " << en.subEntities(codim)
                << std::endl
                << "refElem.size(codim) = " << refElem.size(0,0,codim)
                << std::endl;
      DUNE_THROW(GridError,
                 "wrong number of subEntities of codim " << codim);
    }

    for( int subEntity = 0; subEntity < refElem.size( 0, 0, codim ); ++subEntity )
    {

      // Number of vertices of the current subentity
      int numVertices = refElem.size( subEntity, codim, dim );

      // every entity must have at least one vertex
      assert( numVertices > 0 );

      // The local vertex numbers
      std::vector<int> local (numVertices);

      for( int j = 0; j < numVertices; ++j )
        local[ j ] = refElem.subEntity ( subEntity, codim, j, dim );

      sout << numVertices << " vertices on subEntity< codim = " << codim << " >" << std::endl;
      sout << "check subentity [" << local[ 0 ];
      for( int j = 1; j < numVertices; ++j )
        sout << ", " << local[ j ];
      sout << "]" << std::endl;

      // The global vertex numbers
      std::vector<int> global(numVertices);
      for( int j = 0; j < numVertices; ++j )
        global[ j ] = lset.subIndex( en, local[ j ], dim );

      sout << "Found global numbers of entity [ ";
      for( int j = 0; j < numVertices; ++j )
        sout << global[ j ] << " ";
      sout << "]" << std::endl;

      // Check whether we get the same index if we
      // -- ask for the subEntity itself and then its index
      // -- ask for the subIndex directly
      typedef typename GridType::template Codim< codim >::Entity SubE;
      SubE subE = en.template subEntity< codim >( subEntity );

      if( lset.subIndex( en, subEntity, codim ) != lset.index( subE) )
      {
        std::cerr << "Index for subEntity does not match." << std::endl;
        assert( lset.subIndex( en, subEntity, codim ) == lset.index( subE) );
      }

      // Make a unique identifier for the subEntity.  Since indices are unique only per GeometryType,
      // we need a (index,GeometryType)-pair.
      std::pair<int, GeometryType> globalSubEntity( lset.index( subE ), subE.type() );
      assert( globalSubEntity.first >= 0 );
      sout << "local subentity " << subEntity << " consider subentity with global key (" << globalSubEntity.first << ","
           << globalSubEntity.second << ") on en = " << lset.index(en) << std::endl;

      if( subE.type() != subE.geometry().type() )
      {
        std::cerr << "Geometry types for subEntity don't match." << std::endl;
        assert( subE.type() == subE.geometry().type() );
      }

      // assert that all sub entities have the same level
      assert( subE.level() == en.level() );

      // Loop over all vertices
      for( int j = 0; j < numVertices; ++j )
      {

        // get entity pointer to subEntity vertex
        typedef typename GridType::template Codim< dim >::Entity VertexE;
        VertexE vxE = en.template subEntity< dim >( local[ j ] );

        // Find the global coordinate of the vertex by its index
        if(vertexCoordsMap.find(global[j]) != vertexCoordsMap.end())
        {
          // Check whether index and coordinate match
          FieldVector<coordType,dimworld> vxcheck ( vertexCoordsMap.find(global[j])->second );
          FieldVector< coordType, dimworld > vx1 = vxE.geometry().corner( 0 );
          if( ! compareVec( vxcheck, vx1 ) )
          {
            std::cerr << "ERROR map global vertex [" << global[j] << "] vx " << vxcheck << " is not " << vx1 << "\n";
            assert( compareVec( vxcheck, vx1 ) );
          }
        }
        sout << "vx[" << global[j] << "] = "  <<  subE.geometry().corner( j ) << "\n";
      }
      sout << "sort vector of global vertex\n";

      // sort vector of global vertex number for storage in map
      // the smallest entry is the first entry
      std::sort( global.begin(), global.end() );

      // check whether vertex key is already stored in map
      if(subEntityPerSetOfVertices.find(global) == subEntityPerSetOfVertices.end())
      {
        subEntityPerSetOfVertices[global] = globalSubEntity;
      }
      else
      {
        assert( globalSubEntity == subEntityPerSetOfVertices[global] );
      }

      // check whether subEntity is already stored in map
      if(setOfVerticesPerSubEntity.find(globalSubEntity) == setOfVerticesPerSubEntity.end() )
      {
        setOfVerticesPerSubEntity[globalSubEntity] = global;
      }
      else
      {
        std::vector<int> globalcheck = setOfVerticesPerSubEntity[globalSubEntity];
        if(! (global == globalcheck ))
        {
          std::cerr << "For subEntity key (" << globalSubEntity.first << "," << globalSubEntity.second << ") \n";
          std::cerr << "Got ";
          for(int j=0 ; j<numVertices; j++ )
          {
            std::cerr << global[j] << " ";
          }
          std::cerr << "\n";
          std::cerr << "Found ";
          for(int j=0 ; j<numVertices; j++ )
          {
            std::cerr << globalcheck [j] << " ";
          }
          std::cerr << "\n";
          DUNE_THROW(Dune::GridError, "global != globalcheck");
        }
      }

    }   // end check sub entities
    sout << "end check sub entities\n";
  }


  // check some functionality of grid
  template< int codim, class Grid, class GridView, class OutputStream >
  void checkIndexSetForCodim ( const Grid &grid, const GridView &view,
                               OutputStream &sout, bool levelIndex )
  {
    enum { dim = Grid :: dimension };
    enum { dimworld = Grid :: dimensionworld };

    typedef typename Grid :: ctype coordType;

    //typedef typename GridView :: template Codim< 0 > :: Entity EntityCodim0Type;
    typedef typename GridView :: IndexSet IndexSetType;
    typedef typename GridView :: template Codim< codim > :: Iterator IteratorType;

    const IndexSetType &lset = view.indexSet();

    sout <<"\n\nStart consistency check of index set \n\n";

    // ////////////////////////////////////////////////////////////////
    //   Check whether geomTypes() returns correct result
    // ////////////////////////////////////////////////////////////////

    std :: set< GeometryType > geometryTypes;

    const IteratorType endit = view.template end< codim >();
    IteratorType it = view.template begin< codim >();

    if (it == endit) return;

    for (; it!=endit; ++it) {
      // while we're here: check whether the GeometryTypes returned by the entity
      // and the Geometry match
      assert(it->type()==it->type());
      geometryTypes.insert(it->type());
    }

    const typename IndexSetType::Types types = lset.types( codim );
    bool geomTypesError = false;

    // Check whether all entries in the official geometry types list are contained in our self-computed one
    for( auto it = types.begin(); it != types.end(); ++it )
      if( geometryTypes.find( *it ) == geometryTypes.end() )
        geomTypesError = true;


    // And vice versa
    for( std::set<GeometryType>::iterator it = geometryTypes.begin(); it != geometryTypes.end(); ++it )
    {
      if( std::find( types.begin(), types.end(), *it ) == types.end() )
        geomTypesError = true;
    }



    if (geomTypesError) {

      std::cerr << "There is a mismatch in the list of geometry types of codim " << codim << "." << std::endl;
      std::cerr << "Geometry types present in the grid are:" << std::endl;
      for (std::set<GeometryType>::iterator it = geometryTypes.begin(); it!=geometryTypes.end(); ++it)
        std::cerr << "  " << *it << std::endl;

      std::cerr << std::endl << "but the method types() returned:" << std::endl;
      for( auto it = types.begin(); it != types.end(); ++it )
        std::cerr << "  " << *it << std::endl;

      DUNE_THROW(GridError, "!");
    }

    //*****************************************************************
    // check size of index set
    int gridsize = 0;
    {
      typedef typename GridView :: template Codim< codim > :: Iterator IteratorType;

      int count = 0;
      const IteratorType endit = view.template end< codim >();
      for( IteratorType it = view.template begin< codim >(); it != endit; ++it )
        ++count;

      int lsetsize = lset.size(codim);
      if( count != lsetsize)
      {
        derr << "WARNING: walk = "<< count << " entities | set = "
             << lsetsize << " for codim " << codim << std::endl;
      }
      gridsize = count;
      // lsetsize should be at least the size of iterated entities
      assert( count <= gridsize );
    }

    {
      typedef typename GridView :: template Codim< 0 > :: Iterator Iterator;
      typedef typename Grid :: Traits :: LocalIdSet LocalIdSetType;
      typedef typename LocalIdSetType :: IdType IdType;

      std::set< IdType > entityfound;

      const Iterator endit = view.template end< 0 >();
      Iterator it = view.template begin< 0 >();
      if( it == endit )
        return;

      const LocalIdSetType &localIdSet = grid.localIdSet();
      for( ; it != endit; ++it )
      {
        const typename std::iterator_traits<Iterator>::value_type &entity = *it;
        if( !lset.contains( entity ) )
        {
          std::cerr << "Error: IndexSet does not contain all entities." << std::endl;
          assert( false );
        }
        const int subcount = entity.subEntities(codim);
        for( int i = 0; i < subcount; ++i )
        {
          const IdType id = localIdSet.id( entity.template subEntity< codim >( i ) );
          entityfound.insert( id );
        }
      }

      if( (size_t)gridsize != entityfound.size() )
      {
        std::cerr << "Warning: gridsize = " << gridsize << " entities"
                  << ", set of entities = " << entityfound.size()
                  << " [codim " << codim << "]" << std::endl;
      }

      // gridsize should be at least the size of found entities
      //assert( gridsize <= (int) entityfound.size() );
    }

    //******************************************************************

    typedef std::pair < int , GeometryType > SubEntityKeyType;
    //typedef std::map < int , std::pair<int,int> > subEntitymapType;
    std::map < SubEntityKeyType , std::vector<int> > setOfVerticesPerSubEntity;
    std::map < std::vector<int> , SubEntityKeyType > subEntityPerSetOfVertices;

    std::map < int , FieldVector<coordType,dimworld> > vertexCoordsMap;
    // setup vertex map , store vertex coords for vertex number
    {
      unsigned int count = 0;
      typedef typename GridView :: template Codim< dim > :: Iterator VxIterator;
      const VxIterator end = view.template end< dim >();
      for( VxIterator it = view.template begin< dim >(); it != end; ++it )
      {
        ++count;
        // get coordinates of vertex
        FieldVector< coordType, dimworld > vx ( it->geometry().corner( 0 ) );

        // get index of vertex
        sout << "Vertex " << vx << "\n";
        assert( lset.contains ( *it ) );
        int idx = lset.index( *it );

        sout << "Vertex " << idx << " = [" << vx << "]\n";

        // if vertex not in map insert it
        if( vertexCoordsMap.find(idx) == vertexCoordsMap.end())
          vertexCoordsMap[idx] = vx;
      }
      sout << "Found " << vertexCoordsMap.size() << " vertices for that index set!\n\n";

      // check whether size of map equals all found vertices
      assert( vertexCoordsMap.size() == count );

      // check whether size of vertices of set equals all found vertices
      sout << "Checking size of vertices "
           << count
           << " equals all found vertices "
           << (unsigned int)lset.size(Dune::GeometryTypes::vertex)
           << "\n";
      // assertion goes wrong:
      // - for parallel grid since no iteration over ghost subentities
      // - if there are vertices with geometry type 'none'
      //assert( count == (unsigned int)lset.size(Dune::GeometryTypes::vertex) );
    }

    {
      typedef typename GridView :: template Codim< 0 > :: Iterator Iterator;
      // choose the right reference element
    #ifndef NDEBUG
      const Iterator refend = view.template end< 0 >();
    #endif
      Iterator refit = view.template begin< 0 >();
      assert( refit != refend );

      auto refElem = referenceElement( refit->geometry() );

      // print dune reference element
      sout << "Dune reference element provides: \n";
      for(int i = 0; i < refElem.size( codim ); ++i )
      {
        sout << i << " subEntity [";
        int s = refElem.size(i,codim,dim);
        for(int j=0; j<s; j++)
        {
          sout << refElem.subEntity(i , codim , j , dim );
          if(j != s-1) sout << ",";
        }
        sout << "]\n";
      }
    }

    {
      typedef typename GridView :: template Codim< 0 > :: Iterator Iterator;
      const Iterator endit = view.template end< 0 >();
      for( Iterator it = view.template begin< 0 >(); it != endit; ++it )
      {
        // if (it->partitionType()==4) continue;
        sout << "****************************************\n";
        sout << "Element = " << lset.index(*it) << " on level " << it->level () << "\n";
        sout << "Vertices      = [";
        int svx = it->subEntities(dim);

        // print all vertex numbers
        for( int i = 0; i < svx; ++i )
        {
          const typename IndexSetType::IndexType idx = lset.subIndex( *it, i, dim );
          sout << idx << (i < svx-1 ? ", " : "]\n");
        }

        // print all vertex coordinates
        sout << "Vertex Coords = [";
        for( int i = 0; i < svx; ++i )
        {
          // get entity pointer of sub entity codim=dim (Vertex)
          typedef typename GridView::template Codim< dim >::Entity VertexE;
          VertexE vxE = it->template subEntity< dim >( i );

          // get coordinates of entity pointer
          FieldVector< coordType, dimworld > vx( vxE.geometry().corner( 0 ) );

          // output vertex coordinates
          sout << vx << (i < svx-1 ? ", " : "]\n");

          const typename IndexSetType::IndexType vxidx = lset.subIndex( *it, i, dim );

          // the subIndex and the index for subEntity must be the same
          if( vxidx != lset.index( vxE ) )
          {
            std::cerr << "Error: index( *subEntity< dim >( i ) ) != subIndex( entity, i, dim )" << std::endl;
            assert( vxidx == lset.index( vxE ) );
          }

          // check whether the coordinates are the same
          assert(vertexCoordsMap.find(vxidx)!=vertexCoordsMap.end());
          FieldVector<coordType,dimworld> vxcheck ( vertexCoordsMap[vxidx] );
          if( !compareVec( vxcheck, vx ) )
          {
            std::cerr << "Error: Inconsistent map of global vertex " << vxidx
                      << ": " << vxcheck << " != " << vx
                      << "  (type = " << it->partitionType() << ")" << std::endl;
            assert( compareVec( vxcheck, vx ) );
          }
        }

        ////////////////////////////////////////////////////////////
        // check sub entities
        ////////////////////////////////////////////////////////////
        checkSubEntity< codim >( grid, *it, lset, sout,
                                 setOfVerticesPerSubEntity, subEntityPerSetOfVertices, vertexCoordsMap );

        // check neighbors
        if( codim == 1 )
        {
          typedef typename GridView :: IntersectionIterator IntersectionIterator;

          if( !levelIndex || EnableLevelIntersectionIteratorCheck< typename GridView::Grid >::v )
          {
            const IntersectionIterator endnit = view.iend( *it );
            for( IntersectionIterator nit = view.ibegin( *it ); nit != endnit; ++nit )
            {
              if( !nit->neighbor() )
                continue;

              checkSubEntity< codim >( grid, nit->outside(), lset, sout,
                                       setOfVerticesPerSubEntity, subEntityPerSetOfVertices, vertexCoordsMap );
            }
          }
          else
          {
            static int called = 0;
            if( called++ == 0 )
              std::cerr << "Warning: Skipping index test using LevelIntersectionIterator." << std::endl;
          }
        }
      }
    }
  }


  template< class Grid, class GridView, class OutputStream, int codim, bool hasCodim >
  struct CheckIndexSet
  {
    static void checkIndexSet ( const Grid &grid, const GridView &view,
                                OutputStream &sout, bool levelIndex )
    {
      checkIndexSetForCodim< codim >( grid, view, sout, levelIndex );
      typedef Dune :: Capabilities :: hasEntity< Grid, codim-1 > hasNextCodim;
      CheckIndexSet< Grid, GridView, OutputStream, codim-1, hasNextCodim :: v >
      :: checkIndexSet( grid, view, sout, levelIndex );
    }
  };

  template< class Grid, class GridView, class OutputStream, int codim >
  struct CheckIndexSet< Grid, GridView, OutputStream, codim, false >
  {
    static void checkIndexSet ( const Grid &grid, const GridView &view,
                                OutputStream &sout, bool levelIndex )
    {
      derr << "WARNING: Entities for codim " << codim
           << " are not being tested!" << std::endl;
      typedef Dune :: Capabilities :: hasEntity< Grid, codim-1 > hasNextCodim;
      CheckIndexSet< Grid, GridView, OutputStream, codim-1, hasNextCodim :: v >
      :: checkIndexSet( grid, view, sout, levelIndex );
    }
  };

  template< class Grid, class GridView, class OutputStream >
  struct CheckIndexSet< Grid, GridView, OutputStream, 0, true >
  {
    static void checkIndexSet ( const Grid &grid, const GridView &view,
                                OutputStream &sout, bool levelIndex )
    {
      checkIndexSetForCodim< 0 >( grid, view, sout, levelIndex );
    }
  };

  template< class Grid, class GridView, class OutputStream >
  void checkIndexSet ( const Grid &grid, const GridView &view,
                       OutputStream &sout,  bool levelIndex = false )
  {
    CheckIndexSet< Grid, GridView, OutputStream, Grid :: dimension, true >
    :: checkIndexSet ( grid, view, sout, levelIndex );
  }

} // end namespace Dune

#endif // DUNE_GRID_TEST_CHECKINDEXSET_HH
