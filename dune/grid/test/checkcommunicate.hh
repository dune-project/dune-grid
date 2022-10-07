// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_TEST_CHECKCOMMUNICATE_HH
#define DUNE_GRID_TEST_CHECKCOMMUNICATE_HH

#include <iostream>
#include <fstream>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/gridview.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/datahandleif.hh>

/*  Communication Test for Parallel Grids
 *  -------------------------------------
 *
 *  For a fixed codimension c and a fixed upwind direction u, the test works
 *  as follows:
 *    1) In the center of all upwind codim c subentities of the interioir codim 0
 *       leaf entites a function is stored. Also a flag is set to 1.
 *       The computation is also performed on the subentities of the physical
 *       boundary.
 *    -> For all leaf subentities of codim c the flag should be set to 1
 *       with the exception of the border subentities on the inflow
 *       processor boundary and in the ghost elements - on these
 *       subentities the flag is zero.
 *    2) Exchange both the data and the flags.
 *    3) Test if the flag for all leaf subentities of codim c is set to 1.
 *
 *  Note: This test requires the normals on both sides of an intersection to
 *        sum up to zero, i.e., there is exactly one tangent plane to the grid
 *        in every point of the intersection (actually, the barycenter would be
 *        sufficient).
 */


/*****
   The exchange is done using the ExampleDataHandle given below.
   Together with the function value and the flag the coordinates
   of all corners of the subenties are transmitted, giving
   the possibility for additional testing in the scatter/set
   methods.
 ******/
/*******************************************************************/
namespace
{

  template< class Grid, int codim = Grid::dimension+1 >
  struct NextCodim
  {
    static const bool canCommunicate = Dune::Capabilities::canCommunicate< Grid, codim-1 >::v;

    static const int v = (canCommunicate ? codim-1 : NextCodim< Grid, codim-1 >::v);
  };

  template< class Grid >
  struct NextCodim< Grid, 0 >
  {
    static const int v = -1;
  };

}


template <class IndexSetImp,
    class GlobalIdSetImp,
    class DataVectorType >
class ExampleDataHandle
  : public Dune::CommDataHandleIF< ExampleDataHandle< IndexSetImp, GlobalIdSetImp, DataVectorType >, typename DataVectorType::value_type >
{
  const IndexSetImp & iset_;
  const GlobalIdSetImp & ids_;
  int cdim_;
  DataVectorType & data1_;
  DataVectorType & data2_;
public:
  typedef typename DataVectorType :: value_type DataType;
  ExampleDataHandle(const IndexSetImp & iset,
                    const GlobalIdSetImp & ids,
                    int cdim,
                    DataVectorType & d1, DataVectorType & d2) :
    iset_(iset), ids_(ids) , cdim_(cdim), data1_(d1) , data2_(d2)
  {}

  //! returns true if data for this codim should be communicated
  bool contains ([[maybe_unused]] int dim, int codim) const
  {
    return (codim==cdim_);
  }

  //! returns true if size per entity of given dim and codim is a constant
  bool fixedSize ([[maybe_unused]] int dim, [[maybe_unused]] int codim) const
  {
    // this problem is a fixed size problem,
    // but to simulate also non-fixed size problems
    // we set this to false, should work anyway
    return false;
  }

  /*! how many objects of type DataType have to be sent for a given entity
     Note: Only the sender side needs to know this size.
   */
  template<class EntityType>
  size_t size (EntityType& e) const
  {
    // flag+data+coordinates
    typedef typename EntityType::Geometry Geometry;
    return 2+e.geometry().corners() * Geometry::coorddimension;
  }

  //! pack data from user to message buffer
  template<class MessageBuffer, class EntityType>
  void gather (MessageBuffer& buff, const EntityType& e) const
  {
    int idx = iset_.index(e);

    //typename GlobalIdSetImp :: IdType id = ids_.id( e );
    //buff.write( id );
    buff.write(data2_[idx]);   // flag
    buff.write(data1_[idx]);   // data

    // all corner coordinates
    typedef typename EntityType::Geometry Geometry;
    const Geometry &geometry = e.geometry();
    for( int i = 0; i < geometry.corners(); ++i )
    {
      typedef Dune::FieldVector< typename Geometry::ctype, Geometry::coorddimension > Vector;
      const Vector corner = geometry.corner( i );
      for( int j = 0; j < Geometry::coorddimension; ++j )
        buff.write( corner[ j ] );
    }
  }

  /*! unpack data from message buffer to user
     n is the number of objects sent by the sender
   */
  template<class MessageBuffer, class EntityType>
  void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
  {
    using std::sqrt;
    using Geometry = typename EntityType::Geometry;
    using ctype = typename Geometry::ctype;

    // define a tolerance for floating-point checks
    const ctype tolerance = sqrt(std::numeric_limits< ctype >::epsilon());

    // as this problem is a fixed size problem we can check the sizes
    assert( n == size(e) );

    // make sure that global id on all processors are the same
    // here compare the id of the entity that sended the data with my id
    //typename GlobalIdSetImp :: IdType id;
    //buff.read( id );
    //typename GlobalIdSetImp :: IdType myId = ids_.id( e );
    //std::cout << id << " id | my Id = " << myId << "\n";
    //assert( id == myId );

    // else do normal scatter
    int idx = iset_.index(e);
    DataType x=0.0;
    buff.read(x); // flag

    // for ghost entities x > 0 must be true
    assert( ( e.partitionType() == Dune::GhostEntity ) ? (x>=0.0) : 1);

    if (x>=0)
    { // only overwrite existing data if flag = 1, i.e.,
      // the sending processor acctually computed the value
      data2_[idx] = x;
      x=0.;
      buff.read(x);  // correct function value
      data1_[idx] = x;
    }
    else
    {
      x=0.;
      buff.read(x);
    }

    // test if the sending/receiving entities are geometrically the same
    const Geometry &geometry = e.geometry();
    for( int i = 0; i < geometry.corners(); ++i )
    {
      typedef Dune::FieldVector< ctype, Geometry::coorddimension > Vector;
      const Vector corner = geometry.corner( i );
      for( int j = 0; j < Geometry::coorddimension; ++j )
      {
        buff.read(x);
        if( std::abs( corner[ j ] - x ) > tolerance )
        {
          std::cerr << "ERROR in scatter: Vertex <" << i << "," << j << ">: "
                    << " this : (" << corner[ j ] << ")"
                    << " other : (" << x << ")"
                    << std::endl;
        }
      }
    }
  }

};


/*******************************************************************/
/*******************************************************************/
template< class GridView, int cdim, class OutputStream >
class CheckCommunication
{
  typedef typename GridView :: Grid Grid;
  typedef typename GridView :: IndexSet IndexSet;
  typedef typename GridView :: IntersectionIterator IntersectionIterator;

  typedef typename IntersectionIterator :: Intersection Intersection;

  typedef typename GridView :: template Codim< 0 > :: Entity Entity;
  typedef typename GridView :: template Codim< 0 > :: Iterator Iterator;

  typedef typename GridView :: template Codim< cdim > :: Entity SubEntity;

  constexpr static int dimworld = Grid :: dimensionworld;
  constexpr static int dim = Grid :: dimension;

  typedef typename Grid :: ctype ctype;

  typedef Dune::FieldVector< ctype, dimworld > CoordinateVector;
  typedef std :: vector< ctype > ArrayType;

  CoordinateVector upwind_;
  OutputStream &sout_;

  const GridView &gridView_;
  const IndexSet &indexSet_;
  const int level_;

  // the function
  ctype f ( const CoordinateVector &x )
  {
    CoordinateVector a( 1.0 );
    a[0] = -0.5;
    return a*x+1.5; //+cos(x*x);
  }

  // compute the data on the upwind entities
  void project ( int dataSize, ArrayType &data, ArrayType &weight, int /* rank */ )
  {
    // set initial data
    for(int i=0 ; i<dataSize; ++i)
    {
      data[i] = 0.0;
      weight[i] = -1.0;
    }

    Iterator end = gridView_.template end< 0 >();
    for( Iterator it = gridView_.template begin< 0 >(); it != end ; ++it )
    {
      const Entity &entity = *it;

      if( cdim == 0 )
      {
        CoordinateVector mid( 0.0 );
        const int numVertices = entity.subEntities(dim);
        for( int i = 0; i < numVertices; ++i )
          mid += entity.geometry().corner( i );
        mid /= ctype( numVertices );

        int index = indexSet_.index( entity );
        data[ index ]  = f( mid );
        weight[ index ] = 1.0;
      }
      else
      {
        const IntersectionIterator nend = gridView_.iend( entity );
        for( IntersectionIterator nit = gridView_.ibegin( entity ); nit != nend; ++nit )
        {
          const Intersection &intersection = *nit;

          const typename Intersection::LocalGeometry &geoInSelf = intersection.geometryInInside();

          auto faceRefElement = referenceElement( geoInSelf );
          const Dune::FieldVector< ctype, dim-1 > &bary = faceRefElement.position( 0, 0 );

          const CoordinateVector normal = intersection.integrationOuterNormal( bary );
          ctype calc = normal * upwind_;

          // if testing level, then on non-conform grid also set values on
          // intersections that are not boundary, but has no level
          // neighbor
          const bool proceedAnyway = (level_ < 0 ? false : !intersection.neighbor());
          if( (calc > -1e-8) || intersection.boundary() || proceedAnyway )
          {
            auto insideRefElem = referenceElement( entity.geometry() );

            const int indexInInside = intersection.indexInInside();
            for( int i = 0; i < insideRefElem.size( indexInInside, 1, cdim ); ++i )
            {
              const int e = insideRefElem.subEntity( indexInInside, 1, i, cdim );
              const int idx = indexSet_.subIndex( entity, e, cdim );
              CoordinateVector cmid( 0.0 );
              SubEntity subE = entity.template subEntity< cdim >( e );
              const int c = subE.geometry().corners();
              for( int j = 0; j < c; ++j )
                cmid += subE.geometry().corner( j );
              cmid /= ctype( c );

              data[ idx ] = f( cmid );
              weight[ idx ] = 1.0;
            }

            // on non-conforming grids the neighbor entities might not
            // be the same as those on *it, therefore set data on neighbor
            // as well
            if( intersection.neighbor() )
            {
              Entity neigh = intersection.outside();

              assert( (level_ < 0) ? (neigh.isLeaf()) : 1);
              assert( (level_ < 0) ? 1 : (neigh.level() == level_) );

              auto outsideRefElem = referenceElement( neigh.geometry() );

              const int indexInOutside = intersection.indexInOutside();
              for( int i = 0; i < outsideRefElem.size( indexInOutside, 1, cdim ); ++i )
              {
                const int e = outsideRefElem.subEntity( indexInOutside, 1, i, cdim );
                const int idx = indexSet_.subIndex( neigh, e, cdim );
                CoordinateVector cmid( 0.0 );
                SubEntity subE = neigh.template subEntity< cdim >( e );
                const int c = subE.geometry().corners();
                for( int j = 0; j < c; ++j )
                  cmid += subE.geometry().corner( j );
                cmid /= ctype( c );

                data[ idx ] = f( cmid );
                weight[ idx ] = 1.0;
              }
            }
          }
        }
      }
    }
  }

  // test if all flags are 1 and return the
  // difference in the function values.
  // if testweight is true an error is printed for each
  // flag not equal to 1
  ctype test ( int /* dataSize */, ArrayType &data, ArrayType &weight, bool testweight )
  {
    const int rank = gridView_.comm().rank();
    const int size = gridView_.comm().size();

    // check whether there is an overlap or ghost region of
    // cells for the current grid view.
    if (size > 1 && cdim == 0 &&
        gridView_.overlapSize(0) == 0 &&
        gridView_.ghostSize(0) == 0)
    {
      sout_ << "<" << rank << "/test> Error in communication test.";
      sout_ << " overlapSize+ghostSize:" << gridView_.overlapSize(0) + gridView_.ghostSize(0) << " (should not be 0)";
      sout_ << " communicator size is:" << size;
      sout_ << std :: endl;
    }

    //Variante MIT Geisterzellen
    //typedef typename IndexSet :: template Codim<0> :: template Partition<All_Partition> :: Iterator IteratorType;

    ctype maxerr = 0.0;
    Iterator end = gridView_.template end< 0 >();
    for( Iterator it = gridView_.template begin< 0 >(); it != end ; ++it )
    {
      const Entity &entity = *it;

      CoordinateVector mid( 0.0 );
      const int numVertices = entity.subEntities(dim);
      for( int i = 0; i < numVertices; ++i )
        mid += entity.geometry().corner( i );
      mid /= ctype(numVertices);

      if( cdim == 0 )
      {
        int index = indexSet_.index( entity );
        ctype lerr = std::abs( f( mid ) - data[ index ] );
        maxerr = std :: max( maxerr, lerr );
        if( testweight && (weight[ index ] < 0) )
        {
          sout_ << "<" << rank << "/test> Error in communication test.";
          sout_ << " weight:" << weight[ index ] << " (should be 0)";
          sout_ << " value is : " << data[ index ];
          sout_ << " index is: " << index;
          sout_ << " level:" << entity.level();
          sout_ << std :: endl;
        }
      }
      else
      {
        const int numSubEntities = entity.subEntities(cdim);
        for( int i=0; i < numSubEntities; ++i )
        {
          SubEntity subE = entity.template subEntity< cdim >( i );

          const int index = indexSet_.index( subE );
          CoordinateVector cmid( 0.0 );

          const int numVertices = subE.geometry().corners();
          for( int j = 0; j< numVertices; ++j )
            cmid += subE.geometry().corner( j );
          cmid /= ctype( numVertices );

          ctype lerr = std::abs( f( cmid ) - data[ index ] );
          maxerr = std::max( maxerr, lerr );
          if( testweight && (weight[ index ] < 0) )
          {
            sout_ << "<" << rank << "/test> Error in communication test.";
            sout_ << " weight:" << weight[ index ] << " should be zero!";
            sout_ << " value is : " << data[ index ];
            sout_ << " index is:" << index;
            sout_ << " level: " << entity.level() ;
            sout_ << std :: endl;

            for( int j = 0; j < numVertices; )
            {
              auto refElem = referenceElement ( entity.geometry() );
              const int vx = refElem.subEntity( i, cdim, j, dim );

              sout_ << "index: " << indexSet_.subIndex( entity, vx, dim )
                    << " " << subE.geometry().corner( j );
              (++j < numVertices ? sout_ <<  "/" : sout_ << std :: endl);
            }
          }
        }
      }
    }
    return maxerr;
  }

  // The main ''algorithm''
  bool checkCommunication ()
  {
    using std::sqrt;
    // define a tolerance for floating-point checks
    const ctype tolerance = sqrt(std::numeric_limits< ctype >::epsilon());

    upwind_[ 0 ] = -0.1113;
    int myrank = gridView_.comm().rank();

    if( myrank == 0 )
    {
      std::cout << "TEST ";
      (level_ < 0 ? std :: cout << "Leaf" : std :: cout << "Level<" << level_ << ">");
      std :: cout << " communication for codim " << cdim << std :: endl;
    }

    const int dataSize = indexSet_.size( cdim );
    ArrayType data(dataSize, 0.0);
    ArrayType weight(dataSize, 0.0);
    project( dataSize, data, weight, myrank );

    ctype preresult = test( dataSize, data, weight, false );
    sout_ << "Test before Communication on <" << myrank << "> " << preresult << std :: endl;
    // Communicate
    typedef typename Grid :: Traits :: GlobalIdSet GlobalIdSet;
    ExampleDataHandle< IndexSet, GlobalIdSet, ArrayType >
    dh( indexSet_, gridView_.grid().globalIdSet(), cdim, data, weight );

    // call communication of grid
    try
    {
      // call forward and backward communication
      auto obj1 = gridView_.communicate( dh, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication );
      if( ! obj1.ready() )
        obj1.wait();

      // make sure backward communication does the same, this should change nothing
      auto obj2 = gridView_.communicate( dh, Dune::InteriorBorder_All_Interface, Dune::BackwardCommunication );
    }
    catch( const Dune::NotImplemented &exception )
    {
      if( myrank == 0 )
      {
        sout_ << "Error: Communication for codimension " << cdim
              << " not implemented." << std :: endl;
        sout_ << "       (" << exception << ")" << std :: endl;
      }
      return false;
    }

    ctype result = test( dataSize, data, weight, true );
    sout_ << "Test after Communication on <" << myrank << "> " << result << std :: endl;
    return (std::abs(result) < tolerance);
  }

public:
  CheckCommunication ( const GridView &gridView, OutputStream &sout, int level )
    : upwind_( -1.0 ),
      sout_( sout ),
      gridView_( gridView ),
      indexSet_( gridView_.indexSet() ),
      level_( level )
  {
    // if no overlap and ghost is available we skip the check
    const bool skipCheck = ( cdim == 0 ) ? (gridView_.overlapSize(0) == 0 && gridView_.ghostSize(0) == 0) : false ;

    if( skipCheck )
    {
      std :: cerr << "Codim " << cdim << ": Test skiped because of empty set of overlap and ghosts !" << std :: endl;
    }
    else if ( ! checkCommunication() )
    {
      std :: cerr << "Error in communication test for codim "
                  << cdim << "!" << std :: endl;
    }

    // for automatic testing of all codims
    CheckCommunication< GridView, NextCodim< Grid, cdim >::v, OutputStream >
    test( gridView_, sout_, level_ );
  }
};


template< class GridView, class OutputStream >
class CheckCommunication< GridView, -1, OutputStream >
{
public:
  CheckCommunication ( const GridView &, OutputStream &, int)
  {}
};


template< class Grid, class OutputStream >
void checkCommunication( const Grid &grid, int level, OutputStream &sout )
{
  if( level < 0 )
  {
    typedef typename Grid::LeafGridView GridView;
    GridView gridView = grid.leafGridView();
    CheckCommunication< GridView, NextCodim< Grid >::v, OutputStream >
    test( gridView, sout, level );
  }
  else
  {
    typedef typename Grid::LevelGridView GridView;
    GridView gridView = grid.levelGridView( level );
    CheckCommunication< GridView, NextCodim< Grid >::v, OutputStream >
    test( gridView, sout, level );
  }
}

#endif // DUNE_GRID_TEST_CHECKCOMMUNICATE_HH
