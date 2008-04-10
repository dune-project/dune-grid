// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <iostream>
#include <fstream>
#include <vector>
#include <config.h>


#include <dune/common/fvector.hh>
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/gridpart.hh>
using namespace Dune;

/*******
   Tests the communication on a parallel grid.
   The tests works as follows for a fixed codim c and
   upwind direction:
   1) In the center of all upwind codim c subentities of the interioir codim 0
     leaf entites a function is stored. Also a flag is set to 1.
     The computation is also performed on the subentities of the physical
     boundary.
   -> For all leaf subentities of codim c the flag should be set to 1
     with the exception of the border subentities on the inflow
     processor boundary and in the ghost elements - on these
     subentities the flag is zero.
   2) Exchange both the data and the flags.
   3) Test if the flag for all leaf subentities of codim c is set to 1.
 ******/
/*****
   The exchange is done using the ExampleDataHandle given below.
   Together with the function value and the flag the coordinates
   of all corners of the subenties are transmitted, giving
   the possibility for additional testing in the scatter/set
   methods.
 ******/
/*******************************************************************/
namespace {
  template <class GridType,int c> struct NextCodim {
    static const bool next = Dune::Capabilities::hasEntity<GridType, c-1>::v;
    static const int calc = ((next) ? c-1 : NextCodim<GridType,c-1>::calc);
  };
  template <class GridType> struct NextCodim<GridType,0> {
    static const int calc = -1;
  };
}
template <class IndexSetImp,
    class GlobalIdSetImp,
    class DataVectorType >
class ExampleDataHandle
  : public CommDataHandleIF <
        ExampleDataHandle<IndexSetImp,GlobalIdSetImp,DataVectorType> ,
        typename DataVectorType :: value_type >
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
  bool contains (int dim, int codim) const
  {
    return (codim==cdim_);
  }

  //! returns true if size per entity of given dim and codim is a constant
  bool fixedsize (int dim, int codim) const
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
    return 2+e.geometry().corners()*e.geometry().dimensionworld;
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
    for (int i=0; i<e.geometry().corners(); ++i)
    {
      for (int j=0; j<e.geometry().dimensionworld; ++j)
        buff.write(e.geometry()[i][j]);
    }
  }

  /*! unpack data from message buffer to user
     n is the number of objects sent by the sender
   */
  template<class MessageBuffer, class EntityType>
  void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
  {
    // as this problem is a fixed size problem we can check the sizes
    assert( n == size(e) );

    // make sure that global id on all processors are the same
    // here compare the id of the entity that sended the data with my id
    typename GlobalIdSetImp :: IdType id;
    //buff.read( id );
    //typename GlobalIdSetImp :: IdType myId = ids_.id( e );
    //std::cout << id << " id | my Id = " << myId << "\n";
    //assert( id == myId );

    // else do normal scatter
    int idx = iset_.index(e);
    DataType x=0.0;
    buff.read(x); // flag

    // for ghost entities x > 0 must be true
    assert( ( e.partitionType() == GhostEntity ) ? (x>=0.0) : 1);

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

    // test if the sending/receving entities are geometrically the same
    for (int i=0; i<e.geometry().corners(); ++i)
    {
      for (int j=0; j<e.geometry().dimensionworld; ++j)
      {
        buff.read(x);
        if (fabs(e.geometry()[i][j]-x)>1e-8)
        {
          std::cerr << "ERROR in scatter: Vertex <" << i << "," << j << ">: "
                    << " this : (" << e.geometry()[i][j] << ")"
                    << " other : (" << x << ")"
                    << std::endl;
        }
      }
    }
  }

};

/*******************************************************************/
/*******************************************************************/
template <class GridPartType, int cdim, class OutputStreamImp>
class CheckCommunication {
  typedef typename GridPartType :: GridType GridType;
  typedef typename GridPartType :: IndexSetType IndexSet;

  enum {dimworld = GridType::dimensionworld,
        dim = GridType::dimension};
  FieldVector<double,dimworld> upwind;

  const IndexSet & set;
  const int level_;

  // the function
  double f(const FieldVector<double,dimworld>& x) {
    FieldVector<double,dimworld> a(1.);
    a[0] = -0.5;
    return a*x+1.5; //+cos(x*x);
  }
  // compute the data on the upwind entities
  template <class VectorType >
  void project (const GridType & grid ,
                int dataSize,
                VectorType & data,
                VectorType & weight, int rank)
  {
    // set initial data
    for(int i=0 ; i<dataSize; ++i)
    {
      data[i] = 0.0;
      weight[i] = -1.0;
    }

    enum { dim = GridType :: dimension };
    typedef typename IndexSet :: template Codim<0> ::
    template Partition<Interior_Partition> :: Iterator IteratorType;
    typedef typename GridPartType::IntersectionIteratorType IntersectionIterator;

    IteratorType endit = set.template end<0,Interior_Partition>   ();
    for(IteratorType it = set.template begin<0,Interior_Partition> ();
        it != endit ; ++it )
    {
      if (cdim==0)
      {
        FieldVector<double,dimworld> mid(0.);
        for(int i=0; i<it->template count<dimworld>(); ++i) {
          mid += it->geometry()[i];
        }
        mid /= double(it->template count<dimworld>());
        int idx = set.index(*it);
        data[idx]  = f(mid);
        weight[idx] = 1.0;
      }
      else
      {
        IntersectionIterator endnit = gridPart_.iend(*it);
        for (IntersectionIterator nit = gridPart_.ibegin(*it); nit != endnit; ++nit)
        {
          const typename IntersectionIterator :: Intersection &intersection = *nit;

          const FieldVector< double,dimworld > normal
            = intersection.integrationOuterNormal( FieldVector< double, dim-1 >( 0 ) );
          double calc = normal*upwind;

          // if testing level, then on non-conform grid also set values on
          // intersections that are not boundary, but has no level
          // neighbor
          bool proceedAnyways = ((level_ < 0) ? false : !intersection.neighbor());
          if( ((calc > -1e-8) || intersection.boundary()) || proceedAnyways )
          {
            const ReferenceElement< double, dimworld > &insideRefElem
              =  ReferenceElements< double, dimworld > :: general( it->type() );

            const int numberInSelf = intersection.numberInSelf();
            for( int i = 0; i < insideRefElem.size( numberInSelf, 1, cdim ); ++i )
            {
              int e = insideRefElem.subEntity( numberInSelf, 1, i, cdim );
              int idx = set.template subIndex<cdim>(*it,e);
              FieldVector<double,dimworld> cmid(0.);
              int c = (it->template entity<cdim>(e))->geometry().corners();
              for (int j=0; j<c; j++)
              {
                cmid += it->template entity<cdim>(e)->geometry()[j];
              }

              cmid /= double(c);
              data[idx] = f(cmid);
              weight[idx] = 1.0;
            }

            // on non-conforming grids the neighbor entities might not
            // be the same as those on *it, therefore set data on neighbor
            // as well
            if( intersection.neighbor() )
            {
              typedef typename GridType :: template Codim<0> :: EntityPointer EntityPointerType;
              typedef typename GridType :: template Codim<0> :: Entity EntityType;
              EntityPointerType ep = intersection.outside();
              const EntityType & neigh = *ep;

              assert( (level_ < 0) ? (neigh.isLeaf()) : 1);
              assert( (level_ < 0) ? 1 : (neigh.level() == level_) );

              const ReferenceElement<double, dimworld > & outsideRefElem =
                ReferenceElements<double, dimworld >::general(ep->type());

              const int numberInNeighbor = intersection.numberInNeighbor();
              for( int i = 0; i < outsideRefElem.size(numberInNeighbor, 1, cdim); ++i )
              {
                int e = outsideRefElem.subEntity( numberInNeighbor, 1, i, cdim );
                int idx = set.template subIndex<cdim>(neigh, e);
                FieldVector<double,dimworld> cmid(0.);
                typedef typename GridType:: template Codim<cdim> :: EntityPointer SubEntityPointer;
                SubEntityPointer subEp = neigh. template entity<cdim>(e);
                int c = subEp->geometry().corners();
                for (int j=0; j<c; j++)
                {
                  cmid += subEp->geometry()[j];
                }
                cmid /= double(c);
                data[idx] = f(cmid);
                weight[idx] = 1.0;
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
  template <class VectorType  >
  double test (const GridType & grid ,
               int dataSize,
               VectorType & data,
               VectorType & weight, int rank,
               bool testweight)
  {
    enum { dim = GridType :: dimension };
    //Variante MIT Geisterzellen
    typedef typename IndexSet :: template Codim<0> :: template Partition<All_Partition> :: Iterator IteratorType;
    double maxerr = 0.;
    IteratorType endit  = set.template end<0,All_Partition>   ();
    for(IteratorType it = set.template begin<0,All_Partition> ();
        it != endit ; ++it )
    {
      FieldVector<double,dimworld> mid(0.);
      for(int i=0; i<it->template count<dim>(); ++i)
      {
        mid += it->geometry()[i];
      }
      mid /= double(it->template count<dim>());
      if (cdim==0)
      {
        int idx = set.index(*it);
        double lerr = fabs(f(mid)-data[idx]);
        maxerr = std::max(maxerr,lerr);
        if (testweight)
        {
          if (weight[idx] < 0)
          {
            sout << "<" << rank << "/test> ERROR in communication test. ";
            sout << "weight:" << weight[idx] << " should be zero!";
            sout << " value is : " << data[idx]
                 << "  index is: ";
            sout << idx << "  lvl:" << it->level() << "  " ;
            sout << "\n";
          }
        }
      }
      else
      {
        for(int i=0; i<it->template count<cdim>(); ++i)
        {
          int idx = set.template subIndex<cdim>(*it,i);
          FieldVector<double,dimworld> cmid(0.);
          int c = (it->template entity<cdim>(i))->geometry().corners();
          for (int j=0; j<c; j++)
          {
            cmid += it->template entity<cdim>(i)->geometry()[j];
          }

          cmid /= double(c);
          double lerr = fabs(f(cmid)-data[idx]);
          if (testweight)
          {
            if (weight[idx] < 0) {
              sout << "<" << rank << "/test> ERROR in communication test. ";
              sout << "weight:" << weight[idx] << " should be zero!";
              sout << " value is : " << data[idx]
                   << "  index is:";
              sout << idx << "  lvl:" << it->level() << "  " ;
              sout << "\n";
              int c = (it->template entity<cdim>(i))->geometry().corners();
              for (int j=0; j<c; j++)
              {
                GeometryType type = it->type();
                const ReferenceElement<double, dimworld > & refElem =
                  ReferenceElements<double, dimworld >::general(type);
                int vx = refElem.subEntity(i,cdim,j,3);
                sout << "index:" << set.template subIndex<3>(*it,vx) << " ";
                sout << it->template entity<cdim>(i)->geometry()[j] << "/" ;
              }
              sout << "\n";
            }
          }
          maxerr = std::max(maxerr,lerr);
        }
      }
    }
    return maxerr;
  }
  // The main ''algorithn''
  bool checkCommunication(const GridType &grid)
  {
    upwind[0] = -0.1113;
    int myrank = grid.comm().rank();

    if (myrank == 0)
      std::cout << "TEST ";
    if( level_ < 0)
      std::cout << "Leaf";
    else
      std::cout << "Level<" << level_ << ">";
    std::cout << " Communication for codim: " << cdim << "\n";

    int dataSize = set.size(cdim);

    typedef std::vector<double> ArrayType;
    ArrayType data(dataSize, 0.0);
    ArrayType weight(dataSize, 0.0);

    project(grid,dataSize,data,weight, myrank);
    double preresult = test(grid,dataSize,data,weight,myrank,false);
    sout << "Test before Communication on <" << myrank << "> "
         << preresult << std::endl;
    // Communicate
    typedef typename GridType :: Traits :: GlobalIdSet GlobalIdSetType;
    ExampleDataHandle<IndexSet,GlobalIdSetType,ArrayType>
    dh (set,grid.globalIdSet(),cdim,data,weight);
    // std::cout << "STARTING COMMUNICATION ... " << std::flush;

    // call communication of grid
    // if level < 0 then check leaf level
    gridPart_.communicate(dh,InteriorBorder_All_Interface,ForwardCommunication);
    // make sure backward communication does the same, this should change nothing
    gridPart_.communicate(dh,InteriorBorder_All_Interface,BackwardCommunication);
    //gridPart_.communicate(dh,All_All_Interface,ForwardCommunication);

    double result = test(grid,dataSize,data,weight,myrank,true);
    sout << "Test after Communication on <" << myrank << "> "  << std::flush
         << result << std::endl;
    return (fabs(result)<1e-8);
  }
  OutputStreamImp& sout;
  const GridPartType & gridPart_;
public:
  // --constructor
  CheckCommunication(const GridPartType &gridPart,
                     OutputStreamImp & out,
                     const int level

                     ) : upwind(-1.)
                         , set(gridPart.indexSet())
                         , level_(level)
                         , sout(out)
                         , gridPart_(gridPart)
  {
    if (!checkCommunication(gridPart.grid())) {
      std::cerr << "ERROR in communication test for codim " << cdim << "!!!"
                << std::endl;
    }
    // for automatic testing of all codims
    CheckCommunication<GridPartType,NextCodim<GridType,cdim>::calc,
        OutputStreamImp> test(gridPart_,sout,level_);
  }
};

template <class GridPartType, class OutputStreamImp>
class CheckCommunication<GridPartType,-1,OutputStreamImp> {
public:
  CheckCommunication(const GridPartType &grid,OutputStreamImp & out, const int level) {}
};

template <class GridType, class OutputStreamImp>
void checkCommunication(GridType &grid, int level , OutputStreamImp & out)
{
  {
    if(level < 0)
    {
      typedef LeafGridPart<GridType> GridPartType;
      GridPartType gridPart(grid);
      CheckCommunication<GridPartType,GridType::dimension,OutputStreamImp>
      test(gridPart,out,level);
    }
    else
    {
      typedef LevelGridPart<GridType> GridPartType;
      GridPartType gridPart(grid,level);
      CheckCommunication<GridPartType,GridType::dimension,OutputStreamImp>
      test(gridPart,out,level);
    }
  }
}
