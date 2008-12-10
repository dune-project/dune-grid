// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id: test-ug.cc 4424 2008-09-29 07:46:41Z sander $

#include <config.h>

#include <iostream>

#include <dune/grid/uggrid.hh>
#include <dune/grid/common/mcmgmapper.hh>

// Test parallel interface if a parallel UG is used
#ifdef ModelP
#include <mpi.h>
#endif

using namespace Dune;

//! Parameter for mapper class
template<int dim>
struct P0Layout
{
  bool contains (Dune::GeometryType gt)
  {
    return gt.dim()==dim;
  }
};

// A DataHandle class to exchange entries of a vector
template<class M, class V> // mapper type and vector type
class VectorExchange
  : public Dune::CommDataHandleIF<VectorExchange<M,V>,
        typename V::value_type>
{
public:
  //! export type of data for message buffer
  typedef typename V::value_type DataType;

  //! returns true if data for this codim should be communicated
  bool contains (int dim, int codim) const
  {
    return (codim==0);
  }

  //! returns true if size per entity of given dim and codim is a constant
  bool fixedsize (int dim, int codim) const
  {
    return true;
  }

  /*! how many objects of type DataType have to be sent for a given entity

     Note: Only the sender side needs to know this size.
   */
  template<class EntityType>
  size_t size (EntityType& e) const
  {
    return 1;
  }

  //! pack data from user to message buffer
  template<class MessageBuffer, class EntityType>
  void gather (MessageBuffer& buff, const EntityType& e) const
  {
    buff.write(c[mapper.map(e)]);
  }

  /*! unpack data from message buffer to user

     n is the number of objects sent by the sender
   */
  template<class MessageBuffer, class EntityType>
  void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
  {
    DataType x;
    buff.read(x);
    c[mapper.map(e)]=x;
  }

  //! constructor
  VectorExchange (const M& mapper_, V& c_)
    : mapper(mapper_), c(c_)
  {}

private:
  const M& mapper;
  V& c;
};


int main (int argc , char **argv) try
{

#ifdef ModelP
  // initialize MPI
  MPI_Init(&argc,&argv);
#endif


  int rank;
  int size;

  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  std::cout << "This is process " << rank << " of " << size << "." << std::endl;

  // //////////////////////////////////////////////////////////
  //   Make a uniform grid on rank 0 for testing
  // //////////////////////////////////////////////////////////
  int n = 11;
  typedef UGGrid<2> GridType;
  GridType grid;
  GridFactory<GridType> factory(&grid);

  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      FieldVector<double,2> pos;
      pos[0] = double(i)/(n-1);
      pos[1] = double(j)/(n-1);
      factory.insertVertex(pos);
    }
  }

  for (int i=0; i<n-1; i++) {
    for (int j=0; j<n-1; j++) {
      std::vector<unsigned int> v(4);
      v[0] = i*n + j;
      v[1] = i*n + j+1;
      v[2] = (i+1)*n + j;
      v[3] = (i+1)*n + j+1;
      factory.insertElement(GeometryType(GeometryType::cube,2), v);
    }
  }

  //std::auto_ptr<Dune::UGGrid<2> > grid(factory.createGrid());
  factory.createGrid();

  grid.loadBalance(0,     // strategy
                   1,      // minlevel
                   2,      // depth
                   32,      // maxlevel
                   1       // minelement
                   );

  std::cout << "Process " << rank << " has " << grid.size(0) << " elements." << std::endl;

  // make sure each process has some elements
  assert(grid.size(0) > 0);

  // ////////////////////////////////////////////////////
  //  Test element communication
  // ////////////////////////////////////////////////////

  GridType::LevelGridView gridView = grid.levelView(0);

  // make a mapper for codim 0 entities in the level grid
  typedef LevelMultipleCodimMultipleGeomTypeMapper<GridType,P0Layout> MapperType;
  MapperType mapper(grid, 0);

  typedef std::vector<double> VectorType;
  VectorType data(grid.size(0));
  data.assign(data.size(), rank);

  VectorExchange<MapperType,VectorType> datahandle(mapper,data);
  grid.communicate<VectorExchange<MapperType,VectorType> >(datahandle,Dune::InteriorBorder_All_Interface,
                                                           Dune::ForwardCommunication);

#ifdef ModelP
  // Terminate MPI
  MPI_Finalize();
#endif

  return 0;
}
catch (Dune::Exception& e) {
  std::cerr << e << std::endl;
  return 1;
} catch (...) {
  std::cerr << "Generic exception!" << std::endl;
  return 2;
}
