// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>
#include <dune/grid/test/gridcheck.cc>
#include <dune/grid/common/gridpart.hh>
#include "../dgfgridtype.hh"
#include <dune/grid/common/gridpart.hh>
using namespace Dune;
#if HAVE_GRAPE
#include <dune/grid/io/visual/grapegriddisplay.hh>
#include <dune/grid/io/visual/grapedatadisplay.hh>
template <class GridPart>
void test(GridPart& part,std::vector<double>& dat,int nofParams,int cdim) {
  gridcheck(part.grid());
  return;
  if (nofParams>0) {
    GrapeDataDisplay<typename GridPart::GridType> disp(part);
    if (cdim == 0)
      disp.displayVector("el. Paramters",dat,part.indexSet(),0,
                         nofParams,false);
    else
      disp.displayVector("vtx. Paramters",dat,part.indexSet(),1,
                         nofParams,true);
  } else {
    GrapeGridDisplay<typename GridPart::GridType> disp(part.grid());
    disp.display();
  }
}
#else
template <class GridPart>
void test(GridPart& part,std::vector<double>& dat,int nofParams,int cdim) {
  gridcheck(part.grid());
}
#endif
template <class GridPart>
void test(GridPart& part) {
  gridcheck(part.grid());
}

int main(int argc, char ** argv, char ** envp)
try {
  // this method calls MPI_Init, if MPI is enabled
  MPIHelper & mpiHelper = MPIHelper::instance(argc,argv);
  int myrank = mpiHelper.rank();

  if (argc<2) {
    std::cerr << "supply grid file as parameter!" << std::endl;
    return 1;
  }

  std::cout << "tester: start grid reading" << std::endl;

  // create Grid from DGF parser
  GridPtr<GridType> grid;
  {
    grid = GridPtr<GridType>(argv[1], mpiHelper.getCommunicator() );
  }
  typedef LeafGridPart<GridType> GridPartType;
  typedef GridPartType::IndexSetType IndexSetType;
  GridPartType gridPart(*grid);
  const IndexSetType& index = gridPart.indexSet();
  std::vector<double> eldat(0),vtxdat(0);
  if (grid.nofParameters(0)>0)
  {
    std::cout << "Reading Element Parameters:" << std::endl;
    eldat.resize(index.size(0)*grid.nofParameters(0));
    typedef GridPartType::Codim<0>::IteratorType IteratorType;
    IteratorType endit = gridPart.end<0>();
    for (IteratorType iter=gridPart.begin<0>(); iter!=endit; ++iter) {
      std::vector<double>& param = grid.parameters(*iter);
      assert( (int) param.size() == grid.nofParameters(0));
      for (size_t i=0; i<param.size(); i++) {
        // std::cout << param[i] << " ";
        eldat[index.index(*iter)*param.size()+i] = param[i];
      }
      // std::cout << std::endl;
    }
  }
  if (grid.nofParameters(GridType::dimension)>0)
  {
    std::cout << "Reading Vertex Parameters:" << std::endl;
    vtxdat.resize(index.size(GridType::dimension)*grid.nofParameters(GridType::dimension));
    typedef GridPartType::Codim<GridType::dimension>::IteratorType IteratorType;
    IteratorType endit = gridPart.end<GridType::dimension>();
    for (IteratorType iter=gridPart.begin<GridType::dimension>();
         iter!=endit; ++iter) {
      std::vector<double>& param = grid.parameters(*iter);
      assert( (int) param.size() == grid.nofParameters(GridType::dimension));
      // std::cout << (*iter).geometry()[0] << " -\t ";
      for (size_t i=0; i<param.size(); i++) {
        // std::cout << param[i] << " ";
        vtxdat[index.index(*iter)*param.size()+i] = param[i];
      }
      // std::cout << std::endl;
    }
  }

  // display
  if (myrank <= 0)
  {
    if (0 && grid.nofParameters(0)>0)
      test(gridPart,eldat,grid.nofParameters(0),0);
    else
      test(gridPart,vtxdat,grid.nofParameters(GridType::dimension),GridType::dimension);
  }
  // refine
  std::cout << "tester: refine grid" << std::endl;
  grid->globalRefine(Dune::DGFGridInfo<GridType>::refineStepsForHalf());
  test(gridPart);
  return 0;
}
catch (Dune::Exception &e) {
  std::cerr << e << std::endl;
  return 1;
}
catch (...) {
  std::cerr << "Generic exception!" << std::endl;
  return 1;
}
