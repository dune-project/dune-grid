// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
 #include "config.h"
 #include <dune/grid/io/file/dgfparser/dgfparser.hh>
 #include <dune/grid/io/file/dgfparser/dgfug.hh>



int main(int argc, char** argv)
{
  try{
    // define the problem dimensions
    const int dim=2;

    // create a grid object
    typedef Dune::UGGrid<dim> GridType;
    typedef GridType::ctype DT;
    typedef double NumberType;

    std::stringstream dgfFileName;
    dgfFileName << "test2ug.dgf";

    // create grid pointer
    Dune::GridPtr<GridType> gridPtr1( dgfFileName.str() );

    Dune::GridPtr<GridType> gridPtr2( dgfFileName.str() );

    return 0;
  }
  catch (Dune::Exception &e) {
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
