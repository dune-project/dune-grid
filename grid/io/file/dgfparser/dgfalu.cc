// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
namespace Dune {

  //***********************************************
  template <>
  inline ALUSimplexGrid<3,3>*
  MacroGrid :: Impl<ALUSimplexGrid<3,3> >::generate
    (MacroGrid& mg,const char* filename, int MPICOMM) {
    mg.element=Simplex;
    std::string str(filename);
    MacroGrid::Impl<ALUSimplexGrid<3,3> >().
    generateAlu3d(mg,filename,str,MPICOMM);
  #if defined _ALU3DGRID_PARALLEL_
    return new ALUSimplexGrid<3,3>(str.c_str(),MPICOMM);
  #else
    return new ALUSimplexGrid<3,3>(str.c_str());
  #endif
  }
  template <>
  inline ALUCubeGrid<3,3>*
  MacroGrid :: Impl<ALUCubeGrid<3,3> >::generate
    (MacroGrid& mg,const char* filename, int MPICOMM) {
    mg.element=Cube;
    std::string str(filename);
    MacroGrid::Impl<ALUCubeGrid<3,3> >().
    generateAlu3d(mg,filename,str,MPICOMM);
  #if defined _ALU3DGRID_PARALLEL_
    return new ALUCubeGrid<3,3>(str.c_str(),MPICOMM);
  #else
    return new ALUCubeGrid<3,3>(str.c_str());
  #endif
  }

  template <>
  inline ALUSimplexGrid<2,2>*
  MacroGrid :: Impl<ALUSimplexGrid<2,2> >::generate
    (MacroGrid& mg,const char* filename, int MPICOMM) {
    mg.element=Simplex;
    std::string str(filename);
    MacroGrid::Impl<ALUSimplexGrid<2,2> >().
    generateAlu3d(mg,filename,str,MPICOMM);
    //#if defined _ALU3DGRID_PARALLEL_
    //return new ALUSimplexGrid<2,2>(str.c_str(),MPICOMM);
    //#else
    return new ALUSimplexGrid<2,2>(str.c_str());
    //#endif
  }
  template <int dim,int dimworld>
  inline void
  MacroGrid :: Impl<ALUSimplexGrid<dim,dimworld> > :: generateAlu3d
    (MacroGrid& mg,const char* filename, std::string& str,int MPICOMM) {
    int myrank=-1;
  #if defined _ALU3DGRID_PARALLEL_
    MPI_Comm_rank(MPICOMM,&myrank);
  #endif
    if (myrank<=0) {
      std::ifstream gridin(filename);
      if (mg.readDuneGrid(gridin) == 1) {
        if (mg.dimw != dimworld) {
          std::cerr << "ERROR: "
                    << "Macrofile " << filename << " is for dimension " << mg.dimw
                    << " and connot be used to initialize an ALUGrid of dimension "
                    << dimworld << std::endl;
          abort();
        }
        mg.setOrientation(0);
        str+=".ALUgrid";
        std::string str1(str);
      #if defined _ALU3DGRID_PARALLEL_
        str1+=".0";
      #endif
        std::ofstream out(str1.c_str());
        mg.writeAlu(out);
      }
    }
  }
  template <int dim,int dimworld>
  inline void
  MacroGrid :: Impl<ALUCubeGrid<dim,dimworld> > :: generateAlu3d
    (MacroGrid& mg,const char* filename, std::string& str,int MPICOMM) {
    int myrank=-1;
  #if defined _ALU3DGRID_PARALLEL_
    MPI_Comm_rank(MPICOMM,&myrank);
  #endif
    if (myrank<=0) {
      std::ifstream gridin(filename);
      if (mg.readDuneGrid(gridin) == 1) {
        if (mg.dimw != dimworld) {
          std::cerr << "ERROR: "
                    << "Macrofile " << filename << " is for dimension " << mg.dimw
                    << " and connot be used to initialize an ALUGrid of dimension "
                    << dimworld << std::endl;
          abort();
        }
        mg.setOrientation(0);
        str+=".ALUgrid";
        std::string str1(str);
      #if defined _ALU3DGRID_PARALLEL_
        str1+=".0";
      #endif
        std::ofstream out(str1.c_str());
        mg.writeAlu(out);
      }
    }
  }

}
