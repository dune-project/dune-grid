// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGFPARSERALU_HH
#define DUNE_DGFPARSERALU_HH
#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#include "dgfparser.hh"
namespace Dune {
  /*
     template <int dim,int dimworld,ALU3dGridElementType elType>
     class MacroGrid::Impl<ALU3dGrid<dim,dimworld,elType> > {
     public:
     static ALU3dGrid<dim,dimworld,elType>*
     generate(MacroGrid& mg,
             const char* filename,int MPICOMM=-1);
     private:
     inline void
     generateAlu3d(MacroGrid& mg,
                  const char* filename, std::string& str,int MPICOMM);
     };
   */
  //*********************************
  template <int dim,int dimworld>
  class MacroGrid::Impl<ALUCubeGrid<dim,dimworld> > {
  public:
    static ALUCubeGrid<dim,dimworld>*
    generate(MacroGrid& mg,
             const char* filename,int MPICOMM=-1);
  private:
    inline void
    generateAlu3d(MacroGrid& mg,
                  const char* filename, std::string& str,int MPICOMM);
  };
  template <int dim,int dimworld>
  class MacroGrid::Impl<ALUSimplexGrid<dim,dimworld> > {
  public:
    static ALUSimplexGrid<dim,dimworld>*
    generate(MacroGrid& mg,
             const char* filename,int MPICOMM=-1);
  private:
    inline void
    generateAlu3d(MacroGrid& mg,
                  const char* filename, std::string& str,int MPICOMM);
  };
}
#include "dgfalu.cc"
#endif
#endif
