// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_UTILITY_GRIDNUMBERING_TO_VTK_HH
#define DUNE_GRID_UTILITY_GRIDNUMBERING_TO_VTK_HH

#include <string>

#include <dune/common/fvector.hh>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/common.hh>
#include <dune/grid/io/file/vtk/function.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

template<class GV>
class IndexSetVTKFunction :
  public Dune::VTKFunction<GV>
{
  typedef Dune::VTKFunction<GV> Base;

public:
  typedef typename Base::Entity Entity;
  typedef typename Base::ctype ctype;
  enum { dim = Base::dim };

  IndexSetVTKFunction(const GV &gv, const std::string &name) :
    gv_(gv), name_(name)
  { }

  int ncomps() const
  {
    return 1;
  }

  double evaluate(int comp, const Entity &e,
                  const Dune::FieldVector<ctype,dim> &) const
  {
    return gv_.indexSet().index(e);
  }

  std::string name() const
  {
    return name_;
  }

private:
  GV gv_;
  std::string name_;
};

template<class GV, class Mapper>
class ElementMapperVTKFunction :
  public Dune::VTKFunction<GV>
{
  typedef Dune::VTKFunction<GV> Base;

public:
  typedef typename Base::Entity Entity;
  typedef typename Base::ctype ctype;
  enum { dim = Base::dim };

  ElementMapperVTKFunction(const Mapper &mapper, const std::string &name) :
    mapper_(mapper), name_(name)
  { }

  int ncomps() const
  {
    return 1;
  }

  double evaluate(int comp, const Entity &e,
                  const Dune::FieldVector<ctype,dim> &xi) const
  {
    return mapper_.map(e);
  }

  std::string name() const
  {
    return name_;
  }

private:
  Mapper mapper_;
  std::string name_;
};

template<class GV>
void numberingToVTK(const GV &gv, const std::string &vtkPrefix)
{
  typedef Dune::MultipleCodimMultipleGeomTypeMapper
    <GV, Dune::MCMGElementLayout> ElementMapper;
  ElementMapper elementMapper(gv);

  Dune::VTKWriter<GV> writer(gv);
  writer.addCellData(new IndexSetVTKFunction<GV>(gv, "index"));
  writer.addCellData(new ElementMapperVTKFunction<GV, ElementMapper>
                     (elementMapper, "elementmap"));
  writer.write(vtkPrefix, Dune::VTK::appendedraw);
}


#endif // DUNE_GRID_UTILITY_GRIDNUMBERING_TO_VTK_HH
