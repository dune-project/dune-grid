// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/** \file
    \brief A unit test for the VertexOrder classes
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <algorithm>
#include <array>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <memory>
#include <ostream>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/generalvertexorder.hh>

#include <dune/grid/onedgrid.hh>
#include <dune/grid/uggrid.hh>
#include "../structuredgridfactory.hh"
#include "../vertexorderfactory.hh"

void fail(int &result) {
  result = 1;
}
void pass(int &result) {
  if(result == 77) result = 0;
}

//! test consistency on one dimension and element
template<std::size_t mydim, class VertexOrder>
void testElementDim(const std::integral_constant<std::size_t, mydim>&,
                    const VertexOrder &vo)
{
  static const std::size_t dim = VertexOrder::dimension;
  static const std::size_t codim = dim - mydim;
  auto refelem = Dune::referenceElement<double, dim>(vo.type());

  std::vector<typename VertexOrder::Index> subOrder, tmp;
  for(int subentity = 0; subentity < refelem.size(codim); ++subentity)
  {
    std::size_t voSize = std::distance(vo.begin(codim, subentity),
                                       vo.end(codim, subentity));
    std::size_t refSize = refelem.size(subentity, codim, dim);
    if(voSize != refSize)
      DUNE_THROW(Dune::Exception, "Codim " << codim << ": Subentity " <<
                 subentity << ": Number of vertices: " << voSize << " in "
                 "VertexOrder but " << refSize << " in reference element");
    subOrder.resize(voSize);
    Dune::reduceOrder(vo.begin(codim, subentity), vo.end(codim, subentity),
                      subOrder.begin());
    vo.getReduced(codim, subentity, tmp);
    if(!(subOrder == tmp))
      DUNE_THROW(Dune::Exception, "Codim " << codim << ": Subentity " <<
                 subentity << ": Reduced vertex numbering from VertexOrder "
                 "does not match manually reduced vertex numbering");
  }
}

// test inter-dimensional consistency
template<std::size_t mydim, class VertexOrder>
void testElementInterdim(const std::integral_constant<std::size_t, mydim>&,
                         const VertexOrder &vo)
{
  testElementDim(std::integral_constant<std::size_t, mydim>(), vo);

  static const std::size_t dim = VertexOrder::dimension;
  static const std::size_t codim = dim - mydim;
  auto refelem = Dune::referenceElement<double, dim>(vo.type());

  std::vector<typename VertexOrder::Index> subOrder, tmp, subsubOrder;
  for(int subentity = 0; subentity < refelem.size(codim); ++subentity)
  {
    auto subrefelem = Dune::referenceElement<double, mydim>(refelem.type(subentity, codim));

    subOrder.assign(vo.begin(codim, subentity), vo.end(codim, subentity));

    for(int subsubentity = 0; subsubentity < subrefelem.size(1);
        ++subsubentity)
    {
      tmp.resize(subrefelem.size(subsubentity, 1, mydim));
      for(std::size_t vertex = 0; vertex < tmp.size(); ++vertex)
        tmp[vertex] = subOrder[subrefelem.subEntity(subsubentity, 1,
                                                    vertex, mydim)];
      subsubOrder.resize(tmp.size());
      Dune::reduceOrder(tmp.begin(), tmp.end(), subsubOrder.begin());
      vo.getReduced(codim+1, refelem.subEntity(subentity, codim, subsubentity,
                                               codim+1), tmp);
      if(!(subsubOrder == tmp))
        DUNE_THROW(Dune::Exception, "Codim " << codim << ": Subentity " <<
                   subentity << ": Sub-Codim " << codim+1 <<": Sub-Subentity "
                                             << subsubentity << ": Vertex numbering from subentity does "
                   "not match vertex numbering from sub-subentity.");
    }
  }
}
template<class VertexOrder>
void testElementInterdim(const std::integral_constant<std::size_t, 0>&,
                         const VertexOrder &vo)
{
  testElementDim(std::integral_constant<std::size_t, 0>(), vo);
}

template<std::size_t mydim, class VertexOrder, class Intersection>
void testNeighborDim(const std::integral_constant<std::size_t, mydim>&,
                     const VertexOrder &vo_s, const VertexOrder &vo_n,
                     const Intersection &is)
{
  static const std::size_t dim = VertexOrder::dimension;
  static const std::size_t codim = dim - mydim;
  auto refelem_s = Dune::referenceElement<double, dim>(vo_s.type());
  auto refelem_n = Dune::referenceElement<double, dim>(vo_n.type());

  std::size_t index_s = is.indexInInside();
  std::size_t index_n = is.indexInOutside();

  typedef typename Intersection::Entity Entity;
  Entity inside = is.inside();
  Entity outside = is.outside();

  typedef typename Entity::Geometry Geometry;
  const Geometry &geo_s = inside.geometry();
  const Geometry &geo_n = outside.geometry();

  typedef typename Intersection::ctype DF;
  typedef Dune::FieldVector<DF, Intersection::dimensionworld> DomainW;
  DF eps = std::pow(geo_s.volume(), 1.0/dim)*1e-6;

  std::vector<typename VertexOrder::Index> order_s;
  std::vector<typename VertexOrder::Index> order_n;

  for(int subentity_s = 0; subentity_s < refelem_s.size(index_s, 1, codim);
      ++subentity_s)
  {
    std::size_t subindex_s =
      refelem_s.subEntity(index_s, 1, subentity_s, codim);
    DomainW subpos = geo_s.global(refelem_s.position(subindex_s, codim));
    int subentity_n = 0;
    std::size_t subindex_n = 0;
    for(; subentity_n < refelem_n.size(index_n, 1, codim); ++subentity_n) {
      subindex_n = refelem_n.subEntity(index_n, 1, subentity_n, codim);
      if((subpos - geo_n.global(refelem_s.position(subindex_n, codim))
          ).two_norm() < eps)
        break;
    }
    if(subentity_n == refelem_n.size(index_n, 1, codim))
      DUNE_THROW(Dune::Exception, "Codim " << codim << ": Subentity " <<
                 subindex_s << ": Can't find corresponding subentity in "
                 "neighbor");
    vo_s.getReduced(codim, subindex_s, order_s);
    vo_n.getReduced(codim, subindex_n, order_n);
    if(order_s.size() != order_n.size())
      DUNE_THROW(Dune::Exception, "Codim " << codim << ": Subentity " <<
                 subindex_s << ": Number of vertices is " << order_s.size() <<
                 " in inside but " << order_n.size() << " in outside");
    for(std::size_t i_s = 0; i_s < order_s.size(); ++i_s) {
      std::size_t i_n = 0;
      for(; i_n < order_n.size(); ++i_n)
        if(order_s[i_s] == order_n[i_n])
          break;
      if(i_n == order_n.size())
        DUNE_THROW(Dune::Exception, "Codim " << codim << ": Subentity " <<
                   subindex_s << ": Vertex " << i_s << ": Can't find "
                   "corresponding vertex in neighbor");
      DomainW pos_s = geo_s.corner(refelem_s.subEntity(subindex_s, codim, i_s,
                                                       dim));
      DomainW pos_n = geo_n.corner(refelem_n.subEntity(subindex_n, codim, i_n,
                                                       dim));
      if((pos_s - pos_n).two_norm() > eps)
        DUNE_THROW(Dune::Exception, "Codim " << codim << ": Subentity " <<
                   subindex_s << " at " << subpos << ": Vertex " << i_s << " "
                   "is at " << pos_s << "; Neighbor at " <<
                   geo_n.center() << ": Subentity " << subindex_n << ": "
                   "Vertex " << i_n << " is at " << pos_n);
    }
  }
}
template<class VertexOrder, class Intersection>
void testNeighbor(const std::integral_constant<std::size_t, 0>&,
                  const VertexOrder &vo_s, const VertexOrder &vo_n,
                  const Intersection &is)
{
  testNeighborDim(std::integral_constant<std::size_t, 0>(),
                  vo_s, vo_n, is);
}
template<std::size_t mydim, class VertexOrder, class Intersection>
void testNeighbor(const std::integral_constant<std::size_t, mydim>&,
                  const VertexOrder &vo_s, const VertexOrder &vo_n,
                  const Intersection &is)
{
  testNeighbor(std::integral_constant<std::size_t, mydim-1>(),
               vo_s, vo_n, is);
  testNeighborDim(std::integral_constant<std::size_t, mydim>(),
                  vo_s, vo_n, is);
}

template<std::size_t codim, class GV, class VertexOrderFactory>
void testVertexOrder(const GV& gv, const VertexOrderFactory &voFactory,
                     int& result)
{
  pass(result);

  static const std::size_t dim = GV::dimension;
  typedef typename VertexOrderFactory::template VertexOrder<dim>::type
  VertexOrder;
  typedef typename GV::template Codim<0>::Iterator EIterator;
  typedef typename GV::IntersectionIterator IIterator;

  const EIterator &eend = gv.template end<0>();
  for(EIterator eit = gv.template begin<0>(); eit != eend; ++eit)
    try {
      VertexOrder vo = voFactory.make(*eit);
      testElementInterdim(std::integral_constant<std::size_t, dim-codim>(),
                          vo);
      const IIterator &iend = gv.iend(*eit);
      for(IIterator iit = gv.ibegin(*eit); iit != iend; ++iit)
        if(iit->neighbor()) {
          VertexOrder vo_n = voFactory.make(iit->outside());
          testNeighbor(std::integral_constant<std::size_t,
                           dim-(codim==0 ? 1 : codim)>(),
                       vo, vo_n, *iit);
        }
    }
    catch(const Dune::Exception &e) {
      std::cout << "Element at " << eit->geometry().center() << ": "
                << e.what() << std::endl;
      fail(result);
    }
}

template<class Grid>
void testVertexOrderByIdSimplices(int &result) {
  static const std::size_t dim = Grid::dimension;
  static const std::size_t dimworld = Grid::dimensionworld;
  typedef typename Grid::ctype DF;
  typedef Dune::FieldVector<DF, dimworld> Domain;

  std::array<unsigned int, dim> elements;
  std::fill(elements.begin(), elements.end(), 4);

  std::shared_ptr<Grid> gridp = Dune::StructuredGridFactory<Grid>::
                                 createSimplexGrid(Domain(0), Domain(1), elements);

  typedef typename Grid::GlobalIdSet IdSet;
  typedef Dune::VertexOrderByIdFactory<IdSet> VOFactory;
  VOFactory voFactory(gridp->globalIdSet());

  testVertexOrder<0>(gridp->leafGridView(), voFactory, result);
}

template<class Grid>
void testVertexOrderByIdCubes(int &result) {
  static const std::size_t dim = Grid::dimension;
  static const std::size_t dimworld = Grid::dimensionworld;
  typedef typename Grid::ctype DF;
  typedef Dune::FieldVector<DF, dimworld> Domain;

  std::array<unsigned int, dim> elements;
  std::fill(elements.begin(), elements.end(), 4);

  std::shared_ptr<Grid> gridp = Dune::StructuredGridFactory<Grid>::
                                 createCubeGrid(Domain(0), Domain(1), elements);

  typedef typename Grid::GlobalIdSet IdSet;
  typedef Dune::VertexOrderByIdFactory<IdSet> VOFactory;
  VOFactory voFactory(gridp->globalIdSet());

  testVertexOrder<0>(gridp->leafGridView(), voFactory, result);
}

int main (int argc , char **argv)
try {

  // this method calls MPI_Init, if MPI is enabled
  Dune::MPIHelper::instance(argc,argv);

  int result = 77;

  //////////////////////////////////////////////////////////////////////
  //   Test 1d grids
  //////////////////////////////////////////////////////////////////////

  std::cout << "= Testing 1D" << std::endl;

  std::cout << "== Testing OneDGrid" << std::endl;
  testVertexOrderByIdCubes<Dune::OneDGrid>(result);

  //////////////////////////////////////////////////////////////////////
  //   Test 2d grids
  //////////////////////////////////////////////////////////////////////

  std::cout << "= Testing 2D" << std::endl;

#if HAVE_DUNE_UGGRID
  std::cout << "== Testing UGGrid<2> with simplices" << std::endl;
  testVertexOrderByIdCubes<Dune::UGGrid<2> >(result);
  std::cout << "== Testing UGGrid<2> with cubes" << std::endl;
  testVertexOrderByIdSimplices<Dune::UGGrid<2> >(result);
#endif

  //////////////////////////////////////////////////////////////////////
  //   Test 3d grids
  //////////////////////////////////////////////////////////////////////

  std::cout << "= Testing 3D" << std::endl;

#if HAVE_DUNE_UGGRID
  std::cout << "== Testing UGGrid<3> with simplices" << std::endl;
  testVertexOrderByIdCubes<Dune::UGGrid<3> >(result);
  std::cout << "== Testing UGGrid<3> with cubes" << std::endl;
  testVertexOrderByIdSimplices<Dune::UGGrid<3> >(result);
#endif

  return result;
}
catch (const Dune::Exception &e) {
  std::cerr << e << std::endl;
  throw;
} catch (...) {
  std::cerr << "Generic exception!" << std::endl;
  throw;
}
