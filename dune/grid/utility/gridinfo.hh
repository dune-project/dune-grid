// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_UTILITY_GRIDINFO_HH
#define DUNE_GRID_UTILITY_GRIDINFO_HH

#include <algorithm>
#include <cstddef>
#include <functional>
#include <limits>
#include <map>
#include <ostream>
#include <string>
#include <vector>

#include <dune/common/classname.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/hybridutilities.hh>

#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/grid/common/mcmgmapper.hh>

namespace Dune {

  //! Structure to hold statistical information about one type of entity
  template<class ctype>
  struct EntityInfo {
    //! number of entities in the set
    std::size_t count;
    //! minimum volume of all entities in the set.
    /**
     * If the set is empty, this is \f$+\infty\f$.  If the volume of the
     * entities cannot be determined (some instance of the none GeometryType)
     * this is NaN.
     */
    ctype volumeMin;
    //! maximum volume of all entities in the set.
    /**
     * If the set is empty, this is \f$-\infty\f$.  If the volume of the
     * entities cannot be determined (some instance of the none GeometryType)
     * this is NaN.
     */
    ctype volumeMax;

    //! sum of volumes of all entities in the set.
    /**
     * If the set is empty, this is 0.  If the volume of the entities cannot
     * be determined (some instance of the none GeometryType) this is NaN.
     */
    ctype volumeSum;

    //! initialize the structure
    /**
     * This assumes an empty set of entities so that information can be added
     * later: \c count is set to 0, \c volumeMin to +infinity, and \c
     * volumeMax to -infinity.
     */
    EntityInfo() :
      count(0), volumeMin(std::numeric_limits<ctype>::infinity()),
      volumeMax(-std::numeric_limits<ctype>::infinity()), volumeSum(0)
    { }
  };

  //! Comparison object to sort GeometryType by majorly dimension
  /**
   * This differs from the standard GeometryType::operator<() in that it puts
   * more emphasis on the dimension than the isNone() property.  If it is used
   * as the comparison object in a map where the key is of type GeometryType,
   * all entries of one dimension will be lumped together.
   */
  struct GridViewInfoGTCompare :
    public std::binary_function<GeometryType, GeometryType, bool>
  {
    //! compare two GeometryTypes
    inline bool operator()(const GeometryType &a, const GeometryType &b) const
    {
      return a.dim() < b.dim() ||
             (a.dim() == b.dim() && (a.isNone() < b.isNone() ||
                                     (a.isNone() == b.isNone() && (a.id() >> 1) < (b.id() >> 1))));
      // topologyId is set to 0 for None, so no harm im comparing them even if
      // isNone()==true
    }
  };

  //! structure to hold information about a certain GridView.
  /**
   * This is a map from GeometryType to EntityInfo structures.  The entries in
   * the map are sorted in dimension-first, isNone()-second, and
   * topologyId-last order.
   */
  template<class ctype>
  struct GridViewInfo :
    public std::map<GeometryType, EntityInfo<ctype>, GridViewInfoGTCompare>
  {
    //! name of the grid class this information was extracted from
    std::string gridName;
    //! name of the class of the GridView this information was extracted from
    std::string gridViewName;
    //! name of the partition this information was extracted from
    /**
     * May be empty if not applicable (serial grids, for instance) or may be a
     * combination of partitions such as "interior+border".
     */
    std::string partitionName;

    //! print the information contained in this object
    /**
     * \param stream Stream object to print to.
     * \param prefix Prefix to print in front of each line.
     *
     * Sample output:
     * \verbatim
     * prefix>
     * \endverbatim
     *
     * If \c gridName, \c gridViewName, or \c partitionName is emtpy, the
     * corresponding line is not printed and no extra indentation is added for
     * the subsequent lines.
     */
    void print(std::ostream &stream, std::string prefix) const {
      if(!gridName.empty()) {
        stream << prefix << gridName << ":\n";
        prefix += "  ";
      }
      if(!gridViewName.empty()) {
        stream << prefix << gridViewName << ":\n";
        prefix += "  ";
      }
      if(!partitionName.empty()) {
        stream << prefix << partitionName << ":\n";
        prefix += "  ";
      }

      typedef typename GridViewInfo::const_iterator Iterator;
      std::size_t dim = ~0;
      const Iterator &end = this->end();
      for(Iterator it = this->begin(); it != end; ++it) {
        if(it->first.dim() != dim) {
          dim = it->first.dim();
          stream << prefix << "Dim = " << dim << ":\n";
        }
        stream << prefix << "  " << it->first << ": Count = "
               << it->second.count << ", Volume range = "
               << "(" << it->second.volumeMin << ".."
               << it->second.volumeMax << "), Total volume = "
               << it->second.volumeSum << "\n";
      }
    }
  };

  //! write a GridViewInfo object
  /**
   * \relates GridViewInfo
   *
   * This is equivalent to callinf info.print(stream, "").
   */
  template<class ctype>
  std::ostream &operator<<(std::ostream &stream,
                           const GridViewInfo<ctype> &info)
  {
    info.print(stream, "");
    return stream;
  }

#ifndef DOXYGEN
  //! operation for Hybrid::forEach Internally used by fillGridViewInfoSerial
  template<int codim>
  struct FillGridInfoOperation {
    template<class Entity, class Mapper, class Visited, class RefElem>
    static void apply(const Entity &e, const Mapper &mapper, Visited &visited,
                      const typename Entity::Geometry &geo,
                      RefElem refelem,
                      GridViewInfo<typename Entity::Geometry::ctype> &gridViewInfo)
    {
      typedef typename Entity::Geometry::ctype ctype;
      static const std::size_t dimw = Entity::Geometry::coorddimension;
      static const std::size_t dim = Entity::dimension;
      std::vector<FieldVector<ctype, dimw> > coords;
      for(int i = 0; i < refelem.size(codim); ++i) {
        int index = mapper.map(e, i, codim);
        if(visited[index])
          continue;
        visited[index] = true;

        GeometryType gt = refelem.type(i, codim);
        coords.clear();
        coords.resize( refelem.size(i, codim, dim) );
        for(std::size_t corner = 0; corner < coords.size(); ++corner)
          coords[ corner ] = geo.corner( refelem.subEntity( i, codim, corner, dim ) );
        MultiLinearGeometry<ctype, dim-codim, dimw> mygeo(gt, coords);

        ctype volume = mygeo.volume();
        EntityInfo<ctype> &ei = gridViewInfo[mygeo.type()];
        ei.volumeMin = std::min(ei.volumeMin, volume);
        ei.volumeMax = std::max(ei.volumeMax, volume);
        ei.volumeSum += volume;
      }
    }
  };
#endif // !DOXYGEN

  //! fill a GridViewInfo structure from a serial grid
  /**
   * If used on a parallel grid, it will gather information for entities of
   * all partitions on each rank locally.
   */
  template<class GV>
  void fillGridViewInfoSerial(const GV &gv,
                              GridViewInfo<typename GV::ctype> &gridViewInfo)
  {
    typedef typename GV::ctype ctype;
    static const std::size_t dim = GV::dimension;
    typedef typename GV::template Codim<0>::Iterator EIterator;
    typedef typename GV::template Codim<0>::Geometry EGeometry;
    typedef typename GV::IndexSet IndexSet;

    typedef typename GridViewInfo<ctype>::iterator InfoIterator;

    typedef ReferenceElements<ctype, dim> RefElems;

    MultipleCodimMultipleGeomTypeMapper<GV>
          mapper(gv,
                 [](GeometryType gt, int) { return gt.dim() < GV::dimension; }
                );
    std::vector<bool> visited(mapper.size(), false);

    gridViewInfo.gridName = className<typename GV::Grid>();
    gridViewInfo.gridViewName = className<GV>();
    gridViewInfo.partitionName = "";
    gridViewInfo.clear();

    const EIterator &eend = gv.template end<0>();
    for(EIterator eit = gv.template begin<0>(); eit != eend; ++eit) {
      ctype volume = eit->geometry().volume();
      EntityInfo<ctype> &ei = gridViewInfo[eit->type()];
      ei.volumeMin = std::min(ei.volumeMin, volume);
      ei.volumeMax = std::max(ei.volumeMax, volume);
      ei.volumeSum += volume;

      if(!eit->type().isNone()) {
        const EGeometry &geo = eit->geometry();
        Hybrid::forEach(std::make_index_sequence< dim >{},
          [ & ](auto i){ FillGridInfoOperation< i+1 >::apply(*eit, mapper, visited, geo, RefElems::general(eit->type()), gridViewInfo); } );
      }
    }

    GeometryType gt = GeometryTypes::none(dim);
    if(gridViewInfo.count(gt) > 0) {
      for(std::size_t codim = 0; codim < dim; ++codim)
      {
        gt = GeometryTypes::none(dim-codim);
        EntityInfo<ctype> & ei = gridViewInfo[gt];
        ei.volumeMin = ei.volumeMax = ei.volumeSum =
                                        std::numeric_limits<ctype>::quiet_NaN();
      }
      gt = GeometryTypes::none(0);
      EntityInfo<ctype> & ei = gridViewInfo[gt];
      ei.volumeMin = ei.volumeMax = ei.volumeSum = 0;
    }

    const InfoIterator &end = gridViewInfo.end();
    const IndexSet &is = gv.indexSet();
    for(InfoIterator it = gridViewInfo.begin(); it != end; ++it) {
      it->second.count = is.size(it->first);
      if(it->second.count == 0)
        DUNE_THROW(Exception, "Found Entities of geomentry type " <<
                   it->first << " while iterating through the grid, but "
                   "indexSet.size() == 0 for that geometry type");
    }

  }

} // namespace Dune


#endif // DUNE_GRID_UTILITY_GRIDINFO_HH
