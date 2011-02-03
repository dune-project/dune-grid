// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_UTILITY_GRIDINFO_HH
#define DUNE_GRID_UTILITY_GRIDINFO_HH

#include <cstddef>
#include <functional>
#include <limits>
#include <map>
#include <ostream>
#include <string>

#include <dune/common/geometrytype.hh>

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

    //! initialize the structure
    /**
     * This assumes an empty set of entities so that information can be added
     * later: \c count is set to 0, \c volumeMin to +infinity, and \c
     * volumeMax to -infinity.
     */
    EntityInfo() :
      count(0), volumeMin(std::numeric_limits<ctype>::infinity()),
      volumeMax(-std::numeric_limits<ctype>::infinity())
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
       prefix>
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
               << it->second.volumeMax << ")\n";
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

} // namespace Dune


#endif // DUNE_GRID_UTILITY_GRIDINFO_HH
