#ifndef DUNE_GRID_COMMON_SUBINDEXIDSET_HH
#define DUNE_GRID_COMMON_SUBINDEXIDSET_HH

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/geometry/typeindex.hh>

#include <dune/common/hybridutilities.hh>

#include <set>
#include <map>
#include <unordered_map>
#include <vector>
#include <bitset>
#include <numeric>
#include <utility>

namespace Dune {

  // import range generators to make sure they work with PartitionView
  using Dune::entities;
  // using Dune::elements;
  // using Dune::facets;
  // using Dune::edges;
  // using Dune::vertices;
  // using Dune::descendantElements;
  // using Dune::intersections;

/*
 *
 * Marked entities                   [  x x x     x x   x]
 * Marked gt flips marked entities
 * Base index                        [0 1 2 3 4 5 6 7 8 9]
 * Unused offset                     [0 1 2 3 3 3 4 5 5 6]
 * Red. unused offset                [  3         5     6]
 * New index                         [0       1 2     3  ]
 */

  // TODO support UG grid
  template<typename GV>
  class SubIndexSet
  {
    using single_gt_helper = Dune::Capabilities::hasSingleGeometryType<typename GV::Grid>;

    static constexpr bool has_single_gt = single_gt_helper::v;
    static constexpr unsigned int id_single_gt = single_gt_helper::topologyId;

    static constexpr std::size_t max_gt_types = (has_single_gt
                                      ? GlobalGeometryTypeIndex::size(GV::dimension)
                                      : 1
                                    );
  public:
    /** \brief Export the type of the entity used as parameter in the index(...) method */
    template <int cc>
    struct Codim
    {
      using Entity = typename GV::Traits::template Codim<cc>::Entity;
    };

    /** \brief The type used for the indices */
    using IndexType = typename GV::Traits::IndexSet::IndexType;

    /** \brief iterator range for geometry types in domain */
    using Types = typename GV::Traits::IndexSet::Types;

    /** \brief dimension of the grid (maximum allowed codimension) */
    static constexpr int dimension = GV::Traits::Grid::dimension;

  private:

    // template<class Grid, int dim, int... Args>
    // static constexpr auto iterable_codimensions(std::integer_sequence<int,Args...> i)
    // {
    //   if (dim == 0)
    //     return std::integer_sequence<int,Args...>{};
    //   else if (Dune::Capabilities::hasEntityIterator<Grid,dim>::v)
    //     return iterable_codimensions<Grid,dim-1>(std::integer_sequence<int,Args...,dim>{});
    //   else
    //     return iterable_codimensions<Grid,dim-1>(std::integer_sequence<int,Args...>{});
    // }

    // template<class Grid>
    // static constexpr auto iterable_codimensions()
    // {
    //   return iterable_codimensions<Grid,0>(std::integer_sequence<std::size_t>{});
    // }

  public:

    //! Update flags
    bool _needs_update_base, _needs_update;

    //! Geometry offsets w.r.t. base index set
    std::array<std::size_t,max_gt_types+1> _gt_off_base;

    //! Geometry offsets w.r.t. this index set
    std::array<std::size_t,max_gt_types+1>    _gt_off;

    std::unordered_map<std::size_t,std::map<IndexType,std::size_t>>   _index_off;

    std::bitset<max_gt_types>               _active_gt, _marked_gt;
    std::set<IndexType>                     _marked_entity;
    std::vector<GeometryType>               _gt;
    GV                                      _gv;

  public:
    SubIndexSet(GV grid_view)
      : _needs_update_base(true)
      , _needs_update(true)
      , _gv(grid_view)
    {
      mark();
      update_base_offsets();
    }

  private:
    inline static constexpr std::size_t gt_index(const GeometryType& gt)
    {
      return has_single_gt ? GlobalGeometryTypeIndex::index(gt) : 0;
    }

    bool update_base_offsets(bool force = false)
    {
      if (!(_needs_update_base || force))
        return false;

      std::fill(_gt_off_base.begin(),_gt_off_base.end(),0);

      for (int codim = 0; codim <= dimension; ++codim)
        for (const auto& gt : _gv.indexSet().types(codim))
          _gt_off_base[gt_index(gt)+1] = _gv.indexSet().size(gt);

      std::partial_sum(_gt_off_base.begin(),_gt_off_base.end(),_gt_off_base.begin());
      _needs_update_base = false;
      _needs_update = true;
      return true;
    }

    template<class E>
    IndexType base_index(const E& e) const
    {
      return _gv.indexSet().index(e);
    }

    template<class E>
    IndexType unique_base_index(const E& e) const
    {
      assert(not _needs_update_base);
      return _gv.indexSet().index(e) + _gt_off_base[gt_index(e.type())];
    }

  public:
    bool mark()
    {
      bool modified = false;
      for (int codim = 0; codim < dimension + 1; ++codim)
        for (const auto& gt : _gv.indexSet().types(codim))
          modified |= mark(gt);
      return modified;
    }

    bool unmark()
    {
      bool modified = false;
      for (int codim = 0; codim < dimension + 1; ++codim)
        for (const auto& gt : _gv.indexSet().types(codim))
          modified |= unmark(gt);
      return modified;
    }

    bool mark(GeometryType type)
    {
      const auto gt_id = gt_index(type);
      bool modified = not (_marked_entity.empty() or _marked_gt.test(gt_id));
      _marked_entity.clear();
      _marked_gt.set(gt_id,true);
      _active_gt.set(gt_id,true);
      _needs_update |= modified;
      return modified;
    }

    bool unmark(GeometryType type)
    {
      const auto gt_id = gt_index(type);
      bool modified = (not _marked_entity.empty() or _marked_gt.test(gt_id));
      _marked_entity.clear();
      _marked_gt.set(gt_id,false);
      _active_gt.set(gt_id,false);
      _needs_update |= modified;
      return modified;
    }

    bool mark(int codim)
    {
      bool modified = false;
      for (const auto& gt : _gv.indexSet().types(codim))
        modified |= mark(gt);
      return modified;
    }

    bool unmark(int codim)
    {
      bool modified = false;
      for (const auto& gt : _gv.indexSet().types(codim))
        modified |= unmark(gt);
      return modified;
    }

    template<class E>
    bool mark(const E& e)
    {
      const auto gt_id = gt_index(e.type());
      if (_marked_gt.test(gt_id)) {
        return false;
      } else {
        _active_gt.set(gt_id,true);
        const auto v = _marked_entity.insert(unique_base_index(e));
        _needs_update |= v.second;
        return v.second;
      }
    }

    template<class E>
    bool unmark(const E& e)
    {
      const auto gt_id = gt_index(e.type());
      if (_marked_gt.test(gt_id)) {
        _active_gt.set(gt_id,true);
        const auto v = _marked_entity.insert(unique_base_index(e));
        _needs_update |= v.second;
        return v.second;
      } else {
        return false;
      }
    }

    template<class E>
    bool getMark(const E& e) const
    {
      if (!_gv.indexSet().contains(e))
        return false;

      const auto gt_id = gt_index(e.type());
      if (not _active_gt.test(gt_id))
        return false;

      const auto index = unique_base_index(e);
      if (_marked_gt.test(gt_id)) // masked entities => unmark
        return not _marked_entity.count(index);
      else                           // masked entities => mark
        return _marked_entity.count(index);
    }

    template<PartitionIteratorType partition_iterator = InteriorBorder_Partition>
    bool update(bool force = false)
    {
      update_base_offsets(force);
      if (!(_needs_update || force))
        return false;

      update_offsets<partition_iterator>();
      // update_all_offsets<partition_iterator>(); // test
      update_gt_maps();

      _needs_update = false;
      return true;
    }

    template<PartitionIteratorType partition_iterator>
    void update_offsets()
    {
      // initialize variables
      std::fill(_gt_off.begin(),_gt_off.end(),0);
      _gt.resize(0);

      // static loop over codim
      const auto indices = std::make_index_sequence<dimension+1>{};
      Dune::Hybrid::forEach(indices, [&](auto codim) {
        for (const auto& gt : _gv.indexSet().types(codim))
          if (_active_gt.test(gt_index(gt)))
            _gt.push_back(gt);
        update_codim_offsets<partition_iterator>(Codim<codim>{});
      });

      std::partial_sum(_gt_off.begin(),_gt_off.end(),_gt_off.begin());
    }

    template<PartitionIteratorType partition_iterator, int codim>
    void update_codim_offsets(Codim<codim>)
    {
      constexpr auto partition = partitionSet<partition_iterator>();

      for (const auto& e : entities(_gv,Dune::Codim<codim>{},partition)) {

        if (not partition.contains(e.partitionType()))
          unmark(e);

        auto gt_id = gt_index(e.type());
        if (getMark(e)) {
          _gt_off[gt_id+1] += 1;
        } else {
          _index_off[gt_id][base_index(e)] = 1;
        }
      }
    }

    // template<PartitionIteratorType partition_iterator>
    // void update_all_offsets()
    // {
    //   auto indices = iterable_dimensions<typename GV::Grid>();
    // }

    void update_gt_maps()
    {
      for (const auto& gt : _gt) {
        const auto gt_id = gt_index(gt);
        auto& gt_index_off = _index_off[gt_id];

        // reduce offset map to contain only the last of consecutive unmarked entities
        auto it_i = gt_index_off.rbegin();
        while(it_i != gt_index_off.rend() and std::next(it_i) != gt_index_off.rend())
        {
          auto it_j = std::next(it_i);
          auto id_i = it_i->first;
          auto id_j = it_j->first;
          if (id_i == id_j + 1)
          {
            it_j->second += it_i->second;
            gt_index_off.erase(it_j.base()); // equivalent to erase it_i
            it_i = it_j;
          } else {
            it_i++;
          }
        }

        // partial sum in the maps
        int partial_sum = 0;
        for (auto& [key, value] : gt_index_off) {
          partial_sum += value;
          gt_index_off[key] = partial_sum;
        }
      }
    }


    template<typename E>
    IndexType index(const E& e) const
    {
      assert(not _needs_update);
      // assert(Partitions::contains(e.partitionType()));
      // assert(contains(e.type()));

      const auto id = base_index(e);
      const auto gt_id = gt_index(e.type());

      const auto& gt_index_off = _index_off.find(gt_id)->second;

      auto unused_indices = 0;
      // calculate number of unused indices
      const auto unused_indices_map = gt_index_off.upper_bound(id);
      if (unused_indices_map != gt_index_off.begin())
        unused_indices = std::prev(unused_indices_map)->second;
      return id - unused_indices;
    }

    template<typename E>
    IndexType unique_index(const E& e) const
    {
      // assert(Partitions::contains(e.partitionType()));
      // assert(contains(e.type()));

      assert(not _needs_update);
      return index(e) + _gt_off[gt_index(e.type())];
    }

    template<class Entity>
    bool contains(const Entity& e) const
    {
      // todo: take into account partition used for update
      assert(not _needs_update);
      bool contained = true;

      const auto gt_id = gt_index(e.type());

      if (not _active_gt.test(gt_id)) {
        contained = false;
      } else {
        const auto id = base_index(e);
        const auto& gt_index_off = _index_off.find(gt_id)->second;
        if (not gt_index_off.empty()) {
          // calculate number of unused indices
          auto unused_indices_map = gt_index_off.upper_bound(id);
          if (unused_indices_map != gt_index_off.begin())
          {
            unused_indices_map--;
            int off_unused_indices = unused_indices_map->second;
            if (unused_indices_map != gt_index_off.begin())
              off_unused_indices -= std::prev(unused_indices_map)->second;
            contained = (id - unused_indices_map->first) >= off_unused_indices;
          }
        }
      }

      assert(contained == getMark(e));
      return contained;
    }

    IndexType size(GeometryType type) const
    {
      if constexpr (has_single_gt and (type.id() != id_single_gt))
        return 0;
      else {
        auto gt_id = gt_index(type);
        return _gt_off[gt_id+1]-_gt_off[gt_id];
      }
    }

  };

} // namespace Dune


#endif // DUNE_GRID_COMMON_SUBINDEXIDSET_HH
