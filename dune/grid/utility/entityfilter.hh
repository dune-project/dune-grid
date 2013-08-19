// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UTILITY_ENTITYFILTER_HH
#define DUNE_GRID_UTILITY_ENTITYFILTER_HH

namespace Dune {

  class EntityFilterInterface
  {
  public:
    template<class Entity>
    bool contains(const Entity &e) const;
  };

  template<class IndexSet>
  class StridedEntityFilter
  {
    const IndexSet *is_;
    typename IndexSet::IndexType stride_;
    typename IndexSet::IndexType offset_;
  public:
    StridedEntityFilter(const IndexSet &is,
                        typename IndexSet::IndexType stride,
                        typename IndexSet::IndexType offset) :
      is_(&is), stride_(stride), offset_(offset)
    { }

    template<class Entity>
    bool contains(const Entity &e) const
    {
      return is_->index(e) % stride_ == offset_;
    }
  };

} // namespace Dune

#endif // DUNE_GRID_UTILITY_ENTITYFILTER_HH
