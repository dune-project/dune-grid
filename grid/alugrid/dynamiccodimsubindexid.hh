// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DYNAMIC_CODIM_SUB_INDEX_ID_HH
#define DUNE_DYNAMIC_CODIM_SUB_INDEX_ID_HH

namespace Dune {

  // simple meta program to get subentity ids with dynamic codim
  template<class IdSet, class Element, int d>
  class DynamicCodimSubId
  {
  public:
    static typename IdSet::IdType get(const IdSet& idSet, const Element& element, int subEntity, unsigned int codim)
    {
      if (codim==d)
        return idSet.template subId<d>(element, subEntity);

      return DynamicCodimSubId<IdSet, Element, d-1>::get(idSet, element, subEntity, codim);
    }
  };

  template<class IdSet, class Element>
  class DynamicCodimSubId<IdSet,Element,-1>
  {
  public:
    static typename IdSet::IdType get(const IdSet& idSet, const Element& element, int subEntity, unsigned int codim)
    {
      DUNE_THROW(Dune::GridError, "The requested codim is not in {0,..,dim}!");
      return 0;
    }
  };



  // simple meta program to get subentity indices with dynamic codim
  template<class IndexSet, class Element, int d>
  class DynamicCodimSubIndex
  {
  public:
    static typename IndexSet::IndexType get(const IndexSet& indexSet, const Element& element, int subEntity, unsigned int codim)
    {
      if (codim==d)
        return indexSet.template subIndex<d>(element, subEntity);

      return DynamicCodimSubIndex<IndexSet, Element, d-1>::get(indexSet, element, subEntity, codim);
    }
  };

  template<class IndexSet, class Element>
  class DynamicCodimSubIndex<IndexSet,Element,-1>
  {
  public:
    static typename IndexSet::IndexType get(const IndexSet& indexSet, const Element& element, int subEntity, unsigned int codim)
    {
      DUNE_THROW(Dune::GridError, "The requested codim is not in {0,..,dim}!");
      return 0;
    }
  };

} // end namespace Dune

#endif
