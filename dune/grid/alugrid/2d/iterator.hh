// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU2DGRID_ITERATOR_HH
#define DUNE_ALU2DGRID_ITERATOR_HH

// System includes
#include <stack>
#include <utility>

// Dune includes
#include <dune/grid/common/grid.hh>

// Local includes
#include "entity.hh"

namespace Dune {
  // Forward declarations
  template<int cd, int dim, class GridImp>
  class ALU2dGridEntity;
  template<int cd, PartitionIteratorType pitype, class GridImp >
  class ALU2dGridLevelIterator;
  template<int cd, class GridImp >
  class ALU2dGridEntityPointer;
  template<int mydim, int coorddim, class GridImp>
  class ALU2dGridGeometry;
  template<class GridImp>
  class ALU2dGridHierarchicIterator;
  template<int codim, PartitionIteratorType pitype, class GridImp>
  class ALU2dGridLeafIterator;
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  class ALU2dGrid;


  //********************************************************************
  //  --ALU2dGridLeafIterator
  //  --LeafIterator
  //  --for codim = 0,2
  //
  //********************************************************************

  template<int cdim, PartitionIteratorType pitype, class GridImp>
  class ALU2dGridLeafIterator
    : public ALU2dGridEntityPointer<cdim,GridImp>
      // public LeafIteratorDefaultImplementation<cdim, pitype, GridImp, ALU2dGridLeafIterator>
  {
    static const int dim = GridImp::dimension;
    static const int dimworld  = GridImp::dimensionworld;
    static const ALU2DSPACE ElementType eltype = GridImp::elementType;
    static const int codim = cdim;

    friend class ALU2dGridEntity<0,dimworld,GridImp>;
    friend class ALU2dGridEntity<1,dimworld,GridImp>;
    friend class ALU2dGridEntity<dim,dimworld,GridImp>;
    friend class ALU2dGrid< dim, dimworld, eltype >;

    typedef ALU2dGridEntityPointer<cdim,GridImp> EntityPointerType;
    typedef ALU2dGridEntity<cdim,dim,GridImp> EntityImp;

    typedef ALU2dGridLeafIterator<cdim, pitype, GridImp> ThisType;

    typedef typename GridImp :: ALU2dGridLeafMarkerVectorType LeafMarkerVectorType;

    // default impl for elements
    template <class ElementImp, class MarkerVectorImp, int codim>
    struct GetLevel
    {
      // return level of element
      static int level(const ElementImp & elem, const MarkerVectorImp& marker)
      {
        return elem.level();
      }
    };

    // specialization for vertices
    template <class ElementImp, class MarkerVectorImp>
    struct GetLevel<ElementImp,MarkerVectorImp,2>
    {
      // return level of leaf vertex
      static int level(const ElementImp & elem, const MarkerVectorImp& marker)
      {
        return marker.levelOfVertex(elem.getIndex());
      }
    };

  public:
    typedef typename GridImp :: GridObjectFactoryType FactoryType;

    //! type of entity we iterate (interface)
    typedef typename GridImp::template Codim<cdim>::Entity Entity;
    typedef typename Dune::ALU2dImplTraits< dimworld, eltype >::template Codim<cdim>::InterfaceType ElementType;

    //! Constructor called by LeafIterator
    ALU2dGridLeafIterator(const FactoryType& factory, bool end);

    //! copy Constructor
    ALU2dGridLeafIterator(const ThisType & org);

    //! prefix increment
    void increment ();

    //! assigment of iterator
    ThisType & operator = (const ThisType & org);

  private:
    //! true if iterator is end iterator
    bool endIter_;
    //! actual level
    int level_;
    ElementType * elem_;
    typedef ALU2DSPACE Listwalkptr< ElementType > IteratorType;
    // Listwalkptr, behaves like a proxy for Leafwalk and Levelwalk Ptrs
    IteratorType iter_;

    // for the codim 2 case
    LeafMarkerVectorType & marker_;
  }; // end ALU2dGridLeafIterator


  //********************************************************************
  //  --ALU2dGridLeafIterator
  //  --LeafIterator
  //  --specialized for codim = 1
  //
  //********************************************************************

  template<PartitionIteratorType pitype, class GridImp>
  class ALU2dGridLeafIterator<1,pitype,GridImp>
    : public ALU2dGridEntityPointer<1,GridImp>
      // public LeafIteratorDefaultImplementation<1, pitype, GridImp, ALU2dGridLeafIterator>
  {
    static const int dim = GridImp::dimension;
    static const int dimworld  = GridImp::dimensionworld;
    static const ALU2DSPACE ElementType eltype = GridImp::elementType;
    static const int codim = 1;

    friend class ALU2dGridEntity<0,dimworld,GridImp>;
    friend class ALU2dGridEntity<1,dimworld,GridImp>;
    friend class ALU2dGridEntity<dim,dimworld,GridImp>;
    friend class ALU2dGrid< dim, dimworld, eltype >;

    typedef ALU2dGridEntityPointer<1,GridImp> EntityPointerType;
    typedef ALU2dGridEntity<1,dim,GridImp> EntityImp;

    typedef ALU2dGridLeafIterator<1, pitype, GridImp> ThisType;

    typedef typename GridImp :: ALU2dGridLeafMarkerVectorType LeafMarkerVectorType;

  public:
    typedef typename GridImp :: GridObjectFactoryType FactoryType;

    //! type of entity we iterate (interface)
    typedef typename GridImp::template Codim<1>::Entity Entity;
    typedef typename Dune::ALU2dImplTraits< dimworld, eltype >::template Codim<1>::InterfaceType ElementType;

    //! Constructor called by LeafIterator
    ALU2dGridLeafIterator(const FactoryType& factory, bool end);

    //! copy Constructor
    ALU2dGridLeafIterator(const ThisType & org);

    //! prefix increment
    void increment ();

    //! assigment of iterator
    ThisType & operator = (const ThisType & org);

  private:
    int goNextElement();

    //! true if iterator is end iterator
    bool endIter_;
    //! actual level
    int level_;
    //! information for edges
    int face_;

    //! pointer to element
    ElementType * elem_;

    typedef ALU2DSPACE Listwalkptr< ElementType > IteratorType;

    // Listwalkptr, behaves like a proxy for Leafwalk and Levelwalk Ptrs
    IteratorType iter_;

    // for the codim 1 case
    LeafMarkerVectorType & marker_;

  }; // end ALU2dGridLeafIterator

  //**********************************************************************
  //
  // --ALU2dGridLevelIterator
  // --LevelIterator, specialized for cd=0
  //**********************************************************************

  template<PartitionIteratorType pitype, class GridImp>
  class ALU2dGridLevelIterator<0, pitype, GridImp>
    : public ALU2dGridEntityPointer<0,GridImp>
      // public LevelIteratorDefaultImplementation <0, pitype, GridImp, ALU2dGridLevelIterator>
  {
    static const int dim = GridImp::dimension;
    static const int dimworld  = GridImp::dimensionworld;
    static const ALU2DSPACE ElementType eltype = GridImp::elementType;
    static const int codim = 0;

    friend class ALU2dGridEntity<dim,dimworld,GridImp>;
    friend class ALU2dGridEntity<1,dimworld,GridImp>;
    friend class ALU2dGridEntity<0,dimworld,GridImp>;
    friend class ALU2dGrid< dim, dimworld, eltype >;

    typedef ALU2dGridEntityPointer<codim,GridImp> EntityPointerType;
    typedef ALU2dGridEntity<codim,dim,GridImp> EntityImp;

    typedef typename ALU2dImplTraits< dimworld, eltype >::HElementType HElementType;
    typedef ALU2dGridLevelIterator<0,pitype,GridImp> ThisType;

  public:
    typedef typename GridImp :: GridObjectFactoryType FactoryType;

    typedef typename GridImp::template Codim<codim>::Entity Entity;

    //! Constructor
    ALU2dGridLevelIterator(const FactoryType& factory, int level, bool end);

    //! copy constructor
    ALU2dGridLevelIterator(const ThisType & org);

    //! prefix increment
    void increment ();

    //! do not allow assigment
    ThisType & operator = (const ThisType & org);

  private:
    //! true if iterator is end iterator
    bool endIter_;
    //! actual level
    int level_;

    //! pointer to element
    HElementType * item_;

    //! type of entity we iterate (interface)
    typedef typename Dune::ALU2dImplTraits< dimworld, eltype  >::template Codim<0>::InterfaceType ElementType;
    typedef ALU2DSPACE Listwalkptr< ElementType > IteratorType;
    IteratorType iter_;

  };

  //**********************************************************************
  //
  // --ALU2dGridLevelIterator
  // --LevelIterator, specialized for cd=1
  //**********************************************************************

  template<PartitionIteratorType pitype, class GridImp>
  class ALU2dGridLevelIterator<1, pitype, GridImp>
    : public ALU2dGridEntityPointer<1,GridImp>
      // public LevelIteratorDefaultImplementation <1, pitype, GridImp, ALU2dGridLevelIterator>
  {
    static const int dim = GridImp::dimension;
    static const int dimworld  = GridImp::dimensionworld;
    static const ALU2DSPACE ElementType eltype = GridImp::elementType;
    static const int codim = 1;

    friend class ALU2dGridEntity<dim,dimworld,GridImp>;
    friend class ALU2dGridEntity<1,dimworld,GridImp>;
    friend class ALU2dGridEntity<0,dimworld,GridImp>;
    friend class ALU2dGrid< dim, dimworld, eltype >;

    typedef ALU2dGridEntityPointer<codim,GridImp> EntityPointerType;
    typedef ALU2dGridEntity<codim,dim,GridImp> EntityImp;

    typedef typename ALU2dImplTraits< dimworld, eltype >::HElementType HElementType;

    typedef ALU2dGridLevelIterator<1,pitype,GridImp> ThisType;
  public:
    typedef typename GridImp :: GridObjectFactoryType FactoryType;

    typedef typename GridImp::template Codim<codim>::Entity Entity;

    //! Constructor
    ALU2dGridLevelIterator(const FactoryType& factroy, int level, bool end);

    //! copy constructor
    ALU2dGridLevelIterator(const ThisType & org);

    ~ALU2dGridLevelIterator();

    //! prefix increment
    void increment ();

    //! assigment of iterator
    ThisType & operator = (const ThisType & org);

  private:
    //! true if iterator is end iterator
    bool endIter_;
    //! actual level
    int level_;
    //! information for edges
    int myFace_;

    // current item
    HElementType * item_;
    //! type of entity we iterate (interface)
    typedef typename Dune::ALU2dImplTraits< dimworld, eltype >::template Codim<1>::InterfaceType ElementType;
    ElementType * elem_;
    // Listwalkptr is a proxy for iterator pointers
    typedef ALU2DSPACE Listwalkptr< ElementType > IteratorType;

    IteratorType iter_;

    ALU2dGridMarkerVector * marker_;

    ALU2dGridMarkerVector & marker()
    {
      assert( marker_ );
      return *marker_;
    }
  };

  //**********************************************************************
  //
  // --ALU2dGridLevelIterator
  // --LevelIterator, specialized for cd=2
  //**********************************************************************

  template<PartitionIteratorType pitype, class GridImp>
  class ALU2dGridLevelIterator<2, pitype, GridImp>
    : public ALU2dGridEntityPointer<2,GridImp>
      // public LevelIteratorDefaultImplementation <2, pitype, GridImp, ALU2dGridLevelIterator>
  {
    static const int dim = GridImp::dimension;
    static const int dimworld  = GridImp::dimensionworld;
    static const ALU2DSPACE ElementType eltype = GridImp::elementType;
    static const int codim = 2;

    friend class ALU2dGridEntity<dim,dimworld,GridImp>;
    friend class ALU2dGridEntity<1,dimworld,GridImp>;
    friend class ALU2dGridEntity<0,dimworld,GridImp>;
    friend class ALU2dGrid< dim, dimworld, eltype >;

    typedef ALU2dGridEntityPointer<codim,GridImp> EntityPointerType;
    typedef ALU2dGridEntity<codim,dim,GridImp> EntityImp;

    typedef typename ALU2dImplTraits< dimworld, eltype >::HElementType HElementType ;
    typedef ALU2dGridLevelIterator<2,pitype,GridImp> ThisType;

  public:
    typedef typename GridImp :: GridObjectFactoryType FactoryType;

    typedef typename GridImp::template Codim<codim>::Entity Entity;

    //! Constructor
    ALU2dGridLevelIterator(const FactoryType& factory, int level, bool end);

    //! copy constructor
    ALU2dGridLevelIterator(const ThisType & org);

    ~ALU2dGridLevelIterator();

    //! prefix increment
    void increment ();

    //! assigment of iterator
    ThisType & operator = (const ThisType & org);

  private:
    //! true if iterator is end iterator
    bool endIter_;
    //! actual level
    int level_;
    //! information for edges
    int myFace_;
    //! true if iterator is already a copy

    //! type of entity we iterate (interface)
    typedef typename Dune::ALU2dImplTraits< dimworld, eltype >::template Codim<0>::InterfaceType ElementType;
    typedef typename Dune::ALU2dImplTraits< dimworld, eltype >::template Codim<2>::InterfaceType VertexType;

    HElementType * item_;
    VertexType * vertex_;


    typedef ALU2DSPACE Listwalkptr< ElementType > IteratorType;
    IteratorType iter_;

    // marker vector to tell on which element vertex is visited
    ALU2dGridMarkerVector * marker_;
    ALU2dGridMarkerVector & marker()
    {
      assert( marker_ );
      return *marker_;
    }
  };

  //***************************************************************
  //
  // - HierarchicIteraror
  // --HierarchicIterator
  //***************************************************************

  //! Hierarchic Iterator of ALU2dGrid
  template<class GridImp>
  class ALU2dGridHierarchicIterator
    : public ALU2dGridEntityPointer<0,GridImp>
      // public HierarchicIteratorDefaultImplementation <GridImp,ALU2dGridHierarchicIterator>
  {
    typedef ALU2dGridHierarchicIterator<GridImp> ThisType;

    static const int dim = GridImp::dimension;
    static const int dimworld  = GridImp::dimensionworld;
    static const ALU2DSPACE ElementType eltype = GridImp::elementType;
    static const int codim = 2;

    typedef typename ALU2dImplTraits< dimworld, eltype >::HElementType HElementType ;

  public:
    typedef typename GridImp :: GridObjectFactoryType FactoryType;

    //! type of entities we iterate
    typedef typename GridImp::template Codim<0>::Entity Entity;
    //! type of coordinates, i.e. double
    typedef typename GridImp::ctype ctype;
    //! tpye of entity implementation
    typedef ALU2dGridEntity<0,dim,GridImp> EntityImp;

    //! the normal Constructor
    ALU2dGridHierarchicIterator(const FactoryType& factory,
                                const HElementType & elem, int maxlevel, bool end=false);

    //! the normal Constructor
    ALU2dGridHierarchicIterator(const ALU2dGridHierarchicIterator<GridImp> &org);

    //! increment, go to next entity
    void increment();

    //! the assignment operator
    ThisType & operator = (const ALU2dGridHierarchicIterator<GridImp> &org)
    {
      ALU2dGridEntityPointer<0,GridImp> :: operator = (org);
      elem_ = org.elem_;
      maxlevel_= org.maxlevel_;
      endIter_ = org.endIter_;
      return *this;
    }

  private:

    //! go to next valid element
    HElementType * goNextElement (HElementType * oldEl);

    //! element from where we started
    const HElementType * elem_;

    //! maximal level to go down
    int maxlevel_;
    //! true if iterator is end iterator
    bool endIter_;
  }; // end ALU2dHierarchicIterator

} // end namespace Dune

#include "iterator_imp.cc"

#endif // #ifndef DUNE_ALU2DGRID_ITERATOR_HH
