// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_INTERSECTIONITERATORWRAPPER_HH
#define DUNE_INTERSECTIONITERATORWRAPPER_HH

#include "intersectioniterator.hh"

/** @file
   @author Robert Kloefkorn
   @brief Provides proxy classes for IntersectionsIterators
 */

namespace Dune {

  //! \brief Class that wraps IntersectionIteratorImp of a grid and gets it's
  //! internal object from a object stack hold by the grid
  template <class GridImp, class IntersectionIteratorImpl>
  class IntersectionIteratorWrapper
  {
    enum { dim = GridImp :: dimension };
    enum { dimworld = GridImp :: dimensionworld };

    typedef IntersectionIteratorWrapper<GridImp,IntersectionIteratorImpl> ThisType;

    typedef IntersectionIteratorImpl IntersectionIteratorImp;

    typedef typename IntersectionIteratorImp :: StorageType IntersectionIteratorProviderType;

  public:
    //! dimension
    enum { dimension      = dim };
    //! dimensionworld
    enum { dimensionworld = dimworld };

    //! define type used for coordinates in grid module
    typedef typename GridImp :: ctype ctype;

    //! Entity type
    typedef typename GridImp::template Codim<0>::Entity Entity;
    //! type of EntityPointer
    typedef typename GridImp::template Codim<0>::EntityPointer EntityPointer;

    //! type of intersectionGlobal
    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    //! type of intersection*Local
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;

    //! type of normal vector
    typedef FieldVector<ctype , dimworld> NormalType;

    //! constructor called from the ibegin and iend method
    template <class EntityImp>
    IntersectionIteratorWrapper(const GridImp & grid , const EntityImp & en, int wLevel , bool end ,
                                IntersectionIteratorProviderType & storage)
      : storage_(storage)
        , it_(*(storage_.getObject(grid,wLevel)))
    {
      if(end)
        it().done();
      else
        it().first(en,wLevel);
    }

    //! The copy constructor
    IntersectionIteratorWrapper(const ThisType & org)
      : storage_(org.storage_) , it_(*(storage_.getObjectCopy(org.it_)))
    {
      it().assign(org.it_);
    }

    //! the f*cking assignment operator
    ThisType & operator = (const ThisType & org)
    {
      it().assign(org.it_);
      return *this;
    }

    //! The Destructor puts internal object back to stack
    ~IntersectionIteratorWrapper()
    {
      storage_.freeObject( &it_ );
    }

    //! the equality method
    bool equals (const ThisType & i) const { return it().equals(i.it()); }

    //! increment iterator
    void increment () { it().increment(); }

    //! access neighbor
    EntityPointer outside() const { return it().outside(); }

    //! access entity where iteration started
    EntityPointer inside() const { return it().inside(); }

    //! return true if intersection is with boundary. \todo connection with
    //! boundary information, processor/outer boundary
    bool boundary () const { return it().boundary(); }

    //! return true if across the intersection a neighbor on this level exists
    bool neighbor () const { return it().neighbor(); }

    //! return information about the Boundary
    int boundaryId () const { return it().boundaryId(); }

    //! intersection of codimension 1 of this neighbor with element where
    //! iteration started.
    //! Here returned element is in LOCAL coordinates of the element
    //! where iteration started.
    const LocalGeometry & intersectionSelfLocal () const { return it().intersectionSelfLocal(); }

    //! intersection of codimension 1 of this neighbor with element where
    //!  iteration started.
    //! Here returned element is in GLOBAL coordinates of the element where
    //! iteration started.
    const Geometry & intersectionGlobal () const { return it().intersectionGlobal(); }

    /** \brief obtain the type of reference element for this intersection */
    GeometryType type () const
    {
      return it().type();
    }

    //! local number of codim 1 entity in self where intersection is contained
    //!  in
    int numberInSelf () const { return it().numberInSelf(); }

    //! intersection of codimension 1 of this neighbor with element where
    //! iteration started.
    //! Here returned element is in LOCAL coordinates of neighbor
    const LocalGeometry & intersectionNeighborLocal () const { return it().intersectionNeighborLocal(); }

    //! local number of codim 1 entity in neighbor where intersection is
    //! contained
    int numberInNeighbor () const { return it().numberInNeighbor(); }

    //! twist of the face seen from the inner element
    int twistInSelf() const { return it().twistInSelf(); }

    //! twist of the face seen from the outer element
    int twistInNeighbor() const { return it().twistInNeighbor(); }

    //! return unit outer normal, this should be dependent on local
    //! coordinates for higher order boundary
    const NormalType unitOuterNormal ( const FieldVector< ctype, dim-1 > &local ) const
    {
      return it().unitOuterNormal( local );
    }

    //! return outer normal, this should be dependent on local
    //! coordinates for higher order boundary
    const NormalType outerNormal ( const FieldVector< ctype, dim-1 > &local ) const
    {
      return it().outerNormal( local );
    }

    //! return outer normal, this should be dependent on local
    //! coordinates for higher order boundary
    const NormalType integrationOuterNormal ( const FieldVector< ctype, dim-1 > &local ) const
    {
      return it().integrationOuterNormal( local );
    }

    //! return level of iterator
    int level () const { return it().level(); }

    //! return true if intersection is conform (i.e. only one neighbor)
    bool conforming () const { return it().conforming(); }

  private:
    //! returns reference to underlying intersection iterator implementation
    IntersectionIteratorImp & it() { return it_; }
    const IntersectionIteratorImp & it() const { return it_; }

    IntersectionIteratorProviderType & storage_;
    IntersectionIteratorImp & it_;
  }; // end class IntersectionIteratorWrapper

  //! \brief Class that wraps IntersectionIteratorImp of a grid and gets it's
  //! internal object from a object stack hold by the grid
  template <class GridImp>
  class LeafIntersectionIteratorWrapper
    : public IntersectionIteratorWrapper<GridImp,typename GridImp::LeafIntersectionIteratorImp>
  {
    typedef LeafIntersectionIteratorWrapper<GridImp> ThisType;
    typedef IntersectionIteratorWrapper<GridImp,typename GridImp::LeafIntersectionIteratorImp> BaseType;
  public:
    typedef Dune :: Intersection
    < const GridImp, Dune :: LeafIntersectionIteratorWrapper >
    Intersection;

    //! dimension
    enum { dimension      = GridImp :: dimension  };
    //! dimensionworld
    enum { dimensionworld = GridImp :: dimensionworld };

    //! define type used for coordinates in grid module
    typedef typename GridImp :: ctype ctype;

    //! Entity type
    typedef typename GridImp::template Codim<0>::Entity Entity;
    //! type of EntityPointer
    typedef typename GridImp::template Codim<0>::EntityPointer EntityPointer;

    //! type of intersectionGlobal
    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    //! type of intersection*Local
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;

    //! type of normal vector
    typedef FieldVector<ctype , dimensionworld> NormalType;

    //! constructor called from the ibegin and iend method
    template <class EntityImp>
    LeafIntersectionIteratorWrapper(const GridImp & grid , const EntityImp & en, int wLevel , bool end )
      : BaseType(grid,en,wLevel,end,grid.leafIntersetionIteratorProvider())
    {}

    //! The copy constructor
    LeafIntersectionIteratorWrapper(const ThisType & org)
      : BaseType(org)
    {}

    //! the f*cking assignment operator
    ThisType & operator = (const ThisType & org)
    {
      BaseType::operator = (org);
      return *this;
    }

    const Intersection &dereference () const
    {
      return reinterpret_cast< const Intersection & >( *this );
    }
  }; // end class IntersectionIteratorWrapper

  //! \brief Class that wraps IntersectionIteratorImp of a grid and gets it's
  //! internal object from a object stack hold by the grid
  template <class GridImp>
  class LevelIntersectionIteratorWrapper
    : public IntersectionIteratorWrapper<GridImp,typename GridImp::LevelIntersectionIteratorImp>
  {
    typedef LevelIntersectionIteratorWrapper<GridImp> ThisType;
    typedef IntersectionIteratorWrapper<GridImp,typename GridImp::LevelIntersectionIteratorImp> BaseType;
  public:
    typedef Dune :: Intersection
    < const GridImp, Dune :: LevelIntersectionIteratorWrapper >
    Intersection;

    //! dimension
    enum { dimension      = GridImp :: dimension  };
    //! dimensionworld
    enum { dimensionworld = GridImp :: dimensionworld };

    //! define type used for coordinates in grid module
    typedef typename GridImp :: ctype ctype;

    //! Entity type
    typedef typename GridImp::template Codim<0>::Entity Entity;
    //! type of EntityPointer
    typedef typename GridImp::template Codim<0>::EntityPointer EntityPointer;

    //! type of intersectionGlobal
    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    //! type of intersection*Local
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;

    //! type of normal vector
    typedef FieldVector<ctype , dimensionworld> NormalType;

    //! constructor called from the ibegin and iend method
    template <class EntityImp>
    LevelIntersectionIteratorWrapper(const GridImp & grid , const EntityImp & en, int wLevel , bool end )
      : BaseType(grid,en,wLevel,end,grid.levelIntersetionIteratorProvider())
    {}

    //! The copy constructor
    LevelIntersectionIteratorWrapper(const ThisType & org)
      : BaseType(org)
    {}

    //! the f*cking assignment operator
    ThisType & operator = (const ThisType & org)
    {
      BaseType::operator = (org);
      return *this;
    }

    const Intersection &dereference () const
    {
      return reinterpret_cast< const Intersection & >( *this );
    }
  }; // end class IntersectionIteratorWrapper

} // end namespace Dune
#endif
