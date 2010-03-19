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
    IntersectionIteratorWrapper(const GridImp & grid , const EntityImp & en, int wLevel , bool end)
      : grid_( grid )
        , it_( grid_.getIntersection(wLevel, (IntersectionIteratorImpl *) 0) )
    {
      if(end)
        it().done( en );
      else
        it().first(en,wLevel);
    }

    //! The copy constructor
    IntersectionIteratorWrapper(const ThisType & org)
      : grid_( org.grid_ )
        , it_( grid_.getIntersectionCopy( org.it() ) )
    {
      it().assign(org.it() );
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
      grid_.freeIntersection( it() );
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

    //! return the boundary segment index
    size_t boundarySegmentIndex() const { return it().boundarySegmentIndex(); }

    //! intersection of codimension 1 of this neighbor with element where
    //! iteration started.
    //! Here returned element is in LOCAL coordinates of the element
    //! where iteration started.
    const LocalGeometry &geometryInInside () const
    {
      return it().geometryInInside();
    }

    //! intersection of codimension 1 of this neighbor with element where
    //!  iteration started.
    //! Here returned element is in GLOBAL coordinates of the element where
    //! iteration started.
    const Geometry &geometry () const
    {
      return it().geometry();
    }

    /** \brief obtain the type of reference element for this intersection */
    GeometryType type () const
    {
      return it().type();
    }

    //! local index of codim 1 entity in self where intersection is contained
    //!  in
    int indexInInside () const
    {
      return it().indexInInside();
    }

    //! intersection of codimension 1 of this neighbor with element where
    //! iteration started.
    //! Here returned element is in LOCAL coordinates of neighbor
    const LocalGeometry &geometryInOutside () const
    {
      return it().geometryInOutside();
    }

    //! local index of codim 1 entity in neighbor where intersection is
    //! contained
    int indexInOutside () const
    {
      return it().indexInOutside();
    }

    //! twist of the face seen from the inner element
    int twistInSelf() const { return it().twistInSelf(); }

    //! twist of the face seen from the inner element
    int twistInInside() const { return it().twistInInside(); }

    //! twist of the face seen from the outer element
    int twistInNeighbor() const { return it().twistInNeighbor(); }

    //! twist of the face seen from the outer element
    int twistInOutside() const { return it().twistInOutside(); }

    //! return unit outer normal, this should be dependent on local
    //! coordinates for higher order boundary
    const NormalType unitOuterNormal ( const FieldVector< ctype, dim-1 > &local ) const
    {
      return it().unitOuterNormal( local );
    }

    //! return unit outer normal, this should be dependent on local
    //! coordinates for higher order boundary
    const NormalType centerUnitOuterNormal ( ) const
    {
      GeometryType type = geometry().type();
      const GenericReferenceElement<ctype, dim-1> & refElement =
        GenericReferenceElements<ctype, dim-1>::general(type);
      return unitOuterNormal(refElement.position(0,0));
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

    const GridImp& grid_;
    IntersectionIteratorImp & it_;
  }; // end class IntersectionIteratorWrapper

  template <class GridImp>
  class LeafIntersectionWrapper
    : public IntersectionIteratorWrapper<GridImp,typename GridImp::LeafIntersectionIteratorImp>
  {
    typedef LeafIntersectionWrapper<GridImp> ThisType;
    typedef IntersectionIteratorWrapper<GridImp,typename GridImp::LeafIntersectionIteratorImp> BaseType;
  public:
    //! constructor called from the ibegin and iend method
    template <class EntityImp>
    LeafIntersectionWrapper(const GridImp & grid , const EntityImp & en, int wLevel , bool end )
      : BaseType(grid,en,wLevel,end)
    {}

    //! The copy constructor
    LeafIntersectionWrapper(const ThisType & org)
      : BaseType(org)
    {}

  };

  //! \brief Class that wraps IntersectionIteratorImp of a grid and gets it's
  //! internal object from a object stack hold by the grid
  template <class GridImp>
  class LeafIntersectionIteratorWrapper
  {
    typedef LeafIntersectionIteratorWrapper<GridImp> ThisType;
    typedef LeafIntersectionWrapper<GridImp> IntersectionImp;
  public:
    typedef Dune :: Intersection
    < const GridImp, Dune :: LeafIntersectionWrapper > Intersection;

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
      : intersection_( IntersectionImp(grid,en,wLevel,end) )
    {}

    //! The copy constructor
    LeafIntersectionIteratorWrapper(const ThisType & org)
      : intersection_( org.intersection_ )
    {}

    //! the f*cking assignment operator
    ThisType & operator = (const ThisType & org)
    {
      impl() = org.impl();
      return *this;
    }

    //! return reference to intersection
    const Intersection &dereference () const
    {
      return intersection_;
    }

    //! the equality method
    bool equals (const ThisType & i) const { return impl().equals( i.impl() ); }

    //! increment iterator
    void increment()
    {
      impl().increment();
    }
  protected:
    // intersection object
    Intersection intersection_;

    // return reference to real implementation
    IntersectionImp& impl() { return GridImp :: getRealImplementation( intersection_ ); }
    // return reference to real implementation
    const IntersectionImp& impl() const { return GridImp :: getRealImplementation( intersection_ ); }
  }; // end class IntersectionIteratorWrapper

  //! \brief Class that wraps IntersectionIteratorImp of a grid and gets it's
  //! internal object from a object stack hold by the grid
  template <class GridImp>
  class LevelIntersectionWrapper
    : public IntersectionIteratorWrapper<GridImp,typename GridImp::LevelIntersectionIteratorImp>
  {
    typedef LevelIntersectionWrapper<GridImp> ThisType;
    typedef IntersectionIteratorWrapper<GridImp,typename GridImp::LevelIntersectionIteratorImp> BaseType;
  public:
    //! constructor called from the ibegin and iend method
    template <class EntityImp>
    LevelIntersectionWrapper(const GridImp & grid , const EntityImp & en, int wLevel , bool end )
      : BaseType(grid,en,wLevel,end)
    {}

    //! The copy constructor
    LevelIntersectionWrapper(const ThisType & org)
      : BaseType(org)
    {}
  };

  //! \brief Class that wraps IntersectionIteratorImp of a grid and gets it's
  //! internal object from a object stack hold by the grid
  template <class GridImp>
  class LevelIntersectionIteratorWrapper
  {
    typedef LevelIntersectionIteratorWrapper<GridImp> ThisType;
    typedef LevelIntersectionWrapper<GridImp> IntersectionImp;
  public:
    typedef Dune :: Intersection
    < const GridImp, Dune :: LevelIntersectionWrapper >
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
      : intersection_( IntersectionImp(grid,en,wLevel,end) )
    {}

    //! The copy constructor
    LevelIntersectionIteratorWrapper(const ThisType & org)
      : intersection_( org.intersection_ )
    {}

    //! the f*cking assignment operator
    ThisType & operator = (const ThisType & org)
    {
      impl() = org.impl();
      return *this;
    }

    //! return reference to intersection
    const Intersection &dereference () const
    {
      return intersection_;
    }

    //! the equality method
    bool equals (const ThisType & i) const { return impl().equals( i.impl() ); }

    //! increment iterator
    void increment()
    {
      impl().increment();
    }
  protected:
    // intersection object
    Intersection intersection_;

    // return reference to real implementation
    IntersectionImp& impl() { return GridImp :: getRealImplementation( intersection_ ); }
    // return reference to real implementation
    const IntersectionImp& impl() const { return GridImp :: getRealImplementation( intersection_ ); }
  }; // end class IntersectionIteratorWrapper

} // end namespace Dune
#endif
