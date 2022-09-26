// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_COMMON_GRIDVIEW_HH
#define DUNE_GRID_COMMON_GRIDVIEW_HH

#include <typeinfo>

#include <dune/common/std/type_traits.hh>
#include <dune/common/iteratorrange.hh>
#include <dune/common/parallel/future.hh>

#include <dune/geometry/type.hh>

#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/rangegenerators.hh>

namespace Dune
{

  template< int, int, class, class >
  class GridDefaultImplementation;



  /** \addtogroup GIGridView
   *  @{
   *  Though a DUNE grid is hierarchic, one often only needs access to
   *  a certain subset of the entities in the grid, e.g., all entities
   *  on a given level or the leaf entities in the hierarchy.
   *  These views are provided by an
   *  implementation of GridView. Each grid exports a LevelGridView and
   *  a LeafGridView, corresponding to the two different subsets (views)
   *  described above.
   *
   *  A grid view provides the following functionality:
   *  - The index set for the required subset can be accessed by the indexSet()
   *    method.
   *  - A pair of begin() / end() methods provide iterators for each
   *    codimension.
   *  - A pair of ibegin() / iend() methods return suitable intersection
   *    iterators for a given entity of codimension 0 in the subset.
   *  - For parallel computations, a suitable communicate() method is provided.
   *  - The underlying grid can be accessed through the grid() method.
   *  .
   *
   *  The default GridViews can be obtained from the grid by calling one of the
   *  levelGridView() or leafGridView() methods which are provided as both
   *  free functions and member functions of the Grid interface class.
   */


  /** \brief Grid view abstract base class
   *
   *  Interface class for a view on grids. Grids return two types of view,
   *  a view of the leaf grid and of a level grid, which both satisfy
   *  the same interface. Through the view the user has access to the
   *  iterators, the intersections and the index set.
   *
   *  The interface is implemented using the engine concept.
   *
   */
  template< class ViewTraits >
  class GridView
  {
    typedef GridView< ViewTraits > ThisType;

  public:
    /**
     * \brief type of underlying implementation
     *
     * \warning Implementation details may change without prior notification.
     **/
    typedef typename ViewTraits :: GridViewImp Implementation;

    typedef typename ViewTraits :: GridViewImp GridViewImp;

    /** \brief Traits class */
    typedef ViewTraits Traits;

    /** \brief type of the grid */
    typedef typename Traits :: Grid Grid;

    /** \brief type of the index set */
    typedef typename Traits :: IndexSet IndexSet;

    /** \brief type of the intersection */
    typedef typename Traits :: Intersection Intersection;

    /** \brief type of the intersection iterator */
    typedef typename Traits :: IntersectionIterator IntersectionIterator;

  protected:
    template <class T>
    using Communication_t = typename T :: Communication;
    template <class T>
    using DeprecatedCollectiveCommunication_t = typename T :: CollectiveCommunication;

  public:
    /*! \brief A type that is a model of Dune::Communication.
       It provides a portable way for communication on the set
       of processes used by the grid.
     */
    // if this line produces a warning then the Communication typedef is missing
    // in the Traits
    using Communication = detected_or_fallback_t<DeprecatedCollectiveCommunication_t,
                                                 Communication_t, Traits>;

    /** \deprecated Use Communication instead! Will be removed after Dune 2.9.
     */
    using CollectiveCommunication [[deprecated("CollectiveCommunication is deprecated, use Communication instead!")]] = Communication;

    /** \brief A struct that collects all associated types of one implementation
               from the Traits class.
     */
    template< int cd >
    struct Codim {
      /** \brief type of iterator returned by the grid view */
      typedef typename Traits :: template Codim<cd> :: Iterator Iterator;

      /** \brief type of corresponding entity */
      typedef typename Traits :: template Codim<cd> :: Entity Entity;

      /** \brief type of the geometry implementation */
      typedef typename Traits :: template Codim<cd> :: Geometry Geometry;

      /** \brief type of the implementation for local geometries */
      typedef typename Traits :: template Codim<cd> :: LocalGeometry LocalGeometry;

      /** \brief Define types needed to iterate over entities of a given partition type */
      template< PartitionIteratorType pit >
      struct Partition
      {
        /** \brief iterator over a given codim and partition type */
        typedef typename Traits :: template Codim< cd >
        :: template Partition< pit > :: Iterator Iterator;
      };
    }; //: public Traits :: template Codim<cd> {};

    /** \brief Export if this grid view is guaranteed conforming */
    constexpr static bool conforming = Traits :: conforming;

    /** \brief type used for coordinates in grid */
    typedef typename Grid::ctype ctype;

    //! \brief The dimension of the grid
    constexpr static int dimension = Grid :: dimension;

    //! \brief The dimension of the world the grid lives in
    constexpr static int dimensionworld = Grid :: dimensionworld;

  public:

    //===========================================================
    /** @name Interface for the implementor
     */
    //@{
    //===========================================================
    /** \brief constructor (engine concept) */
    GridView ( const Implementation &imp )
      : impl_( imp )
    {}
    //@}

    /** \brief Copy constructor */
    GridView ( const ThisType &other )
      : impl_( other.impl_ )
    {}

    /** \brief assignment operator */
    ThisType &operator= ( const ThisType &other )
    {
      impl_ = other.impl_;
      return *this;
    }

  public:
    /** \brief obtain a const reference to the underlying hierarchic grid */
    const Grid &grid () const
    {
      return impl().grid();
    }

    /** \brief obtain the index set
     *
     * The lifetime of the returned index set is bound to the lifetime of the
     * grid view. Keep a copy of the grid view to prevent the index set from
     * becoming a dangling reference.
     */
    const IndexSet &indexSet () const
    {
      return impl().indexSet();
    }

    /** \brief obtain number of entities in a given codimension */
    int size ( int codim ) const
    {
      return impl().size( codim );
    }

    /** \brief obtain number of entities with a given geometry type */
    int size ( const GeometryType &type ) const
    {
      return impl().size( type );
    }

    /** \brief return true if current state of grid view represents a conforming grid */
    bool isConforming () const
    {
      // if implementation provides a method isConforming, call it
      if constexpr ( CheckIsConformingImpl< Implementation >::type::value )
      {
        return impl().isConforming();
      }
      else
      {
        // otherwise default to static conforming flag
        return this->isConformingDefaultImplementation();
      }
    }

    /** @brief Return true if the given entity is contained in this grid view
     * @todo Currently we call the implementation on the IndexSet.  This may lead to suboptimal efficiency.
     *
     * \note If the input element e is not an element of the grid, then
     *       the result of contains() is undefined.
     */
    template<class EntityType>
    bool contains (const EntityType& e) const
    {
      return impl().indexSet().contains(e);
    }

    /** \brief obtain begin iterator for this view */
    template< int cd >
    typename Codim< cd > :: Iterator begin () const
    {
      return impl().template begin<cd>();
    }

    /** \brief obtain end iterator for this view */
    template< int cd >
    typename Codim< cd > :: Iterator end () const
    {
      return impl().template end<cd>();
    }

    /** \brief obtain begin iterator for this view */
    template< int cd , PartitionIteratorType pitype >
    typename Codim< cd > :: template Partition< pitype > :: Iterator
    begin () const
    {
      return impl().template begin<cd,pitype>();
    }

    /** \brief obtain end iterator for this view */
    template< int cd, PartitionIteratorType pitype >
    typename Codim< cd > :: template Partition< pitype > :: Iterator
    end () const
    {
      return impl().template end<cd,pitype>();
    }

    /** \brief obtain begin intersection iterator with respect to this view */
    IntersectionIterator
    ibegin ( const typename Codim< 0 > :: Entity &entity ) const
    {
      return impl().ibegin(entity);
    }

    /** \brief obtain end intersection iterator with respect to this view */
    IntersectionIterator
    iend ( const typename Codim< 0 > :: Entity &entity ) const
    {
      return impl().iend(entity);
    }

    /** \brief obtain communication object */
    const Communication &comm () const
    {
      return impl().comm();
    }

    /** \brief Return size of the overlap region for a given codim on the grid view.  */
    int overlapSize(int codim) const
    {
      return impl().overlapSize(codim);
    }

    /** \brief Return size of the ghost region for a given codim on the grid view.  */
    int ghostSize(int codim) const
    {
      return impl().ghostSize(codim);
    }

    /** \brief Communicate data on this view */
    template< class DataHandleImp, class DataType >
    auto communicate ( CommDataHandleIF< DataHandleImp, DataType > &data,
                       InterfaceType iftype,
                       CommunicationDirection dir ) const
    {
      typedef decltype( impl().communicate(data,iftype,dir) ) CommFuture;
      return communicate( data,iftype, dir,
            std::integral_constant< bool, std::is_same< CommFuture, void > :: value >() );
    }

    /**
     * \brief access to the underlying implementation
     *
     * \warning Implementation details may change without prior notification.
     **/
    Implementation &impl () { return impl_; }

    /**
     * \brief access to the underlying implementation
     *
     * \warning Implementation details may change without prior notification.
     **/
    const Implementation &impl () const { return impl_; }

  protected:
    /** \brief Communicate data on this view */
    template< class DataHandleImp, class DataType >
    auto communicate ( CommDataHandleIF< DataHandleImp, DataType > &data,
                       InterfaceType iftype,
                       CommunicationDirection dir, std::integral_constant< bool, false > ) const
    {
      return impl().communicate(data,iftype,dir);
    }

    struct DeprecatedMethodEmptyFuture : public Future<void>
    {
      void printMessage() const
      {
        std::cerr << "WARNING: GridView::communicate of '" <<
          typeid( Implementation ).name() << "' still returns void. Please update implementation to new interface returning a future object!" << std::endl;
      }

      bool ready () {
        printMessage();
        return true;
      }
      void wait () { printMessage(); }
      bool valid () const { printMessage(); return true; }
    };

    /** \brief Communicate data on this view */
    template< class DataHandleImp, class DataType >
    auto communicate ( CommDataHandleIF< DataHandleImp, DataType > &data,
                       InterfaceType iftype,
                       CommunicationDirection dir, std::integral_constant< bool, true > ) const
    {
      impl().communicate(data,iftype,dir);
      return DeprecatedMethodEmptyFuture();
    }

    template <class M>
    class CheckIsConformingImpl
    {
      // check for 'bool isConforming () const'
      template <class T, class R> static std::true_type testSignature(R (T::*)() const);

      template <class T>
      static decltype(testSignature(&T::isConforming)) test(std::nullptr_t);

      template <class T>
      static std::false_type test(...);

    public:
      using type = decltype(test<M>(nullptr));
    };

    [[deprecated("GridView implementation is missing a method 'bool isConforming() const'")]]
    bool isConformingDefaultImplementation() const { return bool(conforming); }

    Implementation impl_;
  };

} // namespace Dune

#endif // #ifndef DUNE_GRID_COMMON_GRIDVIEW_HH
