// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRIDVIEW_HH
#define DUNE_GRIDVIEW_HH

#include <dune/grid/common/gridenums.hh>

namespace Dune
{

  /** \addtogroup GIGridView
   *
   *  Though a DUNE grid is hierarchic, one often only needs access to
   *  a certain level of the grid. These views are provided by an
   *  implementation of GridView. Each grid exports a LevelGridView and
   *  a LeafGridView, corresponding to the two different types of levels.
   *
   *  A grid view provides the following functionality:
   *  - The index set for the grid level can be accessed by the indexSet()
   *    method.
   *  - A pair of begin() / end() methods provide iterators for each
   *    codimension.
   *  - A pair of ibegin() / iend() methods return suitable intersection
   *    iterators for a given entity of codimension 0.
   *  - For parallel computations, a suitable communicate() method is provided.
   *  - The underlying grid can be accessed through the grid() method.
   *  .
   *
   *  A GridView is obtained from the grid by calling one of the levelView() or
   *  leafView() methods.
   */


  /** \brief Grid view abstract base class
   *  \ingroup GIGridView
   *
   *  Interface class for view on grids. Grids return two types of view,
   *  a view of the leaf grid and of a level grid, which both satisfy
   *  the same interface. Through the view the user has access to the
   *  iterators, the intersections and the index set.
   *
   *  The interface is implemented using the engine concept.
   */
  template< class ViewTraits >
  class GridView
  {
    typedef GridView< ViewTraits > ThisType;

  public:
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

    /** \brief type of the collective communication */
    typedef typename Traits :: CollectiveCommunication CollectiveCommunication;

    /** \brief Codim Structure */
    template< int cd >
    struct Codim : public Traits :: template Codim<cd> {};

    enum { conforming = Traits :: conforming };

    enum { dimension = Grid :: dimension };
    enum { dimensionworld = Grid :: dimensionworld };

  public:
    GridView ( const GridViewImp& imp) : imp_(imp) {}
    GridView ( const ThisType &other ) : imp_(other.imp_) {}

  private:
    // prohibit assignment
    ThisType &operator= ( const ThisType & );

  public:
    /** \brief obtain a const reference to the underlying hierarchic grid */
    const Grid &grid () const
    {
      return asImp().grid();
    }

    /** \brief obtain the index set */
    const IndexSet &indexSet () const
    {
      return asImp().indexSet();
    }

    /** \brief obtain begin iterator for this view */
    template< int cd >
    typename Codim< cd > :: Iterator begin () const
    {
      return asImp().template begin<cd>();
    }

    /** \brief obtain end iterator for this view */
    template< int cd >
    typename Codim< cd > :: Iterator end () const
    {
      return asImp().template end<cd>();
    }

    /** \brief obtain begin intersection iterator with respect to this view */
    IntersectionIterator
    ibegin ( const typename Codim< 0 > :: Entity &entity ) const
    {
      return asImp().ibegin(entity);
    }

    /** \brief obtain end intersection iterator with respect to this view */
    IntersectionIterator
    iend ( const typename Codim< 0 > :: Entity &entity ) const
    {
      return asImp().iend(entity);
    }

    /** \brief level, this view belongs to */
    int level () const
    {
      return asImp().level();
    }

    /** \brief obtain collective communication object */
    const CollectiveCommunication &comm () const
    {
      return asImp().comm();
    }

    /** communicate data on this view */
    template< class DataHandleImp, class DataType >
    void communicate ( CommDataHandleIF< DataHandleImp, DataType > &data,
                       InterfaceType iftype,
                       CommunicationDirection dir ) const
    {
      asImp().communicate(data,iftype,dir);
    }

  protected:
    GridViewImp& asImp () {
      return imp_;
    }
    const GridViewImp& asImp () const {
      return imp_;
    }

  private:
    GridViewImp imp_;
  };

}

#endif
