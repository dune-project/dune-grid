// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_TREEITERATOR_HH
#define DUNE_ALBERTA_TREEITERATOR_HH

#include <dune/grid/albertagrid/meshpointer.hh>
#include <dune/grid/albertagrid/entitypointer.hh>

namespace Dune
{

  // AlbertaMarkerVector
  // -------------------

  /** \class   AlbertaMarkerVector
   *  \ingroup AlbertaGrid
   *  \brief   marker assigning subentities to one element containing them
   *
   *  This Helper class is used for the level and leaf iterators of higher
   *  codimension to visit each entity only once (on the element assigned to
   *  it by this marker)
   */
  template< int dim, int dimworld >
  class AlbertaMarkerVector
  {
    typedef AlbertaGrid< dim, dimworld > Grid;

    friend class AlbertaGrid< dim, dimworld >;

    static const int dimension = Grid::dimension;

    typedef typename Grid::HierarchicIndexSet HierarchicIndexSet;
    typedef Alberta::ElementInfo< dimension > ElementInfo;

    template< int codim >
    class MarkSubEntities;

  public:
    //! create AlbertaMarkerVector for Level or Leaf Iterator, true == LevelIterator
    //! the vectors stored inside are empty first
    explicit AlbertaMarkerVector ( const HierarchicIndexSet &hIndexSet )
      : hIndexSet_( &hIndexSet ),
        up2Date_( false )
    {}

    //! visit subentity on this element?
    template< int codim >
    bool subEntityOnElement ( const ElementInfo &elementInfo, int subEntity ) const;

    template< int firstCodim, class Iterator >
    void markSubEntities ( const Iterator &begin, const Iterator &end );

    //! return true if marking is up to date
    bool up2Date () const
    {
      return up2Date_;
    }

    //! unset up2date flag
    void unsetUp2Date () { up2Date_ = false; }

    //! print for debugin' only
    void print ( std::ostream &out = std::cout ) const;

  private:
    const HierarchicIndexSet *hIndexSet_;
    std::vector< int > marker_[ dimension+1 ];

    // true is vertex marker is up to date
    bool up2Date_;
  };



  namespace AlbertaTreeIteratorHelp
  {
    template< class IteratorImp, int dim, int codim >
    struct GoNextEntity;
  }


  // AlbertaGridTreeIterator
  // -----------------------

  /*!
     Enables iteration over all entities of a given codimension and level of a grid.
   */
  template< int codim, class GridImp, bool leafIterator >
  class AlbertaGridTreeIterator
    : public AlbertaGridEntityPointer< codim, GridImp >
  {
    typedef AlbertaGridTreeIterator< codim, GridImp, leafIterator > This;
    typedef AlbertaGridEntityPointer< codim, GridImp > Base;

  public:
    static const int dimension = GridImp::dimension;
    static const int codimension = codim;
    static const int dimensionworld = GridImp::dimensionworld;

  private:
    enum { dim = GridImp::dimension };
    friend class AlbertaGridEntity< 2, dimension, GridImp >;
    friend class AlbertaGridEntity< 1, dimension, GridImp >;
    friend class AlbertaGridEntity< 0, dimension, GridImp >;
    friend class AlbertaGrid< dimension, dimensionworld >;

    friend class AlbertaTreeIteratorHelp::GoNextEntity< This, dim, codim >;

    static const int numSubEntities
      = Alberta::NumSubEntities< dimension, codimension >::value;

    typedef typename GridImp::HierarchicIndexSet HierarchicIndexSet;

  public:
    typedef typename Base::ElementInfo ElementInfo;
    typedef Alberta::MeshPointer< dimension > MeshPointer;
    typedef typename MeshPointer::MacroIterator MacroIterator;

    typedef typename GridImp::template Codim< codim >::Entity Entity;
    typedef MakeableInterfaceObject< Entity > EntityObject;
    typedef typename EntityObject::ImplementationType EntityImp;

    typedef AlbertaMarkerVector< dimension, dimensionworld > MarkerVector;

    //! Constructor making end iterator
    AlbertaGridTreeIterator ( const This &other );

    //! Constructor making end iterator
    This &operator= ( const This &other );

    //! Constructor making end iterator
    AlbertaGridTreeIterator ( const GridImp &grid, int travLevel );

    //! Constructor making begin iterator
    AlbertaGridTreeIterator ( const GridImp &grid,
                              const MarkerVector *marker,
                              int travLevel );

    //! increment
    void increment();

  protected:
    using Base::entityImp;
    using Base::grid;

  private:
    // private Methods
    void makeIterator();

    void nextElement ( ElementInfo &elementInfo );
    void nextElementStop (ElementInfo &elementInfo );
    bool stopAtElement ( const ElementInfo &elementInfo ) const;

    void goNextEntity ( ElementInfo &elementInfo );

    // the real go next methods
    void goNextElement ( ElementInfo &elementInfo );
    void goNextFace ( ElementInfo &elementInfo );
    void goNextEdge ( ElementInfo &elementInfo );
    void goNextVertex ( ElementInfo &elementInfo );

    //! current level
    int level_;

    //! Number of the subentity within the element
    int subEntity_;

    MacroIterator macroIterator_;

    // knows on which element a point,edge,face is viewed
    const MarkerVector *marker_;
  };

}

#endif
