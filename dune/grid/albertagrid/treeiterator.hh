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
    typedef AlbertaMarkerVector< dim, dimworld > This;

    typedef AlbertaGrid< dim, dimworld > Grid;

    //friend class AlbertaGrid< dim, dimworld >;

    static const int dimension = Grid::dimension;

    typedef Alberta::HierarchyDofNumbering< dimension > DofNumbering;
    typedef Alberta::ElementInfo< dimension > ElementInfo;

    template< bool >
    struct NoMarkSubEntities;
    template< bool >
    struct MarkSubEntities;

  public:
    //! create AlbertaMarkerVector with empty vectors
    explicit AlbertaMarkerVector ( const DofNumbering &dofNumbering )
      : dofNumbering_( dofNumbering )
    {
      for( int codim = 0; codim <= dimension; ++codim )
        marker_[ codim ] = 0;
    }

    AlbertaMarkerVector ( const This &other )
      : dofNumbering_( other.dofNumbering_ )
    {
      for( int codim = 0; codim <= dimension; ++codim )
        marker_[ codim ] = 0;
    }

    ~AlbertaMarkerVector ()
    {
      clear();
    }

  private:
    This &operator= ( const This & );

  public:
    //! visit subentity on this element?
    template< int codim >
    bool subEntityOnElement ( const ElementInfo &elementInfo, int subEntity ) const;

    template< int firstCodim, class Iterator >
    void markSubEntities ( const Iterator &begin, const Iterator &end );

    void clear ()
    {
      for( int codim = 0; codim <= dimension; ++codim )
      {
        if( marker_[ codim ] != 0 )
          delete[] marker_[ codim ];
        marker_[ codim ] = 0;
      }
    }

    //! return true if marking is up to date
    bool up2Date () const
    {
      return (marker_[ dimension ] != 0);
    }

    //! print for debugin' only
    void print ( std::ostream &out = std::cout ) const;

  private:
    const DofNumbering &dofNumbering_;
    int *marker_[ dimension+1 ];
  };



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
    friend class AlbertaGrid< dimension, dimensionworld >;

    static const int numSubEntities
      = Alberta::NumSubEntities< dimension, codimension >::value;

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
    void nextElement ( ElementInfo &elementInfo );
    void nextElementStop (ElementInfo &elementInfo );
    bool stopAtElement ( const ElementInfo &elementInfo ) const;

    void goNext ( ElementInfo &elementInfo );
    void goNext ( const Int2Type< 0 > cdVariable, ElementInfo &elementInfo );
    void goNext ( const Int2Type< 1 > cdVariable, ElementInfo &elementInfo );
    template< int cd >
    void goNext ( const Int2Type< cd > cdVariable, ElementInfo &elementInfo );

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
