// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_ALBERTA_TREEITERATOR_HH
#define DUNE_ALBERTA_TREEITERATOR_HH

#include <utility>

#include <dune/common/hybridutilities.hh>
#include <dune/common/typetraits.hh>

#include <dune/grid/albertagrid/elementinfo.hh>
#include <dune/grid/albertagrid/meshpointer.hh>

#if HAVE_ALBERTA

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



  // AlbertaMarkerVector::NoMarkSubEntities
  // --------------------------------------

  template< int dim, int dimworld >
  template< bool >
  struct AlbertaMarkerVector< dim, dimworld >::NoMarkSubEntities
  {
    template< int firstCodim, class Iterator >
    static void mark ( [[maybe_unused]] const DofNumbering & dofNumbering,
                       [[maybe_unused]] int *(&marker)[ dimension + 1 ],
                       [[maybe_unused]] const Iterator &begin,
                       [[maybe_unused]] const Iterator &end )
    {}
  };



  // AlbertaMarkerVector::MarkSubEntities
  // ------------------------------------

  template< int dim, int dimworld >
  template< bool >
  struct AlbertaMarkerVector< dim, dimworld >::MarkSubEntities
  {
    template< int codim >
    struct Codim
    {
      static const int numSubEntities = Alberta::NumSubEntities< dimension, codim >::value;

      typedef Alberta::ElementInfo< dimension > ElementInfo;

      static void apply ( const DofNumbering &dofNumbering,
                          int *(&marker)[ dimension + 1 ],
                          const ElementInfo &elementInfo )
      {
        int *array = marker[ codim ];

        const int index = dofNumbering( elementInfo, 0, 0 );
        for( int i = 0; i < numSubEntities; ++i )
        {
          int &mark = array[ dofNumbering( elementInfo, codim, i ) ];
          mark = std::max( index, mark );
        }
      }
    };

    template< int firstCodim, class Iterator >
    static void mark ( const DofNumbering &dofNumbering, int *(&marker)[ dimension + 1 ],
                       const Iterator &begin, const Iterator &end )
    {
      for( int codim = firstCodim; codim <= dimension; ++codim )
      {
        const int size = dofNumbering.size( codim );
        marker[ codim ] = new int[ size ];

        int *array = marker[ codim ];
        for( int i = 0; i < size; ++i )
          array[ i ] = -1;
      }

      for( Iterator it = begin; it != end; ++it )
      {
        const ElementInfo &elementInfo = it->impl().elementInfo();
        Hybrid::forEach( std::make_index_sequence< dimension+1-firstCodim >{},
          [ & ]( auto i ){ Codim< i+firstCodim >::apply( dofNumbering, marker, elementInfo ); } );
      }
    }
  };



  // AlbertaGridTreeIterator
  // -----------------------

  /*!
     Enables iteration over all entities of a given codimension and level of a grid.
   */
  template< int codim, class GridImp, bool leafIterator >
  class AlbertaGridTreeIterator
  {
    typedef AlbertaGridTreeIterator< codim, GridImp, leafIterator > This;

  public:
    static const int dimension = GridImp::dimension;
    static const int codimension = codim;
    static const int dimensionworld = GridImp::dimensionworld;

  private:
    friend class AlbertaGrid< dimension, dimensionworld >;

    static const int numSubEntities
      = Alberta::NumSubEntities< dimension, codimension >::value;

  public:
    typedef Alberta::MeshPointer< dimension > MeshPointer;
    typedef typename MeshPointer::MacroIterator MacroIterator;

    typedef typename GridImp::template Codim< codim >::Entity Entity;
    typedef MakeableInterfaceObject< Entity > EntityObject;
    typedef typename EntityObject::ImplementationType EntityImp;
    typedef typename EntityImp::ElementInfo ElementInfo;

    typedef AlbertaMarkerVector< dimension, dimensionworld > MarkerVector;

    AlbertaGridTreeIterator ();

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

    //! equality
    bool equals ( const This &other ) const
    {
      return entity_.impl().equals( other.entity_.impl() );
    }

    //! dereferencing
    Entity &dereference () const
    {
      return entity_;
    }

    //! ask for level of entities
    int level () const
    {
      return entity_.impl().level();
    }

    //! increment
    void increment();

  protected:
    //! obtain a reference to the grid
    const GridImp &grid () const
    {
      return entity_.impl().grid();
    }

  private:
    void nextElement ( ElementInfo &elementInfo );
    void nextElementStop (ElementInfo &elementInfo );
    bool stopAtElement ( const ElementInfo &elementInfo ) const;

    void goNext ( ElementInfo &elementInfo );
    void goNext ( const std::integral_constant< int, 0 > cdVariable,
                  ElementInfo &elementInfo );
    void goNext ( const std::integral_constant< int, 1 > cdVariable,
                  ElementInfo &elementInfo );
    template< int cd >
    void goNext ( const std::integral_constant< int, cd > cdVariable,
                  ElementInfo &elementInfo );

    mutable Entity entity_;

    //! current level
    int level_;

    //! Number of the subentity within the element
    int subEntity_;

    MacroIterator macroIterator_;

    // knows on which element a point,edge,face is viewed
    const MarkerVector *marker_;
  };



  // Implementation of AlbertaMarkerVector
  // -------------------------------------

  template< int dim, int dimworld >
  template< int codim >
  inline bool AlbertaMarkerVector< dim, dimworld >
  ::subEntityOnElement ( const ElementInfo &elementInfo, int subEntity ) const
  {
    assert( marker_[ codim ] != 0 );

    const int subIndex = dofNumbering_( elementInfo, codim, subEntity );
    const int markIndex = marker_[ codim ][ subIndex ];
    assert( (markIndex >= 0) );

    const int index = dofNumbering_( elementInfo, 0, 0 );
    return (markIndex == index);
  }


  template< int dim, int dimworld >
  template< int firstCodim, class Iterator >
  inline void AlbertaMarkerVector< dim, dimworld >
  ::markSubEntities ( const Iterator &begin, const Iterator &end )
  {
    clear();
    std::conditional< (firstCodim <= dimension), MarkSubEntities<true>, NoMarkSubEntities<false> >::type
    ::template mark< firstCodim, Iterator >( dofNumbering_, marker_, begin, end );
  }


  template< int dim, int dimworld >
  inline void AlbertaMarkerVector< dim, dimworld >::print ( std::ostream &out ) const
  {
    for( int codim = 1; codim <= dimension; ++codim )
    {
      int *marker = marker_[ codim ];
      if( marker != 0 )
      {
        const int size = dofNumbering_.size( codim );
        out << std::endl;
        out << "Codimension " << codim << " (" << size << " entries)" << std::endl;
        for( int i = 0; i < size; ++i )
          out << "subentity " << i << " visited on Element " << marker[ i ] << std::endl;
      }
    }
  }



  // Implementation of AlbertaGridTreeIterator
  // -----------------------------------------

  template< int codim, class GridImp, bool leafIterator >
  inline AlbertaGridTreeIterator< codim, GridImp, leafIterator >
  ::AlbertaGridTreeIterator ()
    : entity_(),
      level_( -1 ),
      subEntity_( -1 ),
      macroIterator_(),
      marker_( NULL )
  {}

  template< int codim, class GridImp, bool leafIterator >
  inline AlbertaGridTreeIterator< codim, GridImp, leafIterator >
  ::AlbertaGridTreeIterator ( const GridImp &grid,
                              const MarkerVector *marker,
                              int travLevel )
    : entity_( EntityImp( grid ) ),
      level_( travLevel ),
      subEntity_( (codim == 0 ? 0 : -1) ),
      macroIterator_( grid.meshPointer().begin() ),
      marker_( marker )
  {
    ElementInfo elementInfo = *macroIterator_;
    nextElementStop( elementInfo );
    if( codim > 0 )
      goNext( elementInfo );
    // it is ok to set the invalid ElementInfo
    entity_.impl().setElement( elementInfo, subEntity_ );
  }


  // Make LevelIterator with point to element from previous iterations
  template< int codim, class GridImp, bool leafIterator >
  inline AlbertaGridTreeIterator< codim, GridImp, leafIterator >
  ::AlbertaGridTreeIterator ( const GridImp &grid,
                              int travLevel )
    : entity_( EntityImp( grid ) ),
      level_( travLevel ),
      subEntity_( -1 ),
      macroIterator_( grid.meshPointer().end() ),
      marker_( 0 )
  {}


  // Make LevelIterator with point to element from previous iterations
  template< int codim, class GridImp, bool leafIterator >
  inline AlbertaGridTreeIterator< codim, GridImp, leafIterator >
  ::AlbertaGridTreeIterator( const This &other )
    : entity_( other.entity_ ),
      level_( other.level_ ),
      subEntity_( other.subEntity_ ),
      macroIterator_( other.macroIterator_ ),
      marker_( other.marker_ )
  {}


  // Make LevelIterator with point to element from previous iterations
  template< int codim, class GridImp, bool leafIterator >
  inline typename AlbertaGridTreeIterator< codim, GridImp, leafIterator >::This &
  AlbertaGridTreeIterator< codim, GridImp, leafIterator >::operator= ( const This &other )
  {
    entity_ = other.entity_;
    level_ = other.level_;
    subEntity_ =  other.subEntity_;
    macroIterator_ = other.macroIterator_;
    marker_ = other.marker_;

    return *this;
  }


  template< int codim, class GridImp, bool leafIterator >
  inline void AlbertaGridTreeIterator< codim, GridImp, leafIterator >::increment ()
  {
    ElementInfo elementInfo = entity_.impl().elementInfo_;
    goNext ( elementInfo );
    // it is ok to set the invalid ElementInfo
    entity_.impl().setElement( elementInfo, subEntity_ );
  }


  template< int codim, class GridImp, bool leafIterator >
  inline void AlbertaGridTreeIterator< codim, GridImp, leafIterator >
  ::nextElement ( ElementInfo &elementInfo )
  {
    if( elementInfo.isLeaf() || (elementInfo.level() >= level_) )
    {
      while( (elementInfo.level() > 0) && (elementInfo.indexInFather() == 1) )
        elementInfo = elementInfo.father();
      if( elementInfo.level() == 0 )
      {
        ++macroIterator_;
        elementInfo = *macroIterator_;
      }
      else
        elementInfo = elementInfo.father().child( 1 );
    }
    else
      elementInfo = elementInfo.child( 0 );
  }


  template< int codim, class GridImp, bool leafIterator >
  inline void AlbertaGridTreeIterator< codim, GridImp, leafIterator >
  ::nextElementStop ( ElementInfo &elementInfo )
  {
    while( !(!elementInfo || stopAtElement( elementInfo )) )
      nextElement( elementInfo );
  }


  template< int codim, class GridImp, bool leafIterator >
  inline bool AlbertaGridTreeIterator< codim, GridImp, leafIterator >
  ::stopAtElement ( const ElementInfo &elementInfo ) const
  {
    if( !elementInfo )
      return true;
    return (leafIterator ? elementInfo.isLeaf() : (level_ == elementInfo.level()));
  }


  template< int codim, class GridImp, bool leafIterator >
  inline void AlbertaGridTreeIterator< codim, GridImp, leafIterator >
  ::goNext ( ElementInfo &elementInfo )
  {
    std::integral_constant< int, codim > codimVariable;
    goNext( codimVariable, elementInfo );
  }

  template< int codim, class GridImp, bool leafIterator >
  inline void AlbertaGridTreeIterator< codim, GridImp, leafIterator >
  ::goNext ( const std::integral_constant< int, 0 > /* cdVariable */,
             ElementInfo &elementInfo )
  {
    assert( stopAtElement( elementInfo ) );

    nextElement( elementInfo );
    nextElementStop( elementInfo );
  }

  template< int codim, class GridImp, bool leafIterator >
  inline void AlbertaGridTreeIterator< codim, GridImp, leafIterator >
  ::goNext ( const std::integral_constant< int, 1 > cdVariable,
             ElementInfo &elementInfo )
  {
    assert( stopAtElement( elementInfo ) );

    ++subEntity_;
    if( subEntity_ >= numSubEntities )
    {
      subEntity_ = 0;
      nextElement( elementInfo );
      nextElementStop( elementInfo );
      if( !elementInfo )
        return;
    }

    if( leafIterator )
    {
      const int face = (dimension == 1 ? (numSubEntities-1)-subEntity_ : subEntity_);

      const ALBERTA EL *neighbor = elementInfo.elInfo().neigh[ face ];
      if( (neighbor != NULL) && !elementInfo.isBoundary( face ) )
      {
        // face is reached from element with largest number
        const int elIndex = grid().dofNumbering() ( elementInfo, 0, 0 );
        const int nbIndex = grid().dofNumbering() ( neighbor, 0, 0 );
        if( elIndex < nbIndex )
          goNext( cdVariable, elementInfo );
      }
      // uncomment this assertion only if codimension 1 entities are marked
      // assert( marker_->template subEntityOnElement< 1 >( elementInfo, subEntity_ ) );
    }
    else
    {
      assert( marker_ != 0 );
      if( !marker_->template subEntityOnElement< 1 >( elementInfo, subEntity_ ) )
        goNext( cdVariable, elementInfo );
    }
  }

  template< int codim, class GridImp, bool leafIterator >
  template< int cd >
  inline void AlbertaGridTreeIterator< codim, GridImp, leafIterator >
  ::goNext ( const std::integral_constant< int, cd > cdVariable,
             ElementInfo &elementInfo )
  {
    assert( stopAtElement( elementInfo ) );

    ++subEntity_;
    if( subEntity_ >= numSubEntities )
    {
      subEntity_ = 0;
      nextElement( elementInfo );
      nextElementStop( elementInfo );
      if( !elementInfo )
        return;
    }

    assert( marker_ != 0 );
    if( !marker_->template subEntityOnElement< cd >( elementInfo, subEntity_ ) )
      goNext( cdVariable, elementInfo );
  }

}

#endif // #if HAVE_ALBERTA

#endif // #ifndef DUNE_ALBERTA_TREEITERATOR_HH
