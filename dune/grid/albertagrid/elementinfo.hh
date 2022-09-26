// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_ELEMENTINFO_HH
#define DUNE_ALBERTA_ELEMENTINFO_HH

/** \file
 *  \author Martin Nolte
 *  \brief  provides a wrapper for ALBERTA's el_info structure
 */

#include <cassert>
#include <vector>
#include <utility>

#include <dune/grid/albertagrid/geometrycache.hh>
#include <dune/grid/albertagrid/macroelement.hh>

#if HAVE_ALBERTA

namespace Dune
{

  namespace Alberta
  {

    // External Forward Declarations
    // -----------------------------

    template< int dim >
    class MeshPointer;

    struct BasicNodeProjection;



    // ElementInfo
    // -----------

    template< int dim >
    class ElementInfo
    {
      struct Instance;
      class Stack;

      template< int >
      struct Library;

      typedef Instance *InstancePtr;

    public:
      static const int dimension = dim;

      static const int numVertices = NumSubEntities< dimension, dimension >::value;
      static const int numFaces = NumSubEntities< dimension, 1 >::value;

      typedef Alberta::MacroElement< dimension > MacroElement;
      typedef Alberta::MeshPointer< dimension > MeshPointer;
      typedef Alberta::FillFlags< dimension > FillFlags;

      static const int maxNeighbors = N_NEIGH_MAX;

      static const int maxLevelNeighbors = Library< dimWorld >::maxLevelNeighbors;

#if !DUNE_ALBERTA_CACHE_COORDINATES
      typedef GeometryCacheProxy< dim > GeometryCache;
#endif

      struct Seed;

    private:
      explicit ElementInfo ( const InstancePtr &instance );

    public:
      ElementInfo ();
      ElementInfo ( const MeshPointer &mesh, const MacroElement &macroElement,
                    typename FillFlags::Flags fillFlags = FillFlags::standard );
      ElementInfo ( const MeshPointer &mesh, const Seed &seed,
                    typename FillFlags::Flags fillFlags = FillFlags::standard );
      ElementInfo ( const ElementInfo &other );
      ElementInfo ( ElementInfo&& other );

      ~ElementInfo ();

      ElementInfo &operator= ( const ElementInfo &other );
      ElementInfo &operator= ( ElementInfo &&other );

      explicit operator bool () const { return (instance_ != null()); }

      bool operator== ( const ElementInfo &other ) const;
      bool operator!= ( const ElementInfo &other ) const;

      const MacroElement &macroElement () const;
      ElementInfo father () const;
      int indexInFather () const;
      ElementInfo child ( int i ) const;
      bool isLeaf () const;

      Seed seed () const;

      MeshPointer mesh () const;

      bool mightVanish () const;

      int level () const;
      // see ALBERTA documentation for definition of element type
      // values are 0, 1, 2
      int type () const;

      int getMark () const;
      void setMark ( int refCount ) const;

      bool hasLeafNeighbor ( const int face ) const;
      ElementInfo leafNeighbor ( const int face ) const;

      /* obtain all level neighbors of a face
       *
       * param[in]  face            face for which the neighbors are desired
       * param[out] neighbor        array storing the neighbors
       * param[out] faceInNeighbor  array storing the faces in neighbor
       *                            (-1, if this neighbor does not exist)
       *
       * returns (potential) number of neighbors (i.e., the number of valid
       *         entries in the output arrays
       */
      int levelNeighbors ( const int face, ElementInfo (&neighbor)[ maxLevelNeighbors ], int (&faceInNeighbor)[ maxLevelNeighbors ] ) const;

      template< int codim >
      int twist ( int subEntity ) const;
      int twistInNeighbor ( int face ) const;
      bool isBoundary ( int face ) const;
      int boundaryId ( int face ) const;
      AffineTransformation *transformation ( int face ) const;
      BasicNodeProjection *boundaryProjection ( int face ) const;

      bool hasCoordinates () const;
      const GlobalVector &coordinate ( int vertex ) const;
#if !DUNE_ALBERTA_CACHE_COORDINATES
      GeometryCache geometryCache () const
      {
        return GeometryCache( instance_->geometryCache, instance_->elInfo );
      }
#endif

      template< class Functor >
      void hierarchicTraverse ( Functor &functor ) const;

      template< class Functor >
      void leafTraverse ( Functor &functor ) const;

      const Element *element () const;
      const Element *neighbor ( int face ) const;
      Element *el () const;
      ALBERTA EL_INFO &elInfo () const;

      static ElementInfo
      createFake ( const MeshPointer &mesh,
                   const Element *element, int level, int type = 0 );
      static ElementInfo createFake ( const ALBERTA EL_INFO &elInfo );

    private:
      static bool isLeaf ( Element *element );
      static bool mightVanish ( Element *element, int depth );

      static void fill ( Mesh *mesh, const ALBERTA MACRO_EL *mel, ALBERTA EL_INFO &elInfo );
      static void fill ( int ichild, const ALBERTA EL_INFO &parentInfo, ALBERTA EL_INFO &elInfo );

      void addReference () const;
      void removeReference () const;

      static InstancePtr null ();
      static Stack &stack ();

      InstancePtr instance_;
    };



    // ElementInfo::Instance
    // ---------------------

    template< int dim >
    struct ElementInfo< dim >::Instance
    {
      ALBERTA EL_INFO elInfo;
      unsigned int refCount;

      InstancePtr &parent ()
      {
        return parent_;
      }

    private:
      InstancePtr parent_;

#if !DUNE_ALBERTA_CACHE_COORDINATES
    public:
      Alberta::GeometryCache< dim > geometryCache;
#endif
    };



    // ElementInfo::Stack
    // ------------------

    template< int dim >
    class ElementInfo< dim >::Stack
    {
      InstancePtr top_;
      Instance null_;

    public:
      Stack ();
      ~Stack ();

      InstancePtr allocate ();
      void release ( InstancePtr &p );
      InstancePtr null ();
    };



    // ElementInfo::Library
    // --------------------

    template< int dim >
    template< int >
    struct ElementInfo< dim >::Library
    {
      typedef Alberta::ElementInfo< dim > ElementInfo;

      static const int maxLevelNeighbors = (1 << (dim-1));

      static int
      leafNeighbor ( const ElementInfo &element, const int face, ElementInfo &neighbor );

      static int
      levelNeighbors ( const ElementInfo &element, const int face,
                       ElementInfo (&neighbor)[ maxLevelNeighbors ], int (&faceInNeighbor)[ maxLevelNeighbors ] );

    private:
      static int
      macroNeighbor ( const ElementInfo &element, const int face, ElementInfo &neighbor );
    };



    // ElementInfo::Seed
    // -----------------

    template< int dim >
    struct ElementInfo< dim >::Seed
    {
      Seed ()
        : macroIndex_( -1 ), level_( 0 ), path_( 0 )
      {}

      Seed ( const int macroIndex, const int level, const unsigned long path )
        : macroIndex_( macroIndex ), level_( level ), path_( path )
      {}

      bool operator== ( const Seed &other ) const
      {
        return (macroIndex() == other.macroIndex()) && (level() == other.level()) && (path() == other.path());
      }

      bool operator< ( const Seed &other ) const
      {
        const bool ml = (macroIndex() < other.macroIndex());
        const bool me = (macroIndex() == other.macroIndex());
        const bool ll = (level() < other.level());
        const bool le = (level() == other.level());
        const bool pl = (path() < other.path());
        return ml | (me & (ll | (le & pl)));
      }

      bool operator!= ( const Seed &other ) const { return !(*this == other); }
      bool operator<= ( const Seed &other ) const { return !(other < *this); }
      bool operator> ( const Seed &other ) const { return (other < *this); }
      bool operator>= ( const Seed &other ) const { return !(*this < other); }

      bool isValid ( ) const { return macroIndex_ != -1; }

      int macroIndex () const { return macroIndex_; }
      int level () const { return level_; }
      unsigned long path () const { return path_; }

    private:
      int macroIndex_;
      int level_;
      unsigned long path_;
    };



    // Implementation of ElementInfo
    // -----------------------------

    template< int dim >
    inline ElementInfo< dim >::ElementInfo ( const InstancePtr &instance )
      : instance_( instance )
    {
      addReference();
    }


    template< int dim >
    inline ElementInfo< dim >::ElementInfo ()
      : instance_( null() )
    {
      addReference();
    }


    template< int dim >
    inline ElementInfo< dim >
    ::ElementInfo ( const MeshPointer &mesh, const MacroElement &macroElement,
                    typename FillFlags::Flags fillFlags )
    {
      instance_ = stack().allocate();
      instance_->parent() = null();
      ++(instance_->parent()->refCount);

      addReference();

      elInfo().fill_flag = fillFlags;

      // Alberta fills opp_vertex only if there is a neighbor
      for( int k = 0; k < maxNeighbors; ++k )
        elInfo().opp_vertex[ k ] = -1;

      fill( mesh, &macroElement, elInfo() );
    }


    template< int dim >
    inline ElementInfo< dim >
    ::ElementInfo ( const MeshPointer &mesh, const Seed &seed,
                    typename FillFlags::Flags fillFlags )
    {
      instance_ = stack().allocate();
      instance_->parent() = null();
      ++(instance_->parent()->refCount);

      addReference();

      // fill in macro element info
      elInfo().fill_flag = fillFlags;

      // Alberta fills opp_vertex only if there is a neighbor
      for( int k = 0; k < maxNeighbors; ++k )
        elInfo().opp_vertex[ k ] = -1;

      fill( mesh, ((Mesh *)mesh)->macro_els + seed.macroIndex(), elInfo() );

      // traverse the seed's path
      unsigned long path = seed.path();
      for( int i = 0; i < seed.level(); ++i )
      {
        InstancePtr child = stack().allocate();
        child->parent() = instance_;

        // Alberta fills opp_vertex only if there is a neighbor
        for( int k = 0; k < maxNeighbors; ++k )
          child->elInfo.opp_vertex[ k ] = -2;

        fill( path & 1, elInfo(), child->elInfo );

        instance_ = child;
        addReference();

        path = path >> 1;
      }

      assert( this->seed() == seed );
    }


    template< int dim >
    inline ElementInfo< dim >::ElementInfo ( const ElementInfo &other )
      : instance_( other.instance_ )
    {
      addReference();
    }

    template< int dim >
    inline ElementInfo< dim >::ElementInfo ( ElementInfo &&other )
      : instance_( NULL )
    {
      using std::swap;
      swap( instance_, other.instance_ );
    }

    template< int dim >
    inline ElementInfo< dim >::~ElementInfo ()
    {
      removeReference();
    }


    template< int dim >
    inline ElementInfo< dim > &
    ElementInfo< dim >::operator= ( const ElementInfo< dim > &other )
    {
      other.addReference();
      removeReference();
      instance_ = other.instance_;
      return *this;
    }

    template< int dim >
    inline ElementInfo< dim > &
    ElementInfo< dim >::operator= ( ElementInfo< dim > &&other )
    {
      using std::swap;
      swap( instance_, other.instance_ );
      return *this;
    }

    template< int dim >
    inline bool
    ElementInfo< dim >::operator== ( const ElementInfo< dim > &other ) const
    {
      return (instance_->elInfo.el == other.instance_->elInfo.el);
    }


    template< int dim >
    inline bool
    ElementInfo< dim >::operator!= ( const ElementInfo< dim > &other ) const
    {
      return (instance_->elInfo.el != other.instance_->elInfo.el);
    }


    template< int dim >
    inline const typename ElementInfo< dim >::MacroElement &
    ElementInfo< dim >::macroElement () const
    {
      assert( !!(*this) );
      assert( elInfo().macro_el != NULL );
      return static_cast< const MacroElement & >( *(elInfo().macro_el) );
    }


    template< int dim >
    inline ElementInfo< dim > ElementInfo< dim >::father () const
    {
      assert( !!(*this) );
      return ElementInfo< dim >( instance_->parent() );
    }


    template< int dim >
    inline int ElementInfo< dim >::indexInFather () const
    {
      const Element *element = elInfo().el;
      const Element *father = elInfo().parent->el;
      assert( father != NULL );

      const int index = (father->child[ 0 ] == element ? 0 : 1);
      assert( father->child[ index ] == element );
      return index;
    }


    template< int dim >
    inline ElementInfo< dim > ElementInfo< dim >::child ( int i ) const
    {
      assert( !isLeaf() );

      InstancePtr child = stack().allocate();
      child->parent() = instance_;
      addReference();

      // Alberta fills opp_vertex only if there is a neighbor
      for( int k = 0; k < maxNeighbors; ++k )
        child->elInfo.opp_vertex[ k ] = -2;

      fill( i, elInfo(), child->elInfo );
      return ElementInfo< dim >( child );
    }


    template< int dim >
    inline bool ElementInfo< dim >::isLeaf () const
    {
      assert( !(*this) == false );
      return isLeaf( el() );
    }


    template< int dim >
    inline typename ElementInfo< dim >::Seed ElementInfo< dim >::seed () const
    {
      assert( !!(*this) );

      int level = 0;
      unsigned long path = 0;
      for( InstancePtr p = instance_; p->parent() != null(); p = p->parent() )
      {
        const Element *element = p->elInfo.el;
        const Element *father = p->parent()->elInfo.el;
        const unsigned long child = static_cast< unsigned long >( father->child[ 1 ] == element );
        path = (path << 1) | child;
        ++level;
      }

      if( level != elInfo().level )
        DUNE_THROW( NotImplemented, "Seed for fake elements not implemented." );

      return Seed( macroElement().index, level, path );
    }


    template< int dim >
    inline typename ElementInfo< dim >::MeshPointer ElementInfo< dim >::mesh () const
    {
      return MeshPointer( elInfo().mesh );
    }


    template< int dim >
    inline bool ElementInfo< dim >::mightVanish () const
    {
      return mightVanish( el(), 0 );
    }


    template< int dim >
    inline int ElementInfo< dim >::level () const
    {
      return elInfo().level;
    }


    template< int dim >
    inline int ElementInfo< dim >::type () const
    {
      return 0;
    }


    template<>
    inline int ElementInfo< 3 >::type () const
    {
      return instance_->elInfo.el_type;
    }


    template< int dim >
    inline int ElementInfo< dim >::getMark () const
    {
      return el()->mark;
    }


    template< int dim >
    inline void ElementInfo< dim >::setMark ( int refCount ) const
    {
      assert( isLeaf() );
      assert( (refCount >= -128) && (refCount < 127) );
      el()->mark = refCount;
    }


    template< int dim >
    inline bool ElementInfo< dim >::hasLeafNeighbor ( const int face ) const
    {
      assert( !!(*this) );
      assert( (face >= 0) && (face < maxNeighbors) );

      assert( (elInfo().fill_flag & FillFlags::boundaryId) != 0 );
      const int macroFace = elInfo().macro_wall[ face ];
      if( macroFace >= 0 )
        return (macroElement().neighbor( macroFace ) != NULL);
      else
        return true;
    }


    template< int dim >
    inline ElementInfo< dim > ElementInfo< dim >::leafNeighbor ( const int face ) const
    {
      assert( (face >= 0) && (face < numFaces) );
      ElementInfo neighbor;
      Library< dimWorld >::leafNeighbor( *this, face, neighbor );
      return neighbor;
    }


    template< int dim >
    inline int ElementInfo< dim >
    ::levelNeighbors ( const int face, ElementInfo (&neighbor)[ maxLevelNeighbors ], int (&faceInNeighbor)[ maxLevelNeighbors ] ) const
    {
      assert( (face >= 0) && (face < numFaces) );
      return Library< dimWorld >::levelNeighbors( *this, face, neighbor, faceInNeighbor );
    }


    template< int dim >
    template< int codim >
    inline int ElementInfo< dim >::twist ( int subEntity ) const
    {
      return Twist< dim, dim-codim >::twist( element(), subEntity );
    }


    template< int dim >
    inline int ElementInfo< dim >::twistInNeighbor ( const int face ) const
    {
      assert( neighbor( face ) != NULL );
      return Twist< dim, dim-1 >::twist( neighbor( face ), elInfo().opp_vertex[ face ] );
    }


    template< int dim >
    inline bool ElementInfo< dim >::isBoundary ( int face ) const
    {
      assert( !!(*this) );
      assert( (face >= 0) && (face < maxNeighbors) );

      assert( (elInfo().fill_flag & FillFlags::boundaryId) != 0 );
      const int macroFace = elInfo().macro_wall[ face ];
      if( macroFace >= 0 )
        return macroElement().isBoundary( macroFace );
      else
        return false;
    }


    template< int dim >
    inline int ElementInfo< dim >::boundaryId ( int face ) const
    {
      assert( !!(*this) );
      assert( (face >= 0) && (face < N_WALLS_MAX) );

      assert( (elInfo().fill_flag & FillFlags::boundaryId) != 0 );
      const int macroFace = elInfo().macro_wall[ face ];
      const int id = macroElement().boundaryId( macroFace );
      // this assertion is only allowed, if FILL_BOUND is set
      // assert( id == elInfo().wall_bound[ face ] );
      return id;
    }


    template< int dim >
    inline AffineTransformation *
    ElementInfo< dim >::transformation ( int face ) const
    {
      assert( !!(*this) );
      assert( (face >= 0) && (face < N_WALLS_MAX) );

      assert( (elInfo().fill_flag & FillFlags::boundaryId) != 0 );
      const int macroFace = elInfo().macro_wall[ face ];
      return (macroFace < 0 ? NULL : macroElement().wall_trafo[ macroFace ]);
    }


    template< int dim >
    inline BasicNodeProjection *
    ElementInfo< dim >::boundaryProjection ( int face ) const
    {
      assert( !!(*this) );
      assert( (face >= 0) && (face < N_WALLS_MAX) );

      assert( (elInfo().fill_flag & FillFlags::boundaryId) != 0 );
      const int macroFace = elInfo().macro_wall[ face ];
      if( macroFace >= 0 )
        return static_cast< BasicNodeProjection * >( macroElement().projection[ macroFace+1 ] );
      else
        return 0;
    }


    template< int dim >
    inline bool ElementInfo< dim >::hasCoordinates () const
    {
      return ((elInfo().fill_flag & FillFlags::coords) != 0);
    }

    template< int dim >
    inline const GlobalVector &ElementInfo< dim >::coordinate ( int vertex ) const
    {
      assert( hasCoordinates() );
      assert( (vertex >= 0) && (vertex < numVertices) );
      return elInfo().coord[ vertex ];
    }


    template< int dim >
    template< class Functor >
    inline void ElementInfo< dim >::hierarchicTraverse ( Functor &functor ) const
    {
      functor( *this );
      if( !isLeaf() )
      {
        child( 0 ).hierarchicTraverse( functor );
        child( 1 ).hierarchicTraverse( functor );
      }
    }


    template< int dim >
    template< class Functor >
    inline void ElementInfo< dim >::leafTraverse ( Functor &functor ) const
    {
      if( !isLeaf() )
      {
        child( 0 ).leafTraverse( functor );
        child( 1 ).leafTraverse( functor );
      }
      else
        functor( *this );
    }


    template< int dim >
    inline const Element *ElementInfo< dim >::element () const
    {
      return elInfo().el;
    }


    template< int dim >
    inline const Element *ElementInfo< dim >::neighbor ( int face ) const
    {
      assert( (face >= 0) && (face < numFaces) );
      assert( (elInfo().fill_flag & FillFlags::neighbor) != 0 );
      return elInfo().neigh[ face ];
    }


    template< int dim >
    inline Element *ElementInfo< dim >::el () const
    {
      return elInfo().el;
    }


    template< int dim >
    inline ALBERTA EL_INFO &ElementInfo< dim >::elInfo () const
    {
      return (instance_->elInfo);
    }


    template< int dim >
    inline ElementInfo< dim >
    ElementInfo< dim >::createFake ( const MeshPointer &mesh,
                                     const Element *element, int level, int type )
    {
      InstancePtr instance = stack().allocate();
      instance->parent() = null();
      ++(instance->parent()->refCount);

      instance->elInfo.mesh = mesh;
      instance->elInfo.macro_el = NULL;
      instance->elInfo.el = const_cast< Element * >( element );
      instance->elInfo.parent = NULL;
      instance->elInfo.fill_flag = FillFlags::nothing;
      instance->elInfo.level = level;
      instance->elInfo.el_type = type;

      return ElementInfo< dim >( instance );
    }


    template< int dim >
    inline ElementInfo< dim >
    ElementInfo< dim >::createFake ( const ALBERTA EL_INFO &elInfo )
    {
      InstancePtr instance = stack().allocate();
      instance->parent() = null();
      ++(instance->parent()->refCount);

      instance->elInfo = elInfo;
      return ElementInfo< dim >( instance );
    }


    template< int dim >
    inline bool ElementInfo< dim >::isLeaf ( Element *element )
    {
      return IS_LEAF_EL( element );
    }


    template< int dim >
    inline bool ElementInfo< dim >::mightVanish ( Alberta::Element *element, int depth )
    {
      if( isLeaf( element ) )
        return (element->mark < depth);
      else
        return (mightVanish( element->child[ 0 ], depth-1 ) && mightVanish( element->child[ 1 ], depth-1 ));
    }


    template< int dim >
    inline void ElementInfo< dim >
    ::fill ( Mesh *mesh, const ALBERTA MACRO_EL *mel, ALBERTA EL_INFO &elInfo )
    {
      ALBERTA fill_macro_info( mesh, mel, &elInfo );
    }

    template< int dim >
    inline void ElementInfo< dim >
    ::fill ( int ichild, const ALBERTA EL_INFO &parentInfo, ALBERTA EL_INFO &elInfo )
    {
      ALBERTA fill_elinfo( ichild, FILL_ANY, &parentInfo, &elInfo );
    }


    template< int dim >
    inline void ElementInfo< dim >::addReference () const
    {
      ++(instance_->refCount);
    }


    template< int dim >
    inline void ElementInfo< dim >::removeReference () const
    {
      // short-circuit for rvalues that have been drained as argument to a move operation
      if ( !instance_ )
        return;
      // this loop breaks when instance becomes null()
      for( InstancePtr instance = instance_; --(instance->refCount) == 0; )
      {
        const InstancePtr parent = instance->parent();
        stack().release( instance );
        instance = parent;
      }
    }


    template< int dim >
    inline typename ElementInfo< dim >::InstancePtr
    ElementInfo< dim >::null ()
    {
      return stack().null();
    }


    template< int dim >
    inline typename ElementInfo< dim >::Stack &
    ElementInfo< dim >::stack ()
    {
      static Stack s;
      return s;
    }



    // Implementation of ElementInfo::Stack
    // ------------------------------------

    template< int dim >
    inline ElementInfo< dim >::Stack::Stack ()
      : top_( 0 )
    {
      null_.elInfo.el = NULL;
      null_.refCount = 1;
      null_.parent() = 0;
    }


    template< int dim >
    inline ElementInfo< dim >::Stack::~Stack ()
    {
      while( top_ != 0 )
      {
        InstancePtr p = top_;
        top_ = p->parent();
        delete p;
      }
    }


    template< int dim >
    inline typename ElementInfo< dim >::InstancePtr
    ElementInfo< dim >::Stack::allocate ()
    {
      InstancePtr p = top_;
      if( p != 0 )
        top_ = p->parent();
      else
        p = new Instance;
      p->refCount = 0;
      return p;
    }


    template< int dim >
    inline void ElementInfo< dim >::Stack::release ( InstancePtr &p )
    {
      assert( (p != null()) && (p->refCount == 0) );
      p->parent() = top_;
      top_ = p;
    }


    template< int dim >
    inline typename ElementInfo< dim >::InstancePtr
    ElementInfo< dim >::Stack::null ()
    {
      return &null_;
    }

  } // namespace Alberta

} // namespace Dune

#endif // #if HAVE_ALBERTA

#endif // #ifndef DUNE_ALBERTA_ELEMENTINFO_HH
