// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRAPEGRIDDISPLAY_HH
#define DUNE_GRAPEGRIDDISPLAY_HH

//- system includes
#include <list>
#include <set>
#include <stack>

//- Dune includes
#include <dune/common/stdstreams.hh>
#include <dune/grid/common/grid.hh>

#if HAVE_GRAPE
//-local includes
#include "grape/grapeinclude.hh"
#endif

/** @file
   @author Robert Kloefkorn
   @brief Provides a GridDisplay class using the GRAPE-HMesh Interface.
 */

namespace Dune
{

  //! the internal grape partition iterator types
  enum GrapePartitionIteratorType
  {
    g_Interior_Partition       = Interior_Partition,
    g_InteriorBorder_Partition = InteriorBorder_Partition,
    g_Overlap_Partition        = Overlap_Partition,
    g_OverlapFront_Partition   = OverlapFront_Partition,
    g_All_Partition            = All_Partition,
    g_Ghost_Partition          = Ghost_Partition
  };

  /** \brief the internal grape partition iterator types
      need to be these exact values to associate with combo button value. */
  enum GrapeIteratorType
  {
    g_LeafIterator       = 0,
    g_LevelIterator      = 1,
    g_HierarchicIterator = 2,
    g_GridPart           = 3
  };

  /** \todo Please doc me!
      \ingroup Grape
   */
  template<class GridType>
  class GrapeGridDisplay
  {
    typedef GrapeGridDisplay < GridType > MyDisplayType;
    typedef  MyDisplayType ThisType;
    enum { dim = GridType::dimension };
    enum { dimworld = GridType::dimensionworld };

#if HAVE_GRAPE
    // defined in griddisplay.hh
    typedef typename GrapeInterface<dim,dimworld>::DUNE_ELEM DUNE_ELEM;
    typedef typename GrapeInterface<dim,dimworld>::DUNE_DAT DUNE_DAT;
    typedef typename GrapeInterface<dim,dimworld>::DUNE_FDATA DUNE_FDATA;
    typedef typename GrapeInterface<dim,dimworld>::F_DATA F_DATA;
    typedef typename GrapeInterface<dim,dimworld>::STACKENTRY STACKENTRY;

    typedef typename std::stack < STACKENTRY * > StackEntryType;
    typedef void setGridPartIterators_t (DUNE_DAT * , void * gridPart);
  protected:
    //! store the actual element pointer
    DUNE_ELEM hel_;
    DUNE_DAT dune_;
    setGridPartIterators_t * setGridPartIter_;

    StackEntryType stackEntry_;
#endif

  public:
    typedef typename GridType::template Codim<0>:: HierarchicIterator
    HierarchicIteratorType;

    typedef typename GridType::Traits::LocalIdSet LocalIdSetType;
    typedef typename GridType::Traits::LeafIndexSet LeafIndexSetType;

  protected:
    //! the grid we want to display
    const GridType &grid_;

    //! true if we can use LevelIntersectionIterator
    const bool hasLevelIntersections_;

    void * gridPart_;

    //! leaf index set of the grid
    void * indexSet_;

    //! leaf index set of the grid
    const LocalIdSetType & lid_;

    //! my process number
    const int myRank_;

    // no better way than this canot export HMESH structure to here
    //! pointer to hmesh
    void * hmesh_;

    typedef std::list<HierarchicIteratorType *> HierarchicIteratorList;
    typedef typename HierarchicIteratorList::iterator ListIteratorType;
    HierarchicIteratorList hierList_;

  private:
    //! copy Constructor
    GrapeGridDisplay(const GrapeGridDisplay &);
  public:
    //! Constructor, make a GrapeGridDisplay for given grid
    inline GrapeGridDisplay(const GridType &grid, const int myrank = -1);

    //! Constructor, make a GrapeGridDisplay for given grid
    template <class GridPartType>
    inline GrapeGridDisplay(const GridPartType &gridPart, const int myrank = -1);

    template< class VT >
    inline GrapeGridDisplay ( const GridView< VT > &gridView, const int myrank = -1 );

    //! Destructor for GrapeGridDisplay
    inline ~GrapeGridDisplay();

    //! Calls the display of the grid and draws the discrete function
    //! if discretefunction is NULL, then only the grid is displayed
    inline void display();

    //! return rank of this display, for visualisation of parallel grid
    int myRank () const { return myRank_; }

    //! return reference to Dune Grid
    inline const GridType& getGrid() const ;

#if HAVE_GRAPE
    //! return pointer to Grape Hmesh
    inline void * getHmesh();

    /** \todo Please doc me! */
    DUNE_DAT & getDuneDat () { return dune_; }

    /** \todo Please doc me! */
    inline void addMyMeshToTimeScene(void * timescene, double time, int proc);

    bool hasData () { return false; }

  protected:
    // generate hmesh
    inline void * setupHmesh();

    inline void deleteHmesh();

    typedef typename GridType::template Codim<0>::Entity EntityCodim0Type;

    // type of index method
    typedef int EntityIndexFuncType (void * iset, const EntityCodim0Type & en);
    // type of vertex method
    typedef int VertexIndexFuncType (void * iset, const EntityCodim0Type & en, int vx);

    // pointer to index method
    const EntityIndexFuncType * entityIndex;
    // pointer to vertex method
    const VertexIndexFuncType * vertexIndex;

    // return element index from given index set
    template <class IndexSetType>
    static int getEntityIndex(void * iset, const EntityCodim0Type & en)
    {
      assert( iset );
      const IndexSetType * set = ((const IndexSetType *) iset);
      return (en.isLeaf()) ? set->index(en) : -1;
    }

    // return vertex index from given index set
    template< class IndexSet >
    static int getVertexIndex ( void *iset, const EntityCodim0Type &entity, int vx )
    {
      assert( iset != 0 );
      const IndexSet *indexSet = (const IndexSet *)iset;
      //return set->template subIndex< dim >( entity, vx );
      return indexSet->subIndex( entity, vx, dim );
    }

  public:
    //****************************************************************
    //
    // --GrapeGridDisplay, Some Subroutines needed for display with GRAPE
    //
    //****************************************************************
    // update element from entity
    template <class IntersectionIteratorType>
    inline void checkNeighbors(IntersectionIteratorType&,
                               const IntersectionIteratorType&, DUNE_ELEM *) ;

    // update element from entity
    template <class Entity>
    inline void el_update_base (Entity& en , DUNE_ELEM *) ;

    // update element from entity
    template <class EntityPointerType>
    inline int el_update (EntityPointerType *, DUNE_ELEM *) ;

    // update element from entity
    template <class EntityPointerType, class GridPartType>
    inline int el_update (EntityPointerType *, DUNE_ELEM *, GridPartType& );

    template< class EntityPointer, class VT >
    int el_update ( EntityPointer *, DUNE_ELEM *, const GridView< VT > & );

    // update child element
    template <class EntityPointerType>
    inline int child_update (EntityPointerType * , DUNE_ELEM *) ;

    template <class EntityPointerType>
    inline int child_n_update (EntityPointerType *, DUNE_ELEM *) ;

    // first and next macro element via LevelIterator level 0
    template <PartitionIteratorType pitype>
    inline int first_leaf (DUNE_ELEM * he) ;
    template <PartitionIteratorType pitype>
    inline int next_leaf (DUNE_ELEM * he) ;

    // first and next macro element via LevelIterator level 0
    template <class GridPartImp>
    inline int first_item (DUNE_ELEM * he) ;
    template <class GridPartImp>
    inline int next_item (DUNE_ELEM * he) ;

    // first and next macro element via LevelIterator level
    template <PartitionIteratorType pitype>
    inline int first_level (DUNE_ELEM * he, int level) ;

    template <PartitionIteratorType pitype>
    inline int next_level (DUNE_ELEM * he) ;

    // methods to call for combined display
    inline int firstMacro (DUNE_ELEM * elem) { return dune_.first_macro(elem); }
    inline int nextMacro  (DUNE_ELEM * elem) { return dune_.next_macro(elem);  }
    inline int firstChild (DUNE_ELEM * elem) { return (dune_.first_child) ? dune_.first_child(elem) : 0; }
    inline int nextChild  (DUNE_ELEM * elem) { return (dune_.next_child) ? dune_.next_child(elem) : 0; }

    // first and next child via HierarchicIterator with given maxlevel in Grape
    inline int first_child (DUNE_ELEM * he) ;
    inline int next_child (DUNE_ELEM * he) ;

  public:
    // delete leaf iterators
    template <PartitionIteratorType pitype>
    inline void delete_leaf (DUNE_ELEM * he) ;
    // delete level iterators
    template <PartitionIteratorType pitype>
    inline void delete_level (DUNE_ELEM * he) ;
    // delete level and hierarchical iterators
    template <PartitionIteratorType pitype>
    inline void delete_hier (DUNE_ELEM * he) ;

    // delete iterators
    template <class IteratorType>
    inline void delete_iterators(DUNE_ELEM * he) ;
  public:

    // fake function for copy iterator
    inline static void * copy_iterator (const void * i) ;

    // local to world
    inline void local2world (DUNE_ELEM * he, const double * c, double * w);

    // world to local
    inline int world2local (DUNE_ELEM * he, const double * w, double * c);

    // check inside reference element
    inline int checkWhetherInside (DUNE_ELEM * he, const double * w);

    //*********************************
    //  wrapper functions
    //*********************************
    // local to world
    inline static void ctow (DUNE_ELEM * he, const double * c, double * w);

    // world to local
    inline static int wtoc (DUNE_ELEM * he, const double * w, double * c);

    // check inside reference element
    inline static int check_inside (DUNE_ELEM * he, const double * w);

    // dito
    template< class Entity >
    int checkInside ( const Entity &entity, const double *w );

    // dito
    template< class Entity >
    int world_to_local ( const Entity &entity, const double *w, double *c );

    // dito
    template <class EntityType>
    inline void local_to_world(EntityType &en, const double * c, double * w);

    template <PartitionIteratorType pitype>
    inline void selectIterators(DUNE_DAT *, void *, setGridPartIterators_t *) const;

    inline void setIterationMethods(DUNE_DAT *, DUNE_FDATA * ) const;

    inline void changeIterationMethods(int iterType, int partType, DUNE_FDATA *);

    template <PartitionIteratorType pitype>
    struct IterationMethods
    {
      // wrapper methods for first_child and next_child
      inline static int first_mac (DUNE_ELEM * he)
      {
        MyDisplayType & disp = *((MyDisplayType *) he->display);
        return disp.template first_level<pitype>(he,0);
      }

      // wrapper methods for first_child and next_child
      inline static int first_lev (DUNE_ELEM * he)
      {
        MyDisplayType & disp = *((MyDisplayType *) he->display);
        return disp.template first_level<pitype>(he,he->level_of_interest);
      }

      inline static int next_lev  (DUNE_ELEM * he)
      {
        MyDisplayType & disp = *((MyDisplayType *) he->display);
        return disp.template next_level<pitype>(he);
      }

      // wrapper methods for first_child and next_child
      inline static int fst_leaf (DUNE_ELEM * he)
      {
        MyDisplayType & disp = *((MyDisplayType *) he->display);
        return disp.template first_leaf<pitype>(he);
      }
      inline static int nxt_leaf (DUNE_ELEM * he)
      {
        MyDisplayType & disp = *((MyDisplayType *) he->display);
        return disp.template next_leaf<pitype>(he);
      }

      // wrapper methods for first_child and next_child
      inline static int fst_child (DUNE_ELEM * he)
      {
        MyDisplayType & disp = *((MyDisplayType *) he->display);
        return disp.first_child(he);
      }
      inline static int nxt_child (DUNE_ELEM * he)
      {
        MyDisplayType & disp = *((MyDisplayType *) he->display);
        return disp.next_child(he);
      }

      // wrapper methods for deleting iterators
      inline static void del_leaf (DUNE_ELEM * he)
      {
        MyDisplayType & disp = *((MyDisplayType *) he->display);
        disp.template delete_leaf<pitype>(he);
      }

      // wrapper methods for deleting iterators
      inline static void del_level (DUNE_ELEM * he)
      {
        MyDisplayType & disp = *((MyDisplayType *) he->display);
        disp.template delete_level<pitype>(he);
      }

      // wrapper methods for deleting iterators
      inline static void del_hier (DUNE_ELEM * he)
      {
        MyDisplayType & disp = *((MyDisplayType *) he->display);
        disp.template delete_hier<pitype>(he);
      }

    };

  protected:
    template <class GridPartType>
    struct IterationMethodsGP
    {
      // wrapper methods for first_item and next_item
      inline static int fst_item (DUNE_ELEM * he)
      {
        assert( he->display );
        MyDisplayType & disp = *((MyDisplayType *) he->display);
        return disp.template first_item<GridPartType>(he);
      }
      inline static int nxt_item (DUNE_ELEM * he)
      {
        assert( he->display );
        MyDisplayType & disp = *((MyDisplayType *) he->display);
        return disp.template next_item<GridPartType>(he);
      }

      // delete iterators
      inline static void del_iter (DUNE_ELEM * he)
      {
        assert( he->display );
        MyDisplayType & disp = *((MyDisplayType *) he->display);
        typedef typename GridPartType :: template Codim<0> :: IteratorType IteratorType;
        disp.template delete_iterators<IteratorType> (he);
      }
    };

    template <class GridPartImp>
    struct SetIter
    {
      static void setGPIterator (DUNE_DAT * dune ,void * gridPart)
      {
        assert( gridPart );
        dune->gridPart = gridPart;
        dune->first_macro = &IterationMethodsGP<GridPartImp>::fst_item;
        dune->next_macro  = &IterationMethodsGP<GridPartImp>::nxt_item;
        dune->delete_iter = &IterationMethodsGP<GridPartImp>::del_iter;

        dune->first_child = 0;
        dune->next_child = 0;
      }
    };

    template< class ViewTraits >
    struct GridViewIterators
    {
      typedef Dune::GridView< ViewTraits > GridView;
      typedef typename GridView::template Codim< 0 >::Iterator Iterator;

      static int first_macro ( DUNE_ELEM *he )
      {
        assert( he->display != 0 );
        MyDisplayType &display = *static_cast< MyDisplayType * >( he->display );

        if( he->liter != 0 )
          display.template delete_iterators< Iterator >( he );

        assert( he->gridPart != 0 );
        const GridView &gridView = *static_cast< const GridView * >( he->gridPart );

        assert( he->liter   == 0 );
        assert( he->enditer == 0 );

        Iterator *it  = new Iterator( gridView.template begin< 0 >() );
        Iterator *end = new Iterator( gridView.template end  < 0 >() );

        he->liter   = it;
        he->enditer = end;

        if( *it == *end )
        {
          display.template delete_iterators< Iterator >( he );
          return 0;
        }

        return display.el_update( it, he, gridView );
      }

      static int next_macro ( DUNE_ELEM *he )
      {
        assert( he->display != 0 );
        MyDisplayType &display = *static_cast< MyDisplayType * >( he->display );

        assert( he->gridPart != 0 );
        const GridView &gridView = *static_cast< const GridView * >( he->gridPart );

        Iterator *it  = static_cast< Iterator * >( he->liter );
        Iterator *end = static_cast< Iterator * >( he->enditer );
        assert( (it != 0) && (end != 0) );

        ++(*it);
        if( *it == *end )
        {
          display.template delete_iterators< Iterator >( he );
          return 0;
        }

        return display.el_update( it, he, gridView );
      }

      static void delete_iter ( DUNE_ELEM *he )
      {
        assert( he->display );
        MyDisplayType &display = *static_cast< MyDisplayType * >( he->display );
        display.template delete_iterators< Iterator >( he );
      }

      static void set ( DUNE_DAT *dune, void *gridView )
      {
        assert( gridView );
        dune->gridPart = gridView;
        dune->first_macro = &first_macro;
        dune->next_macro = &next_macro;
        dune->delete_iter = &delete_iter;

        dune->first_child = 0;
        dune->next_child = 0;
      }
    };

    inline static void setIterationModus(DUNE_DAT * , DUNE_FDATA *);

  public:
    // create STACKENTRY or get from stack
    inline static void * getStackEntry(StackEntryType & stackEntry);

    // get StackEntry Wrapper
    inline static void * getStackEn(DUNE_DAT * dune);

    // free StackEntry Wrapper
    inline static void freeStackEn(DUNE_DAT * dune, void * entry);

    inline static void deleteStackEntry(StackEntryType &);

    // push STACKENTRY to stack
    inline static void freeStackEntry(StackEntryType & stackEntry, void * entry);
#endif
  }; // end class GrapeGridDisplay

#if HAVE_GRAPE
  /**************************************************************************/
  //  element types, see dune/grid/common/grid.hh
  // and also geldesc.hh for GR_ElementTypes
  enum GRAPE_ElementType
  {  g_vertex         = GrapeInterface_three_three::gr_vertex
     ,  g_line           = GrapeInterface_three_three::gr_line
     ,  g_triangle       = GrapeInterface_three_three::gr_triangle
     ,  g_quadrilateral  = GrapeInterface_three_three::gr_quadrilateral
     ,  g_tetrahedron    = GrapeInterface_three_three::gr_tetrahedron
     ,  g_pyramid        = GrapeInterface_three_three::gr_pyramid
     ,  g_prism          = GrapeInterface_three_three::gr_prism
     ,  g_hexahedron     = GrapeInterface_three_three::gr_hexahedron
     ,  g_iso_triangle   = GrapeInterface_three_three::gr_iso_triangle
     ,  g_iso_quadrilateral  = GrapeInterface_three_three::gr_iso_quadrilateral
     ,  g_unknown            = GrapeInterface_three_three::gr_unknown};

  //! convert dune geometry types to grape geometry types with numbers
  static inline GRAPE_ElementType convertToGrapeType ( GeometryType type , int dim )
  {
    if(dim < 3)
    {
      if(type.isTriangle()) return g_triangle;
      if(type.isQuadrilateral()) return g_quadrilateral;
      if(type.isVertex()) return g_vertex;
      if(type.isLine()) return g_line;
    }
    else
    {
      if(type.isTetrahedron()) return g_tetrahedron;
      if(type.isHexahedron()) return g_hexahedron;
      if(type.isPyramid()) return g_pyramid;
      if(type.isPrism()) return g_prism;
    }

    std::cerr << "No requested conversion for GeometryType " << type << "!\n";
    return g_unknown;
  }

  // see geldesc.hh for definition of this mapping
  // this is the same for all namespaces (two_two , and two_three, ...)
  static const int * const * vxMap = GrapeInterface_three_three::dune2GrapeVertex;
  static inline int mapDune2GrapeVertex( int geomType , int vx )
  {
    enum { usedTypes = GrapeInterface_three_three::numberOfUsedGrapeElementTypes };
    assert( geomType >= 0 );
    assert( geomType <  usedTypes ); // at the moment only defined from 2 to 7
    return vxMap[geomType][vx];
  }

  // see geldesc.hh for definition of this mapping
  // this is the same for all namespaces (two_two , and two_three, ...)
  static const int * const * faceMap = GrapeInterface_three_three::dune2GrapeFace;
  static inline int mapDune2GrapeFace( int geomType , int duneFace )
  {
    enum { usedTypes = GrapeInterface_three_three::numberOfUsedGrapeElementTypes };
    assert( geomType >= 0 );
    assert( geomType <  usedTypes ); // at the moment only defined from 2 to 7
    return faceMap[geomType][ duneFace ];
  }
#endif

} // end namespace Dune

#include "grape/grapegriddisplay.cc"

// undefs all defines
#include "grape/grape_undefs.hh"
#endif
