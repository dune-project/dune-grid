// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_TREEITERATOR_HH
#define DUNE_ALBERTA_TREEITERATOR_HH

#include <dune/grid/albertagrid/entitypointer.hh>

namespace Dune
{

  // AlbertaMarkerVector
  // -------------------

  //! Class to mark the Vertices on the leaf level
  //! to visit every vertex only once
  //! for the LevelIterator codim == dim
  class AlbertaMarkerVector
  {
    friend class AlbertaGrid<2,2>;
    friend class AlbertaGrid<2,3>;
    friend class AlbertaGrid<3,3>;

    enum { vxBufferSize_ = 10000 };
  public:
    //! create AlbertaMarkerVector for Level or Leaf Iterator, true == LevelIterator
    //! the vectors stored inside are empty first
    AlbertaMarkerVector (bool meLevel=true) : up2Date_(false), meLevel_(meLevel) {} ;

    //! return true if vertex is not watched on this element
    bool vertexNotOnElement(const int elIndex, const int vertex) const;

    //! return true if edge is not watched on this element
    bool edgeNotOnElement(const int elIndex, const int edge) const;

    //! return true if edge is not watched on this element
    bool faceNotOnElement(const int elIndex, const int face) const;

    //! mark vertices for LevelIterator and given level
    template <class GridType>
    void markNewVertices(GridType &grid, int level);

    //! mark vertices for LeafIterator , uses leaf level
    template <class GridType>
    void markNewLeafVertices(GridType &grid);

    //! return true if marking is up to date
    bool up2Date () const { return up2Date_; }

    //! unset up2date flag
    void unsetUp2Date () { up2Date_ = false; }

    //! print for debugin' only
    void print() const;

  private:
    typedef std::vector<int> ArrayType;

    // built in array to mark on which element a vertex is reached
    ArrayType vec_;
    ArrayType edgevec_;
    ArrayType facevec_;

    // number of vertices
    int numVertex_;

    // true is vertex marker is up to date
    bool up2Date_;
    bool meLevel_;
  };



  namespace AlbertaTreeIteratorHelp
  {
    template <class IteratorImp, int dim, int codim>
    struct GoNextEntity;
  }


  // AlbertaGridTreeIterator
  // -----------------------

  /*!
     Enables iteration over all entities of a given codimension and level of a grid.
   */
  template< int codim, PartitionIteratorType pitype, class GridImp >
  class AlbertaGridTreeIterator
    : public AlbertaGridEntityPointer< codim, GridImp >
  {
    typedef AlbertaGridTreeIterator< codim, pitype, GridImp > This;
    typedef AlbertaGridEntityPointer< codim, GridImp > Base;

    enum { dim = GridImp::dimension };
    friend class AlbertaGridEntity<2,dim,GridImp>;
    friend class AlbertaGridEntity<1,dim,GridImp>;
    friend class AlbertaGridEntity<0,dim,GridImp>;
    friend class AlbertaGrid < dim , GridImp::dimensionworld >;


    typedef AlbertaGridTreeIterator< codim, pitype, GridImp > AlbertaGridTreeIteratorType;
    friend class AlbertaTreeIteratorHelp::GoNextEntity< This, dim, codim >;

  public:
    typedef typename GridImp::template Codim< codim >::Entity Entity;
    typedef MakeableInterfaceObject< Entity > EntityObject;
    typedef typename EntityObject::ImplementationType EntityImp;

    //! Constructor making end iterator
    AlbertaGridTreeIterator ( const This &other );

    //! Constructor making end iterator
    This &operator= ( const This &other );

    //! Constructor making end iterator
    AlbertaGridTreeIterator ( const GridImp &grid,
                              int travLevel,
                              int proc,
                              bool leafIt=false );

    //! Constructor making begin iterator
    AlbertaGridTreeIterator ( const GridImp & grid,
                              const AlbertaMarkerVector * vec,
                              int travLevel,
                              int proc,
                              bool leafIt=false);

    //! increment
    void increment();
    //! equality

  private:
    // private Methods
    void makeIterator();

    ALBERTA EL_INFO * goFirstElement(ALBERTA TRAVERSE_STACK *stack,
                                     ALBERTA MESH *mesh,
                                     int level, ALBERTA FLAGS fill_flag);
    ALBERTA EL_INFO * traverseElLevel(ALBERTA TRAVERSE_STACK * stack);
    ALBERTA EL_INFO * traverseElLevelInteriorBorder(ALBERTA TRAVERSE_STACK * stack);
    ALBERTA EL_INFO * traverseElLevelGhosts(ALBERTA TRAVERSE_STACK * stack);

    // the default is, go to next elInfo
    //template <int cc>
    ALBERTA EL_INFO * goNextEntity(ALBERTA TRAVERSE_STACK *stack,ALBERTA EL_INFO *elinfo_old);

    // the real go next methods
    ALBERTA EL_INFO * goNextElInfo(ALBERTA TRAVERSE_STACK *stack,ALBERTA EL_INFO *elinfo_old);
    ALBERTA EL_INFO * goNextFace(ALBERTA TRAVERSE_STACK *stack,
                                 ALBERTA EL_INFO *elInfo);
    ALBERTA EL_INFO * goNextEdge(ALBERTA TRAVERSE_STACK *stack,
                                 ALBERTA EL_INFO *elInfo);
    ALBERTA EL_INFO * goNextVertex(ALBERTA TRAVERSE_STACK *stack,
                                   ALBERTA EL_INFO *elInfo);

    // search next macro el
    ALBERTA MACRO_EL * nextGhostMacro(ALBERTA MACRO_EL *mel);

    //! level :)
    int level_;

    //! level :)
    int enLevel_;

    //! reference to entity of entity pointer class
    EntityImp & virtualEntity_;

    // contains ALBERTA traversal stack
    ALBERTA ManageTravStack manageStack_;

    //! Number of the subentity within the element
    int subEntity_;

    // knows on which element a point,edge,face is viewed
    const AlbertaMarkerVector * vertexMarker_;

    // variable for operator++
    bool okReturn_;

    // store processor number of elements
    // for ghost walktrough, i.e. walk over ghosts which belong
    // tp processor 2
    const int proc_;
  };

}

#endif
