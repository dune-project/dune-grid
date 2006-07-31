// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU2DGRIDGRID_HH
#define DUNE_ALU2DGRIDGRID_HH

//- System includes
#include "alu2dinclude.hh"
#include <vector>

//- Dune includes
#include <dune/grid/common/capabilities.hh>
#include <dune/common/interfaces.hh>

#include <dune/grid/common/grid.hh>
#include <dune/grid/common/referenceelements.hh>
#include <dune/grid/common/defaultindexsets.hh>
#include <dune/grid/common/sizecache.hh>
#include <dune/common/collectivecommunication.hh>
#include <dune/grid/common/intersectioniteratorwrapper.hh>

//- Local includes
#include "indexsets.hh"
#include "memory.hh"

namespace Dune {

  typedef double alu2d_ctype;

  // Forward declarations
  template<int cd, int dim, class GridImp>
  class ALU2dGridEntity;
  template<int cd, PartitionIteratorType pitype, class GridImp >
  class ALU2dGridLevelIterator;
  template<int cd, class GridImp >
  class ALU2dGridEntityPointer;
  template<int mydim, int coorddim, class GridImp>
  class ALU2dGridMakeableGeometry;
  template<int mydim, int cdim, class GridImp>
  class ALU2dGridGeometry;
  template<class GridImp>
  class ALU2dGridHierarchicIterator;
  //template<class GridImp>
  //class ALU2dGridIntersectionIterator;
  template<class GridImp>
  class ALU2dGridIntersectionBase;
  template<class GridImp>
  class ALU2dGridLevelIntersectionIterator;
  template<class GridImp>
  class ALU2dGridLeafIntersectionIterator;
  template<int codim, PartitionIteratorType pitype, class GridImp>
  class ALU2dGridLeafIterator;
  template <int mydim, int coorddim, class GridImp>
  class ALU2dGridMakeableEntity;
  template <class GridImp>
  class ALU2dGridFaceGeometryInfo;
  template<int dim, int dimworld>
  class ALU2dGridGlobalIdSet;
  template<int dim, int dimworld>
  class ALU2dGridLocalIdSet;
  template<int dim, int dimworld>
  class ALU2dGridHierarchicIndexSet;
  template <class EntityImp>
  class ALUMemoryProvider;
  template<int dim, int dimworld>
  class ALU2dGrid;
  template <class GeometryImp, int nChild>
  class ALU2DLocalGeometryStorage;
  template <class GridImp, int codim>
  struct ALU2dGridEntityFactory;



  //**********************************************************************
  //
  // --ALU2dGrid
  // --Grid
  //
  //**********************************************************************
  template <int dim, int dimworld>
  struct ALU2dGridFamily
  {
    //! Type of the global id set
    typedef ALU2dGridLocalIdSet<dim,dimworld> GlobalIdSetImp;

    //! Type of the local id set
    typedef ALU2dGridLocalIdSet<dim,dimworld> LocalIdSetImp;

    //! Type of the level index set
    typedef DefaultLevelIndexSet<ALU2dGrid < dim, dimworld > >  LevelIndexSetImp;
    //! Type of the leaf index set
    typedef DefaultLeafIndexSet<ALU2dGrid < dim, dimworld > >   LeafIndexSetImp;

    typedef int GlobalIdType;
    typedef int LocalIdType;

    typedef ALU2dGrid<dim,dimworld> GridImp;

    struct Traits
    {
      typedef ALU2dGrid<dim,dimworld> Grid;

      typedef Dune::IntersectionIterator<const GridImp, LeafIntersectionIteratorWrapper> IntersectionIterator;
      typedef Dune::IntersectionIterator<const GridImp, LeafIntersectionIteratorWrapper> LeafIntersectionIterator;
      typedef Dune::IntersectionIterator<const GridImp, LevelIntersectionIteratorWrapper> LevelIntersectionIterator;

      typedef Dune::HierarchicIterator<const GridImp, ALU2dGridHierarchicIterator> HierarchicIterator;

      template <int cd>
      struct Codim
      {
        // IMPORTANT: Codim<codim>::Geometry == Geometry<dim-codim,dimw>
        typedef Dune::Geometry<dim-cd, dimworld, const GridImp, ALU2dGridGeometry> Geometry;
        typedef Dune::Geometry<dim-cd, dim, const GridImp, ALU2dGridGeometry> LocalGeometry;
        // we could - if needed - introduce an other struct for dimglobal of Geometry
        typedef Dune::Entity<cd, dim, const GridImp, ALU2dGridEntity> Entity;

        typedef Dune::LevelIterator<cd,All_Partition,const GridImp,ALU2dGridLevelIterator> LevelIterator;

        typedef Dune::LeafIterator<cd,All_Partition,const GridImp,ALU2dGridLeafIterator> LeafIterator;

        typedef Dune::EntityPointer<const GridImp,ALU2dGridEntityPointer<cd,const GridImp> > EntityPointer;

        template <PartitionIteratorType pitype>
        struct Partition
        {
          typedef Dune::LevelIterator<cd,pitype,const GridImp,ALU2dGridLevelIterator> LevelIterator;
          typedef Dune::LeafIterator<cd,pitype,const GridImp,ALU2dGridLeafIterator> LeafIterator;
        };

      };

      typedef IndexSet<GridImp,LevelIndexSetImp,DefaultLevelIteratorTypes<GridImp> > LevelIndexSet;
      typedef LeafIndexSetImp LeafIndexSet;
      typedef IdSet<GridImp,GlobalIdSetImp,GlobalIdType> GlobalIdSet;
      typedef IdSet<GridImp,LocalIdSetImp,LocalIdType> LocalIdSet;

      typedef CollectiveCommunication<GridImp> CollectiveCommunication;
    };
  }; // end of ALU2dGridFamily


  /**
     \brief [<em> provides \ref Dune::Grid </em>]
     \brief 2D grid, will provide non conform grids
     @ingroup GridImplementations
     The ALU2dGrid implements the Dune GridInterface for 2d meshes.
     This grid can be locally adapted and will provide non conform grids.

     @note
     Adaptive grid, written mainly by Bernard Schupp.
     This grid supports non conform grids.

     (see ALUGrid homepage: http://www.mathematik.uni-freiburg.de/IAM/Research/alugrid/)

   */

  template <int dim, int dimworld>
  class ALU2dGrid :
    public GridDefaultImplementation<dim,dimworld,alu2d_ctype,ALU2dGridFamily<dim,dimworld> >,
    public HasObjectStream
  {
    CompileTimeChecker<(dim      == 2)> ALU2dGrid_only_implemented_for_2dp;
    CompileTimeChecker<(dimworld == 2)> ALU2dGrid_only_implemented_for_2dw;

    typedef ALU2dGrid<dim,dimworld> ThisType;

    friend class ALU2dGridEntity<0,dim,const ThisType>;
    friend class ALU2dGridEntity<1,dim,const ThisType>;
    friend class ALU2dGridEntity<dim,dim,const ThisType>;

    friend class ALU2dGridGeometry<0,dimworld,const ThisType>;
    friend class ALU2dGridGeometry<1,dimworld,const ThisType>;
    friend class ALU2dGridGeometry<dim,dimworld,const ThisType>;

    friend class ALU2dGridEntityPointer<0,const ThisType>;
    friend class ALU2dGridEntityPointer<1,const ThisType>;
    friend class ALU2dGridEntityPointer<dim,const ThisType>;

    friend class ALU2dGridHierarchicIndexSet<dim,dimworld>;
    //friend class ALU2dGridGlobalIdSet<dim,dimworld>;
    //friend class ALU2dGridLocalIdSet<dim,dimworld>;

    friend class ALU2dGridIntersectionBase < const ThisType > ;
    friend class ALU2dGridLevelIntersectionIterator< const ThisType > ;
    friend class ALU2dGridLeafIntersectionIterator< const ThisType > ;

    //**********************************************************
    // The Interface Methods
    //**********************************************************
  public:
    //! my Traits class
    typedef typename ALU2dGridFamily < dim , dimworld > :: Traits Traits;

  protected:

    typedef MakeableInterfaceObject<typename Traits::template
        Codim<0>::Geometry> GeometryObject;
    friend class ALU2DLocalGeometryStorage<GeometryObject, 4 >;

  public:
    //! my Traits class
    typedef ALU2dGridFamily < dim , dimworld > GridFamily;
    //! Type of the hierarchic index set
    typedef ALU2dGridHierarchicIndexSet<dim,dimworld> HierarchicIndexSet;

    //! Type of the local id set
    typedef ALU2dGridLocalIdSet<dim,dimworld> LocalIdSetImp;
    typedef LocalIdSetImp GlobalIdSetImp;

    //! Type of the global id set
    typedef typename Traits :: GlobalIdSet GlobalIdSet;

    //! Type of the local id set
    typedef typename Traits :: LocalIdSet LocalIdSet;


    //! Type of the level index set
    typedef typename GridFamily :: LevelIndexSetImp LevelIndexSetImp;
    //! Type of the leaf index set
    typedef typename GridFamily :: LeafIndexSetImp LeafIndexSetImp;

    //! a standard leaf iterator
    typedef ALU2dGridLeafIterator<0, All_Partition, const ThisType> LeafIteratorImp;
    typedef typename Traits::template Codim<0>::LeafIterator LeafIteratorType;
    typedef typename Traits::template Codim<0>::LeafIterator LeafIterator;

    //! a standard leaf iterator
    typedef ALU2dGridLevelIterator<0, All_Partition, const ThisType> LevelIteratorImp;
    typedef typename Traits::template Codim<0>::LevelIterator LevelIteratorType;
    typedef typename Traits::template Codim<0>::LevelIterator LevelIterator;

    typedef ALU2dGridHierarchicIterator<ThisType> HierarchicIteratorImp;

    typedef typename Traits::CollectiveCommunication CollectiveCommunicationType;

    //! maximal number of levels
    enum { MAXL = 64 };

  protected:
    //! Constructor which reads an ALU2dGrid Macro Triang file
    //! or given GridFile
    ALU2dGrid(std::string macroTriangFilename );

  public:
    //! Desctructor
    ~ALU2dGrid();

    //! for grid identification
    std::string name () const;

    //! Return maximum level defined in this grid. Levels are numbered
    //! 0 ... maxLevel with 0 the coarsest level.
    int maxLevel() const;

    //! --Leveliterator
    //! Iterator to first entity of given codim on level
    template<int cd, PartitionIteratorType pitype>
    typename Traits::template Codim<cd>::template Partition<pitype>::LevelIterator
    lbegin (int level) const;

    //! one past the end on this level
    template<int cd, PartitionIteratorType pitype>
    typename Traits::template Codim<cd>::template Partition<pitype>::LevelIterator
    lend (int level) const;

    //! Iterator to first entity of given codim on level
    template<int cd>
    typename Traits::template Codim<cd>::
    template Partition<All_Partition>::LevelIterator
    lbegin (int level) const;

    //! one past the end on this level
    template<int cd>
    typename Traits::template Codim<cd>::
    template Partition<All_Partition>::LevelIterator
    lend (int level) const;

    //! Iterator to first entity of codim 0 on level
    LevelIteratorType lbegin (int level) const;

    //! last entity of codim 0 on level
    LevelIteratorType lend (int level) const;

    //! General definiton for a leaf iterator
    template <int codim, PartitionIteratorType pitype>
    typename Traits::template Codim<codim>::template Partition<pitype>::LeafIterator
    leafbegin() const;

    //! General definition for an end iterator on leaf level
    template <int codim, PartitionIteratorType pitype>
    typename Traits::template Codim<codim>::template Partition<pitype>::LeafIterator
    leafend() const;

    //! General definiton for a leaf iterator
    template <int codim>
    typename Traits::template Codim<codim>::LeafIterator
    leafbegin() const;

    //! General definition for an end iterator on leaf level
    template <int codim>
    typename Traits::template Codim<codim>::LeafIterator
    leafend() const;

    //! Iterator to first entity of codim 0 on leaf level (All_Partition)
    LeafIteratorType leafbegin () const;

    //! one past the end on this leaf level (codim 0 and All_Partition)
    LeafIteratorType leafend () const;

    //! number of grid entities per level and codim
    int size (int level, int cd) const;

    //! number of leaf entities per codim in this process
    int size (int codim) const;

    //! number of entities per level, codim and geometry type in this process
    int size (int level, GeometryType type) const;

    //! number of leaf entities per codim and geometry type in this process
    int size (GeometryType type) const;

    //! deliver all geometry types used in this grid
    const std::vector<GeometryType>& geomTypes (int codim) const { return geomTypes_[codim]; }


    //****************************************************************
    // index and id sets
    //****************************************************************

    //! get global id set of grid
    const GlobalIdSet & globalIdSet () const {
      return localIdSet();
    }

    //! get global id set of grid
    const LocalIdSet & localIdSet () const {
      return localIdSet_;
    }

    //! number of grid entities in the entire grid for given codim
    int hierSetSize (int cd) const;

    //! get hierarchic index set of the grid
    const HierarchicIndexSet & hierarchicIndexSet () const { return hIndexSet_; }

    //! Not yet implemented...
    //! get leaf index set of the grid
    const typename Traits :: LeafIndexSet & leafIndexSet () const;

    //! Not yet implemented...
    //! get level index set of the grid
    const typename Traits :: LevelIndexSet & levelIndexSet (int level) const;


    //**********************************************************
    // End of Interface Methods
    //**********************************************************

    //! return reference to org ALU2dGrid
    //! private method, but otherwise we have to friend class all possible
    //! types of LevelIterator ==> later
    ALU2DSPACE Hmesh & myGrid() { return mesh_; }
    ALU2DSPACE Hmesh & myGrid() const { return mesh_; }

    //! refine grid refCount times
    bool globalRefine(int refCount);

    //! returns if a least one entity was marked for coarsening
    bool preAdapt ( );

    //! clear all entity new markers
    void postAdapt ( );

    /**! refine all positive marked leaf entities,
       return true if a least one entity was refined
     */
    bool adapt ( );

    //! refine grid
    bool refineGrid();

    //! mark entities for refinement or coarsening, refCount < 0 will mark
    //! the entity for one coarsen step and refCount > 0 will mark for one
    //! refinement, one refinement will create 8 children per element
    bool mark( int refCount , const typename Traits::template Codim<0>::EntityPointer & ep );

    //! return dummy communication
    const CollectiveCommunicationType & comm() const { return comm_; }
  private:
    CollectiveCommunicationType comm_;

    bool mark( int refCount , const typename Traits::template Codim<0>::Entity & en );

    void updateStatus();

    void calcMaxlevel();

    void calcExtras();

    // fake dummy
    typedef int EntityProviderType;
    int entityProvider_;

    template <int codim>
    MakeableInterfaceObject<typename Traits::template Codim<codim>:: Entity> * getNewEntity ( int level ) const
    {
      return ALU2dGridEntityFactory<ThisType,codim>::getNewEntity(*this,entityProvider_,level);
    }

    template <int codim>
    void freeEntity (MakeableInterfaceObject<typename Traits::template Codim<codim>:: Entity> * en) const
    {
      return ALU2dGridEntityFactory<ThisType,codim>::freeEntity(entityProvider_, en);
    }

    // create GeomTypes
    void makeGeomTypes ();
    friend class Conversion<ALU2dGrid<dim, dimworld>, HasObjectStream>;
    friend class Conversion<const ALU2dGrid<dim, dimworld>, HasObjectStream>;

    //! Copy constructor should not be used
    ALU2dGrid( const ThisType & g );

    //! assignment operator should not be used
    ThisType & operator = (const ThisType & g);

    //! the real grid
    mutable ALU2DSPACE Hmesh mesh_;

    //! the hierarchic index set
    HierarchicIndexSet hIndexSet_;

    //! out global id set
    LocalIdSetImp localIdSet_;

    //! the level index set ( default type )
    mutable std::vector < LevelIndexSetImp * > levelIndexVec_;

    // at the moment the number of different geom types is 1
    enum { numberOfGeomTypes = 1 };
    std::vector< std::vector<GeometryType> > geomTypes_;

    //! the leaf index set
    mutable LeafIndexSetImp * leafIndexSet_;

    int maxLevel_;
    int refineMarked_ , coarsenMarked_;

    //typedef ALU2dGridVertexList VertexListType;
    //mutable VertexListType vertexList_[MAXL];

    //! the type of our size cache
    typedef SingleTypeSizeCache<ThisType> SizeCacheType;
    SizeCacheType * sizeCache_;

    // new intersection iterator is a wrapper which get itersectioniteratoimp as pointers
  public:
    typedef ALU2dGridLeafIntersectionIterator <const ThisType>  LeafIntersectionIteratorImp;
    typedef ALU2dGridLevelIntersectionIterator<const ThisType> LevelIntersectionIteratorImp;

    typedef ALUMemoryProvider< LeafIntersectionIteratorImp > LeafIntersectionIteratorProviderType;
    typedef ALUMemoryProvider< LevelIntersectionIteratorImp > LevelIntersectionIteratorProviderType;

  private:
    friend class LeafIntersectionIteratorWrapper< const ThisType > ;
    friend class LevelIntersectionIteratorWrapper< const ThisType > ;
    // return reference to intersectioniterator storage
    LeafIntersectionIteratorProviderType & leafIntersetionIteratorProvider() const
    {
      return leafInterItProvider_;
    }
    mutable LeafIntersectionIteratorProviderType leafInterItProvider_;

    LevelIntersectionIteratorProviderType & levelIntersetionIteratorProvider() const
    {
      return levelInterItProvider_;
    }
    mutable LevelIntersectionIteratorProviderType levelInterItProvider_;

    mutable ALU2dGridMarkerVector marker_[MAXL];
  public:
    ALU2dGridMarkerVector & getMarkerVector(int level) const
    {
      assert( level >= 0);
      assert( level <= MAXL);
      return marker_[level];
    }

  }; // end class ALU2dGrid

  template <class GridImp, int codim>
  struct ALU2dGridEntityFactory
  {
    typedef typename GridImp:: template Codim<codim>::Entity Entity;
    typedef MakeableInterfaceObject<Entity> EntityObj;
    typedef ALU2dGridEntity<codim,GridImp::dimension,GridImp> EntityImp;

    template <class EntityProviderType>
    static EntityObj *
    getNewEntity (const GridImp & grid, EntityProviderType & ep, int level)
    {
      return new EntityObj(EntityImp(grid, level));
    }

    template <class EntityProviderType>
    static void freeEntity( EntityProviderType & ep, EntityObj * e )
    {
      //if( e )
      delete e;
    }
  };

  template <class GridImp>
  struct ALU2dGridEntityFactory<GridImp,0>
  {
    enum {codim =0};
    typedef typename GridImp:: template Codim<codim>::Entity Entity;
    typedef MakeableInterfaceObject<Entity> EntityObj;
    //typedef typename Entity :: ImplementationType EntityImp;

    template <class EntityProviderType>
    static EntityObj *
    getNewEntity (const GridImp & grid, EntityProviderType & ep, int level)
    {
      //return ep.template getEntityObject( grid, level, (EntityImp *) 0)
      return ep.getObject( grid, level);
    }

    template <class EntityProviderType>
    static void freeEntity( EntityProviderType & ep, EntityObj * e )
    {
      ep.freeObject( e );
    }
  };


  namespace Capabilities
  {
    template<int dim,int dimw>
    struct hasLeafIterator< ALU2dGrid<dim,dimw> >
    {
      static const bool v = true;
    };

    template<int dim, int dimw, int cdim>
    struct hasEntity<ALU2dGrid<dim,dimw>, cdim >
    {
      static const bool v = true;
    };

    template<int dim, int dimw>
    struct isLevelwiseConforming< ALU2dGrid<dim,dimw> >
    {
      static const bool v = false;
    };

    template<int dim, int dimw>
    struct hasHangingNodes< ALU2dGrid<dim,dimw> >
    {
      static const bool v = false;
    };
  } // end namespace Capabilities

} // end namespace Dune

#include "entity.hh"
#include "geometry.hh"
#include "iterator.hh"

#include "grid_imp.cc"

#endif
