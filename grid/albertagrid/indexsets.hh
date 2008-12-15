// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTAGRIDINDEXSETS_HH
#define DUNE_ALBERTAGRIDINDEXSETS_HH

#if HAVE_ALBERTA

#include <dune/common/stdstreams.hh>

#include <dune/grid/common/grid.hh>
#include <dune/grid/common/indexidset.hh>

#include <dune/grid/genericgeometry/misc.hh>

#include <dune/grid/albertagrid/albertaheader.hh>
#include <dune/grid/albertagrid/exceptions.hh>
#include <dune/grid/albertagrid/referencetopo.hh>
#include <dune/grid/albertagrid/dofadmin.hh>
#include <dune/grid/albertagrid/elementinfo.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< int dim, int dimworld >
  class AlbertaGrid;

  template< int codim, int dim, class GridImp >
  class AlbertaGridEntity;

  template< class GridType, int dim >
  struct MarkEdges;



  //! HierarchicIndexSet uses LeafIterator types for all codims and partition types
  template <class GridImp>
  struct AlbertaGridHierarchicIteratorTypes
  {
    //! The types of the iterator
    template<int cd>
    struct Codim
    {
      template<PartitionIteratorType pitype>
      struct Partition
      {
        /*
           We use the remove_const to extract the Type from the mutable class,
           because the const class is not instantiated yet.
         */
        typedef typename remove_const<GridImp>::type::Traits::template Codim<cd>::template Partition<pitype>::LeafIterator Iterator;
      };
    };
  };


  template< int dim, int dimworld >
  class AlbertaGridHierarchicIndexSet
    : public IndexSetDefaultImplementation
      < AlbertaGrid< dim, dimworld >,
          AlbertaGridHierarchicIndexSet< dim,dimworld >,
          AlbertaGridHierarchicIteratorTypes< AlbertaGrid< dim, dimworld > > >
  {
    typedef AlbertaGridHierarchicIndexSet< dim, dimworld > This;

    typedef AlbertaGrid< dim, dimworld > Grid;

    typedef typename remove_const< Grid >::type::Traits Traits;

    typedef typename Traits::template Codim< 0 >::Entity EntityCodim0Type;
    enum { numVecs  = AlbertHelp::numOfElNumVec };
    enum { numCodim = dim + 1 };

    friend class AlbertaGrid< dim, dimworld >;
    friend class MarkEdges< Grid, 3 >;
    friend class MarkEdges< const Grid, 3 >;

    AlbertaGridHierarchicIndexSet ( const Grid &grid )
      : grid_( grid )
    {}

  public:
    enum { ncodim = numCodim };

    //! return true if entity is contained in set
    template <class EntityType>
    bool contains(const EntityType &) const
    {
      // always true for this set
      return true;
    }

    //! return index of entity
    template< class Entity >
    int index ( const Entity &entity ) const
    {
      const int codim = Entity::codimension;
      const AlbertaGridEntity< codim, dim, const Grid > &entityImp
        = Grid::getRealImplementation( entity );
      Int2Type< dim-codim > dimVariable;
      return getIndex( entityImp.elementInfo().el(), entityImp.getFEVnum(), dimVariable );
    }

    //! return subIndex of given enitiy's sub entity
    template< int codim >
    int subIndex ( const EntityCodim0Type &entity, int i ) const
    {
      const AlbertaGridEntity< 0, dim, const Grid > &entityImp
        = Grid::getRealImplementation( entity );
      Int2Type< dim-codim > dimVariable;
      return getIndex( entityImp.elementInfo().el(), i, dimVariable );
    }

    //! return size of set for given GeometryType
    int size (GeometryType type) const
    {
      if( !type.isSimplex() ) return 0;
      return this->size(Grid::dimension-type.dim());
    }

    //! return size of set
    int size (int codim) const
    {
      return grid_.global_size(codim);
    }

    //! return geometry types this set has indices for
    const std::vector< GeometryType > & geomTypes(int codim) const
    {
      // returns all simplex
      return grid_.geomTypes(codim);
    }

#ifdef INDEXSET_HAS_ITERATORS
    /** @brief Iterator to one past the last entity of given codim for partition type
     */
    template<int cd, PartitionIteratorType pitype>
    typename AlbertaGridHierarchicIteratorTypes<Grid>::template Codim<cd>::
    template Partition<pitype>::Iterator end () const
    {
      return grid_.template leafend<cd,pitype> ();
    }

    /** @brief Iterator to first entity of given codimension and partition type.
     */
    template<int cd, PartitionIteratorType pitype>
    typename AlbertaGridHierarchicIteratorTypes<Grid>::template Codim<cd>::
    template Partition<pitype>::Iterator begin () const
    {
      return grid_.template leafbegin<cd,pitype> ();
    }
#endif

  private:
    // out grid
    const Grid &grid_;
    // constains the mapping from dune to alberta numbers
    const ALBERTA AlbertHelp :: AlbertaGridReferenceTopology<dim> refTopo_;
    // the vectors containing the numbers
    const int * elNumVec_[numVecs];

    // stores offset of dof numbers on given elements
    int nv_[numVecs];
    int dof_[numVecs];

    template< int codim >
    struct SetDofIdentifier
    {
      static void apply ( This &indexSet,
                          ALBERTA AlbertHelp::DOFVEC_STACK &dofvecs )
      {
        if( codim < numVecs )
        {
          const ALBERTA DOF_ADMIN *elAdmin_
            = dofvecs.elNumbers[ codim ]->fe_space->admin;
          // see Albert Doc. , should stay the same

          const int codimtype = Alberta::CodimType< dim, codim >::value;
          indexSet.nv_[ codim ] = elAdmin_->n0_dof[ codimtype ];
          assert( indexSet.nv_[ codim ] == 0 );
          indexSet.dof_[ codim ] = elAdmin_->mesh->node[ codimtype ];
        }
      }
    };

    // update vec pointer of the DOF_INT_VECs, which can change during resize
    void updatePointers(ALBERTA AlbertHelp::DOFVEC_STACK & dofvecs)
    {
      for(int i=0; i<numVecs; i++)
      {
        elNumVec_[i] = (dofvecs.elNumbers[i])->vec;
        assert(elNumVec_[i]);
      }
      GenericGeometry::ForLoop< SetDofIdentifier, 0, dim >::apply( *this, dofvecs );
    }

    // codim = 0 means we get from dim-cd = dim
    // this is the method for the element numbers
    // --element
    int getIndex ( const ALBERTA EL * el, int i , Int2Type<dim> fake ) const
    {
      enum { cd = 0 };
      assert(el);
      return elNumVec_[cd][ el->dof[ dof_[cd] ][nv_[cd]] ];
    }

    enum { cd1 = (dim == 2) ? 1 : 2 };
    // method for face numbers
    // codim = 0 means we get from dim-cd = dim
    // --face
    int getIndex ( const ALBERTA EL * el, int i , Int2Type<cd1> fake ) const
    {
      enum { cd = 1 };
      assert(el);
      assert( (dim == 2) || (dim == 3));
      // dof_[cd] marks the insertion point form which this dofs start
      // then i is the i-th dof
      return elNumVec_[cd][ el->dof[ dof_[cd] + i ][ nv_[cd] ] ];
    }

    enum { cd2 = (dim > 2) ? 1 : 6 };
    // codim = 0 means we get from dim-cd = dim
    // this method we have only in 3d, for edges
    // --edges
    int getIndex ( const ALBERTA EL * el, int i , Int2Type<cd2> fake ) const
    {
      enum { cd = 2 };
      assert(el);

      // dof_[cd] marks the insertion point form which this dofs start
      // then i is the i-th dof, here we addionally have to use the edge
      // mapping
      return elNumVec_[cd][ el->dof[ dof_[cd] + refTopo_.dune2albertaEdge(i) ][ nv_[cd] ] ];
    }

    // codim = dim  means we get from dim-cd = 0
    // return index of vertices
    // --vertex
    int getIndex ( const ALBERTA EL * el, int i , Int2Type<0> fake ) const
    {
      assert(el);
      return (el->dof[i][0]);
    }

    int getIndex ( const ALBERTA EL * el, int i , Int2Type<-1> fake ) const
    {
      assert(false);
      DUNE_THROW(AlbertaError,"Error, wrong codimension!\n");
      return -1;
    }
  }; // end class AlbertaGridHierarchicIndexSet


  //! hierarchic index set of AlbertaGrid
  template <int dim, int dimworld>
  class AlbertaGridIdSet :
    public IdSetDefaultImplementation < AlbertaGrid<dim,dimworld> ,
        AlbertaGridIdSet<dim,dimworld> , int >
  {
    typedef AlbertaGrid<dim,dimworld> GridType;
    typedef typename GridType :: HierarchicIndexSet HierarchicIndexSetType;

    // this means that only up to 300000000 entities are allowed
    enum { codimMultiplier = 300000000 };
    /*
       We use the remove_const to extract the Type from the mutable class,
       because the const class is not instantiated yet.
     */
    typedef typename remove_const<GridType>::type::Traits::template Codim<0>::Entity EntityCodim0Type;

    //! create id set, only allowed for AlbertaGrid
    AlbertaGridIdSet(const GridType & grid) : hset_(grid.hierarchicIndexSet())
    {
      for(int i=0; i<dim+1; i++)
        codimStart_[i] = i*codimMultiplier;
    }

    friend class AlbertaGrid<dim,dimworld>;
  public:
    //! export type of id
    typedef int IdType;

    //! return global id of given entity
    template <class EntityType>
    int id (const EntityType & ep) const
    {
      enum { cd = EntityType :: codimension };
      assert( hset_.size(cd) < codimMultiplier );
      return codimStart_[cd] + hset_.index(ep);
    }

    //! return global id of given entity
    template <int codim>
    int id (const typename GridType::template Codim<codim>::Entity& ep) const
    {
      //enum { cd = EntityType :: codimension };
      assert( hset_.size(codim) < codimMultiplier );
      return codimStart_[codim] + hset_.index(ep);
    }

    //! return subId of given entity
    template <int cd>
    int subId (const EntityCodim0Type & ep, int i) const
    {
      assert( hset_.size(cd) < codimMultiplier );
      return codimStart_[cd] + hset_.template subIndex<cd>(ep,i);
    }

  private:
    // our Grid
    const HierarchicIndexSetType & hset_;

    // store start of each codim numbers
    int codimStart_[dim+1];
  };
} // namespace Dune

#endif // HAVE_ALBERTA

#endif
