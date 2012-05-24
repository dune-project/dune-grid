// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALUGRID_ENTITY_CC
#define DUNE_ALUGRID_ENTITY_CC

#if COMPILE_ALUGRID_INLINE == 0
#include <config.h>
#endif

#include "alu3dinclude.hh"
#include <dune/grid/alugrid/3d/alugrid.hh>
#include "entity.hh"

#if COMPILE_ALUGRID_INLINE == 0
#include <dune/grid/alugrid/3d/alu3dgridfactory.hh>
#endif
#include <dune/grid/alugrid/common/geostorage.hh>

#if COMPILE_ALUGRID_INLINE
#define alu_inline inline
#else
#define alu_inline
#endif

namespace Dune {

  // --Entity
  template <int cd, int dim, class GridImp>
  ALU3dGridEntity<cd,dim,GridImp> ::
  ALU3dGridEntity(const FactoryType& factory, int level) :
    factory_( factory ),
    item_(0),
    level_(0),
    gIndex_(-1),
    twist_(0),
    face_(-1),
    partitionType_(InteriorEntity)
  {}

  // --Entity
  template <int cd, int dim, class GridImp>
  alu_inline ALU3dGridEntity<cd,dim,GridImp> ::
  ALU3dGridEntity(const ALU3dGridEntity<cd,dim,GridImp> & org) :
    factory_(org.factory_),
    item_(org.item_),
    level_(org.level_),
    gIndex_(org.gIndex_),
    twist_(org.twist_),
    face_(org.face_),
    partitionType_(org.partitionType_)
  {}

  template<int cd, int dim, class GridImp>
  alu_inline void ALU3dGridEntity<cd,dim,GridImp> ::
  setEntity(const ALU3dGridEntity<cd,dim,GridImp> & org)
  {
    item_   = org.item_;
    gIndex_ = org.gIndex_;
    twist_  = org.twist_;
    level_  = org.level_;
    face_   = org.face_;
    partitionType_ = org.partitionType_;
  }

  template<int cd, int dim, class GridImp>
  alu_inline void ALU3dGridEntity<cd,dim,GridImp> ::
  setElement(const HItemType & item)
  {
    setElement(item,GetLevel<GridImp,cd>::getLevel(grid(),item));
  }

  template<int cd, int dim, class GridImp>
  alu_inline void ALU3dGridEntity<cd,dim,GridImp> ::
  setElement(const EntitySeed& seed )
  {
    setElement(*seed.item(), seed.level(), seed.twist(), seed.face());
  }

  template<int cd, int dim, class GridImp>
  alu_inline void ALU3dGridEntity<cd,dim,GridImp> ::
  setElement(const HItemType & item, const int level, int twist , int face )
  {
    // cast interface to implementation
    item_   = static_cast<const ItemType *> (&item);
    gIndex_ = (*item_).getIndex();
    twist_  = twist;
    level_  = level;
    face_   = face;
    partitionType_ = this->convertBndId( *item_ );

    // reset geometry information
    geo_.invalidate();
  }

  template<int cd, int dim, class GridImp>
  alu_inline void ALU3dGridEntity<cd,dim,GridImp> ::
  setGhost(const HBndSegType &ghost)
  {
    // this method only exists, that we don't have to specialize the
    // Iterators for each codim, this method should not be called otherwise
    // error
    DUNE_THROW(GridError,"This method should not be called!");
  }

  template<int cd, int dim, class GridImp>
  alu_inline PartitionType ALU3dGridEntity<cd,dim,GridImp> ::
  convertBndId(const HItemType & item) const
  {
    if(item.isGhost())
    {
      return GhostEntity;
    }
    else if(item.isBorder())
    {
      return BorderEntity;
    }
    else
    {
      assert( item.isInterior() );
      return InteriorEntity;
    }
  }

  template< int cd, int dim, class GridImp >
  alu_inline typename ALU3dGridEntity< cd, dim, GridImp >::Geometry
  ALU3dGridEntity< cd, dim, GridImp >::geometry () const
  {
    //assert( (cd == 1) ? (face_ >= 0) : 1 );
    if( ! geo_.valid() )
      geo_.buildGeom( *item_, twist_, face_ );
    return Geometry( geo_ );
  }

  /////////////////////////////////////////////////
  //
  //  --Entity0
  //  --Codim0Entity
  //
  /////////////////////////////////////////////////

  template<int dim, class GridImp>
  alu_inline ALU3dGridEntity<0,dim,GridImp> ::
  ALU3dGridEntity(const FactoryType &factory, int wLevel)
    : factory_( factory )
      , item_( 0 )
      , ghost_( 0 )
      , level_(-1)
      , isLeaf_ (false)
  {  }

  template<int dim, class GridImp>
  alu_inline ALU3dGridEntity<0,dim,GridImp> ::
  ALU3dGridEntity(const ALU3dGridEntity<0,dim,GridImp> & org)
    : factory_(org.factory_)
      , item_(org.item_)
      , ghost_( org.ghost_ )
      , level_(org.level_)
      , isLeaf_ (org.isLeaf_)
  {  }

  template< int dim, class GridImp >
  alu_inline typename ALU3dGridEntity< 0, dim, GridImp >::Geometry
  ALU3dGridEntity< 0, dim, GridImp >::geometry () const
  {
    assert(item_ != 0);
    if( ! geo_.valid() )
      geo_.buildGeom( *item_ );
    return Geometry( geo_ );
  }

  template< int dim, class GridImp >
  alu_inline typename ALU3dGridEntity<0,dim,GridImp>::LocalGeometry
  ALU3dGridEntity< 0, dim, GridImp >::geometryInFather () const
  {
    assert( item_ );
    // this method should only be called if a father exists
    assert( item_->up() );

    // get child number
    const int child = item_->nChild();

    // if the rule of the farher is not refine_element, it has to be bisection
    // this can only be true for tetrahedral elements
    if( (GridImp::elementType == tetra) && (item_->up()->getrule() != ImplTraits::refine_element_t) )
    {
      static LocalGeometryImpl geom;
      geom.buildGeomInFather( father()->geometry(), geometry() );
      return LocalGeometry( geom );
    }
    else
    {
      // to be improved, when we using not the refine 8 rule
      // see dune/grid/alugrid/common/geostrage.hh for implementation
      typedef ALULocalGeometryStorage< GridImp, LocalGeometryImpl, 8 > GeometryStroage ;
      return LocalGeometry( GeometryStroage::geom( type(), true, child ) );
    }
  }

  //********* begin method subIndex ********************
  // partial specialisation of subIndex
  template <class IMPLElemType, ALU3dGridElementType type, int codim>
  struct IndexWrapper {};

  // specialisation for vertices
  template <class IMPLElemType, ALU3dGridElementType type>
  struct IndexWrapper<IMPLElemType, type, 3>
  {
    typedef ElementTopologyMapping<type> ElemTopo;

    static int subIndex(const IMPLElemType &elem, int i)
    {
      return elem.myvertex( ElemTopo::dune2aluVertex(i) )->getIndex(); // element topo
    }
  };

  // specialisation for faces
  template <class IMPLElemType, ALU3dGridElementType type>
  struct IndexWrapper<IMPLElemType, type , 1>
  {
    static int subIndex(const IMPLElemType &elem, int i)
    {
      // is specialised for each element type and uses
      // the dune2aluFace mapping
      return (getFace(elem,i))->getIndex();
    }
  };

  // specialisation for edges
  template <class IMPLElemType, ALU3dGridElementType type>
  struct IndexWrapper<IMPLElemType, type, 2>
  {
    typedef ElementTopologyMapping<type> ElemTopo;

    // return subIndex of given edge
    static int subIndex(const IMPLElemType &elem, int i)
    {
      // get hedge1 corresponding to dune reference element and return number
      return elem.myhedge1( ElemTopo::dune2aluEdge(i) )->getIndex();
    }
  };

  // specialisation for elements
  template <class IMPLElemType, ALU3dGridElementType type>
  struct IndexWrapper<IMPLElemType, type, 0>
  {
    static int subIndex(const IMPLElemType &elem, int i) {
      // just return the elements index
      return elem.getIndex();
    }
  };

  template<int dim, class GridImp>
  template<int cc>
  alu_inline int ALU3dGridEntity<0,dim,GridImp> :: getSubIndex (int i) const
  {
    assert(item_ != 0);
    typedef typename  ImplTraits::IMPLElementType IMPLElType;
    return IndexWrapper<IMPLElType,GridImp::elementType,cc>::subIndex ( *item_, i);
  }

  template<int dim, class GridImp>
  alu_inline int ALU3dGridEntity<0,dim,GridImp> :: subIndex (int i, unsigned int codim ) const
  {
    typedef ElementTopologyMapping<GridImp::elementType> ElemTopo;

    assert(item_ != 0);
    switch (codim)
    {
    case 0 :
      return this->getIndex();
    case 1 :
      return (ALU3dGridFaceGetter< Comm >::getFace( *item_, i ))->getIndex();
    case 2 :
      return item_->myhedge1( ElemTopo::dune2aluEdge( i ) )->getIndex();
    case 3 :
      return item_->myvertex( ElemTopo::dune2aluVertex( i ) )->getIndex();
    default :
      assert(false);
      abort();
    }
    return -1;
  }

  //******** begin method entity ******************
  template <class GridImp, int dim, int cd> struct SubEntities {};

  // specialisation for elements
  template <class GridImp, int dim>
  struct SubEntities<GridImp, dim, 0>
  {
    typedef typename GridImp::GridObjectFactoryType FactoryType;
    typedef ALU3dGridEntity<0,dim,GridImp> EntityType;

    static typename ALU3dGridEntity<0,dim,GridImp>::template Codim<0>::EntityPointer
    entity (const FactoryType& factory,
            const int level,
            const EntityType & entity,
            const typename ALU3dImplTraits<GridImp::elementType, typename GridImp::MPICommunicatorType>::IMPLElementType & item,
            int i)
    {
      return ALU3dGridEntityPointer<0, GridImp>( entity );
    }
  };

  // specialisation for faces
  template <class GridImp, int dim>
  struct SubEntities<GridImp,dim,1>
  {
    typedef typename GridImp::GridObjectFactoryType FactoryType;
    typedef ElementTopologyMapping<GridImp::elementType> Topo;
    typedef ALU3dGridEntity<0,dim,GridImp> EntityType;

    static typename ALU3dGridEntity<0,dim,GridImp> :: template Codim<1>:: EntityPointer
    entity (const FactoryType& factory,
            const int level,
            const EntityType & en,
            const typename ALU3dImplTraits<GridImp::elementType, typename GridImp::MPICommunicatorType>::IMPLElementType & item,
            int duneFace)
    {
      int aluFace = Topo::dune2aluFace(duneFace);
      return
        ALU3dGridEntityPointer<1,GridImp>
          (factory,
          level,
          *(ALU3dGridFaceGetter< typename GridImp::MPICommunicatorType >::getFace(item, duneFace)),    // getFace already constains dune2aluFace
          item.twist(aluFace),
          duneFace    // we need the duneFace number here for the buildGeom method
          );
    }
  };

  // specialisation for edges
  template <class GridImp, int dim>
  struct SubEntities<GridImp,dim,2>
  {
    typedef typename GridImp::GridObjectFactoryType FactoryType;
    typedef ElementTopologyMapping<GridImp::elementType> Topo;
    typedef ALU3dGridEntity<0,dim,GridImp> EntityType;
    typedef typename GridImp::ctype coordType;

    typedef typename GridImp :: ReferenceElementType ReferenceElementType;

    static typename ALU3dGridEntity<0,dim,GridImp> :: template Codim<2>:: EntityPointer
    entity (const FactoryType& factory,
            const int level,
            const EntityType & en,
            const typename ALU3dImplTraits<GridImp::elementType, typename GridImp::MPICommunicatorType>::IMPLElementType & item,
            int i)
    {
      // get reference element
      const ReferenceElementType & refElem = factory.grid().referenceElement();

      // get first local vertex number of edge i
      int localNum = refElem.subEntity(i,2,0,dim);

      // get number of first vertex on edge
      int v = en.template getSubIndex<dim> (localNum);

      // get the hedge object
      const typename ALU3dImplTraits<GridImp::elementType, typename GridImp::MPICommunicatorType>::GEOEdgeType &
      edge = *(item.myhedge1(Topo::dune2aluEdge(i)));

      int vx = edge.myvertex(0)->getIndex();

      // check whether vertex numbers are equal, otherwise twist is 1
      int twst = (v != vx) ? 1 : 0;
      return ALU3dGridEntityPointer<2,GridImp> (factory, level, edge, twst );
    }
  };

  // specialisation for vertices
  template <class GridImp, int dim>
  struct SubEntities<GridImp,dim,3>
  {
    typedef typename GridImp::GridObjectFactoryType FactoryType;
    typedef ElementTopologyMapping<GridImp::elementType> Topo;
    typedef ALU3dGridEntity<0,dim,GridImp> EntityType;

    static typename ALU3dGridEntity<0,dim,GridImp> :: template Codim<3>:: EntityPointer
    entity (const FactoryType& factory,
            const int level,
            const EntityType & en,
            const typename ALU3dImplTraits<GridImp::elementType, typename GridImp::MPICommunicatorType>::IMPLElementType & item,
            int i)
    {
      return ALU3dGridEntityPointer<3,GridImp>
               (factory, level, *item.myvertex(Topo::dune2aluVertex(i))); // element topo
    }
  };

  template<int dim, class GridImp>
  template<int cc>
  typename ALU3dGridEntity<0,dim,GridImp>::template Codim<cc>:: EntityPointer
  ALU3dGridEntity<0,dim,GridImp> :: subEntity (int i) const
  {
    return SubEntities<GridImp,dim,cc>::entity(factory_,level(),*this,*item_,i);
  }
  //**** end method entity *********

  template<int dim, class GridImp>
  typename ALU3dGridEntity<0,dim,GridImp> :: EntityPointer
  ALU3dGridEntity<0,dim,GridImp> :: father() const
  {
    HElementType* up = item_->up();
    if( ! up )
    {
      std::cerr << "ALU3dGridEntity<0," << dim << "," << dimworld << "> :: father() : no father of entity globalid = " << getIndex() << "\n";
      return ALU3dGridEntityPointer<0,GridImp> (factory_, static_cast<HElementType &> (*item_));
    }

    if( isGhost () )
    {
      return ALU3dGridEntityPointer<0,GridImp> (factory_, static_cast<const HBndSegType &> (*(getGhost().up())));
    }

    return ALU3dGridEntityPointer<0,GridImp> (factory_, static_cast<HElementType &> ( *up ));
  }

  // Adaptation methods
  template<int dim, class GridImp>
  bool ALU3dGridEntity<0,dim,GridImp> :: mark (int ref) const
  {
    assert(item_ != 0);

    // do not allow to mark ghost cells or non-leaf cells
    // this will lead to unpredictable results errors
    if( isGhost() || ! isLeaf() ) return false ;

    // mark for coarsening
    if(ref < 0)
    {
      // don't mark macro elements for coarsening ;)
      if(level() <= 0) return false;

      item_->request(coarse_element_t);
      return true;
    }

    // mark for refinement
    if(ref > 0)
    {
      // for tetrahedral elements check whether to use bisection
      if( GridImp :: elementType == tetra &&
          grid().conformingRefinement() )
      {
        item_->request( bisect_element_t );
      }
      else
      {
        item_->request( refine_element_t );
      }
      return true;
    }

    // mark for none
    item_->request( nosplit_element_t );
    return true;
  }

  // return mark of entity
  template<int dim, class GridImp>
  alu_inline int ALU3dGridEntity<0,dim,GridImp> :: getMark () const
  {
    assert(item_ != 0);

    const MarkRuleType rule = (*item_).requestrule();

    if(rule == refine_element_t) return 1;
    if(rule == coarse_element_t) return -1;
    assert( rule == nosplit_element_t );
    return 0;
  }


  template<int dim, class GridImp>
  bool ALU3dGridEntity<0,dim,GridImp> :: hasBoundaryIntersections () const
  {
    // on ghost elements return false
    if( isGhost() ) return false;

    enum { numFaces = EntityCount<GridImp::elementType>::numFaces };
    typedef typename ImplTraits::HasFaceType HasFaceType;
    typedef typename ImplTraits::GEOFaceType GEOFaceType;

    assert( item_ );
    for(int i=0; i<numFaces; ++i)
    {
      const GEOFaceType &face = *ALU3dGridFaceGetter< Comm >::getFace( *item_, i );

      // don't count internal boundaries as boundary
      if( face.isBorder() ) continue ;

      // check both
      const HasFaceType * outerElement = face.nb.front().first;
      // if we got our own element, get other side
      if( item_ == outerElement )
      {
        outerElement = face.nb.rear().first;
      }

      assert( outerElement );
      if( outerElement->isboundary() ) return true;
    }
    return false;
  }

#if COMPILE_ALUGRID_LIB
  // Instantiation
  template class ALU3dGrid< hexa, No_Comm >;
  template class ALU3dGrid< tetra, No_Comm >;

  // Instantiation
  template class ALU3dGridEntity<0, 3, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridEntity<0, 3, const ALU3dGrid< hexa, No_Comm > >;

  template class ALU3dGridEntity<1, 3, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridEntity<1, 3, const ALU3dGrid< hexa, No_Comm > >;

  template class ALU3dGridEntity<2, 3, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridEntity<2, 3, const ALU3dGrid< hexa, No_Comm > >;

  template class ALU3dGridEntity<3, 3, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridEntity<3, 3, const ALU3dGrid< hexa, No_Comm > >;

  template ALU3dGrid< tetra, No_Comm > :: Traits :: Codim< 0 > :: EntityPointer
  ALU3dGridEntity<0, 3, const ALU3dGrid< tetra, No_Comm > > :: subEntity< 0 >( int ) const;
  template ALU3dGrid< hexa, No_Comm > :: Traits :: Codim< 0 > :: EntityPointer
  ALU3dGridEntity<0, 3, const ALU3dGrid< hexa, No_Comm > > :: subEntity< 0 >( int ) const;

  template ALU3dGrid< tetra, No_Comm > :: Traits :: Codim< 1 > :: EntityPointer
  ALU3dGridEntity<0, 3, const ALU3dGrid< tetra, No_Comm > > :: subEntity< 1 >( int ) const;
  template ALU3dGrid< hexa, No_Comm > :: Traits :: Codim< 1 > :: EntityPointer
  ALU3dGridEntity<0, 3, const ALU3dGrid< hexa, No_Comm > > :: subEntity< 1 >( int ) const;

  template ALU3dGrid< tetra, No_Comm > :: Traits :: Codim< 2 > :: EntityPointer
  ALU3dGridEntity<0, 3, const ALU3dGrid< tetra, No_Comm > > :: subEntity< 2 >( int ) const;
  template ALU3dGrid< hexa, No_Comm > :: Traits :: Codim< 2 > :: EntityPointer
  ALU3dGridEntity<0, 3, const ALU3dGrid< hexa, No_Comm > > :: subEntity< 2 >( int ) const;

  template ALU3dGrid< tetra, No_Comm > :: Traits :: Codim< 3 > :: EntityPointer
  ALU3dGridEntity<0, 3, const ALU3dGrid< tetra, No_Comm > > :: subEntity< 3 >( int ) const;
  template ALU3dGrid< hexa, No_Comm > :: Traits :: Codim< 3 > :: EntityPointer
  ALU3dGridEntity<0, 3, const ALU3dGrid< hexa, No_Comm > > :: subEntity< 3 >( int ) const;

#if ALU3DGRID_PARALLEL
  // Instantiation
  template class ALU3dGrid< hexa, MPI_Comm >;
  template class ALU3dGrid< tetra, MPI_Comm >;

  // Instantiation with MPI
  template class ALU3dGridEntity<0, 3, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridEntity<0, 3, const ALU3dGrid< hexa, MPI_Comm > >;

  template class ALU3dGridEntity<1, 3, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridEntity<1, 3, const ALU3dGrid< hexa, MPI_Comm > >;

  template class ALU3dGridEntity<2, 3, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridEntity<2, 3, const ALU3dGrid< hexa, MPI_Comm > >;

  template class ALU3dGridEntity<3, 3, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridEntity<3, 3, const ALU3dGrid< hexa, MPI_Comm > >;

  template ALU3dGrid< tetra, MPI_Comm > :: Traits :: Codim< 0 > :: EntityPointer
  ALU3dGridEntity<0, 3, const ALU3dGrid< tetra, MPI_Comm > > :: subEntity< 0 >( int ) const;
  template ALU3dGrid< hexa, MPI_Comm > :: Traits :: Codim< 0 > :: EntityPointer
  ALU3dGridEntity<0, 3, const ALU3dGrid< hexa, MPI_Comm > > :: subEntity< 0 >( int ) const;

  template ALU3dGrid< tetra, MPI_Comm > :: Traits :: Codim< 1 > :: EntityPointer
  ALU3dGridEntity<0, 3, const ALU3dGrid< tetra, MPI_Comm > > :: subEntity< 1 >( int ) const;
  template ALU3dGrid< hexa, MPI_Comm > :: Traits :: Codim< 1 > :: EntityPointer
  ALU3dGridEntity<0, 3, const ALU3dGrid< hexa, MPI_Comm > > :: subEntity< 1 >( int ) const;

  template ALU3dGrid< tetra, MPI_Comm > :: Traits :: Codim< 2 > :: EntityPointer
  ALU3dGridEntity<0, 3, const ALU3dGrid< tetra, MPI_Comm > > :: subEntity< 2 >( int ) const;
  template ALU3dGrid< hexa, MPI_Comm > :: Traits :: Codim< 2 > :: EntityPointer
  ALU3dGridEntity<0, 3, const ALU3dGrid< hexa, MPI_Comm > > :: subEntity< 2 >( int ) const;

  template ALU3dGrid< tetra, MPI_Comm > :: Traits :: Codim< 3 > :: EntityPointer
  ALU3dGridEntity<0, 3, const ALU3dGrid< tetra, MPI_Comm > > :: subEntity< 3 >( int ) const;
  template ALU3dGrid< hexa, MPI_Comm > :: Traits :: Codim< 3 > :: EntityPointer
  ALU3dGridEntity<0, 3, const ALU3dGrid< hexa, MPI_Comm > > :: subEntity< 3 >( int ) const;
#endif // #if ALU3DGRID_PARALLEL

#endif // #if COMPILE_ALUGRID_LIB
}
#undef alu_inline
#endif
