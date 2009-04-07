// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU3DGRIDINDEXSETS_HH
#define DUNE_ALU3DGRIDINDEXSETS_HH

//- System includes
#include <vector>

//- Dune includes
#include <dune/common/stdstreams.hh>
#include <dune/common/bigunsignedint.hh>

#include <dune/grid/common/grid.hh>
#include <dune/grid/common/indexidset.hh>
#include <dune/grid/common/dynamicsubindexid.hh>


//- Local includes
#include "alu3dinclude.hh"
#include "topology.hh"
#include "alu3diterators.hh"

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template <int dim, int dimworld, ALU3dGridElementType elType>
  class ALU3dGrid;

  template<int cd, int dim, class GridImp>
  class ALU3dGridEntity;



  // ALU3dGridHierarchicIndexSet
  // ---------------------------

  //! hierarchic index set of ALU3dGrid
  template <int dim, int dimworld, ALU3dGridElementType elType>
  class ALU3dGridHierarchicIndexSet :
    public IndexSet <ALU3dGrid<dim,dimworld,elType>,
        ALU3dGridHierarchicIndexSet<dim,dimworld,elType> >
  {
    typedef ALU3dGridHierarchicIndexSet< dim, dimworld, elType > This;

    typedef ALU3dGrid<dim,dimworld,elType> GridType;
    enum { numCodim = dim+1 }; // i.e. 4

    friend class ALU3dGrid<dim,dimworld,elType>;

    // constructor
    ALU3dGridHierarchicIndexSet( const GridType &grid )
      : grid_( grid ),
        dynamicSubIndex_( *this )
    {}

  public:
    typedef typename GridType::Traits::template Codim<0>::Entity EntityCodim0Type;

    //! return hierarchic index of given entity
    template <class EntityType>
    int index (const EntityType & ep) const
    {
      enum { cd = EntityType :: codimension };
      return index<cd>(ep);
    }

    //! return hierarchic index of given entity
    template< int codim >
    int index ( const typename GridType::Traits::template Codim< codim >::Entity &entity ) const
    {
      return GridType::getRealImplementation( entity ).getIndex();
    }

    //! return subIndex of given entity
    template< int codim >
    int subIndex ( const EntityCodim0Type &e, int i ) const
    {
      typedef ALU3dGridEntity< 0, dim, const GridType > EntityImpl;
      const EntityImpl &entity = GridType::getRealImplementation( e );
      return entity.template getSubIndex< codim >( i );
    }

    int subIndex ( const EntityCodim0Type &e, int i, unsigned int codim ) const
    {
      return dynamicSubIndex_( e, i, codim );
    }

    //! return size of indexset, i.e. maxindex+1
    //! for given type, if type is not exisiting within grid 0 is returned
    int size ( GeometryType type ) const
    {
      if( elType == tetra && !type.isSimplex() ) return 0;
      if( elType == hexa  && !type.isCube() ) return 0;
      // return size of hierarchic index set
      return this->size(GridType::dimension-type.dim());
    }

    //! return size of indexset, i.e. maxindex+1
    int size ( int codim ) const
    {
      // return size of hierarchic index set
      return grid_.hierSetSize(codim);
    }

    //! deliver all geometry types used in this grid
    const std::vector<GeometryType>& geomTypes (int codim) const
    {
      return grid_.geomTypes(codim);
    }

    //! return true because all entities are contained in this set
    template <class EntityType>
    bool contains (const EntityType &) const { return true; }

  private:
    // our Grid
    const GridType & grid_;

    // dynamic caller for subIndex
    const DynamicSubIndex< GridType, This > dynamicSubIndex_;
  };

  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////

  class ALUMacroKey : public ALU3DSPACE Key4<int>
  {
    typedef int A;
    typedef ALUMacroKey ThisType;
    typedef ALU3DSPACE Key4<A> BaseType;

  public:
    ALUMacroKey() : BaseType(-1,-1,-1,-1) {}
    ALUMacroKey(const A&a,const A&b,const A&c,const A&d) : BaseType(a,b,c,d) {}
    ALUMacroKey(const ALUMacroKey & org ) : BaseType(org) {}
    ALUMacroKey & operator = (const ALUMacroKey & org )
    {
      BaseType::operator = (org);
      return *this;
    }

    bool operator == (const ALUMacroKey & org) const
    {
      return ( (this->_a == org._a) &&
               (this->_b == org._b) &&
               (this->_c == org._c) &&
               (this->_d == org._d) );
    }

    // operator < is already implemented in BaseType
    bool operator > (const ALUMacroKey & org) const
    {
      return ( (!this->operator == (org)) && (!this->operator <(org)) );
    }

    void print(std::ostream & out) const
    {
      out << "[" << this->_a << "," << this->_b << "," << this->_c << "," << this->_d << "]";
    }
  };

  template <class MacroKeyImp>
  class ALUGridId
  {
    MacroKeyImp key_;
    int nChild_;
    int codim_;
  public:
    ALUGridId() : key_()
                  , nChild_(-1)
                  , codim_(-1)
    {}

    ALUGridId(const MacroKeyImp & key, int nChild , int cd)
      : key_(key) , nChild_(nChild)
        , codim_(cd)
    {}

    ALUGridId(const ALUGridId & org )
      : key_(org.key_)
        , nChild_(org.nChild_)
        , codim_(org.codim_)
    {}

    ALUGridId & operator = (const ALUGridId & org )
    {
      key_    = org.key_;
      nChild_ = org.nChild_;
      codim_  = org.codim_;
      return *this;
    }

    bool operator == (const ALUGridId & org) const
    {
      return equals(org);
    }

    bool operator != (const ALUGridId & org) const
    {
      return ! equals(org);
    }

    bool operator <= (const ALUGridId & org) const
    {
      if(equals(org)) return true;
      else return lesser(org);
    }

    bool operator >= (const ALUGridId & org) const
    {
      if(equals(org)) return true;
      else return ! lesser(org);
    }

    bool operator < (const ALUGridId & org) const
    {
      return lesser(org);
    }

    bool operator > (const ALUGridId & org) const
    {
      return (!equals(org) && ! lesser(org));
    }

    const MacroKeyImp & getKey() const { return key_; }
    int nChild() const { return nChild_; }
    int codim() const { return codim_; };

    bool isValid () const
    {
      return ( (nChild_ >= 0) && (codim_  >= 0) );
    }

    void reset()
    {
      nChild_ = -1;
      codim_  = -1;
    }

    void print(std::ostream & out) const
    {
      out << "(" << getKey() << "," << nChild_ << "," << codim_ << ")";
    }

  protected:
    // returns true is the id is lesser then org
    bool lesser(const ALUGridId & org) const
    {
      if(getKey() < org.getKey() ) return true;
      if(getKey() > org.getKey() ) return false;
      if(getKey() == org.getKey() )
      {
        if(nChild_ == org.nChild_)
        {
          return codim_ < org.codim_;
        }
        else
          return nChild_ < org.nChild_;
      }
      assert( equals(org) );
      return false;
    }

    // returns true if this id equals org
    bool equals(const ALUGridId & org) const
    {
      return ( (getKey() == org.getKey() ) && (nChild_ == org.nChild_)
               && (codim_ == org.codim_) );
    }
  };

  inline std::ostream& operator<< (std::ostream& s, const ALUMacroKey & key)
  {
    key.print(s);
    return s;
  }

  template <class KeyImp>
  inline std::ostream& operator<< (std::ostream& s, const ALUGridId<KeyImp> & id)
  {
    id.print(s);
    return s;
  }

  //*****************************************************************
  //
  //  --GlobalIdSet
  //
  //*****************************************************************
  //! global id set for ALU3dGrid
  template <int dim, int dimworld, ALU3dGridElementType elType>
  class ALU3dGridGlobalIdSet :
    public IdSetDefaultImplementation  < ALU3dGrid<dim,dimworld,elType> ,
        ALU3dGridGlobalIdSet<dim,dimworld,elType> ,
        typename  ALU3dGrid<dim,dimworld,elType>::Traits::GlobalIdType >
    , public ALU3DSPACE AdaptRestrictProlongType
  {
    typedef ALU3dGrid<dim,dimworld,elType> GridType;
    typedef typename GridType :: HierarchicIndexSet HierarchicIndexSetType;

    typedef ALU3dImplTraits<elType> ImplTraitsType;
    typedef typename ImplTraitsType::IMPLElementType IMPLElementType;
    typedef typename ImplTraitsType::GEOElementType GEOElementType;
    typedef typename ImplTraitsType::GEOFaceType GEOFaceType;
    typedef typename ImplTraitsType::GEOEdgeType GEOEdgeType;

    typedef ALU3DSPACE HElementType HElementType;
    typedef ALU3DSPACE HFaceType HFaceType;
    typedef ALU3DSPACE HEdgeType HEdgeType;
    typedef ALU3DSPACE VertexType VertexType;
    typedef ALU3DSPACE HBndSegType HBndSegType;

    typedef EntityCount<elType> EntityCountType;

  public:
    //! export type of id
    // id is of type ALUGridId <  MacroKeyType >
    typedef typename ALU3dGrid<dim,dimworld,elType>::Traits::GlobalIdType IdType;

  private:
    typedef ALUMacroKey MacroKeyType;

    //typedef std::pair < MacroKeyType , int > MacroIdType;
    typedef ALUGridId <  MacroKeyType > MacroIdType; // same as IdType
    enum { numCodim = dim+1 };

    // this means that only up to 300000000 entities are allowed
    typedef typename GridType::Traits::template Codim<0>::Entity EntityCodim0Type;
  private:
    mutable std::map< int , IdType > ids_[numCodim];
    //mutable std::map< int , MacroKeyType > macroKeys_[numCodim];

    // our Grid
    const GridType & grid_;

    // the hierarchicIndexSet
    const HierarchicIndexSetType & hset_;

    int vertexKey_[4];

    int chunkSize_ ;

    enum { startOffSet_ = 0 };

  public:
    //! create id set, only allowed for ALU3dGrid
    ALU3dGridGlobalIdSet(const GridType & grid)
      : grid_(grid), hset_(grid.hierarchicIndexSet())
        , chunkSize_(100)
    {
      if(elType == hexa)
      {
        // see ALUGrid/src/serial/gitter_mgb.cc
        // InsertUniqueHexa
        const int vxKey[4] = {0,1,3,4};
        for(int i=0; i<4; i++) vertexKey_[i] = vxKey[i];
      }
      else
      {
        assert( elType == tetra );
        // see ALUGrid/src/serial/gitter_mgb.cc
        // InsertUniqueTetra
        const int vxKey[4] = {0,1,2,3};
        for(int i=0; i<4; i++) vertexKey_[i] = vxKey[i];
      }

      // setup the id set
      buildIdSet();
    }

    virtual ~ALU3dGridGlobalIdSet() {}

    // update id set after adaptation
    void updateIdSet()
    {
      // to be revised
      buildIdSet();
    }

    // print all ids
    void print () const
    {
      for(int i=0 ; i<numCodim; ++i)
      {
        std::cout << "*****************************************************\n";
        std::cout << "Ids for codim " << i << "\n";
        std::cout << "*****************************************************\n";
        for(unsigned int k=0; k<ids_[i].size(); ++k)
        {
          std::cout << "Item[" << i << "," << k <<"] has id " << ids_[i][k] << "\n";
        }
        std::cout << "\n\n\n";
      }
    }

    template <class IterType>
    void checkId(const IdType & macroId, const IterType & idIter) const //int codim , unsigned int num ) const
    {

      IdType id = getId(macroId);
      for(int i=0 ; i<numCodim; ++i)
      {
        typedef typename std::map<int,IdType>::iterator IteratorType;
        IteratorType end = ids_[i].end();
        for(IteratorType it = ids_[i].begin(); it != end; ++it)
        //for(unsigned int k=0; k<ids_[i].size(); ++k)
        {
          if(idIter == it) continue;
          //if((i == codim) && (k == num)) continue;
          const IdType & checkMId = (*it).second; //ids_[i][k];
          IdType checkId = getId(checkMId);
          if( id == checkId )
          {
            //std::cout << "Check(codim,num = " << codim<< "," << num <<") failed for k="<<k << " codim = " << i << "\n";
            std::cout << id << " equals " << checkId << "\n";
            assert( id != checkId );
            DUNE_THROW(GridError," " << id << " equals " << checkId << "\n");
          }
          else
          {
            bool lesser  = (id < checkId);
            bool greater = (id > checkId);
            assert( lesser != greater );
            if( lesser == greater )
            {
              assert( lesser != greater );
              DUNE_THROW(GridError," lesser equals greater of one id ");
            }
          }
        }
      }
    }

    // check id set for uniqueness
    void uniquenessCheck() const
    {
      for(int i=0 ; i<numCodim; i++)
      {
        typedef typename std::map<int,IdType>::iterator IteratorType;
        IteratorType end = ids_[i].end();
        for(IteratorType it = ids_[i].begin(); it != end; ++it)
        //unsigned int k=0; k<ids_[i].size(); ++k)
        {
          const IdType & id = (*it).second; //ids_[i][k];
          if( id.isValid() )
            checkId(id,it); //i,k);
        }
      }
    }

    void setChunkSize( int chunkSize )
    {
      chunkSize_ = chunkSize;
    }

    // creates the id set
    void buildIdSet ()
    {
      for(int i=0; i<numCodim; ++i)
      {
        ids_[i].clear();
      }

      ALU3DSPACE GitterImplType & gitter = const_cast<ALU3DSPACE
                                                      GitterImplType &> (grid_.myGrid());

      // all interior and border vertices
      {
        ALU3DSPACE AccessIterator <VertexType>::Handle fw (gitter.container ()) ;
        for (fw.first () ; ! fw.done () ; fw.next ())
        {
          int idx = fw.item().getIndex();
          ids_[3][idx] = buildMacroVertexId( fw.item() );
        }
      }

      // all ghost vertices
      {
        typedef typename ALU3DSPACE ALU3dGridLevelIteratorWrapper<3,Ghost_Partition> IteratorType;
        IteratorType fw (grid_ , 0 , grid_.nlinks() );
        typedef typename IteratorType :: val_t val_t;
        for (fw.first () ; ! fw.done () ; fw.next ())
        {
          val_t & item = fw.item();
          assert( item.first );
          VertexType & vx = * (item.first);
          int idx = vx.getIndex();
          ids_[3][idx] = buildMacroVertexId( vx );
        }
      }

      // create ids for all macro edges
      {
        ALU3DSPACE AccessIterator <HEdgeType> :: Handle w (gitter.container ()) ;
        for (w.first(); !w.done(); w.next())
        {
          int idx = w.item().getIndex();
          ids_[2][idx] = buildMacroEdgeId( w.item() );
          buildEdgeIds( w.item() , ids_[2][idx] , startOffSet_ );
        }
      }

      // all ghost edges
      {
        typedef typename ALU3DSPACE ALU3dGridLevelIteratorWrapper<2,Ghost_Partition> IteratorType;
        IteratorType fw (grid_ , 0 , grid_.nlinks() );
        typedef typename IteratorType :: val_t val_t;
        for (fw.first () ; ! fw.done () ; fw.next ())
        {
          val_t & item = fw.item();
          assert( item.first );
          HEdgeType & edge = * (item.first);
          int idx = edge.getIndex();

          ids_[2][idx] = buildMacroEdgeId( edge );
          buildEdgeIds( edge , ids_[2][idx] , startOffSet_ );
        }
      }


      // for all macro faces and all children
      {
        ALU3DSPACE AccessIterator <HFaceType>::Handle w (gitter.container ()) ;
        for (w.first () ; ! w.done () ; w.next ())
        {
          int idx = w.item().getIndex();
          ids_[1][idx] = buildMacroFaceId( w.item() );
          buildFaceIds( w.item() , ids_[1][idx] , startOffSet_ );
        }
      }

      // all ghost faces
      {
        typedef typename ALU3DSPACE ALU3dGridLevelIteratorWrapper<1,Ghost_Partition> IteratorType;
        IteratorType fw (grid_ , 0 , grid_.nlinks() );
        typedef typename IteratorType :: val_t val_t;
        for (fw.first () ; ! fw.done () ; fw.next ())
        {
          val_t & item = fw.item();
          assert( item.first );
          HFaceType & face = * (item.first);
          int idx = face.getIndex();
          ids_[1][idx] = buildMacroFaceId( face );
          buildFaceIds( face , ids_[1][idx] , startOffSet_ );
        }
      }

      // for all macro elements and all internal entities
      {
        ALU3DSPACE AccessIterator <HElementType> :: Handle w (gitter.container ()) ;
        for (w.first () ; ! w.done () ; w.next ())
        {
          int idx = w.item().getIndex();
          ids_[0][idx] = buildMacroElementId( w.item() );
          buildElementIds( w.item() , ids_[0][idx] , startOffSet_ );
        }
      }

      // all ghost elements
      {
        typedef typename ALU3DSPACE ALU3dGridLevelIteratorWrapper<0,Ghost_Partition> IteratorType;
        IteratorType fw (grid_ , 0 , grid_.nlinks() );
        typedef typename IteratorType :: val_t val_t;
        for (fw.first () ; ! fw.done () ; fw.next ())
        {
          val_t & item = fw.item();
          assert( item.second );
          HElementType & elem = * ( item.second->getGhost().first );
          int idx = elem.getIndex();
          ids_[0][idx] = buildMacroElementId( elem );
          buildElementIds( elem , ids_[0][idx] , startOffSet_ );
        }
      }

      // check uniqueness of id only in serial, because
      // in parallel some faces and edges of ghost exists more than once
      // but have the same id, but not the same index, there for the check
      // will fail for ghost elements
#if ! ALU3DGRID_PARALLEL
      // be carefull with this check, it's complexity is O(N^2)
      //uniquenessCheck();
#endif
    }

    IdType buildMacroVertexId(const VertexType & item )
    {
      int vx[4] = { item.ident(), -1, -1, -1};
      enum {codim = 3 };
      MacroKeyType key(vx[0],vx[1],vx[2],vx[3]);
      MacroIdType id(key,1, codim + startOffSet_ );
      return id;
    }

    IdType buildMacroEdgeId(const HEdgeType & item )
    {
      const GEOEdgeType & edge = static_cast<const GEOEdgeType &> (item);
      int vx[4] = {-1,-1,-1,-1};
      for(int i=0; i<2; ++i)
      {
        vx[i] = edge.myvertex(i)->ident();
      }

      enum { codim = 2 };
      MacroKeyType key(vx[0],vx[1],vx[2],vx[3]);
      MacroIdType id( key,1,  codim + startOffSet_ );
      return id;
    }

    IdType buildMacroFaceId(const HFaceType & item )
    {
      const GEOFaceType & face = static_cast<const GEOFaceType &> (item);
      int vx[4] = {-1,-1,-1,-1};
      for(int i=0; i<3; ++i)
      {
        vx[i] = face.myvertex(i)->ident();
      }

      enum { codim = 1 };
      MacroKeyType key(vx[0],vx[1],vx[2],vx[3]);
      MacroIdType id(key,1,  codim + startOffSet_ );
      return id;
    }

    IdType buildMacroElementId(const HElementType & item )
    {
      const GEOElementType & elem = static_cast<const GEOElementType &> (item);
      int vx[4] = {-1,-1,-1,-1};
      for(int i=0; i<4; ++i)
      {
        vx[i] = elem.myvertex(vertexKey_[i])->ident();
      }
      enum { codim = 0 };
      MacroKeyType key(vx[0],vx[1],vx[2],vx[3]);
      return MacroIdType(key,1,  codim + startOffSet_ );
    }

    template <int cd>
    IdType createId(const typename ImplTraitsType::
                    template Codim<cd>::InterfaceType & item , const IdType & creatorId , int nChild )
    {
      assert( creatorId.isValid() );

      // we have up to 12 internal hexa faces, therefore need 100 offset
      enum { childOffSet = ((cd == 1) && (elType == hexa)) ? 16 : 8 };
      enum { codimOffSet = 4 };

      assert( nChild < childOffSet );

      int newChild = (creatorId.nChild() * childOffSet ) + nChild;
      int newCodim = (creatorId.codim()  * codimOffSet ) + ( cd + startOffSet_ );

      IdType newId( creatorId.getKey() , newChild , newCodim );
      assert( newId != creatorId );
      return newId;
    }

    // build ids for all children of this element
    void buildElementIds(const HElementType & item , const IdType & macroId , int nChild)
    {
      enum { codim = 0 };
      ids_[codim][item.getIndex()] = createId<codim>(item,macroId,nChild);

      const IdType & itemId = ids_[codim][item.getIndex()];

      buildInteriorElementIds(item,itemId);
    }

    // build ids for all children of this element
    void buildInteriorElementIds(const HElementType & item , const IdType & fatherId)
    {
      assert( fatherId.isValid() );

      // build id for inner vertex
      {
        const VertexType * v = item.innerVertex() ;
        // for tetras there is no inner vertex, therefore check
        if(v) buildVertexIds(*v,fatherId );
      }

      // build edge ids for all inner edges
      {
        int inneredge = startOffSet_;
        for(const HEdgeType * e = item.innerHedge () ; e ; e = e->next ())
        {
          buildEdgeIds(*e,fatherId,inneredge);
          ++inneredge;
        }
      }

      // build face ids for all inner faces
      {
        int innerface = startOffSet_;
        for(const HFaceType * f = item.innerHface () ; f ; f = f->next ())
        {
          buildFaceIds(*f,fatherId,innerface);
          ++innerface;
        }
      }

      // build ids of all children
      {
        int numChild = startOffSet_;
        for(const HElementType * child = item.down(); child; child =child->next() )
        {
          //assert( numChild == child->nChild() );
          buildElementIds(*child, fatherId, numChild);
          ++numChild;
        }
      }
    }

    // build ids for all children of this face
    void buildFaceIds(const HFaceType & face, const IdType & fatherId , int innerFace )
    {
      enum { codim = 1 };
      ids_[codim][face.getIndex()] = createId<codim>(face,fatherId,innerFace);
      const IdType & faceId = ids_[codim][face.getIndex()];

      buildInteriorFaceIds(face,faceId);
    }

    // build ids for all children of this face
    void buildInteriorFaceIds(const HFaceType & face, const IdType & faceId)
    {
      assert( faceId.isValid () );

      // build id for inner vertex
      {
        const VertexType * v = face.innerVertex() ;
        //std::cout << "create inner vertex of face " << face.getIndex() << "\n";
        if(v) buildVertexIds(*v,faceId );
      }

      // build ids for all inner edges
      {
        int inneredge = startOffSet_;
        for (const HEdgeType * e = face.innerHedge () ; e ; e = e->next ())
        {
          buildEdgeIds(*e,faceId ,inneredge );
          ++inneredge;
        }
      }

      // build ids for all child faces
      {
        int child = startOffSet_;
        for(const HFaceType * f = face.down () ; f ; f = f->next ())
        {
          assert( child == f->nChild()+startOffSet_);
          buildFaceIds(*f,faceId,child);
          ++child;
        }
      }
    }

    // build ids for all children of this edge
    void buildEdgeIds(const HEdgeType & edge, const IdType & fatherId , int inneredge)
    {
      enum { codim = 2 };
      ids_[codim][edge.getIndex()] = createId<codim>(edge,fatherId,inneredge);
      const IdType & edgeId = ids_[codim][edge.getIndex()];
      buildInteriorEdgeIds(edge,edgeId);
    }

    void buildInteriorEdgeIds(const HEdgeType & edge, const IdType & edgeId)
    {
      assert( edgeId.isValid() );

      // build id for inner vertex
      {
        const VertexType * v = edge.innerVertex() ;
        if(v) buildVertexIds(*v,edgeId );
      }

      // build ids for all inner edges
      {
        int child = startOffSet_;
        for (const HEdgeType * e = edge.down () ; e ; e = e->next ())
        {
          assert( child == e->nChild()+ startOffSet_ );
          buildEdgeIds(*e,edgeId , child );
          ++child;
        }
      }
    }

    // build id for this vertex
    void buildVertexIds(const VertexType & vertex, const IdType & fatherId )
    {
      enum { codim = 3 };
      // inner vertex number is 1
      ids_[codim][vertex.getIndex()] = createId<codim>(vertex,fatherId,1);
      assert( ids_[codim][vertex.getIndex()].isValid() );
    }

    friend class ALU3dGrid<dim,dimworld,elType>;

    const IdType & getId(const IdType & macroId) const
    {
      return macroId;
    }
  public:
    //! return global id of given entity
    template <class EntityType>
    IdType id (const EntityType & ep) const
    {
      enum { cd = EntityType :: codimension };
      assert( ids_[cd].find( hset_.index(ep) ) != ids_[cd].end() );
      const IdType & macroId = ids_[cd][hset_.index(ep)];
      assert( macroId.isValid() );
      return getId(macroId);
    }

    //! return global id of given entity
    template <int codim>
    IdType id (const typename GridType:: template Codim<codim> :: Entity & ep) const
    {
      assert( ids_[codim].find( hset_.index(ep) ) != ids_[codim].end() );
      const IdType & macroId = ids_[codim][hset_.index(ep)];
      assert( macroId.isValid() );
      return getId(macroId);
    }

    //! return subId of given entity
    template< int codim >
    IdType subId ( const EntityCodim0Type &e, int i ) const
    {
      const int hIndex = hset_.template subIndex< codim >( e, i );
      assert( ids_[ codim ].find( hIndex ) != ids_[ codim ].end() );
      const IdType &macroId = ids_[ codim ][ hIndex ];
      assert( macroId.isValid() );
      return getId( macroId );
    }

    //! return subId of given entity
    IdType subId ( const EntityCodim0Type &e, int i, unsigned int codim ) const
    {
      const int hIndex = hset_.subIndex( e, i, codim );
      assert( ids_[ codim ].find( hIndex ) != ids_[ codim ].end() );
      const IdType &macroId = ids_[ codim ][ hIndex ];
      assert( macroId.isValid() );
      return getId( macroId );
    }

    template <int d, ALU3dGridElementType element_t >
    struct BuildIds;

    template <int d>
    struct BuildIds<d,tetra>
    {
      //static const IdType zero;
      template <class MyIdSet, class IdStorageType>
      static void buildFace(MyIdSet & set, const HElementType & item, int faceNum,
                            IdStorageType & ids )
      {
        const IMPLElementType & elem = static_cast<const IMPLElementType &> (item);
        const HFaceType & face  = *(elem.myhface3(faceNum));
        const IdType & id = ids[face.getIndex()];
        assert( id.isValid() );
        set.buildInteriorFaceIds(face,id);
      }
    };

    template <int d>
    struct BuildIds<d,hexa>
    {
      //static const IdType zero;
      template <class MyIdSet, class IdStorageType>
      static void buildFace(MyIdSet & set, const HElementType & item, int faceNum,
                            IdStorageType & ids )
      {
        const IMPLElementType & elem = static_cast<const IMPLElementType &> (item);
        const HFaceType & face  = *(elem.myhface4(faceNum));
        const IdType & id = ids[face.getIndex()];
        assert( id.isValid() );
        set.buildInteriorFaceIds(face,id);
      }
    };

    // create ids for refined elements
    int postRefinement( HElementType & item )
    {
      {
        enum { elCodim = 0 };
        const IdType & fatherId = ids_[elCodim][item.getIndex()];
        assert( fatherId.isValid() );
        buildInteriorElementIds(item, fatherId );
      }

      for(int i=0; i<EntityCountType::numFaces; ++i)
      {
        enum { faceCodim = 1 };
        BuildIds<dim,elType>::buildFace(*this,item,i,ids_[faceCodim]);
      }

      for(int i=0; i<EntityCountType::numEdges; ++i)
      {
        enum { edgeCodim = 2 };
        const IMPLElementType & elem = static_cast<const IMPLElementType &> (item);
        const HEdgeType & edge  = *( elem.myhedge1(i));
        const IdType & id = ids_[edgeCodim][edge.getIndex()];
        assert( id.isValid() );
        buildInteriorEdgeIds(edge,id);
      }
      return 0;
    }

    // dummy functions
    int preCoarsening( HElementType & elem )
    {
      /*
         const IdType & fatherId = ids_[0][item.getIndex()];

         removeElementIds(item,fatherId,item.nChild());

         for(int i=0; i<EntityCountType::numFaces; ++i)
         BuildIds<dim,elType>::buildFace(*this,item,i,ids_[1]);

         for(int i=0; i<EntityCountType::numEdges; ++i)
         {
         const IMPLElementType & elem = static_cast<const IMPLElementType &> (item);
         const HEdgeType & edge  = *( elem.myhedge1(i));
         const HEdgeType * child = edge.down();
         assert( child );
         if( ids_[2][child->getIndex() ] > zero_ ) continue;
         buildEdgeIds(edge,ids_[2][edge.getIndex()],0);
         }
         #ifndef NDEBUG
         //uniquenessCheck();
         #endif
       */
      return 0;
    }

    // dummy functions
    int preCoarsening ( HBndSegType & el ) { return 0; }

    //! prolong data, elem is the father
    int postRefinement ( HBndSegType & el ) { return 0; }

  };

  //***********************************************************
  //
  //  --LocalIdSet
  //
  //***********************************************************

  //! hierarchic index set of ALU3dGrid
  template <int dim, int dimworld, ALU3dGridElementType elType>
  class ALU3dGridLocalIdSet :
    public IdSetDefaultImplementation < ALU3dGrid<dim,dimworld,elType> ,
        ALU3dGridLocalIdSet<dim,dimworld,elType> ,
        int >
    , public ALU3DSPACE AdaptRestrictProlongType
  {
    typedef ALU3dGridLocalIdSet< dim, dimworld, elType > This;

    typedef ALU3DSPACE HElementType HElementType;
    typedef ALU3DSPACE HBndSegType HBndSegType;

    typedef ALU3dGrid<dim,dimworld,elType> GridType;
    typedef typename GridType :: HierarchicIndexSet HierarchicIndexSetType;

    // this means that only up to 300000000 entities are allowed
    enum { codimMultiplier = 300000000 };
    typedef typename GridType::Traits::template Codim<0>::Entity EntityCodim0Type;

    // create local id set , only for the grid allowed
    ALU3dGridLocalIdSet(const GridType & grid) : hset_(grid.hierarchicIndexSet())
    {
      for(int i=0; i<dim+1; i++)
        codimStart_[i] = i*codimMultiplier;
    }
    friend class ALU3dGrid<dim,dimworld,elType>;

    // fake method to have the same method like GlobalIdSet
    void updateIdSet() {}

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
    int id (const typename GridType:: template Codim<codim> :: Entity & ep) const
    {
      //enum { cd = EntityType :: codimension };
      assert( hset_.size(codim) < codimMultiplier );
      return codimStart_[codim] + hset_.index(ep);
    }

    //! return subId of given entity
    template< int codim >
    int subId ( const EntityCodim0Type &e, int i ) const
    {
      assert( hset_.size( codim ) < codimMultiplier );
      return codimStart_[ codim ] + hset_.template subIndex<codim >( e, i );
    }

    //! return subId of given entity
    IdType subId ( const EntityCodim0Type &e, int i, unsigned int codim ) const
    {
      assert( hset_.size( codim ) < codimMultiplier );
      return codimStart_[ codim ] + hset_.subIndex( e, i, codim );
    }

    // dummy functions
    int preCoarsening( HElementType & elem )  { return 0; }
    // create ids for refined elements
    int postRefinement( HElementType & item )  { return 0; }

    // dummy functions
    int preCoarsening ( HBndSegType & el ) { return 0; }

    //! prolong data, elem is the father
    int postRefinement ( HBndSegType & el ) { return 0; }

    void setChunkSize( int chunkSize ) {}

  private:
    // our HierarchicIndexSet
    const HierarchicIndexSetType & hset_;

    // store start of each codim numbers
    int codimStart_[dim+1];
  };

} // end namespace Dune
#endif
