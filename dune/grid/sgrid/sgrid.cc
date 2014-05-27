// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_SGRID_CC
#define DUNE_SGRID_CC

#include <algorithm>
#include <iostream>
#include <assert.h>

#include <dune/common/stdstreams.hh>
#include <dune/common/typetraits.hh>

#include <dune/grid/sgrid/generic2dune.hh>

namespace Dune {


  //************************************************************************
  // inline methods for SEntityBase

  template<int codim, int dim, class GridImp, template<int,int,class> class EntityImp>
  void SEntityBase<codim,dim,GridImp,EntityImp>::make (GridImp* _grid, int _l, int _id)
  {
    grid = _grid;
    l = _l;
    index = _id;
    z = grid->z(_l,_id,codim);
    builtgeometry = false;
  }

  template<int codim, int dim, class GridImp, template<int,int,class> class EntityImp>
  void SEntityBase<codim,dim,GridImp,EntityImp>::make (int _l, int _id)
  {
    l = _l;
    index = _id;
    z = grid->z(_l,_id,codim);
    builtgeometry = false;
  }

  template<int codim, int dim, class GridImp, template<int,int,class> class EntityImp>
  inline int SEntityBase<codim,dim,GridImp,EntityImp>::globalIndex () const
  {
    int ind = 0;
    for(int i=0; i<l; i++)
      ind += grid->size(i,codim);
    return ind+compressedIndex();
  }

  template<int codim, int dim, class GridImp, template<int,int,class> class EntityImp>
  void SEntityBase<codim,dim,GridImp,EntityImp>::makegeometry () const
  {
    // find dim-codim direction vectors and reference point
    FieldMatrix<ctype,dim-codim,dimworld> A(0);

    // count number of direction vectors found
    int dir=0;
    FieldVector<ctype, dimworld> p1,p2;
    array<int,dim> t=z;

    // check all directions
    for (int i=0; i<dim; i++)
      if (t[i]%2==1)
      {
        // coordinate i is odd => gives one direction vector
        t[i] += 1;                 // direction i => even
        p2 = grid->pos(l,t);
        t[i] -= 2;                 // direction i => even
        p1 = grid->pos(l,t);
        t[i] += 1;                 // revert t to original state
        A[dir] = p2-p1;
        dir++;
      }

    // find reference point, subtract 1 from all odd directions
    for (int i=0; i<dim; i++)
      t[i] -= t[i]%2;

    // make element
    geo.make(grid->pos(l,t),A);
    builtgeometry = true;
  }

  //************************************************************************
  // inline methods for SEntity

  // singleton holding mapper of unit cube
  template<int dim>
  struct SUnitCubeMapper {
    static CubeMapper<dim> mapper;      // one cube per direction
  };

  // initialize static variable with default constructor (which makes reference elements)
  template<int dim>
  CubeMapper<dim> SUnitCubeMapper<dim>::mapper;


  // codim 0
  template<int dim, class GridImp> template<int cc>
  inline int SEntity<0,dim,GridImp>::count () const
  {
    return SUnitCubeMapper<dim>::mapper.elements(cc);
  }

  template<int dim, class GridImp>
  inline unsigned int SEntity<0,dim,GridImp>::count (unsigned int codim) const
  {
    return SUnitCubeMapper<dim>::mapper.elements(codim);
  }

  // subentity construction
  template<int dim, class GridImp> template<int cc>
  inline typename SEntity<0,dim,GridImp>::template Codim<cc>::EntityPointer SEntity<0,dim,GridImp>::subEntity (int i) const
  {
    // make Iterator
    return SLevelIterator<cc,All_Partition,const GridImp>(this->grid,this->l,this->subCompressedIndex(cc,i));
  }

  template<int dim, class GridImp>
  inline typename SEntity<0,dim,GridImp>::IntersectionIterator SEntity<0,dim,GridImp>::ibegin () const
  {
    return IntersectionIterator(SIntersectionIterator<GridImp>(this->grid,this,0));
  }
  template<int dim, class GridImp>
  inline typename SEntity<0,dim,GridImp>::IntersectionIterator SEntity<0,dim,GridImp>::ileafbegin () const
  {
    // only obtain leaf intersections on maxLevel
    if (isLeaf())
      return ibegin();
    else
      return iend();
  }

  template<int dim, class GridImp>
  inline typename SEntity<0,dim,GridImp>::IntersectionIterator SEntity<0,dim,GridImp>::ilevelbegin () const
  {
    return ibegin();
  }

  template<int dim, class GridImp>
  inline typename SEntity<0,dim,GridImp>::IntersectionIterator SEntity<0,dim,GridImp>::iend () const
  {
    return IntersectionIterator(SIntersectionIterator<GridImp>(this->grid,this,count<1>()));
  }
  template<int dim, class GridImp>
  inline typename SEntity<0,dim,GridImp>::IntersectionIterator SEntity<0,dim,GridImp>::ileafend () const
  {
    return iend();
  }

  template<int dim, class GridImp>
  inline typename SEntity<0,dim,GridImp>::IntersectionIterator SEntity<0,dim,GridImp>::ilevelend () const
  {
    return iend();
  }

  template<int dim, class GridImp>
  void SEntity<0,dim,GridImp>::make_father () const
  {
    // check level
    if (this->l<=0)
    {
      father_index = 0;
      built_father = true;
      return;
    }

    // reduced coordinates from expanded coordinates
    array<int,dim> zz = this->grid->compress(this->l,this->z);

    // look for odd coordinates
    FieldVector<ctype, dim> delta;
    for (int i=0; i<dim; i++)
    {
      delta[i] = zz[i] % 2;
      zz[i] = zz[i] / 2;
    }

    // zz is now the reduced coordinate of the father, compute index
    int partition = this->grid->partition(this->l,this->z);
    father_index = this->grid->n((this->l)-1,this->grid->expand((this->l)-1,zz,partition));

    // now make a subcube of size 1/2 in each direction
    FieldMatrix<ctype,dim,dim> A;
    FieldVector<ctype, dim> v;
    for (int i=0; i<dim; i++)
    {
      v = 0.0; v[i] = 0.5;
      A[i] = v;
    }
    for (int i=0; i<dim; i++) v[i] = 0.5*delta[i];

    in_father_local.make(v,A);     // build geometry

    built_father = true;
  }

  template<int dim, class GridImp>
  inline typename SEntity<0,dim,GridImp>::EntityPointer SEntity<0,dim,GridImp>::father () const
  {
    if (!built_father) make_father();
    if (this->l>0)
      return SLevelIterator<0,All_Partition,const GridImp>((this->grid),(this->l)-1,father_index);
    else
      return SLevelIterator<0,All_Partition,const GridImp>((this->grid),this->l,index);
  }

  template<int dim, class GridImp>
  inline typename GridImp::template Codim<0>::LocalGeometry
  SEntity<0,dim,GridImp>::geometryInFather () const
  {
    if (!built_father) make_father();
    return LocalGeometry( in_father_local );
  }

  template<int dim, class GridImp>
  inline typename SEntity<0,dim,GridImp>::HierarchicIterator SEntity<0,dim,GridImp>::hbegin (int maxLevel) const
  {
    return HierarchicIterator(SHierarchicIterator<GridImp>(this->grid,*this,maxLevel,false));
  }

  template<int dim, class GridImp>
  inline typename SEntity<0,dim,GridImp>::HierarchicIterator SEntity<0,dim,GridImp>::hend (int maxLevel) const
  {
    return HierarchicIterator(SHierarchicIterator<GridImp>(this->grid,*this,maxLevel,true));
  }

  //************************************************************************
  // inline methods for HierarchicIterator

  template<class GridImp>
  inline void SHierarchicIterator<GridImp>::push_sons (int level, int fatherid)
  {
    // check level
    if (level+1>maxLevel) return;     // nothing to do

    // compute reduced coordinates of element
    array<int,dim> z =
      grid->z(level,fatherid,0);      // expanded coordinates from index
    array<int,dim> zred =
      grid->compress(level,z);     // reduced coordinates from expaned coordinates

    // refine to first son
    for (int i=0; i<dim; i++) zred[i] = 2*zred[i];

    // generate all \f$2^{dim}\f$ sons
    int partition = grid->partition(level,z);
    for (int b=0; b<(1<<dim); b++)
    {
      array<int,dim> zz = zred;
      for (int i=0; i<dim; i++)
        if (b&(1<<i)) zz[i] += 1;
      // zz is reduced coordinate of a son on level level+1
      int sonid = grid->n(level+1,grid->expand(level+1,zz,partition));

      // push son on stack
      SHierarchicStackElem son(level+1,sonid);
      //stack.push(StackElem(level+1,sonid));
      stack.push(son);
    }
  }

  template<class GridImp>
  inline void SHierarchicIterator<GridImp>::increment ()
  {
    // check empty stack
    if (stack.empty()) return;

    // OK, lets pop
    SHierarchicStackElem newe = stack.top();
    stack.pop();
    l = newe.l;
    index = newe.index;
    realEntity().make(l,index);     // here is our new element

    // push all sons of this element if it is not the original element
    if (newe.l!=orig_l || newe.index!=orig_index)
      push_sons(newe.l,newe.index);
  }

  //************************************************************************
  // inline methods for IntersectionIterator

  template<class GridImp>
  void SIntersectionIterator<GridImp>::make (int _count) const
  {
    // reset cache flags
    built_intersections = false;
    valid_nb = false;
    valid_count = false;

    // start with given neighbor
    count = _count;

    // check if count is valid
    if (count<0 || count>=grid->getRealImplementation(self).entity().template count<1>())
    {
      grid->getRealImplementation(ne).index = -1;
      return;   // done, this is end iterator
    }
    valid_count = true;

    // and compute compressed coordinates of neighbor
    array<int,dim> zrednb = zred;
    zrednb[count/2] += -1+2*(count%2); // (count%2) ? +1 : -1

    // now check if neighbor exists
    is_on_boundary = !grid->exists(grid->getRealImplementation(self).l,zrednb);
    if (is_on_boundary)
    {
      grid->getRealImplementation(ne).index = -1;
      return;   // ok, done it
    }

    // now neighbor is in the grid and must be initialized.
    // First compute its index
    grid->getRealImplementation(ne).index =
      grid->n(grid->getRealImplementation(self).l,
              grid->expand(grid->getRealImplementation(self).l,zrednb,partition));
    grid->getRealImplementation(ne).realEntity().make(
      grid->getRealImplementation(ne).l,
      grid->getRealImplementation(ne).index);
  }

  template<class GridImp>
  inline bool SIntersectionIterator<GridImp>::equals (const SIntersectionIterator<GridImp>& i) const
  {
    return (self == i.self) && (count==i.count);
  }

  template<class GridImp>
  inline typename SIntersectionIterator<GridImp>::EntityPointer SIntersectionIterator<GridImp>::inside () const
  {
    return self;
  }

  template<class GridImp>
  inline typename SIntersectionIterator<GridImp>::EntityPointer SIntersectionIterator<GridImp>::outside () const
  {
    return ne;
  }

  template<class GridImp>
  inline void SIntersectionIterator<GridImp>::increment ()
  {
    count++;
    make(count);
  }

  template<class GridImp>
  inline bool SIntersectionIterator<GridImp>::boundary () const
  {
    return is_on_boundary;
  }

  template<class GridImp>
  inline bool SIntersectionIterator<GridImp>::neighbor () const
  {
    return (!is_on_boundary);
  }

  template<class GridImp>
  void SIntersectionIterator<GridImp>::makeintersections () const
  {
    // compute direction and value in direction
    int dir = count/2;
    int c = count%2;

    // compute expanded coordinates of entity
    array<int,dim> z1 =
      grid->getRealImplementation(grid->getRealImplementation(self).entity()).z;
    if (c==1)
      z1[dir] += 1;   // odd
    else
      z1[dir] -= 1;   // even

    // z1 is even in direction dir, all others must be odd because it is codim 1
    FieldMatrix<ctype,dim-1,dim> __AsLocal;
    FieldVector<ctype, dim> p1Local,p2Local;

    int t;

    // local coordinates in self
    p1Local = 0.0;
    p1Local[dir] = c;    // all points have p[dir]=c in entity
    t = 0;
    for (int i=0; i<dim; ++i) // this loop makes dim-1 direction vectors
      if (i!=dir)
      {
        // each i!=dir gives one direction vector
        p2Local = p1Local;
        p2Local[i] = 1.0;
        __AsLocal[t] = p2Local-p1Local;     // a direction vector
        ++t;
      }
    // update geometry
    is_self_local.make(p1Local,__AsLocal);

    // local coordinates in neighbor
    p1Local = 0.0;
    p1Local[dir] = 1-c;    // all points have p[dir]=1-c in entity
    t = 0;
    for (int i=0; i<dim; ++i) // this loop makes dim-1 direction vectors
      if (i!=dir)
      {
        // each i!=dir gives one direction vector
        p2Local = p1Local;
        p2Local[i] = 1.0;
        __AsLocal[t] = p2Local-p1Local;     // a direction vector
        ++t;
      }
    // update geometry
    is_nb_local.make(p1Local,__AsLocal);

    // global coordinates
    FieldMatrix<ctype,dim-1,dimworld> A;
    FieldVector<ctype, dimworld> p1,p2;
    t = 0;
    for (int i=0; i<dim; i++)
      if (i!=dir)
      {
        // each i!=dir gives one direction vector
        z1[i] += 1;     // direction i => even
        p2 = grid->pos(self->level(),z1);
        z1[i] -= 2;     // direction i => even
        p1 = grid->pos(self->level(),z1);
        z1[i] += 1;     // revert t to original state
        A[t] = p2-p1;
        ++t;
      }
    for (int i=0; i<dim; i++)
      if (i!=dir)
        z1[i] -= 1;

    // update geometry
    is_global.make(grid->pos(self->level(),z1), A);

    built_intersections = true;
  }

  template<class GridImp>
  inline typename SIntersectionIterator< GridImp >::LocalGeometry
  SIntersectionIterator< GridImp >::geometryInInside () const
  {
    assert (valid_count);
    if (!built_intersections) makeintersections();
    return LocalGeometry( is_self_local );
  }

  template<class GridImp>
  inline typename SIntersectionIterator< GridImp >::LocalGeometry
  SIntersectionIterator< GridImp >::geometryInOutside () const
  {
    assert (valid_count);
    if (!built_intersections) makeintersections();
    return LocalGeometry( is_nb_local );
  }

  template<class GridImp>
  inline typename SIntersectionIterator< GridImp >::Geometry
  SIntersectionIterator< GridImp >::geometry () const
  {
    assert (valid_count);
    if (!built_intersections) makeintersections();
    return Geometry( is_global );
  }

  template<class GridImp>
  inline int SIntersectionIterator<GridImp>::indexInInside () const
  {
    return count;
  }

  template<class GridImp>
  inline int SIntersectionIterator<GridImp>::indexInOutside () const
  {
    return (count/2)*2 + (1-count%2);
  }

  //************************************************************************
  // inline methods for SLevelIterator

  template<int codim, PartitionIteratorType pitype, class GridImp>
  inline void SLevelIterator<codim,pitype,GridImp>::increment ()
  {
    ++index;
    realEntity().make(l,index);
  }

  //************************************************************************
  // inline methods for SEntityPointer

  template<int codim, class GridImp>
  inline bool SEntityPointer<codim,GridImp>::equals (const SEntityPointer<codim,GridImp>& i) const
  {
    return (index==i.index)&&(l==i.l)&&(grid==i.grid);
  }

  template<int codim, class GridImp>
  inline typename SEntityPointer<codim,GridImp>::Entity& SEntityPointer<codim,GridImp>::dereference () const
  {
    return entity();
  }

  template<int codim, class GridImp>
  inline int SEntityPointer<codim,GridImp>::level () const
  {
    return l;
  }


  //************************************************************************
  // inline methods for SGrid
  template<int dim, int dimworld, typename ctype>
  inline void SGrid<dim,dimworld,ctype>::makeSGrid (const array<int,dim>& N_,
                                                    const FieldVector<ctype,dim>& L_, const FieldVector<ctype,dim>& H_)
  {
    static_assert(dimworld <= std::numeric_limits<int>::digits,"world dimension too high, must be <= # of bits of int");

#ifndef NDEBUG
    bool correct = true;
    for (int i=0; i<dim; i++)
      if (H_[i] < L_[i])
        correct = false;
    if (!correct)
    {
      DUNE_THROW(GridError, "Orientation of lower left and upper right corner is wrong: lower = "
                 << L_ << " upper = " << H_);
    }
#endif

    N.resize(MAXL);
    h.resize(MAXL);
    mapper = new CubeMapper<dim>[MAXL];

    indexsets.push_back( new SGridLevelIndexSet<const SGrid<dim,dimworld> >(*this,0) );

    L = 1;
    low = L_;
    H = H_;
    N[0] = N_;

    // define coarse mesh
    for (int i=0; i<MAXL; i++)
    {
      for (int d=0; d<dim; d++)
      {
        N[i][d] = N_[d] * (1<<i);
        h[i][d] = (H[d]-low[d])/((ctype)N[i][d]);
      }
      mapper[i].make(N[i]);
    }

    dinfo << "level=" << L-1 << " size=(" << N[L-1][0];
    for (int i=1; i<dim; i++) dinfo << "," <<  N[L-1][i];
    dinfo << ")" << std::endl;

    // initialize boundary segment mappers
    boundarysize = 0;
    for (int d=0; d<dim; d++)
    {
      array<int,dim-1> fsize;
      for (int i=0; i<d; i++)
        fsize[i] = N_[i];
      for (int i=d+1; i<dim; i++)
        fsize[i-1] = N_[i];
      boundarymapper[d].make(fsize);
      boundarysize += 2 * boundarymapper[d].elements(0);
    }
  }

  template<int dim, int dimworld, typename ctype>
  inline SGrid<dim,dimworld,ctype>::SGrid (const int * const N_, const ctype * const H_)
  {
    static_assert(dimworld <= std::numeric_limits<int>::digits,"world dimension too high, must be <= # of bits of int");

    array<int,dim> N;
    FieldVector<ctype,dim> L(0.0);
    FieldVector<ctype,dim> H(0.0);
    for (int i=0; i<dim; i++ ) N[i] = N_[i];
    for (int i=0; i<dim; i++ ) H[i] = H_[i];

    makeSGrid(N, L, H);
  }

  template<int dim, int dimworld, typename ctype>
  inline SGrid<dim,dimworld,ctype>::SGrid (const int * const N_, const ctype * const L_, const ctype * const H_)
  {
    static_assert(dimworld <= std::numeric_limits<int>::digits, "dimworld is too large!");

    array<int,dim> N;
    FieldVector<ctype,dim> L(0.0);
    FieldVector<ctype,dim> H(0.0);
    for (int i=0; i<dim; i++ ) N[i] = N_[i];
    for (int i=0; i<dim; i++ ) L[i] = L_[i];
    for (int i=0; i<dim; i++ ) H[i] = H_[i];

    makeSGrid(N, L, H);
  }

  template<int dim, int dimworld, typename ctype>
  inline SGrid<dim,dimworld,ctype>::SGrid (FieldVector<int,dim> N_, FieldVector<ctype,dim> L_,
                                           FieldVector<ctype,dim> H_)
  {
    static_assert(dimworld <= std::numeric_limits<int>::digits, "dimworld is too large!");

    array<int,dim> N;
    for (int i=0; i<dim; i++ ) N[i] = N_[i];
    makeSGrid(N, L_, H_);
  }


  template<int dim, int dimworld, typename ctype>
  inline SGrid<dim,dimworld,ctype>::SGrid ()
  {
    static_assert(dimworld <= std::numeric_limits<int>::digits, "dimworld is too large!");

    array<int,dim> N_;
    FieldVector<ctype,dim> L_(0.0);
    FieldVector<ctype,dim> H_(1.0);

    for(int i = 0; i < dim; ++i)
      N_[i] = 1;

    makeSGrid(N_, L_, H_);
  }

  template<int dim, int dimworld, typename ctype>
  inline SGrid<dim,dimworld,ctype>::~SGrid ()
  {
    for (size_t i=0; i<indexsets.size(); i++)
      delete indexsets[i];

    delete[] mapper;
  }

  template<int dim, int dimworld, typename ctype>
  inline void SGrid<dim,dimworld,ctype>::globalRefine (int refCount)
  {
    for(int ref=0; ref<refCount; ref++)
    {
      L++;
      indexsets.push_back( new SGridLevelIndexSet<const SGrid<dim,dimworld> >(*this,maxLevel()) );
    }
  }

  template<int dim, int dimworld, typename ctype>
  inline int SGrid<dim,dimworld,ctype>::maxLevel () const
  {
    return L-1;
  }

  template <int dim, int dimworld,class ctype> template <int cd, PartitionIteratorType pitype>
  inline typename SGrid<dim,dimworld,ctype>::Traits::template Codim<cd>::template Partition<pitype>::LevelIterator
  SGrid<dim,dimworld,ctype>::lbegin (int level) const
  {
    if (pitype == Ghost_Partition)
      return lend<cd, pitype> (level);
    return SLevelIterator<cd,pitype,const SGrid<dim,dimworld> > (this,level,0);
  }

  template <int dim, int dimworld,class ctype> template <int cd, PartitionIteratorType pitype>
  inline typename SGrid<dim,dimworld,ctype>::Traits::template Codim<cd>::template Partition<pitype>::LevelIterator
  SGrid<dim,dimworld,ctype>::lend (int level) const
  {
    return SLevelIterator<cd,pitype,const SGrid<dim,dimworld> > (this,level,size(level,cd));
  }

  template <int dim, int dimworld,class ctype> template <int cd, PartitionIteratorType pitype>
  inline typename SGrid<dim,dimworld,ctype>::Traits::template Codim<cd>::template Partition<pitype>::LeafIterator
  SGrid<dim,dimworld,ctype>::leafbegin () const
  {
    if (pitype == Ghost_Partition)
      return leafend<cd, pitype> ();
    return SLevelIterator<cd,pitype,const SGrid<dim,dimworld> > (this,maxLevel(),0);
  }

  template <int dim, int dimworld,class ctype> template <int cd, PartitionIteratorType pitype>
  inline typename SGrid<dim,dimworld,ctype>::Traits::template Codim<cd>::template Partition<pitype>::LeafIterator
  SGrid<dim,dimworld,ctype>::leafend () const
  {
    return SLevelIterator<cd,pitype,const SGrid<dim,dimworld> > (this,maxLevel(),size(maxLevel(),cd));
  }

  template<int dim, int dimworld, typename ctype>
  inline int SGrid<dim,dimworld,ctype>::size (int level, int codim) const
  {
    return mapper[level].elements(codim);
  }

  template<int dim, int dimworld, typename ctype>
  inline int SGrid<dim,dimworld,ctype>::global_size (int codim) const
  {
    int gSize = 0;
    for(int i=0; i <= maxLevel(); i++)
      gSize += size(i,codim);
    return gSize;
  }

  template<int dim, int dimworld, typename ctype>
  inline FieldVector<ctype, dimworld> SGrid<dim,dimworld,ctype>::pos (int level, array<int,dim>& z) const
  {
    FieldVector<ctype, dimworld> x;
    for (int k=0; k<dim; k++)
      x[k] = (z[k]*h[level][k])*0.5 + low[k];

    // Fill up additional coordinates with zero
    for (int k=dim; k<dimworld; k++)
      x[k] = 0;

    return x;
  }

  template<int dim, int dimworld, typename ctype>
  inline int SGrid<dim,dimworld,ctype>::calc_codim (int level, const array<int,dim>& z) const
  {
    return mapper[level].codim(z);
  }

  template<int dim, int dimworld, typename ctype>
  inline int SGrid<dim,dimworld,ctype>::n (int level, const array<int,dim>& z) const
  {
    return mapper[level].n(z);
  }

  template<int dim, int dimworld, typename ctype>
  inline array<int,dim> SGrid<dim,dimworld,ctype>::z (int level, int i, int codim) const
  {
    return mapper[level].z(i,codim);
  }

  template<int dim, int dimworld, typename ctype>
  inline array<int,dim> SGrid<dim,dimworld,ctype>::subz (const array<int,dim> & z, int i, int codim) const
  {
    // map to old numbering
    const int j = SGridInternal::CubeNumberingTable<dim>::generic2dune( i, codim );

    // find expanded coordinates of entity in reference cube
    // has components in {0,1,2}
    array<int,dim> zref = SUnitCubeMapper<dim>::mapper.z(j,codim);

    // compute expanded coordinates of entity in global coordinates
    array<int,dim> zentity;
    for (int k=0; k<dim; k++) zentity[k] = z[k] + zref[k] - 1;

    return zentity;
  }



  template<int dim, int dimworld, typename ctype>
  inline array<int,dim> SGrid<dim,dimworld,ctype>::compress (int level, const array<int,dim>& z) const
  {
    return mapper[level].compress(z);
  }

  template<int dim, int dimworld, typename ctype>
  inline array<int,dim> SGrid<dim,dimworld,ctype>::expand (int level, const array<int,dim>& r, int b) const
  {
    return mapper[level].expand(r,b);
  }

  template<int dim, int dimworld, typename ctype>
  inline int SGrid<dim,dimworld,ctype>::partition (int level, const array<int,dim>& z) const
  {
    return mapper[level].partition(z);
  }

  template<int dim, int dimworld, typename ctype>
  inline bool SGrid<dim,dimworld,ctype>::exists (int level, const array<int,dim>& zred) const
  {
    for (int i=0; i<dim; i++)
    {
      if (zred[i]<0) return false;
      if (zred[i]>=N[level][i]) return false;
    }
    return true;
  }


} // end namespace

#endif
