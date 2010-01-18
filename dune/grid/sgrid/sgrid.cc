// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_SGRID_CC
#define DUNE_SGRID_CC

#include <algorithm>
#include <iostream>
#include <assert.h>

#include <dune/common/stdstreams.hh>
#include <dune/common/typetraits.hh>

namespace Dune {


  //************************************************************************
  // SGeometry

  template<int mydim, int cdim, class GridImp>
  void SGeometry<mydim,cdim,GridImp>::make(FieldMatrix<typename GridImp::ctype,mydim+1,cdim>& __As)
  {
    // clear jacobian
    builtinverse = false;

    // copy arguments
    s = __As[mydim];
    for (int j=0; j<mydim; j++) A[j] = __As[j];

    // make corners
    for (int i=0; i<(1<<mydim); i++)     // there are 2^d corners
    {
      // use binary representation of corner number to assign corner coordinates
      int mask = 1;
      c[i] = s;
      for (int k=0; k<cdim; k++)
      {
        if (i&mask) c[i] = c[i]+A[k];
        mask = mask<<1;
      }
    }

    // compute centroid
    centroid = 0.0;
    for (int i=0; i<(1<<mydim); i++)
      centroid += c[i];
    centroid *= 1.0/(1<<mydim);
  }

  template<int mydim, int cdim, class GridImp>
  inline FieldVector<typename GridImp::ctype, cdim> SGeometry<mydim,cdim,GridImp>::global (const FieldVector<typename GridImp::ctype, mydim>& local) const
  {
    FieldVector<ctype, cdim> global = s;
    // global += A^t * local
    A.umtv(local,global);

    return global;
  }

  template<int mydim, int cdim, class GridImp>
  inline FieldVector<typename GridImp::ctype, mydim> SGeometry<mydim,cdim,GridImp>::local (const FieldVector<typename GridImp::ctype, cdim>& global) const
  {
    FieldVector<ctype, mydim> l;     // result
    FieldVector<ctype, cdim> rhs = global-s;
    for (int k=0; k<mydim; k++)
      l[k] = (rhs*A[k]) / (A[k]*A[k]);
    return l;
  }

  template<int mydim, int cdim, class GridImp>
  inline typename GridImp::ctype SGeometry<mydim,cdim,GridImp>::volume () const
  {
    sgrid_ctype s = 1.0;
    for (int j=0; j<mydim; j++) s *= A[j].one_norm();

    return s;
  }

  template< int mydim, int cdim, class GridImp >
  inline const FieldMatrix< typename GridImp::ctype, mydim, cdim > &
  SGeometry< mydim, cdim, GridImp >::jacobianTransposed ( const FieldVector< typename GridImp::ctype, mydim > &local ) const
  {
    return A;
  }

  template<int mydim, int cdim, class GridImp>
  inline const FieldMatrix<typename GridImp::ctype,cdim,mydim>& SGeometry<mydim,cdim,GridImp>::jacobianInverseTransposed (const FieldVector<typename GridImp::ctype, mydim>& local) const
  {
    if (!builtinverse)
    {
      // transpose A and invert non-zero entries
      for (int j=0; j<cdim; ++j)
      {
        for (int i=0; i<mydim; ++i)
        {
          if (j<i || std::abs(A[i][j]) < 1e-15)
            Jinv[j][i] = 0.0;
          else
            Jinv[j][i] = 1.0/A[i][j];
        }
      }
      builtinverse = true;
    }
    return Jinv;
  }

  template<int mydim, int cdim, class GridImp>
  inline void SGeometry<mydim,cdim,GridImp>::print (std::ostream& ss, int indent) const
  {
    for (int k=0; k<indent; k++) ss << " ";ss << "SGeometry<" << mydim << "," << cdim << ">" << std::endl;
    for (int k=0; k<indent; k++) ss << " ";ss << "{" << std::endl;
    for (int k=0; k<indent+2; k++) ss << " ";ss << "Position: " << s << std::endl;
    for (int j=0; j<mydim; j++)
    {
      for (int k=0; k<indent+2; k++) ss << " ";
      ss << "direction " << j << "  " << A(j) << std::endl;
    }
    for (int j=0; j<1<<mydim; j++)
    {
      for (int k=0; k<indent+2; k++) ss << " ";
      ss << "corner " << j << "  " << c[j] << std::endl;
    }
    if (builtinverse)
    {
      for (int k=0; k<indent+2; k++) ss << " ";ss << "Jinv ";
      Jinv.print(ss,indent+2);
    }
    for (int k=0; k<indent+2; k++) ss << " ";ss << "builtinverse " << builtinverse << std::endl;
    for (int k=0; k<indent; k++) ss << " ";ss << "}";
  }

  template<int cdim, class GridImp>
  inline void SGeometry<0,cdim,GridImp>::make (FieldMatrix<typename GridImp::ctype,1,cdim>& __As)
  {
    s = __As[0];
  }

  template<int cdim, class GridImp>
  inline void SGeometry<0,cdim,GridImp>::print (std::ostream& ss, int indent) const
  {
    for (int i=0; i<indent; i++) ss << " ";
    ss << "SGeometry<0," << cdim << "> at position " << s;
  }

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
    FieldMatrix<ctype,dim-codim+1,dimworld> __As;

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
        __As[dir] = p2-p1;
        dir++;
      }

    // find reference point, subtract 1 from all odd directions
    for (int i=0; i<dim; i++)
      t[i] -= t[i]%2;
    __As[dir] =grid->pos(l,t);     // all components of t are even

    // make element
    grid->getRealImplementation(geo).make(__As);
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

  // subentity construction
  template<int dim, class GridImp> template<int cc>
  inline typename SEntity<0,dim,GridImp>::template Codim<cc>::EntityPointer SEntity<0,dim,GridImp>::subEntity (int i) const
  {
    // make Iterator
    return SLevelIterator<cc,All_Partition,const GridImp>(grid,l,this->subCompressedIndex(cc,i));
  }

  template<int dim, class GridImp>
  inline typename SEntity<0,dim,GridImp>::IntersectionIterator SEntity<0,dim,GridImp>::ibegin () const
  {
    return IntersectionIterator(SIntersectionIterator<GridImp>(grid,this,0));
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
    return IntersectionIterator(SIntersectionIterator<GridImp>(grid,this,count<1>()));
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
    if (l<=0)
    {
      father_index = 0;
      built_father = true;
      return;
    }

    // reduced coordinates from expanded coordinates
    array<int,dim> zz = grid->compress(l,z);

    // look for odd coordinates
    FieldVector<ctype, dim> delta;
    for (int i=0; i<dim; i++)
    {
      delta[i] = zz[i] % 2;
      zz[i] = zz[i] / 2;
    }

    // zz is now the reduced coordinate of the father, compute index
    int partition = grid->partition(l,z);
    father_index = grid->n((l)-1,grid->expand((l)-1,zz,partition));

    // now make a subcube of size 1/2 in each direction
    FieldMatrix<ctype,dim+1,dim> __As;
    FieldVector<ctype, dim> v;
    for (int i=0; i<dim; i++)
    {
      v = 0.0; v[i] = 0.5;
      __As[i] = v;
    }
    for (int i=0; i<dim; i++) v[i] = 0.5*delta[i];
    __As[dim] =v;
    grid->getRealImplementation(in_father_local).make(__As);     // build geometry

    built_father = true;
  }

  template<int dim, class GridImp>
  inline typename SEntity<0,dim,GridImp>::EntityPointer SEntity<0,dim,GridImp>::father () const
  {
    if (!built_father) make_father();
    if (l>0)
      return SLevelIterator<0,All_Partition,const GridImp>((grid),(l)-1,father_index);
    else
      return SLevelIterator<0,All_Partition,const GridImp>((grid),l,index);
  }

  template<int dim, class GridImp>
  inline
  const typename GridImp::template Codim<0>::LocalGeometry&
  SEntity<0,dim,GridImp>::geometryInFather () const
  {
    if (!built_father) make_father();
    return in_father_local;
  }

  template<int dim, class GridImp>
  inline typename SEntity<0,dim,GridImp>::HierarchicIterator SEntity<0,dim,GridImp>::hbegin (int maxLevel) const
  {
    return HierarchicIterator(SHierarchicIterator<GridImp>(grid,*this,maxLevel,false));
  }

  template<int dim, class GridImp>
  inline typename SEntity<0,dim,GridImp>::HierarchicIterator SEntity<0,dim,GridImp>::hend (int maxLevel) const
  {
    return HierarchicIterator(SHierarchicIterator<GridImp>(grid,*this,maxLevel,true));
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
    SHierarchicStackElem newe = stack.pop();
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
      return;     // done, this is end iterator
    }
    valid_count = true;

    // and compute compressed coordinates of neighbor
    array<int,dim> zrednb = zred;
    zrednb[count/2] += -1+2*(count%2);     // (count%2) ? +1 : -1

    // now check if neighbor exists
    is_on_boundary = !grid->exists(grid->getRealImplementation(self).l,zrednb);
    if (is_on_boundary)
    {
      grid->getRealImplementation(ne).index = -1;
      return;     // ok, done it
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
  inline bool SIntersectionIterator<GridImp>::conforming () const
  {
    return true;
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
      z1[dir] += 1;           // odd
    else
      z1[dir] -= 1;           // even

    // z1 is even in direction dir, all others must be odd because it is codim 1
    FieldMatrix<ctype,dim,dim> __AsLocal;
    FieldVector<ctype, dim> p1Local,p2Local;

    int t;

    // local coordinates in self
    p1Local = 0.0;
    p1Local[dir] = c;        // all points have p[dir]=c in entity
    __AsLocal[dim-1] = p1Local;     // position vector
    t = 0;
    for (int i=0; i<dim; ++i)     // this loop makes dim-1 direction vectors
      if (i!=dir)
      {
        // each i!=dir gives one direction vector
        p2Local = p1Local;
        p2Local[i] = 1.0;
        __AsLocal[t] = p2Local-p1Local;                 // a direction vector
        ++t;
      }
    // update geometry
    grid->getRealImplementation(is_self_local).make(__AsLocal);

    // local coordinates in neighbor
    p1Local = 0.0;
    p1Local[dir] = 1-c;        // all points have p[dir]=1-c in entity
    __AsLocal[dim-1] = p1Local;       // position vector
    t = 0;
    for (int i=0; i<dim; ++i)     // this loop makes dim-1 direction vectors
      if (i!=dir)
      {
        // each i!=dir gives one direction vector
        p2Local = p1Local;
        p2Local[i] = 1.0;
        __AsLocal[t] = p2Local-p1Local;                 // a direction vector
        ++t;
      }
    // update geometry
    grid->getRealImplementation(is_nb_local).make(__AsLocal);

    // global coordinates
    FieldMatrix<ctype,dim,dimworld> __As;
    FieldVector<ctype, dimworld> p1,p2;
    t = 0;
    for (int i=0; i<dim; i++)
      if (i!=dir)
      {
        // each i!=dir gives one direction vector
        z1[i] += 1;                 // direction i => even
        p2 = grid->pos(self.level(),z1);
        z1[i] -= 2;                 // direction i => even
        p1 = grid->pos(self.level(),z1);
        z1[i] += 1;                 // revert t to original state
        __As[t] = p2-p1;
        ++t;
      }
    for (int i=0; i<dim; i++)
      if (i!=dir)
        z1[i] -= 1;
    __As[t] = grid->pos(self.level(),z1);
    // update geometry
    grid->getRealImplementation(is_global).make(__As);

    built_intersections = true;
  }

  template<class GridImp>
  inline const typename SIntersectionIterator< GridImp >::LocalGeometry &
  SIntersectionIterator< GridImp >::geometryInInside () const
  {
    assert (valid_count);
    if (!built_intersections) makeintersections();
    return is_self_local;
  }

  template<class GridImp>
  inline const typename SIntersectionIterator< GridImp >::LocalGeometry &
  SIntersectionIterator< GridImp >::geometryInOutside () const
  {
    assert (valid_count);
    if (!built_intersections) makeintersections();
    return is_nb_local;
  }

  template<class GridImp>
  inline const typename SIntersectionIterator< GridImp >::Geometry &
  SIntersectionIterator< GridImp >::geometry () const
  {
    assert (valid_count);
    if (!built_intersections) makeintersections();
    return is_global;
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
  inline void SGrid<dim,dimworld,ctype>::makeSGrid (const int* N_,
                                                    const ctype* L_, const ctype* H_)
  {
    dune_static_assert(dimworld <= std::numeric_limits<int>::digits,"world dimension too high, must be <= # of bits of int");

#ifndef NDEBUG
    bool correct = true;
    for (int i=0; i<dim; i++)
      if (H_[i] < L_[i])
        correct = false;
    if (!correct)
    {
      FieldVector<ctype,dim> L, H;
      for (int i=0; i<dim; ++i)
      {
        L[i] = L_[i]; H[i] = H_[i];
      }
      DUNE_THROW(GridError, "Orientation of lower left and upper right corner is wrong: lower = "
                 << L << " upper = " << H);
    }
#endif

    N = new array<int,dim>[MAXL];
    h = new FieldVector<ctype, dim>[MAXL];
    mapper = new CubeMapper<dim>[MAXL];

    indexsets.push_back( new SGridLevelIndexSet<const SGrid<dim,dimworld> >(*this,0) );
    theleafindexset = new SGridLeafIndexSet<const SGrid<dim,dimworld> >(*this);
    theglobalidset = new SGridGlobalIdSet<const SGrid<dim,dimworld> >(*this);

    L = 1;
    for (int i=0; i<dim; i++) low[i] = L_[i];
    for (int i=0; i<dim; i++) H[i] = H_[i];
    for (int i=0; i<dim; i++) N[0][i] = N_[i];

    // define coarse mesh
    mapper[0].make(N[0]);
    for (int i=0; i<dim; i++)
      h[0][i] = (H[i]-low[i])/((ctype)N[0][i]);

    // define boudary segments
    boundarysize = 0;
    {
      array<int, dim-1> bN;
      for (int dir = 0; dir < dim; dir++)
      {
        // extent of this face
        for (int i=0; i<dir; i++) bN[i] = N[0][i];
        for (int i=dir+1; i<dim; i++) bN[i-1] = N[0][i];
        // setup mapper for this face
        boundarymapper[dir].make(bN);
        // compute size
        boundarysize += 2 * boundarymapper[dir].elements(0);
      }
    }

    dinfo << "level=" << L-1 << " size=(" << N[L-1][0];
    for (int i=1; i<dim; i++) dinfo << "," <<  N[L-1][i];
    dinfo << ")" << std::endl;
  }

  template<int dim, int dimworld, typename ctype>
  inline SGrid<dim,dimworld,ctype>::SGrid (const int * const N_, const ctype * const H_)
  {
    dune_static_assert(dimworld <= std::numeric_limits<int>::digits,"world dimension too high, must be <= # of bits of int");

    ctype L_[dim];
    for (int i=0; i<dim; i++)
      L_[i] = 0;

    makeSGrid(N_,L_, H_);
  }

  template<int dim, int dimworld, typename ctype>
  inline SGrid<dim,dimworld,ctype>::SGrid (const int * const N_, const ctype * const L_, const ctype * const H_)
  {
    dune_static_assert(dimworld <= std::numeric_limits<int>::digits, "dimworld is too large!");

    makeSGrid(N_, L_, H_);
  }

  template<int dim, int dimworld, typename ctype>
  inline SGrid<dim,dimworld,ctype>::SGrid (FieldVector<int,dim> N_, FieldVector<ctype,dim> L_,
                                           FieldVector<ctype,dim> H_)
  {
    dune_static_assert(dimworld <= std::numeric_limits<int>::digits, "dimworld is too large!");

    ctype LL[dim], HH[dim];
    int NN[dim];

    for (int i=0; i<dim; ++i)
    {
      LL[i] = L_[i]; HH[i] = H_[i]; NN[i] = N_[i];
    }

    makeSGrid(NN, LL, HH);
  }


  template<int dim, int dimworld, typename ctype>
  inline SGrid<dim,dimworld,ctype>::SGrid ()
  {
    dune_static_assert(dimworld <= std::numeric_limits<int>::digits, "dimworld is too large!");

    int N_[dim];
    ctype L_[dim];
    ctype H_[dim];

    for(int i = 0; i < dim; ++i) {
      N_[i] = 1;
      L_[i] = 0.0;
      H_[i] = 1.0;
    }

    makeSGrid(N_, L_, H_);
  }

  template<int dim, int dimworld, typename ctype>
  inline SGrid<dim,dimworld,ctype>::~SGrid ()
  {
    for (size_t i=0; i<indexsets.size(); i++)
      delete indexsets[i];

    delete theleafindexset;
    delete theglobalidset;


    delete[] N;
    delete[] h;
    delete[] mapper;
  }

  template<int dim, int dimworld, typename ctype>
  inline void SGrid<dim,dimworld,ctype>::globalRefine (int refCount)
  {
    for(int ref=0; ref<refCount; ref++)
    {

      // refine the mesh
      for (int i=0; i<dim; i++) N[L][i] = 2*N[L-1][i];
      mapper[L].make(N[L]);

      // compute mesh size
      for (int i=0; i<dim; i++)
        h[L][i] = (H[i]-low[i])/((ctype)N[L][i]);
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
    static const GeometryType cubeType(dim);

    // map to old numbering
    typedef GenericGeometry::MapNumberingProvider< dim > Numbering;
    const unsigned int cubeid = GenericGeometry::topologyId( cubeType );
    const int j = Numbering::template generic2dune( cubeid, i, codim );

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
