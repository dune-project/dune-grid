// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_YASPGRIDENTITY_HH
#define DUNE_GRID_YASPGRIDENTITY_HH

/** \file
 * \brief the YaspEntity class and its specializations
 *
   YaspEntity realizes the concept a mesh entity.

   We have specializations for codim==0 (elements) and
   codim=dim (vertices).
   The general version throws a GridError.
 */
//========================================================================




namespace Dune {


  template<int codim, int dim, class GridImp>
  class YaspEntity
    :  public EntityDefaultImplementation <codim,dim,GridImp,YaspEntity>
  {
  public:
    typedef typename GridImp::ctype ctype;

    typedef typename GridImp::template Codim<codim>::Geometry Geometry;

    //! level of this element
    int level () const
    {
      DUNE_THROW(GridError, "YaspEntity not implemented");
    }

    //! index is unique and consecutive per level and codim used for access to degrees of freedom
    int index () const
    {
      DUNE_THROW(GridError, "YaspEntity not implemented");
    }

    //! geometry of this entity
    Geometry geometry () const
    {
      DUNE_THROW(GridError, "YaspEntity not implemented");
    }

    //! return partition type attribute
    PartitionType partitionType () const
    {
      DUNE_THROW(GridError, "YaspEntity not implemented");
    }

    const GridImp * yaspgrid() const
    {
      DUNE_THROW(GridError, "YaspEntity not implemented");
    }

    typedef typename GridImp::YGridLevelIterator YGLI;
    typedef typename SubYGrid<dim,ctype>::TransformingSubIterator TSI;
    YaspEntity (const GridImp* yg, const YGLI& g, const TSI& it)
    {
      DUNE_THROW(GridError, "YaspEntity not implemented");
    }

    // IndexSets needs access to the private index methods
    friend class Dune::YaspIndexSet<GridImp,true>;
    friend class Dune::YaspIndexSet<GridImp,false>;
    friend class Dune::YaspGlobalIdSet<GridImp>;
    typedef typename GridImp::PersistentIndexType PersistentIndexType;

    //! globally unique, persistent index
    PersistentIndexType persistentIndex () const
    {
      DUNE_THROW(GridError, "YaspEntity not implemented");
    }

    //! consecutive, codim-wise, level-wise index
    int compressedIndex () const
    {
      DUNE_THROW(GridError, "YaspEntity not implemented");
    }

    //! consecutive, codim-wise, level-wise index
    int compressedLeafIndex () const
    {
      DUNE_THROW(GridError, "YaspEntity not implemented");
    }

    //! subentity compressed index (not available here)
    int subCompressedIndex (int, unsigned int ) const
    {
      DUNE_THROW(NotImplemented,"subIndex for entities with codimension > 0 is not implemented");
      return -1;
    }

    //! subentity compressed leaf index (not available here)
    int subCompressedLeafIndex (int, unsigned int ) const
    {
      DUNE_THROW(NotImplemented,"subIndex for entities with codimension > 0 is not implemented");
      return -1;
    }
  };


  // specialization for codim=0
  template<int dim, class GridImp>
  class YaspEntity<0,dim,GridImp>
    : public EntityDefaultImplementation <0,dim,GridImp,YaspEntity>
  {
    enum { dimworld = GridImp::dimensionworld };

    typedef typename GridImp::Traits::template Codim< 0 >::GeometryImpl GeometryImpl;

  public:
    typedef typename GridImp::ctype ctype;

    typedef typename GridImp::YGridLevelIterator YGLI;
    typedef typename SubYGrid<dim,ctype>::TransformingSubIterator TSI;

    typedef typename GridImp::template Codim< 0 >::Geometry Geometry;
    typedef typename GridImp::template Codim< 0 >::LocalGeometry LocalGeometry;

    template <int cd>
    struct Codim
    {
      typedef typename GridImp::template Codim<cd>::EntityPointer EntityPointer;
    };

    typedef typename GridImp::template Codim<0>::EntityPointer EntityPointer;
    typedef typename GridImp::template Codim<0>::EntitySeed EntitySeed;
    typedef typename GridImp::LevelIntersectionIterator IntersectionIterator;
    typedef typename GridImp::LevelIntersectionIterator LevelIntersectionIterator;
    typedef typename GridImp::LeafIntersectionIterator LeafIntersectionIterator;
    typedef typename GridImp::HierarchicIterator HierarchicIterator;

    //! define the type used for persisitent indices
    typedef typename GridImp::PersistentIndexType PersistentIndexType;

    //! define type used for coordinates in grid module
    typedef typename YGrid<dim,ctype>::iTupel iTupel;

    // constructor
    YaspEntity (const GridImp * yg, const YGLI& g, const TSI& it)
      : _yg(yg), _it(it), _g(g)
    {}

    //! level of this element
    int level () const { return _g->level(); }

    //! index is unique and consecutive per level
    int index () const { return _it.superindex(); } // superindex works also for iteration over subgrids

    //! globalIndex is unique and consecutive per global level
    int globalIndex () const {
      return _g.cell_global().index(_it.coord());
    }

    /** \brief Return the entity seed which contains sufficient information
     *  to generate the entity again and uses as little memory as possible
     */
    EntitySeed seed () const {
      return EntitySeed(YaspEntitySeed<0,GridImp>(_g->level(), _it.coord()));
    }

    //! return partition type attribute
    PartitionType partitionType () const
    {
      if (_g->cell_interior.inside(_it.coord()))
        return InteriorEntity;
      if (_g->cell_overlap.inside(_it.coord()))
        return OverlapEntity;
      DUNE_THROW(GridError, "Impossible GhostEntity " << _it.coord() << "\t"
                                                      << _g->cell_interior.origin() << "/" << _g->cell_interior.size());
      return GhostEntity;
    }

    //! geometry of this entity
    Geometry geometry () const {
      // the element geometry
      GeometryImpl _geometry(_it.position(),_it.meshsize());
      return Geometry( _geometry );
    }

    /*! Return number of subentities with codimension cc.
     */
    template<int cc> int count () const
    {
      if (cc==dim) return 1<<dim;
      if (cc==1) return 2*dim;
      if (cc==dim-1) return dim*(1<<(dim-1));
      if (cc==0) return 1;
      DUNE_THROW(GridError, "codim " << cc << " (dim=" << dim << ") not (yet) implemented");
    }

    /*! Intra-element access to subentities of codimension cc > codim.
     */
    template<int cc>
    typename Codim<cc>::EntityPointer subEntity (int i) const
    {
      dune_static_assert( cc == dim || cc == 0 ,
                          "YaspGrid only supports Entities with codim=dim and codim=0");
      // coordinates of the cell == coordinates of lower left corner
      if (cc==dim)
      {
        iTupel coord = _it.coord();

        // get corner from there
        for (int k=0; k<dim; k++)
          if (i&(1<<k)) (coord[k])++;

        return YaspEntityPointer<cc,GridImp>(_yg,_g,_g->vertex_overlapfront.tsubbegin(coord));
      }
      if (cc==0)
      {
        return YaspEntityPointer<cc,GridImp>(_yg,_g,_it);
      }
      DUNE_THROW(GridError, "codim " << cc << " (dim=" << dim << ") not (yet) implemented");
    }

    //! Inter-level access to father element on coarser grid. Assumes that meshes are nested.
    EntityPointer father () const
    {
      // check if coarse level exists
      if (_g->level()<=0)
        DUNE_THROW(GridError, "tried to call father on level 0");

      // yes, get iterator to it
      YGLI cg = _g;
      --cg;

      // coordinates of the cell
      iTupel coord = _it.coord();

      // get coordinates on next coarser level
      for (int k=0; k<dim; k++) coord[k] = coord[k]/2;

      return YaspEntityPointer<0,GridImp>(_yg,cg,cg->cell_overlap.tsubbegin(coord));
    }

    //! returns true if father entity exists
    bool hasFather () const
    {
      return (_g->level()>0);
    }

    /*! Location of this element relative to the reference element of its father
     */
    LocalGeometry geometryInFather () const
    {
      // configure one of the 2^dim transformations
      FieldVector<ctype,dim> midpoint;
      FieldVector<ctype,dim> extension(0.5);

      for (int k=0; k<dim; k++)
        midpoint[k] = (_it.coord(k)%2) ? 0.75 : 0.25;

      return LocalGeometry( YaspGeometry<dim,dim,GridImp>(midpoint,extension) );
    }

    const TSI& transformingsubiterator () const
    {
      return _it;
    }

    const YGLI& gridlevel () const
    {
      return _g;
    }

    const GridImp* yaspgrid () const
    {
      return _yg;
    }

    bool isLeaf() const
    {
      return (_g->level() == _yg->maxLevel());
    }

    /**\brief Returns true, if the entity has been created during the last call to adapt()
     */
    bool isNew () const { return _yg->adaptRefCount > 0 && _yg->maxLevel() < _g->level() + _yg->adaptRefCount; }

    /**\brief Returns true, if entity might disappear during the next call to adapt()
     */
    bool mightVanish () const { return false; }
    // { return _yg->adaptRefCount < 0 && _g.mg()->maxlevel() < _g.level() - _yg->adaptRefCount; }

    //! returns intersection iterator for first intersection
    IntersectionIterator ibegin () const
    {
      return YaspIntersectionIterator<GridImp>(*this,false);
    }

    //! returns intersection iterator for first intersection
    LeafIntersectionIterator ileafbegin () const
    {
      // only if entity is leaf this iterator delivers intersections
      return YaspIntersectionIterator<GridImp>(*this, ! isLeaf() );
    }

    //! returns intersection iterator for first intersection
    LevelIntersectionIterator ilevelbegin () const
    {
      return ibegin();
    }

    //! Reference to one past the last neighbor
    IntersectionIterator iend () const
    {
      return YaspIntersectionIterator<GridImp>(*this,true);
    }

    //! Reference to one past the last neighbor
    LeafIntersectionIterator ileafend () const
    {
      return iend();
    }

    //! Reference to one past the last neighbor
    LevelIntersectionIterator ilevelend () const
    {
      return iend();
    }

    /*! Inter-level access to son elements on higher levels<=maxlevel.
          This is provided for sparsely stored nested unstructured meshes.
          Returns iterator to first son.
     */
    HierarchicIterator hbegin (int maxlevel) const
    {
      return YaspHierarchicIterator<GridImp>(_yg,_g,_it,maxlevel);
    }

    //! Returns iterator to one past the last son
    HierarchicIterator hend (int maxlevel) const
    {
      return YaspHierarchicIterator<GridImp>(_yg,_g,_it,_g->level());
    }

  private:
    // IndexSets needs access to the private index methods
    friend class Dune::YaspIndexSet<GridImp,true>;
    friend class Dune::YaspIndexSet<GridImp,false>;
    friend class Dune::YaspGlobalIdSet<GridImp>;

    //! globally unique, persistent index
    PersistentIndexType persistentIndex () const
    {
      // get size of global grid
      const iTupel& size =  _g->cell_global.size();

      // get coordinate correction for periodic boundaries
      int coord[dim];
      for (int i=0; i<dim; i++)
      {
        coord[i] = _it.coord(i);
        if (coord[i]<0) coord[i] += size[i];
        if (coord[i]>=size[i]) coord[i] -= size[i];
      }

      // encode codim
      PersistentIndexType id(0);

      // encode level
      id = id << yaspgrid_level_bits;
      id = id+PersistentIndexType(_g->level());


      // encode coordinates
      for (int i=dim-1; i>=0; i--)
      {
        id = id << yaspgrid_dim_bits;
        id = id+PersistentIndexType(coord[i]);
      }

      return id;
    }

    //! consecutive, codim-wise, level-wise index
    int compressedIndex () const
    {
      return _it.superindex();
    }

    //! consecutive, codim-wise, level-wise index
    int compressedLeafIndex () const
    {
      return _it.superindex();
    }

    //! subentity persistent index
    PersistentIndexType subPersistentIndex (int i, int cc) const
    {
      if (cc==0)
        return persistentIndex();

      // get position of cell, note that global origin is zero
      // adjust for periodic boundaries
      int coord[dim];
      for (int k=0; k<dim; k++)
      {
        coord[k] = _it.coord(k);
        if (coord[k]<0)
          coord[k] += _g->cell_global.size(k);
        if (coord[k]>=_g->cell_global.size(k))
          coord[k] -= _g->cell_global.size(k);
      }

      if (cc==dim)
      {
        // transform to vertex coordinates
        for (int k=0; k<dim; k++)
          if (i&(1<<k)) (coord[k])++;

        // determine min number of trailing zeroes
        int trailing = 1000;
        for (int i=0; i<dim; i++)
        {
          // count trailing zeros
          int zeros = 0;
          for (int j=0; j<_g->level(); j++)
            if (coord[i]&(1<<j))
              break;
            else
              zeros++;
          trailing = std::min(trailing,zeros);
        }

        // determine the level of this vertex
        int level = _g->level()-trailing;

        // encode codim
        PersistentIndexType id(dim);

        // encode level
        id = id << yaspgrid_level_bits;
        id = id+PersistentIndexType(level);

        // encode coordinates
        for (int i=dim-1; i>=0; i--)
        {
          id = id << yaspgrid_dim_bits;
          id = id+PersistentIndexType(coord[i]>>trailing);
        }

        return id;
      }

      if (cc==1) // faces, i.e. for dim=2 codim=1 is treated as a face
      {
        // Idea: Use the doubled grid to assign coordinates to faces

        // ivar is the direction that varies
        int ivar=i/2;

        // compute position from cell position
        for (int k=0; k<dim; k++)
          coord[k] = coord[k]*2 + 1; // the doubled grid
        if (i%2)
          coord[ivar] += 1;
        else
          coord[ivar] -= 1;

        // encode codim
        PersistentIndexType id(1);

        // encode level
        id = id << yaspgrid_level_bits;
        id = id+PersistentIndexType(_g->level());

        // encode coordinates
        for (int i=dim-1; i>=0; i--)
        {
          id = id << yaspgrid_dim_bits;
          id = id+PersistentIndexType(coord[i]);
        }

        return id;
      }

      // map to old numbering
      static unsigned int edge[ 12 ] = { 0, 1, 2, 3, 4, 5, 8, 9, 6, 7, 10, 11 };
      i = edge[i];

      if (cc==dim-1) // edges, exist only for dim>2
      {
        // Idea: direction i is fixed, all others are vary, i.e. 2^(dim-1) possibilities per direction

        // number of entities per direction
        int m=1<<(dim-1);

        // ifix is the direction that is fixed
        int ifix=(dim-1)-(i/m);

        // compute position from cell position
        int bit=1;
        for (int k=0; k<dim; k++)
        {
          coord[k] = coord[k]*2+1;   // cell position in doubled grid
          if (k==ifix) continue;
          if ((i%m)&bit) coord[k] += 1;else coord[k] -= 1;
          bit *= 2;
        }

        // encode codim
        PersistentIndexType id(dim-1);

        // encode level
        id = id << yaspgrid_level_bits;
        id = id+PersistentIndexType(_g->level());

        // encode coordinates
        for (int i=dim-1; i>=0; i--)
        {
          id = id << yaspgrid_dim_bits;
          id = id+PersistentIndexType(coord[i]);
        }

        return id;
      }

      DUNE_THROW(GridError, "codim " << cc << " (dim=" << dim << ") not (yet) implemented");
    }

    //! subentity compressed index
    int subCompressedIndex (int i, int cc) const
    {
      if (cc==0)
        return compressedIndex();

      // get cell position relative to origin of local cell grid
      iTupel coord;
      for (int k=0; k<dim; ++k)
        coord[k] = _it.coord(k)-_g->cell_overlap.origin(k);

      if (cc==dim) // vertices
      {
        // transform cell coordinate to corner coordinate
        for (int k=0; k<dim; k++)
          if (i&(1<<k)) (coord[k])++;

        // do lexicographic numbering
        int index = coord[dim-1];
        for (int k=dim-2; k>=0; --k)
          index = (index*(_g->cell_overlap.size(k)+1))+coord[k];
        return index;
      }

      if (cc==1) // faces, i.e. for dim=2 codim=1 is treated as a face
      {
        // Idea: direction ivar varies, all others are fixed, i.e. 2 possibilities per direction

        // ivar is the direction that varies
        int ivar=i/2;

        // compute position from cell position
        if (i%2) coord[ivar] += 1;

        // do lexicographic numbering
        int index = coord[dim-1];
        for (int k=dim-2; k>=0; --k)
          if (k==ivar)
            index = (index*(_g->cell_overlap.size(k)+1))+coord[k]; // one more
          else
            index = (index*(_g->cell_overlap.size(k)))+coord[k];

        // add size of all subsets for smaller directions
        for (int j=0; j<ivar; j++)
        {
          int n=_g->cell_overlap.size(j)+1;
          for (int l=0; l<dim; l++)
            if (l!=j) n *= _g->cell_overlap.size(l);
          index += n;
        }

        return index;
      }

      // map to old numbering
      static unsigned int edge[ 12 ] = { 0, 1, 2, 3, 4, 5, 8, 9, 6, 7, 10, 11 };
      i = edge[i];

      if (cc==dim-1) // edges, exist only for dim>2
      {
        // Idea: direction i is fixed, all others are vary, i.e. 2^(dim-1) possibilities per direction

        // number of entities per direction
        int m=1<<(dim-1);

        // ifix is the direction that is fixed
        int ifix=(dim-1)-(i/m);

        // compute position from cell position
        int bit=1;
        for (int k=0; k<dim; k++)
        {
          if (k==ifix) continue;
          if ((i%m)&bit) coord[k] += 1;
          bit *= 2;
        }

        // do lexicographic numbering
        int index = coord[dim-1];
        for (int k=dim-2; k>=0; --k)
          if (k!=ifix)
            index = (index*(_g->cell_overlap.size(k)+1))+coord[k]; // one more
          else
            index = (index*(_g->cell_overlap.size(k)))+coord[k];

        // add size of all subsets for smaller directions
        for (int j=dim-1; j>ifix; j--)
        {
          int n=_g->cell_overlap.size(j);
          for (int l=0; l<dim; l++)
            if (l!=j) n *= _g->cell_overlap.size(l)+1;
          index += n;
        }

        return index;
      }

      DUNE_THROW(GridError, "codim " << cc << " (dim=" << dim << ") not (yet) implemented");
    }

    //! subentity compressed index
    int subCompressedLeafIndex (int i, int cc) const
    {
      if (cc==0)
        return compressedIndex();

      // get cell position relative to origin of local cell grid
      iTupel coord;
      for (int k=0; k<dim; ++k)
        coord[k] = _it.coord(k)-_g.cell_overlap().origin(k);

      if (cc==dim) // vertices
      {
        // transform cell coordinate to corner coordinate
        for (int k=0; k<dim; k++)
          if (i&(1<<k)) (coord[k])++;

        // move coordinates up to maxlevel
        for (int k=0; k<dim; k++)
          coord[k] = coord[k]<<(_g.mg()->maxlevel()-_g.level());

        // do lexicographic numbering
        int index = coord[dim-1];
        for (int k=dim-2; k>=0; --k)
          index = (index*(_g.mg()->rbegin().cell_overlap().size(k)+1))+coord[k];
        return index;
      }

      if (cc==1) // faces, i.e. for dim=2 codim=1 is treated as a face
      {
        // Idea: direction ivar varies, all others are fixed, i.e. 2 possibilities per direction

        // ivar is the direction that varies
        int ivar=i/2;

        // compute position from cell position
        if (i%2) coord[ivar] += 1;

        // do lexicographic numbering
        int index = coord[dim-1];
        for (int k=dim-2; k>=0; --k)
          if (k==ivar)
            index = (index*(_g.cell_overlap().size(k)+1))+coord[k]; // one more
          else
            index = (index*(_g.cell_overlap().size(k)))+coord[k];

        // add size of all subsets for smaller directions
        for (int j=0; j<ivar; j++)
        {
          int n=_g.cell_overlap().size(j)+1;
          for (int l=0; l<dim; l++)
            if (l!=j) n *= _g.cell_overlap().size(l);
          index += n;
        }

        return index;
      }

      // map to old numbering
      static unsigned int edge[ 12 ] = { 0, 1, 2, 3, 4, 5, 8, 9, 6, 7, 10, 11 };
      i = edge[i];

      if (cc==dim-1) // edges, exist only for dim>2
      {
        // Idea: direction i is fixed, all others are vary, i.e. 2^(dim-1) possibilities per direction

        // number of entities per direction
        int m=1<<(dim-1);

        // ifix is the direction that is fixed
        int ifix=(dim-1)-(i/m);

        // compute position from cell position
        int bit=1;
        for (int k=0; k<dim; k++)
        {
          if (k==ifix) continue;
          if ((i%m)&bit) coord[k] += 1;
          bit *= 2;
        }

        // do lexicographic numbering
        int index = coord[dim-1];
        for (int k=dim-2; k>=0; --k)
          if (k!=ifix)
            index = (index*(_g.cell_overlap().size(k)+1))+coord[k]; // one more
          else
            index = (index*(_g.cell_overlap().size(k)))+coord[k];

        // add size of all subsets for smaller directions
        for (int j=dim-1; j>ifix; j--)
        {
          int n=_g.cell_overlap().size(j);
          for (int l=0; l<dim; l++)
            if (l!=j) n *= _g.cell_overlap().size(l)+1;
          index += n;
        }

        return index;
      }

      DUNE_THROW(GridError, "codim " << cc << " (dim=" << dim << ") not (yet) implemented");
    }

    const GridImp * _yg;    // access to YaspGrid
    const TSI& _it;         // position in the grid level
    const YGLI& _g;         // access to grid level
  };


  // specialization for codim=dim (vertex)
  template<int dim, class GridImp>
  class YaspEntity<dim,dim,GridImp>
    : public EntityDefaultImplementation <dim,dim,GridImp,YaspEntity>
  {
    enum { dimworld = GridImp::dimensionworld };

    typedef typename GridImp::Traits::template Codim<dim>::GeometryImpl GeometryImpl;

  public:
    typedef typename GridImp::ctype ctype;

    typedef typename GridImp::YGridLevelIterator YGLI;
    typedef typename SubYGrid<dim,ctype>::TransformingSubIterator TSI;

    typedef typename GridImp::template Codim<dim>::Geometry Geometry;

    template <int cd>
    struct Codim
    {
      typedef typename GridImp::template Codim<cd>::EntityPointer EntityPointer;
    };

    typedef typename GridImp::template Codim<dim>::EntityPointer EntityPointer;
    typedef typename GridImp::template Codim<dim>::EntitySeed EntitySeed;

    //! define the type used for persisitent indices
    typedef typename GridImp::PersistentIndexType PersistentIndexType;

    //! define type used for coordinates in grid module
    typedef typename YGrid<dim,ctype>::iTupel iTupel;

    // constructor
    YaspEntity (const GridImp* yg, const YGLI& g, const TSI& it)
      : _yg(yg), _it(it), _g(g)
    {}

    //! level of this element
    int level () const {return _g->level();}

    //! index is unique and consecutive per level
    int index () const {return _it.superindex();}

    //! globally unique, persistent index
    int globalIndex () const { return _g.cell_global().index(_it.coord()); }

    /** \brief Return the entity seed which contains sufficient information
     *  to generate the entity again and uses as little memory as possible
     */
    EntitySeed seed () const {
      return EntitySeed(YaspEntitySeed<dim,GridImp>(_g->level(), _it.coord()));
    }

    //! geometry of this entity
    Geometry geometry () const {
      GeometryImpl _geometry(_it.position());
      return Geometry( _geometry );
    }

    //! return partition type attribute
    PartitionType partitionType () const
    {
      if (_g->vertex_interior.inside(_it.coord()))
        return InteriorEntity;
      if (_g->vertex_interiorborder.inside(_it.coord()))
        return BorderEntity;
      if (_g->vertex_overlap.inside(_it.coord()))
        return OverlapEntity;
      if (_g->vertex_overlapfront.inside(_it.coord()))
        return FrontEntity;
      return GhostEntity;
    }

    //! subentity compressed index simply returns compressedIndex
    int subCompressedIndex (int, unsigned int ) const
    {
      return compressedIndex();
    }

    //! subentity compressed leaf index simply returns compressedLeafIndex
    int subCompressedLeafIndex (int, unsigned int ) const
    {
      return compressedLeafIndex();
    }

  private:
    // IndexSets needs access to the private index methods
    friend class Dune::YaspIndexSet<GridImp,true>;
    friend class Dune::YaspIndexSet<GridImp,false>;
    friend class Dune::YaspGlobalIdSet<GridImp>;

    //! globally unique, persistent index
    PersistentIndexType persistentIndex () const
    {
      // get coordinate and size of global grid
      const iTupel& size =  _g->vertex_global.size();
      int coord[dim];

      // correction for periodic boundaries
      for (int i=0; i<dim; i++)
      {
        coord[i] = _it.coord(i);
        if (coord[i]<0)
          coord[i] += size[i];
        if (coord[i]>=size[i])
          coord[i] -= size[i];
      }

      // determine min number of trailing zeroes
      int trailing = 1000;
      for (int i=0; i<dim; i++)
      {
        // count trailing zeros
        int zeros = 0;
        for (int j=0; j<_g->level(); j++)
          if (coord[i]&(1<<j))
            break;
          else
            zeros++;
        trailing = std::min(trailing,zeros);
      }

      // determine the level of this vertex
      int level = _g->level()-trailing;

      // encode codim
      PersistentIndexType id(dim);

      // encode level
      id = id << yaspgrid_level_bits;
      id = id+PersistentIndexType(level);

      // encode coordinates
      for (int i=dim-1; i>=0; i--)
      {
        id = id << yaspgrid_dim_bits;
        id = id+PersistentIndexType(coord[i]>>trailing);
      }

      return id;
    }

    //! consecutive, codim-wise, level-wise index
    int compressedIndex () const { return _it.superindex();}

    //! consecutive, codim-wise, level-wise index
    int compressedLeafIndex () const
    {
      if (_g.level()==_g.mg()->maxlevel())
        return _it.superindex();

      // not on leaf level, interpolate to finest grid
      int coord[dim];
      for (int i=0; i<dim; i++) coord[i] = _it.coord(i)-(_g).vertex_overlap().origin(i);

      // move coordinates up to maxlevel (multiply by 2 for each level
      for (int k=0; k<dim; k++)
        coord[k] = coord[k]*(1<<(_g.mg()->maxlevel()-_g.level()));

      // do lexicographic numbering
      int index = coord[dim-1];
      for (int k=dim-2; k>=0; --k)
        index = (index*(_g.mg()->rbegin().cell_overlap().size(k)+1))+coord[k];
      return index;
    }

  public:
    const TSI& transformingsubiterator() const { return _it; }
    const YGLI& gridlevel() const { return _g; }
    const GridImp * yaspgrid() const { return _yg; }
  protected:
    const GridImp * _yg;          // access to YaspGrid
    const TSI& _it;               // position in the grid level
    const YGLI& _g;               // access to grid level
  };

}   // namespace Dune

#endif  // DUNE_GRID_YASPGRIDENTITY_HH
