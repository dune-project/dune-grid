// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ONE_D_GRID_ENTITY_HH
#define DUNE_ONE_D_GRID_ENTITY_HH

#include <array>

#include <dune/common/fvector.hh>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/entity.hh>

#include "nulliteratorfactory.hh"

/** \file
 * \brief The OneDGridEntity class and its specializations
 */

namespace Dune {

  // forward declarations
  template <class GridImp>
  class OneDGridLeafIntersectionIterator;
  template <class GridImp>
  class OneDGridLevelIntersectionIterator;
  template<int codim, PartitionIteratorType pitype, class GridImp>
  class OneDGridLevelIterator;
  template<class GridImp>
  class OneDGridHierarchicIterator;

  template <int mydim>
  class OneDEntityImp {};

  /** \brief Specialization for vertices */
  template <>
  class OneDEntityImp<0>
  {
  public:

    OneDEntityImp(int level, double pos)
      : pos_(pos), levelIndex_(0), leafIndex_(0), level_(level),
        son_(OneDGridNullIteratorFactory<0>::null()),
        pred_(OneDGridNullIteratorFactory<0>::null()),
        succ_(OneDGridNullIteratorFactory<0>::null())
    {}

    OneDEntityImp(int level, const FieldVector<double, 1>& pos, unsigned int id)
      : pos_(pos), levelIndex_(0), leafIndex_(0), id_(id), level_(level),
        son_(OneDGridNullIteratorFactory<0>::null()),
        pred_(OneDGridNullIteratorFactory<0>::null()),
        succ_(OneDGridNullIteratorFactory<0>::null())
    {}

    //private:
    bool isLeaf() const {
      return son_==OneDGridNullIteratorFactory<0>::null();
    }

    FieldVector<double, 1> pos_;

    //! entity number
    unsigned int levelIndex_;

    unsigned int leafIndex_;

    unsigned int id_;

    //! level
    int level_;

    //! Son vertex on the next finer grid
    OneDEntityImp<0>* son_;

    //!
    OneDEntityImp<0>* pred_;

    OneDEntityImp<0>* succ_;


  };


  /** \brief Specialization for elements */
  template <>
  class OneDEntityImp<1>
  {
  public:

    /** \brief The different ways to mark an element for grid changes */
    enum MarkState { DO_NOTHING , COARSEN , REFINE };

    OneDEntityImp(int level, unsigned int id, bool reversedBoundarySegmentNumbering)
      : levelIndex_(0), leafIndex_(0), id_(id), level_(level),
        markState_(DO_NOTHING), isNew_(false),
        reversedBoundarySegmentNumbering_(reversedBoundarySegmentNumbering),
        pred_(OneDGridNullIteratorFactory<1>::null()),
        succ_(OneDGridNullIteratorFactory<1>::null())
    {
      father_ = OneDGridNullIteratorFactory<1>::null();
      sons_[0] = sons_[1] = OneDGridNullIteratorFactory<1>::null();
    }

    bool isLeaf() const {
      assert( (sons_[0]==OneDGridNullIteratorFactory<1>::null() && sons_[1]==OneDGridNullIteratorFactory<1>::null())
              || (sons_[0]!=OneDGridNullIteratorFactory<1>::null() && sons_[1]!=OneDGridNullIteratorFactory<1>::null()) );
      return sons_[0]==OneDGridNullIteratorFactory<1>::null() && sons_[1]==OneDGridNullIteratorFactory<1>::null();
    }

    std::array<OneDEntityImp<1>*, 2> sons_;

    OneDEntityImp<1>* father_;

    OneDEntityImp<0>* vertex_[2];

    //! element number
    unsigned int levelIndex_;

    unsigned int leafIndex_;

    /** \brief Unique and persistent id for elements */
    unsigned int id_;

    //! the level of the entity
    int level_;

    /** \brief Stores requests for refinement and coarsening */
    MarkState markState_;

    /** \brief This flag is set by adapt() if this element has been newly created. */
    bool isNew_;

    /** Since a OneDGrid is one-dimensional and connected, there can only be two possible numberings
        of the boundary segments.  Either the left one is '0' and the right one is '1' or the reverse.
        This flag stores which is the case. It has the same value throughout the entire grid.
        This is a waste, because the Intersections class (which hands out the information), can
        access data here much easier than data in the central grid class.
     */
    bool reversedBoundarySegmentNumbering_;

    /** \brief Predecessor in the doubly linked list of elements */
    OneDEntityImp<1>* pred_;

    /** \brief Successor in the doubly linked list of elements */
    OneDEntityImp<1>* succ_;

  };

  // forward declarations
  template<class GridImp> class OneDGridLevelIndexSet;
  template<class GridImp> class OneDGridLeafIndexSet;
  template<class GridImp> class OneDGridIdSet;

  //**********************************************************************
  //
  // --OneDGridEntity
  // --Entity
  //
  /** \brief The implementation of entities in a OneDGrid
     \ingroup OneDGrid

     A Grid is a container of grid entities. An entity is parametrized by the codimension.
     An entity of codimension c in dimension d is a d-c dimensional object.

   */
  template<int cd, int dim, class GridImp>
  class OneDGridEntity :
    public EntityDefaultImplementation <cd,dim,GridImp,OneDGridEntity>
  {

    template <int codim_, PartitionIteratorType PiType_, class GridImp_>
    friend class OneDGridLevelIterator;

    friend class OneDGrid;

    // IndexSets and IdSets need to access indices and ids
    friend class OneDGridLevelIndexSet<GridImp>;
    friend class OneDGridLeafIndexSet<GridImp>;
    friend class OneDGridIdSet<GridImp>;

    typedef typename GridImp::Traits::template Codim< cd >::GeometryImpl GeometryImpl;

  public:
    /** \brief The type of OneDGrid Entity seeds */
    typedef typename GridImp::Traits::template Codim<cd>::EntitySeed EntitySeed;

    //! Default constructor
    OneDGridEntity()
      : target_(OneDGridNullIteratorFactory<0>::null())
    {}

    explicit OneDGridEntity(OneDEntityImp<0>* target)
      : target_(target)
    {}

    typedef typename GridImp::template Codim<cd>::Geometry Geometry;

    bool equals(const OneDGridEntity& other) const
    {
      return target_ == other.target_;
    }

    //! level of this element
    int level () const {return target_->level_;}

    //! only interior entities
    PartitionType partitionType () const { return InteriorEntity; }

  private:
    unsigned int levelIndex() const {return target_->levelIndex_;}

    unsigned int leafIndex() const {return target_->leafIndex_;}

    unsigned int globalId() const {return target_->id_;}

    /** \brief Returns vertex level index, the arguments are ignored
     */
    int subLevelIndex ([[maybe_unused]] int i, [[maybe_unused]] unsigned int codim) const
    {
      return target_->levelIndex_;
    }

    /** \brief Returns vertex leaf index, the arguments are ignored
     */
    int subLeafIndex ([[maybe_unused]] int i, [[maybe_unused]] unsigned int codim) const
    {
      return target_->leafIndex_;
    }

  public:
    //! return the element type identifier (segment)
    GeometryType type () const {return GeometryTypes::vertex;}

    /** \brief Return the number of subentities of codimension codim.
     */
    unsigned int subEntities (unsigned int codim) const
    {
      assert(codim==1);
      return 1;
    }

    //! geometry of this entity
    Geometry geometry () const { return Geometry(GeometryImpl(target_->pos_)); }

    /** \brief Get the seed corresponding to this entity */
    EntitySeed seed () const { return EntitySeed( *this ); }

    void setToTarget(OneDEntityImp<0>* target) {
      target_ = target;
    }

    OneDEntityImp<0>* target_;

  };

  //***********************
  //
  //  --OneDGridEntity
  //  --0Entity
  //
  //***********************



  /** \brief Specialization for codim-0-entities.
   * \ingroup OneDGrid
   *
   * This class embodies the topological parts of elements of the grid.
   * It has an extended interface compared to the general entity class.
   * For example, Entities of codimension 0  allow to visit all neighbors.
   *
   * OneDGrid only implements the case dim==dimworld==1
   */
  template<int dim, class GridImp>
  class OneDGridEntity<0,dim, GridImp> :
    public EntityDefaultImplementation<0,dim,GridImp, OneDGridEntity>
  {
    friend class OneDGrid;
    template <class GridImp_>
    friend class OneDGridLevelIntersectionIterator;
    template <class GridImp_>
    friend class OneDGridLeafIntersectionIterator;
    friend class OneDGridHierarchicIterator <GridImp>;
    friend class OneDGridLevelIterator <0,All_Partition,GridImp>;

    // IndexSets and IdSets need to access indices and ids
    friend class OneDGridLevelIndexSet<GridImp>;
    friend class OneDGridLeafIndexSet<GridImp>;
    friend class OneDGridIdSet<GridImp>;

    typedef typename GridImp::Traits::template Codim< 0 >::GeometryImpl GeometryImpl;
    typedef typename GridImp::Traits::template Codim< 0 >::GeometryImpl LocalGeometryImpl;

  public:
    typedef typename GridImp::template Codim<0>::Geometry Geometry;
    typedef typename GridImp::template Codim<0>::LocalGeometry LocalGeometry;
    typedef typename GridImp::template Codim<0>::LevelIterator LevelIterator;
    typedef typename GridImp::LeafIntersectionIterator LeafIntersectionIterator;
    typedef typename GridImp::LevelIntersectionIterator LevelIntersectionIterator;
    typedef typename GridImp::HierarchicIterator HierarchicIterator;

    /** \brief The type of OneDGrid Entity seeds */
    typedef typename GridImp::Traits::template Codim<0>::EntitySeed EntitySeed;

    template<int codim>
    struct Codim
    {
      typedef typename GridImp::Traits::template Codim<codim>::Entity Entity;
    };

    typedef typename GridImp::Traits::template Codim<0>::Entity Entity;

    //! Default Constructor
    OneDGridEntity ()
      : target_( OneDGridNullIteratorFactory<1>::null() )
    {}

    explicit OneDGridEntity (OneDEntityImp<1>* target)
      : target_( target )
    {}

    bool equals(const OneDGridEntity& other) const
    {
      return target_ == other.target_;
    }

    //! Level of this element
    int level () const {return target_->level_;}

    //! only interior entities
    PartitionType partitionType () const { return InteriorEntity; }

  private:
    //! Level index is unique and consecutive per level and codim
    unsigned int levelIndex() const {return target_->levelIndex_;}

    unsigned int leafIndex() const {return target_->leafIndex_;}

    unsigned int globalId() const {return target_->id_;}

  public:
    //! Geometry of this entity
    Geometry geometry () const { return Geometry( GeometryImpl(target_->vertex_[0]->pos_, target_->vertex_[1]->pos_) ); }

    /** \brief Get the seed corresponding to this entity */
    EntitySeed seed () const { return EntitySeed( *this ); }

    //! return the element type identifier (segment)
    GeometryType type () const {return GeometryTypes::line;}

    /** \brief Return the number of subentities of codimension codim.
     */
    unsigned int subEntities (unsigned int codim) const
    {
      assert(codim==0 || codim==1);
      return (codim==0) ? 1 : 2;
    }

  private:
    /** \brief Return index of sub entity with codim = cc and local number i
     */
    int subLevelIndex (int i,unsigned int codim) const {
      assert(i==0 || i==1);
      return (codim==0)
             ? target_->levelIndex_
             : target_->vertex_[i]->levelIndex_;
    }

    /** \brief Return leaf index of sub entity with codim = cc and local number i
     */
    int subLeafIndex (int i, unsigned int codim) const {
      assert(i==0 || i==1);
      return (codim==0)
             ? target_->leafIndex_
             : target_->vertex_[i]->leafIndex_;
    }

    /** \brief Return leaf index of sub entity with codim = cc and local number i
     */
    int subId (int i, unsigned int codim) const {
      assert(i==0 || i==1);
      return (codim==0)
             ? target_->id_
             : target_->vertex_[i]->id_;
    }

  public:
    /** \brief Access to codim 0 subentities */
    template<int cc>
    typename std::enable_if<
      cc == 0,
      typename Codim<0>::Entity
      >::type
    subEntity (int i) const
    {
      assert(i==0);
      return typename Codim<0>::Entity(OneDGridEntity<0,dim,GridImp>(this->target_));
    }

    /** \brief Access to codim 1 subentities */
    template<int cc>
    typename std::enable_if<
      cc == 1,
      typename Codim<1>::Entity
      >::type
    subEntity (int i) const
    {
      assert(i==0 || i==1);
      return typename Codim<1>::Entity(OneDGridEntity<1,dim,GridImp>(this->target_->vertex_[i]));
    }

    LeafIntersectionIterator ileafbegin () const {
      return OneDGridLeafIntersectionIterator<GridImp>(target_, (isLeaf()) ? 0 : 2);
    }

    LevelIntersectionIterator ilevelbegin () const {
      return OneDGridLevelIntersectionIterator<GridImp>(target_, 0);
    }

    LeafIntersectionIterator ileafend () const {
      return OneDGridLeafIntersectionIterator<GridImp>(target_);
    }

    LevelIntersectionIterator ilevelend () const {
      return OneDGridLevelIntersectionIterator<GridImp>(target_);
    }

    //! returns true if Entity has no children
    bool isLeaf () const {
      return (target_->sons_[0]==OneDGridNullIteratorFactory<1>::null())
             && (target_->sons_[1]==OneDGridNullIteratorFactory<1>::null());
    }

    //! Inter-level access to father element on coarser grid.
    //! Assumes that meshes are nested.
    Entity father () const {
      return Entity(OneDGridEntity(target_->father_));
    }
    //! returns true if father entity exists
    bool hasFather () const
    {
      return (this->level()>0);
    }

    /*! Location of this element relative to the reference element
       of the father. This is sufficient to interpolate all
       dofs in conforming case. Nonconforming may require access to
       neighbors of father and computations with local coordinates.
       On the fly case is somewhat inefficient since dofs  are visited
       several times. If we store interpolation matrices, this is tolerable.
       We assume that on-the-fly implementation of numerical algorithms
       is only done for simple discretizations. Assumes that meshes are nested.
     */
    LocalGeometry geometryInFather () const {
      assert(target_->father_);
      assert(target_->father_->sons_[0] == target_ || target_->father_->sons_[1] == target_);

      if (target_->father_->sons_[0] == target_ && target_->father_->sons_[1] == target_) {
        // Copied element?
        return LocalGeometry(LocalGeometryImpl(FieldVector<double,1>(0), FieldVector<double,1>(1)));
      } else if (target_->father_->sons_[0] == target_) {
        // Left son?
        return LocalGeometry(LocalGeometryImpl(FieldVector<double,1>(0), FieldVector<double,1>(0.5)));
      }

      // Right son!
      return LocalGeometry(LocalGeometryImpl(FieldVector<double,1>(0.5), FieldVector<double,1>(1)));
    }

    /*! Inter-level access to son elements on higher levels<=maxlevel.
       This is provided for sparsely stored nested unstructured meshes.
       Returns iterator to first son.
     */
    OneDGridHierarchicIterator<GridImp> hbegin (int maxlevel) const {

      OneDGridHierarchicIterator<GridImp> it(maxlevel);

      if (level()<=maxlevel) {

        // Load sons of old target onto the iterator stack
        if (!isLeaf()) {

          it.elemStack.push(target_->sons_[0]);
          it.elemStack.push(target_->sons_[1]);

        }

      }

      it.virtualEntity_.impl().setToTarget((it.elemStack.empty())
                                                                    ? OneDGridNullIteratorFactory<1>::null() : it.elemStack.top());

      return it;
    }

    //! Returns iterator to one past the last son
    HierarchicIterator hend (int maxlevel) const {
      return HierarchicIterator(maxlevel);
    }

    // ***************************************************************
    //  Interface for Adaptation
    // ***************************************************************

    /** returns true, if entity might be coarsened during next adaptation cycle */
    bool mightVanish () const { return target_->markState_ == OneDEntityImp<1> :: COARSEN; }

    /** returns true, if entity was refined during last adaptation cycle */
    bool isNew () const { return target_->isNew_; }

    void setToTarget(OneDEntityImp<1>* target) {
      target_ = target;
    }

    OneDEntityImp<1>* target_;

  }; // end of OneDGridEntity codim = 0

} // namespace Dune

#endif
