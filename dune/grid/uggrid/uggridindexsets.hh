// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_UGGRID_INDEXSETS_HH
#define DUNE_UGGRID_INDEXSETS_HH

/** \file
    \brief The index and id sets for the UGGrid class
 */

#include <vector>
#include <set>

#include <dune/common/std/type_traits.hh>
#include <dune/grid/common/grid.hh>

namespace Dune {

  template<class GridImp>
  class UGGridLevelIndexSet : public IndexSet<GridImp,UGGridLevelIndexSet<GridImp>, UG::UINT>
  {
    constexpr static int dim = GridImp::dimension;

  public:

    /** \brief Default constructor

       Unfortunately we can't force the user to init grid_ and level_, because
       UGGridLevelIndexSets are meant to be stored in an array.

       \todo I want to make this constructor private, but I can't, because
       it is called by UGGrid through a std::vector::resize()
     */
    UGGridLevelIndexSet ()
      : level_(0),
        numSimplices_(0),
        numPyramids_(0),
        numPrisms_(0),
        numCubes_(0),
        numVertices_(0),
        numEdges_(0),
        numTriFaces_(0),
        numQuadFaces_(0)
    {}

    //! get index of an entity
    template<int cd>
    unsigned int index (const typename GridImp::Traits::template Codim<cd>::Entity& e) const
    {
      return UG_NS<dim>::levelIndex(e.impl().getTarget());
    }

    /** \brief Get index of subEntity of a codim cc entity
     * \param codim Codimension WITH RESPECT TO THE GRID of the subentity whose index we want
     */
    template<int cc>
    unsigned int subIndex (const typename GridImp::Traits::template Codim<cc>::Entity& entity,
                           int i,
                           unsigned int codim) const
    {
      // The entity is a vertex, so each subentity must be a vertex too (anything else is not supported)
      if (cc==dim)
      {
        assert(codim==dim);
        return UG_NS<dim>::levelIndex(entity.impl().getTarget());
      }

      // The following block is for 2d grids
      if constexpr (dim==2)
      {
        // The entity is an element
        if constexpr (cc==0)
        {
          // Element indices
          if (codim==0)
            return UG_NS<dim>::levelIndex(entity.impl().getTarget());

          // Edge indices
          if (codim==1)
          {
            auto ref_el = referenceElement<double,dim>(entity.type());
            auto a = ref_el.subEntity(i,dim-1,0,dim);
            auto b = ref_el.subEntity(i,dim-1,1,dim);
            return UG_NS<dim>::levelIndex(UG_NS<dim>::GetEdge(UG_NS<dim>::Corner(entity.impl().getTarget(),
                                                                                 UGGridRenumberer<dim>::verticesDUNEtoUG(a,entity.type())),
                                                              UG_NS<dim>::Corner(entity.impl().getTarget(),
                                                                                 UGGridRenumberer<dim>::verticesDUNEtoUG(b,entity.type()))));
          }

          // Vertex indices
          if (codim==dim)
            return UG_NS<dim>::levelIndex(UG_NS<dim>::Corner(entity.impl().getTarget(),
                                                             UGGridRenumberer<dim>::verticesDUNEtoUG(i,entity.type())));
        }

        // The entity is an edge
        if constexpr (cc==1)
          return UG_NS<dim>::levelIndex(entity.impl().getTarget()->links[i].nbnode);
      }

      // The following block is for 3d grids
      if constexpr (dim==3)
      {
        // The entity is an element
        if constexpr (cc==0)
        {
          // Element indices
          if (codim==0)
            return UG_NS<dim>::levelIndex(entity.impl().getTarget());

          // Face indices
          if (codim==1)
            return UG_NS<dim>::levelIndex(UG_NS<dim>::SideVector(entity.impl().getTarget(),
                                                                 UGGridRenumberer<dim>::facesDUNEtoUG(i,entity.type())));

          // Edge indices
          if (codim==2)
          {
            auto ref_el = referenceElement<double,dim>(entity.type());
            auto a = ref_el.subEntity(i,dim-1,0,dim);
            auto b = ref_el.subEntity(i,dim-1,1,dim);
            return UG_NS<dim>::levelIndex(UG_NS<dim>::GetEdge(UG_NS<dim>::Corner(entity.impl().getTarget(),
                                                                                 UGGridRenumberer<dim>::verticesDUNEtoUG(a,entity.type())),
                                                              UG_NS<dim>::Corner(entity.impl().getTarget(),
                                                                                 UGGridRenumberer<dim>::verticesDUNEtoUG(b,entity.type()))));
          }

          // Vertex indices
          if (codim==3)
            return UG_NS<dim>::levelIndex(UG_NS<dim>::Corner(entity.impl().getTarget(),
                                                               UGGridRenumberer<dim>::verticesDUNEtoUG(i,entity.type())));
        }

        // The entity is a face
        if constexpr (cc==1)
        {
          DUNE_THROW(NotImplemented, "Subindices of an element face");
        }

        // The entity is an edge
        if constexpr (cc==2)
          return UG_NS<dim>::levelIndex(entity.impl().getTarget()->links[i].nbnode);
      }

      // Should never happen
      return std::numeric_limits<unsigned int>::max();
    }

    //! get number of entities of given codim, type and on this level
    std::size_t size (int codim) const {
      if (codim==0)
        return numSimplices_+numPyramids_+numPrisms_+numCubes_;
      if (codim==dim)
        return numVertices_;
      if (codim==dim-1)
        return numEdges_;
      if (codim==1)
        return numTriFaces_+numQuadFaces_;
      DUNE_THROW(NotImplemented, "wrong codim!");
    }

    //! get number of entities of given codim, type and on this level
    std::size_t size (GeometryType type) const
    {
      int codim = GridImp::dimension-type.dim();

      if (codim==0) {
        if (type.isSimplex())
          return numSimplices_;
        else if (type.isPyramid())
          return numPyramids_;
        else if (type.isPrism())
          return numPrisms_;
        else if (type.isCube())
          return numCubes_;
        else
          return 0;

      }

      if (codim==dim) {
        return numVertices_;
      }
      if (codim==dim-1) {
        return numEdges_;
      }
      if (codim==1) {
        if (type.isSimplex())
          return numTriFaces_;
        else if (type.isCube())
          return numQuadFaces_;
        else
          return 0;
      }

      DUNE_THROW(NotImplemented, "Wrong codim!");
    }

    std::vector< GeometryType > types ( int codim ) const { return myTypes_[ codim ]; }

    /** \brief Deliver all geometry types used in this grid */
    const std::vector<GeometryType>& geomTypes (int codim) const
    {
      return myTypes_[codim];
    }

    /** \brief Return true if e is contained in the index set.

        \note If 'entity' is from another grid this method may still return 'true'.
        This is acceptable by the Dune grid interface specification.
     */
    template <class EntityType>
    bool contains (const EntityType& entity) const
    {
      return entity.level() == level_;
    }

    /** \brief Update the level indices.  This method is called after each grid change */
    void update(const GridImp& grid, int level, std::vector<unsigned int>* nodePermutation=0);

    const GridImp* grid_;
    int level_;

    int numSimplices_;
    int numPyramids_;
    int numPrisms_;
    int numCubes_;
    int numVertices_;
    int numEdges_;
    int numTriFaces_;
    int numQuadFaces_;

    std::vector<GeometryType> myTypes_[dim+1];
  };

  template<class GridImp>
  class UGGridLeafIndexSet : public IndexSet<GridImp,UGGridLeafIndexSet<GridImp>, UG::UINT>
  {
  public:

    /*
       We use the remove_const to extract the Type from the mutable class,
       because the const class is not instantiated yet.
     */
    constexpr static int dim = std::remove_const<GridImp>::type::dimension;

    //! constructor stores reference to a grid and level
    UGGridLeafIndexSet (const GridImp& g)
      : grid_(g), coarsestLevelWithLeafElements_(0)
    {}

    //! get index of an entity
    /*
       We use the remove_const to extract the Type from the mutable class,
       because the const class is not instantiated yet.
     */
    template<int cd>
    int index (const typename std::remove_const<GridImp>::type::Traits::template Codim<cd>::Entity& e) const
    {
      return UG_NS<dim>::leafIndex(e.impl().getTarget());
    }

    //! get index of subEntity of a codim 0 entity
    /*
       We use the remove_const to extract the Type from the mutable class,
       because the const class is not instantiated yet.
     */
    template<int cc>
    unsigned int subIndex (const typename std::remove_const<GridImp>::type::Traits::template Codim<cc>::Entity& entity,
                           int i,
                           unsigned int codim) const
    {
      // The entity is a vertex, so each subentity must be a vertex too (anything else is not supported)
      if (cc==dim)
      {
        assert(codim==dim);
        return UG_NS<dim>::leafIndex(entity.impl().getTarget());
      }

      // The following block is for 2d grids
      if constexpr (dim==2)
      {
        // The entity is an element
        if constexpr (cc==0)
        {
          // Element indices
          if (codim==0)
            return UG_NS<dim>::leafIndex(entity.impl().getTarget());

          // Edge indices
          if (codim==1)
          {
            auto ref_el = referenceElement<double,dim>(entity.type());
            auto a = ref_el.subEntity(i,dim-1,0,dim);
            auto b = ref_el.subEntity(i,dim-1,1,dim);
            return UG_NS<dim>::leafIndex(UG_NS<dim>::GetEdge(UG_NS<dim>::Corner(entity.impl().getTarget(),
                                                                                 UGGridRenumberer<dim>::verticesDUNEtoUG(a,entity.type())),
                                                              UG_NS<dim>::Corner(entity.impl().getTarget(),
                                                                                 UGGridRenumberer<dim>::verticesDUNEtoUG(b,entity.type()))));
          }

          // Vertex indices
          if (codim==dim)
            return UG_NS<dim>::leafIndex(UG_NS<dim>::Corner(entity.impl().getTarget(),
                                                             UGGridRenumberer<dim>::verticesDUNEtoUG(i,entity.type())));
        }

        // The entity is an edge
        if constexpr (cc==1)
          return UG_NS<dim>::leafIndex(entity.impl().getTarget()->links[i].nbnode);
      }

      // The following block is for 3d grids
      if constexpr (dim==3)
      {
        // The entity is an element
        if constexpr (cc==0)
        {
          // Element indices
          if (codim==0)
            return UG_NS<dim>::leafIndex(entity.impl().getTarget());

          // Face indices
          if (codim==1)
            return UG_NS<dim>::leafIndex(UG_NS<dim>::SideVector(entity.impl().getTarget(),
                                                                 UGGridRenumberer<dim>::facesDUNEtoUG(i,entity.type())));

          // Edge indices
          if (codim==2)
          {
            auto ref_el = referenceElement<double,dim>(entity.type());
            auto a = ref_el.subEntity(i,dim-1,0,dim);
            auto b = ref_el.subEntity(i,dim-1,1,dim);
            return UG_NS<dim>::leafIndex(UG_NS<dim>::GetEdge(UG_NS<dim>::Corner(entity.impl().getTarget(),
                                                                                 UGGridRenumberer<dim>::verticesDUNEtoUG(a,entity.type())),
                                                              UG_NS<dim>::Corner(entity.impl().getTarget(),
                                                                                 UGGridRenumberer<dim>::verticesDUNEtoUG(b,entity.type()))));
          }

          // Vertex indices
          if (codim==3)
            return UG_NS<dim>::leafIndex(UG_NS<dim>::Corner(entity.impl().getTarget(),
                                                              UGGridRenumberer<dim>::verticesDUNEtoUG(i,entity.type())));
        }

        // The entity is a face
        if constexpr (cc==1)
        {
          DUNE_THROW(NotImplemented, "Subindices of an element face");
        }

        // The entity is an edge
        if constexpr (cc==2)
          return UG_NS<dim>::leafIndex(entity.impl().getTarget()->links[i].nbnode);
      }

      return std::numeric_limits<unsigned int>::max();
    }

    //! get number of entities of given codim and type
    std::size_t size (GeometryType type) const
    {
      if (type.dim()==GridImp::dimension) {
        if (type.isSimplex())
          return numSimplices_;
        else if (type.isPyramid())
          return numPyramids_;
        else if (type.isPrism())
          return numPrisms_;
        else if (type.isCube())
          return numCubes_;
        else
          return 0;
      }
      if (type.dim()==0) {
        return numVertices_;
      }
      if (type.dim()==1) {
        return numEdges_;
      }
      if (type.isTriangle())
        return numTriFaces_;
      else if (type.isQuadrilateral())
        return numQuadFaces_;

      return 0;
    }

    //! get number of entities of given codim
    std::size_t size (int codim) const
    {
      int s=0;
      const std::vector<GeometryType>& geomTs = geomTypes(codim);
      for (unsigned int i=0; i<geomTs.size(); ++i)
        s += size(geomTs[i]);
      return s;
    }

    std::vector< GeometryType > types ( int codim ) const { return myTypes_[ codim ]; }

    /** deliver all geometry types used in this grid */
    const std::vector<GeometryType>& geomTypes (int codim) const
    {
      return myTypes_[codim];
    }

    /** \brief Return true if e is contained in the index set.

        \note If 'entity' is from another grid this method may still return 'true'.
        This is acceptable by the Dune grid interface specification.
     */
    template <class EntityType>
    bool contains (const EntityType& entity) const
    {
      return UG_NS<dim>::isLeaf(entity.impl().getTarget());
    }


    /** \brief Update the leaf indices.  This method is called after each grid change. */
    void update(std::vector<unsigned int>* nodePermutation=0);

    const GridImp& grid_;

    /** \brief The lowest level that contains leaf elements

       This corresponds to UG's fullRefineLevel, which is, unfortunately only
       computed if you use some nontrivial UG algebra.  Thus we compute it
       ourselves, and use it to speed up the leaf iterators.
     */
    unsigned int coarsestLevelWithLeafElements_;



    int numSimplices_;
    int numPyramids_;
    int numPrisms_;
    int numCubes_;
    int numVertices_;
    int numEdges_;
    int numTriFaces_;
    int numQuadFaces_;

    std::vector<GeometryType> myTypes_[dim+1];
  };


  /** \brief Implementation class for the UGGrid Id sets

     The UGGridGlobalIdSet and the UGGridLocalIdSet are identical. This
     class implements them both at once.
   */
  template <class GridImp>
  class UGGridIdSet : public IdSet<GridImp,UGGridIdSet<GridImp>,typename UG_NS<std::remove_const<GridImp>::type::dimension>::UG_ID_TYPE>
  {
    constexpr static int dim = std::remove_const<GridImp>::type::dimension;

    typedef typename std::pair<const typename UG_NS<dim>::Element*, int> Face;

    /** \brief Look for copy of a face on the next-lower grid level.

       \todo This method is not implemented very efficiently
     */
    static Face getFatherFace(const Face& face) {

      // set up result object
      Face resultFace;
      resultFace.first = UG_NS<dim>::EFather(face.first);

      // If there is no father element then we know there is no father face
      /** \bug This is not true when doing vertical load balancing. */
      if (resultFace.first == nullptr)
        return resultFace;

      // Get all corners of the face
      std::set<const typename UG_NS<dim>::Vertex*> corners;

      for (int i=0; i<UG_NS<dim>::Corners_Of_Side(face.first, face.second); i++)
        corners.insert(UG_NS<dim>::Corner(face.first, UG_NS<dim>::Corner_Of_Side(face.first, face.second, i))->myvertex);

      // Loop over all faces of the father element and look for a face that has the same vertices
      for (int i=0; i<UG_NS<dim>::Sides_Of_Elem(resultFace.first); i++) {

        // Copy father face into set
        std::set<const typename UG_NS<dim>::Vertex*> fatherFaceCorners;

        for (int j=0; j<UG_NS<dim>::Corners_Of_Side(resultFace.first, i); j++)
          fatherFaceCorners.insert(UG_NS<dim>::Corner(resultFace.first, UG_NS<dim>::Corner_Of_Side(resultFace.first, i, j))->myvertex);

        // Do the vertex sets match?
        if (corners==fatherFaceCorners) {
          resultFace.second = i;
          return resultFace;
        }

      }

      // No father face found
      resultFace.first = nullptr;
      return resultFace;
    }

  public:
    //! constructor stores reference to a grid
    UGGridIdSet (const GridImp& g) : grid_(g) {}

    /** \brief Get id of an entity

       We use the remove_const to extract the Type from the mutable class,
       because the const class is not instantiated yet.

       \bug Since copies of different entities on different levels are supposed to have the
       same id, we look for the ancestor on the coarsest level that is still a copy of
       the entity we are interested in.  However, the current implementation only searches
       on one processor, while with UG's vertical load balancing the ancestors of an entity
       may be distributed across different processors.  This will lead to very-difficult-to-fix
       bugs.  Unfortunately, the proper fix for this is not easy, either.
     */
    template<int cd>
    typename UG_NS<dim>::UG_ID_TYPE id (const typename std::remove_const<GridImp>::type::Traits::template Codim<cd>::Entity& e) const
    {
      if (cd==0) {
        // If we're asked for the id of an element, and that element is a copy of its father, then
        // we return the id of the lowest ancestor that the element is a copy from.  That way copies
        // of elements have the same id
        const typename UG_NS<dim>::Element* ancestor = (typename UG_NS<dim>::Element* const)(e.impl().getTarget());
        /** \todo We should really be using an isCopy() method rather than hasCopy() */
        while (UG_NS<dim>::EFather(ancestor) && UG_NS<dim>::hasCopy(UG_NS<dim>::EFather(ancestor)))
          ancestor = UG_NS<dim>::EFather(ancestor);

        return UG_NS<dim>::id(ancestor);
      }

      if (dim-cd==1) {

        const typename UG_NS<dim>::Edge* edge = (typename UG_NS<dim>::Edge* const)(e.impl().getTarget());

        // If this edge is the copy of an edge on a lower level we return the id of that lower
        // edge, because Dune wants entities which are copies of each other to have the same id.
        // BUG: in the parallel setting, we only search on our own processor, but the lowest
        // copy may actually be on a different processor!
        const typename UG_NS<dim>::Edge* fatherEdge;
        fatherEdge = GetFatherEdge(edge);

        while (fatherEdge   // fatherEdge exists
               // ... and it must be a true copy father
               && ( (fatherEdge->links[0].nbnode->myvertex == edge->links[0].nbnode->myvertex
                     && fatherEdge->links[1].nbnode->myvertex == edge->links[1].nbnode->myvertex)
                    ||
                    (fatherEdge->links[0].nbnode->myvertex == edge->links[1].nbnode->myvertex
                     && fatherEdge->links[1].nbnode->myvertex == edge->links[0].nbnode->myvertex) ) ) {
          edge = fatherEdge;
          fatherEdge = GetFatherEdge(edge);
        }

        return UG_NS<dim>::id(edge);
      }


      if (cd == dim) {
        typename UG_NS<dim>::Node *node =
          reinterpret_cast<typename UG_NS<dim>::Node *>(e.impl().getTarget());

        return UG_NS<dim>::id(node);
      }

      // The entity must be a facet in a 3d grid
      assert(cd==1 && dim==3);

      typename UG_NS<dim>::Vector *facet =
          reinterpret_cast<typename UG_NS<dim>::Vector *>(e.impl().getTarget());

      // If 'facet' is the copy of a facet on a lower level we return the id of that lower-level
      // facet, because Dune wants entities that are copies of each other to have the same id.
      // BUG: in a distributed setting, we only search on our own process, but the lowest
      // copy may actually be on a different process!

      // Going up and down the refinement hierarchy can only be done through an
      // element that has 'facet' as a face.
      // TODO: Simplify this if ever SideVector objects get hierarchical information
      typename UG_NS<dim>::Element* supportingElement = (typename UG_NS<dim>::Element*)facet->object;
      int thisSideVectorUGNumbering = -1;   // Dummy value
      for (int side=0; side<UG_NS<dim>::Sides_Of_Elem(supportingElement); side++)
      {
        auto sideVector = UG_NS<dim>::SideVector(supportingElement,side);
        if (sideVector==facet)
          thisSideVectorUGNumbering = side;
      }

      Face face(supportingElement, thisSideVectorUGNumbering);

      // Go down in the grid hierarchy as long as there is a copy father-facet
      Face fatherFace;
      fatherFace = getFatherFace(face);
      while (fatherFace.first) {
        face = fatherFace;
        fatherFace = getFatherFace(face);
      }

      // Actually return the facet id
      return UG_NS<dim>::id(UG_NS<dim>::SideVector(face.first, face.second));
    }

    /** \brief Get id of subentity

       We use the remove_const to extract the Type from the mutable class,
       because the const class is not instantiated yet.

       \bug Since copies of different entities on different levels are supposed to have the
       same id, we look for the ancestor on the coarsest level that is still a copy of
       the entity we are interested in.  However, the current implementation only searches
       on one processor, while with UG's vertical load balancing the ancestors of an entity
       may be distributed across different processors.  This will lead to very-difficult-to-fix
       bugs.  Unfortunately, the proper fix for this is not easy, either.
     */
    typename UG_NS<dim>::UG_ID_TYPE subId (const typename std::remove_const<GridImp>::type::Traits::template Codim<0>::Entity& e,
                                           int i,
                                           unsigned int codim) const
    {
      if (codim==0)
        return id<0>(e);

      const typename UG_NS<dim>::Element* target = e.impl().getTarget();
      GeometryType type = e.type();

      if (dim-codim==1) {

        auto ref_el = referenceElement<double,dim>(type);
        auto a = ref_el.subEntity(i,dim-1,0,dim);
        auto b = ref_el.subEntity(i,dim-1,1,dim);
        const typename UG_NS<dim>::Edge* edge = UG_NS<dim>::GetEdge(UG_NS<dim>::Corner(target, UGGridRenumberer<dim>::verticesDUNEtoUG(a,type)),
                                                                    UG_NS<dim>::Corner(target, UGGridRenumberer<dim>::verticesDUNEtoUG(b,type)));

        // If this edge is the copy of an edge on a lower level we return the id of that lower
        // edge, because Dune wants entities which are copies of each other to have the same id.
        // BUG: in the parallel setting, we only search on our own processor, but the lowest
        // copy may actually be on a different processor!
        const typename UG_NS<dim>::Edge* fatherEdge;
        fatherEdge = GetFatherEdge(edge);

        while (fatherEdge   // fatherEdge exists
               // ... and it must be a true copy father
               && ( (fatherEdge->links[0].nbnode->myvertex == edge->links[0].nbnode->myvertex
                     && fatherEdge->links[1].nbnode->myvertex == edge->links[1].nbnode->myvertex)
                    ||
                    (fatherEdge->links[0].nbnode->myvertex == edge->links[1].nbnode->myvertex
                     && fatherEdge->links[1].nbnode->myvertex == edge->links[0].nbnode->myvertex) ) ) {
          edge = fatherEdge;
          fatherEdge = GetFatherEdge(edge);
        }

        return UG_NS<dim>::id(edge);
      }

      if (codim==1) {  // Faces

        Face face(target, UGGridRenumberer<dim>::facesDUNEtoUG(i,type));

        // If this face is the copy of a face on a lower level we return the id of that lower
        // face, because Dune wants entities which are copies of each other to have the same id.
        // BUG: in the parallel setting, we only search on our own processor, but the lowest
        // copy may actually be on a different processor!
        Face fatherFace;
        fatherFace = getFatherFace(face);
        while (fatherFace.first) {
          face = fatherFace;
          fatherFace = getFatherFace(face);
        }

        return UG_NS<dim>::id(UG_NS<dim>::SideVector(face.first, face.second));
      }

      if (codim==dim) {
        return UG_NS<dim>::id(UG_NS<dim>::Corner(target, UGGridRenumberer<dim>::verticesDUNEtoUG(i,type)));
      }

      DUNE_THROW(GridError, "UGGrid<" << dim << ">::subId isn't implemented for codim==" << codim );
    }

    //private:

    const GridImp& grid_;
  };

}  // namespace Dune

#endif
