// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_COMMON_REFINEMENT_HCUBETRIANGULATION_CC
#define DUNE_GRID_COMMON_REFINEMENT_HCUBETRIANGULATION_CC

// This file is part of DUNE, a Distributed and Unified Numerics Environment
// This file is copyright (C) 2005 Jorrit Fahlke <jorrit@jorrit.de>
// This file is licensed under version 2 of the GNU General Public License,
// with a special "runtime exception."  See COPYING at the top of the source
// tree for the full licence.

/*! @file

   @brief This file contains the @ref Refinement implementation for
         triangulating hypercubes (quadrilateral -> triangle,
         hexahedron -> tetrahedron)

   See @ref HCubeTriangulation.
   @verbatim
   $Id$
   @endverbatim
 */

/*! @defgroup HCubeTriangulation Refinement implementation for triangulating hypercubes
    @ingroup Refinement

   Most stuff here is explained in @ref SimplexRefinement.

   We simply triangulate the hypercube into Kuhn simplices, and than
   use the @ref SimplexRefinement to do the refinement.

   We explicitly use some of the utilities from the @ref SimplexRefinement.

 */

#include <dune/common/fvector.hh>
#include <dune/common/misc.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>
#include "base.cc"
#include "simplex.cc"

namespace Dune {

  namespace RefinementImp {

    /*! @brief This namespace contains the @ref Refinement
               implementation for triangulating hypercubes
               (GeometryType::cube -> GeometryType::simplex)

       See @ref HCubeTriangulation.
     */
    namespace HCubeTriangulation {

      // ////////////
      //
      //  Utilities
      //

      using Simplex::factorial;
      using Simplex::binomial;
      using Simplex::getPermutation;
      using Simplex::referenceToKuhn;
      using Simplex::kuhnToReference;

      // ////////////////////////////////////
      //
      //  Refine a hypercube with simplices
      //

      template<int mydimension, int coorddimension, class GridImp>
      class Geometry;

      // forward declaration of the iterator base
      template<int dimension, class CoordType, int codimension>
      class RefinementIteratorSpecial;

      template<int dimension_, class CoordType>
      class RefinementImp
      {
        friend class Geometry<dimension_, dimension_, RefinementImp>;

      public:
        enum { dimension = dimension_ };
        // to make Dune::Geometry work:
        struct GridFamily;
        typedef CoordType ctype;
        enum { dimensionworld = dimension };

        template<int codimension>
        struct Codim;
        typedef typename Codim<dimension>::SubEntityIterator VertexIterator;
        typedef FieldVector<CoordType, dimension> CoordVector;
        typedef typename Codim<0>::SubEntityIterator ElementIterator;
        typedef FieldVector<int, dimension+1> IndexVector;

        static int nVertices(int level);
        static VertexIterator vBegin(int level);
        static VertexIterator vEnd(int level);

        static int nElements(int level);
        static ElementIterator eBegin(int level);
        static ElementIterator eEnd(int level);
      private:
        friend class RefinementIteratorSpecial<dimension, CoordType, 0>;
        friend class RefinementIteratorSpecial<dimension, CoordType, dimension>;

        typedef Simplex::RefinementImp<dimension, CoordType> BackendRefinement;
      };

      template<int dimension, class CoordType>
      template<int codimension>
      struct RefinementImp<dimension, CoordType>::Codim
      {
        class SubEntityIterator;
        typedef Dune::Geometry<dimension-codimension, dimension, RefinementImp<dimension, CoordType>, Geometry> Geometry;
      };

      template<int dimension, class CoordType>
      int
      RefinementImp<dimension, CoordType>::
      nVertices(int level)
      {
        return BackendRefinement::nVertices(level) * Factorial<dimension>::factorial;
      }

      template<int dimension, class CoordType>
      typename RefinementImp<dimension, CoordType>::VertexIterator
      RefinementImp<dimension, CoordType>::
      vBegin(int level)
      {
        return VertexIterator(level);
      }

      template<int dimension, class CoordType>
      typename RefinementImp<dimension, CoordType>::VertexIterator
      RefinementImp<dimension, CoordType>::
      vEnd(int level)
      {
        return VertexIterator(level, true);
      }

      template<int dimension, class CoordType>
      int
      RefinementImp<dimension, CoordType>::
      nElements(int level)
      {
        return BackendRefinement::nElements(level) * Factorial<dimension>::factorial;
      }

      template<int dimension, class CoordType>
      typename RefinementImp<dimension, CoordType>::ElementIterator
      RefinementImp<dimension, CoordType>::
      eBegin(int level)
      {
        return ElementIterator(level);
      }

      template<int dimension, class CoordType>
      typename RefinementImp<dimension, CoordType>::ElementIterator
      RefinementImp<dimension, CoordType>::
      eEnd(int level)
      {
        return ElementIterator(level, true);
      }

      // //////////////
      //
      // The iterator
      //

      // vertices

      template<int dimension, class CoordType>
      class RefinementIteratorSpecial<dimension, CoordType, dimension>
      {
      public:
        typedef RefinementImp<dimension, CoordType> Refinement;
        typedef typename Refinement::CoordVector CoordVector;

        RefinementIteratorSpecial(int level, bool end = false);

        void increment();

        CoordVector coords() const;
        int index() const;
      protected:
        typedef typename Refinement::BackendRefinement BackendRefinement;
        typedef typename BackendRefinement::template Codim<dimension>::SubEntityIterator BackendIterator;
        enum { nKuhnSimplices = Factorial<dimension>::factorial };

        int level;

        int kuhnIndex;
        BackendIterator backend;
        const BackendIterator backendEnd;
      };

      template<int dimension, class CoordType>
      RefinementIteratorSpecial<dimension, CoordType, dimension>::
      RefinementIteratorSpecial(int level_, bool end)
        : level(level_), kuhnIndex(0),
          backend(BackendRefinement::vBegin(level)),
          backendEnd(BackendRefinement::vEnd(level))
      {
        if(end)
          kuhnIndex = nKuhnSimplices;
      }

      template<int dimension, class CoordType>
      void
      RefinementIteratorSpecial<dimension, CoordType, dimension>::
      increment()
      {
        ++backend;
        if(backend == backendEnd) {
          backend = BackendRefinement::vBegin(level);
          ++kuhnIndex;
        }
      }

      template<int dimension, class CoordType>
      typename RefinementIteratorSpecial<dimension, CoordType, dimension>::CoordVector
      RefinementIteratorSpecial<dimension, CoordType, dimension>::
      coords() const
      {
        return referenceToKuhn(backend.coords(), getPermutation<dimension>(kuhnIndex));
      }

      template<int dimension, class CoordType>
      int
      RefinementIteratorSpecial<dimension, CoordType, dimension>::
      index() const
      {
        return kuhnIndex*BackendRefinement::nVertices(level) + backend.index();
      }

      // elements

      template<int dimension, class CoordType>
      class RefinementIteratorSpecial<dimension, CoordType, 0>
      {
      public:
        typedef RefinementImp<dimension, CoordType> Refinement;
        typedef typename Refinement::IndexVector IndexVector;
        typedef typename Refinement::CoordVector CoordVector;
        typedef typename Refinement::template Codim<0>::Geometry Geometry;

        RefinementIteratorSpecial(int level, bool end = false);
        RefinementIteratorSpecial(const RefinementIteratorSpecial<dimension, CoordType, 0> &other);

        void increment();

        IndexVector vertexIndices() const;
        int index() const;
        CoordVector coords() const;
        Geometry geometry() const;

      protected:
        typedef typename Refinement::BackendRefinement BackendRefinement;
        typedef typename BackendRefinement::template Codim<0>::SubEntityIterator BackendIterator;
        enum { nKuhnSimplices = Factorial<dimension>::factorial };

        int level;

        int kuhnIndex;
        BackendIterator backend;
        const BackendIterator backendEnd;
      private:
        mutable bool builtGeometry;
        mutable HCubeTriangulation::Geometry< dimension, dimension, RefinementImp< dimension, CoordType > > geometry_;
      };

      template<int dimension, class CoordType>
      RefinementIteratorSpecial<dimension, CoordType, 0>::
      RefinementIteratorSpecial(int level_, bool end)
        : level(level_), kuhnIndex(0),
          backend(BackendRefinement::eBegin(level)),
          backendEnd(BackendRefinement::eEnd(level)),
          builtGeometry(false), geometry_(backend)
      {
        if(end)
          kuhnIndex = nKuhnSimplices;
      }
      template<int dimension, class CoordType>
      RefinementIteratorSpecial<dimension, CoordType, 0>::
      RefinementIteratorSpecial(const RefinementIteratorSpecial<dimension, CoordType, 0> &other)
        : level(other.level), kuhnIndex(other.kuhnIndex),
          backend(other.backend),
          backendEnd(other.backendEnd),
          builtGeometry(false), geometry_(backend)
      {}

      template<int dimension, class CoordType>
      void
      RefinementIteratorSpecial<dimension, CoordType, 0>::
      increment()
      {
        ++backend;
        builtGeometry = false;
        if(backend == backendEnd) {
          backend = BackendRefinement::eBegin(level);
          ++kuhnIndex;
        }
      }

      template<int dimension, class CoordType>
      typename RefinementIteratorSpecial<dimension, CoordType, 0>::IndexVector
      RefinementIteratorSpecial<dimension, CoordType, 0>::
      vertexIndices() const
      {
        IndexVector indices = backend.vertexIndices();

        int base = kuhnIndex * BackendRefinement::nVertices(level);
        indices += base;

        return indices;
      }

      template<int dimension, class CoordType>
      int
      RefinementIteratorSpecial<dimension, CoordType, 0>::
      index() const
      {
        return kuhnIndex*BackendRefinement::nElements(level) + backend.index();
      }

      template<int dimension, class CoordType>
      typename RefinementIteratorSpecial<dimension, CoordType, 0>::CoordVector
      RefinementIteratorSpecial<dimension, CoordType, 0>::
      coords() const
      {
        return geometry()
               .global(GenericReferenceElements<CoordType, dimension>
                       ::simplex().position(0,0));
      }

      template<int dimension, class CoordType>
      typename RefinementIteratorSpecial<dimension, CoordType, 0>::Geometry
      RefinementIteratorSpecial<dimension, CoordType, 0>::geometry () const
      {
        if(!builtGeometry) {
          geometry_.make(kuhnIndex);
          builtGeometry = true;
        }

        return Geometry( geometry_ );
      }

      // common

      template<int dimension, class CoordType>
      template<int codimension>
      class RefinementImp<dimension, CoordType>::Codim<codimension>::SubEntityIterator
        : public ForwardIteratorFacade<typename RefinementImp<dimension, CoordType>::template Codim<codimension>::SubEntityIterator, int>,
          public RefinementIteratorSpecial<dimension, CoordType, codimension>
      {
      public:
        typedef RefinementImp<dimension, CoordType> Refinement;
        typedef SubEntityIterator This;

        SubEntityIterator(int level, bool end = false);

        bool equals(const This &other) const;
      protected:
        using RefinementIteratorSpecial<dimension, CoordType, codimension>::kuhnIndex;
        using RefinementIteratorSpecial<dimension, CoordType, codimension>::backend;
      };

#ifndef DOXYGEN

      template<int dimension, class CoordType>
      template<int codimension>
      RefinementImp<dimension, CoordType>::Codim<codimension>::SubEntityIterator::
      SubEntityIterator(int level, bool end)
        : RefinementIteratorSpecial<dimension, CoordType, codimension>(level, end)
      {}

      template<int dimension, class CoordType>
      template<int codimension>
      bool
      RefinementImp<dimension, CoordType>::Codim<codimension>::SubEntityIterator::
      equals(const This &other) const
      { return kuhnIndex == other.kuhnIndex && backend == other.backend; }

#endif // DOXYGEN

      // ///////////
      //
      //  Geometry
      //

      template<int mydimension, int coorddimension, class GridImp>
      class Geometry : public GeometryDefaultImplementation<mydimension, coorddimension, GridImp, Geometry>
      {
        typedef typename GridImp::ctype ct;
        enum { dimension = GridImp::dimension };

        typedef typename GridImp::BackendRefinement BackendRefinement;
        typedef typename BackendRefinement::template Codim<dimension-mydimension>::SubEntityIterator BackendIterator;
      public:
        GeometryType type() const
        { return GeometryType(GeometryType::simplex, mydimension); }

        int corners() const
        { return mydimension + 1; }

        const FieldVector<ct, coorddimension>& operator[] (int i) const
        {
          if(!builtCoords) {
            for(int i = 0; i < corners(); ++i)
              coords[i] = referenceToKuhn(backend->geometry()[i], getPermutation<dimension>(kuhnIndex));
            builtCoords = true;
          }
          return coords[i];
        }

        FieldVector<ct, coorddimension> global(const FieldVector<ct, mydimension>& local) const
        { return referenceToKuhn(backend.geometry().global(local), getPermutation<dimension>(kuhnIndex)); }

        FieldVector<ct, mydimension> local(const FieldVector<ct, coorddimension>& global) const
        { return backend->geometry().local(kuhnToReference(global, getPermutation<dimension>(kuhnIndex))); }

        bool checkInside (const FieldVector<ct, mydimension>& local) const
        { return backend->geometry().checkInside(local); }

        ct integrationElement(const FieldVector<ct, mydimension>& local) const
        { return backend->geometry().integrationElement(local) / Factorial<dimension>::factorial; }

        const FieldMatrix<ct, mydimension, mydimension>& jacobianInverse(const FieldVector<ct, mydimension>& local) const
        {
          if(!builtJinv) {
            // create unit vectors
            FieldMatrix<int, mydimension, mydimension> M = 0;
            for(int i = 0; i < mydimension; ++i)
              M[i][i] = 1;
            // transform them into local coordinates
            for(int i = 0; i < mydimension; ++i)
              M[i] = kuhnToReference(M[i], getPermutation<mydimension>(kuhnIndex));
            // transpose the matrix
            for(int i = 0; i < mydimension; ++i)
              for(int j = 0; j < mydimension; ++j)
                Jinv[i][j] = M[j][i];
            // take the backends inverse Jacobian into account
            Jinv.leftmultiply(backend->geometry().jacobianInverse(local));

            builtJinv = true;
          }

          return Jinv;
        }

        Geometry(const BackendIterator &backend_)
          : coords(), builtCoords(false), Jinv(), builtJinv(false),
            backend(backend_), kuhnIndex(0)
        {
          dune_static_assert(mydimension == coorddimension, "mydimension != coorddimension");
        }

        void make(int kuhnIndex_)
        {
          kuhnIndex = kuhnIndex_;
          builtCoords = false;
          builtJinv = false;
        }
      private:
        mutable FieldVector<FieldVector<ct, coorddimension>, mydimension+1> coords;
        mutable bool builtCoords;
        mutable FieldMatrix<ct, mydimension, mydimension> Jinv;
        mutable bool builtJinv;
        const BackendIterator &backend;
        int kuhnIndex;
      };

    } // namespace HCubeTriangulation

    // ///////////////////////
    //
    // The refinement traits
    //

#ifndef DOXYGEN
    template<unsigned topologyId, class CoordType, unsigned coerceToId,
        int dim>
    struct Traits<
        topologyId, CoordType, coerceToId, dim,
        typename enable_if<
            (dim >= 2 &&
             (GenericGeometry::CubeTopology<dim>::type::id >> 1) ==
             (topologyId >> 1) &&
             (GenericGeometry::SimplexTopology<dim>::type::id >> 1) ==
             (coerceToId >> 1)
            )>::type
        >
    {
      typedef HCubeTriangulation::RefinementImp<dim, CoordType> Imp;
    };
#endif

  } // namespace RefinementImp

} // namespace Dune

#endif //DUNE_GRID_COMMON_REFINEMENT_HCUBETRIANGULATION_CC