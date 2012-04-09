// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_COMMON_REFINEMENT_HCUBE_CC
#define DUNE_GRID_COMMON_REFINEMENT_HCUBE_CC

// This file is part of DUNE, a Distributed and Unified Numerics Environment
// This file is copyright (C) 2005 Jorrit Fahlke <jorrit@jorrit.de>, 2012 Christian Engwer <christi@dune-project.org>
// This file is licensed under version 2 of the GNU General Public License,
// with a special "runtime exception."  See COPYING at the top of the source
// tree for the full licence.

/*! @file

   @brief This file contains the @ref Refinement implementation for
         hypercubes (quadrilaterals, hexahedrons, etc.).

   See @ref HCubeRefinement.
   @verbatim
   $Id$
   @endverbatim
 */

/*! @defgroup HCubeRefinement Refinement implementation for hypercubes
    @ingroup Refinement

   @section Iterators The Iterators
   <!--=========================-->

   For the iterators we have to hack around a bit.  The problem is as
   follows:
   @code
   template<int A>
   class outer
   {
    template<int B>
    class inner;
   };
   @endcode
   C++ does not allow specialisation of the inner class when the outer
   class is not specialized.

   So I had to create a baseclass for the iterators which is not inside
   another class.  This base class can then be specialized, and the
   real Iterator class inherits from it.  I gave it the somewhat clumsy
   name RefinementSubEntityIteratorSpecial.
 */

#include <cassert>

#include <dune/common/fvector.hh>
#include <dune/common/iteratorfacades.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/genericgeometry/geometry.hh>

#include <dune/grid/common/geometry.hh>
#include "base.cc" // for RefinementTraits

namespace Dune {

  namespace RefinementImp {

    /*! @brief This namespace contains the @ref Refinement implementation
       for hypercubes (GeometryType::cube).

       See @ref HCubeRefinement.
     */
    namespace HCube {

#ifndef DOXYGEN
      template <int p, bool odd = p%2>
      struct PowerImp {};
#endif

      //! Helper class for power computation
      template <int p>
      struct Power
      {
        template <typename T>
        static T eval(const T & a)
        {
          return PowerImp<p>::eval(a);
        }
      };

#ifndef DOXYGEN
      template <int p>
      struct PowerImp<p,false>
      {
        template <typename T>
        static T eval(const T & a)
        {
          T t = Power<p/2>::eval(a);
          return t*t;
        }
      };

      template <int p>
      struct PowerImp<p,true>
      {
        template <typename T>
        static T eval(const T & a)
        {
          return a*Power<p-1>::eval(a);;
        }
      };

      template <>
      struct PowerImp<1,true>
      {
        template <typename T>
        static T eval(const T & a)
        {
          return a;
        }
      };
#endif

      /*! @brief @ref Refinement implementation for hypercubes

          @param dimension_ Dimension of the refined hypercube
          @param CoordType  Coordinate type of the refined hypercube

          We use @ref RefinementGrid as backend to do all the work.

          The interface is the same as for @ref Dune::StaticRefinement (apart
          from the template parameters).
       */
      template<int dimension_, class CoordType>
      class RefinementImp
      {
      public:
        enum { dimension = dimension_ /*!< Know your own dimension @hideinitializer */ };
        //- Know yourself
        typedef RefinementImp<dimension, CoordType> Refinement;

        template<int codimension>
        struct Codim;
        typedef typename Codim<dimension>::SubEntityIterator VertexIterator;
        typedef FieldVector<CoordType, dimension> CoordVector;
        typedef typename Codim<0>::SubEntityIterator ElementIterator;
        typedef FieldVector<int, (1<<dimension)> IndexVector;

        static int nVertices(int level);
        static VertexIterator vBegin(int level);
        static VertexIterator vEnd(int level);

        static int nElements(int level);
        static ElementIterator eBegin(int level);
        static ElementIterator eEnd(int level);
      };

      template<int dimension, class CoordType>
      template<int codimension>
      struct RefinementImp<dimension, CoordType>::Codim
      {
        class SubEntityIterator;
        typedef Dune::Geometry<dimension-codimension, dimension, RefinementImp<dimension, CoordType>, GenericGeometry::Geometry> Geometry;
      };

      template<int dimension, class CoordType>
      int
      RefinementImp<dimension, CoordType>::
      nVertices(int level)
      {
        // return (2^level + 1)^dim
        return Power<dimension>::eval((1<<level)+1);
      }

      template<int dimension, class CoordType>
      typename RefinementImp<dimension, CoordType>::VertexIterator
      RefinementImp<dimension, CoordType>::
      vBegin(int level)
      {
        return VertexIterator(0,level);
      }

      template<int dimension, class CoordType>
      typename RefinementImp<dimension, CoordType>::VertexIterator
      RefinementImp<dimension, CoordType>::
      vEnd(int level)
      {
        return VertexIterator(nVertices(level),level);
      }

      template<int dimension, class CoordType>
      int
      RefinementImp<dimension, CoordType>::
      nElements(int level)
      {
        // return (2^level)^dim
        return 1<<(level*dimension);
      }

      template<int dimension, class CoordType>
      typename RefinementImp<dimension, CoordType>::ElementIterator
      RefinementImp<dimension, CoordType>::
      eBegin(int level)
      {
        return ElementIterator(0,level);
      }

      template<int dimension, class CoordType>
      typename RefinementImp<dimension, CoordType>::ElementIterator
      RefinementImp<dimension, CoordType>::
      eEnd(int level)
      {
        return ElementIterator(nElements(level),level);
      }

      //
      // The iterators
      //

#ifdef DOXYGEN
      /*! @brief SubEntityIterator base class for hypercube refinement

          @tparam dimension   Dimension of the refined element
          @tparam CoordType   Coordinate type of the refined element
          @tparam codimension Codimension of the iterator

          This is the base class for SubEntityIterators.  We have to use
          this construct because RefinementImp<...>::%codim<...> cannot
          be specialized without first specializing RefinementImp.
       */
      template<int dimension, class CoordType, int codimension>
      class RefinementSubEntityIteratorSpecial {};
#else //!DOXYGEN
      template<int dimension, class CoordType, int codimension>
      class RefinementSubEntityIteratorSpecial;
#endif //DOXYGEN

      // for vertices

      template<int dimension, class CoordType>
      class RefinementSubEntityIteratorSpecial<dimension, CoordType, dimension>
      {
      public:
        typedef RefinementImp<dimension, CoordType> Refinement;
        typedef typename Refinement::template Codim<dimension>::SubEntityIterator Common;
        typedef typename Refinement::CoordVector CoordVector;

        CoordVector coords() const;

      private:
        const Common & asCommon() const
        {
          return *static_cast<const Common*>(this);
        };
      };

      template<int dimension, class CoordType>
      typename RefinementSubEntityIteratorSpecial<dimension, CoordType, dimension>::CoordVector
      RefinementSubEntityIteratorSpecial<dimension, CoordType, dimension>::
      coords() const
      {
        array<unsigned int, dimension> v(asCommon().vertexCoord());
        CoordVector c;
        for (int d=0; d<dimension; d++)
        {
          c[d] = v[d]*1.0 / (1<<asCommon()._level);
        }
        return c;
      }

      // for elements

      template<int dimension, class CoordType>
      class RefinementSubEntityIteratorSpecial<dimension, CoordType, 0>
      {
      public:
        typedef RefinementImp<dimension, CoordType> Refinement;
        typedef typename Refinement::template Codim<0>::SubEntityIterator Common;
        typedef typename Refinement::IndexVector IndexVector;
        typedef typename Refinement::CoordVector CoordVector;

        IndexVector vertexIndices() const;
        CoordVector coords() const;

      private:
        const Common & asCommon() const
        {
          return *static_cast<const Common*>(this);
        };
      };

      template<int dimension, class CoordType>
      typename RefinementSubEntityIteratorSpecial<dimension, CoordType, 0>::IndexVector
      RefinementSubEntityIteratorSpecial<dimension, CoordType, 0>::
      vertexIndices() const
      {
        enum { nIndices = (1 << dimension) };

        // cell index tuple
        array<unsigned int, dimension> e(asCommon().cellCoord());

        // vertices
        IndexVector vec;
        array<unsigned int, dimension> v;
        for(int i = 0; i < nIndices; ++i)
        {
          for (int d = 0; d < dimension; d++)
          {
            v[d] = e[d];
            if (i & (1<<d))
              v[d]++;
          }
          // compute vertex index tuple from cell tuple
          vec[nIndices-1-i] = asCommon().vertexIdx(v);
        }
        return vec;
      }

      template<int dimension, class CoordType>
      typename RefinementSubEntityIteratorSpecial<dimension, CoordType, 0>::CoordVector
      RefinementSubEntityIteratorSpecial<dimension, CoordType, 0>::
      coords() const
      {
        array<unsigned int, dimension> v(asCommon().cellCoord());
        CoordVector c;
        for (int d=0; d<dimension; d++)
        {
          c[d] = (v[d]*1.0 + 0.5) / (1<<asCommon()._level);
        }
        return c;
      }



      // common
      template<int dimension, class CoordType>
      template<int codimension>
      class RefinementImp<dimension, CoordType>::Codim<codimension>::SubEntityIterator
        : public ForwardIteratorFacade<typename RefinementImp<dimension, CoordType>::template Codim<codimension>::SubEntityIterator, int>,
          public RefinementSubEntityIteratorSpecial<dimension, CoordType, codimension>
      {
      public:
        typedef RefinementImp<dimension, CoordType> Refinement;
        typedef typename Refinement::template Codim<codimension>::SubEntityIterator This;

        SubEntityIterator(unsigned int index, unsigned int level);

        bool equals(const This &other) const;
        void increment();

        int index() const;
        Geometry geometry () const;
      private:
        friend class RefinementSubEntityIteratorSpecial<dimension, CoordType, codimension>;
        unsigned int _index;
        unsigned int _level;

        array<unsigned int, dimension>
        cellCoord(unsigned int idx) const
        {
          return idx2coord(idx, 1<<_level);
        }

        array<unsigned int, dimension>
        vertexCoord(unsigned int idx) const
        {
          return idx2coord(idx, (1<<_level)+1);
        }

        array<unsigned int, dimension>
        cellCoord() const
        {
          return cellCoord(_index);
        }

        array<unsigned int, dimension>
        vertexCoord() const
        {
          return vertexCoord(_index);
        }

        array<unsigned int, dimension>
        idx2coord(unsigned int idx, unsigned int w) const
        {
          array<unsigned int, dimension> c;
          for (unsigned int d = 0; d < dimension; d++)
          {
            c[d] = idx%w;
            idx = idx/w;
          }
          return c;
        }

        unsigned int
        coord2idx(array<unsigned int, dimension> c, unsigned int w) const
        {
          unsigned int i = 0;
          for (unsigned int d = dimension; d > 0; d--)
          {
            i *= w;
            i += c[d-1];
          }
          return i;
        }

        unsigned int
        vertexIdx(array<unsigned int, dimension> c) const
        {
          return coord2idx(c, (1<<_level)+1);
        }

      };

#ifndef DOXYGEN

      template<int dimension, class CoordType>
      template<int codimension>
      RefinementImp<dimension, CoordType>::Codim<codimension>::SubEntityIterator::
      SubEntityIterator(unsigned int index, unsigned int level)
        : _index(index), _level(level)
      {}

      template<int dimension, class CoordType>
      template<int codimension>
      bool
      RefinementImp<dimension, CoordType>::Codim<codimension>::SubEntityIterator::
      equals(const This &other) const
      {
        return _index == other._index && _level == other._level;
      }

      template<int dimension, class CoordType>
      template<int codimension>
      void
      RefinementImp<dimension, CoordType>::Codim<codimension>::SubEntityIterator::
      increment()
      {
        ++_index;
      }

      template<int dimension, class CoordType>
      template<int codimension>
      int
      RefinementImp<dimension, CoordType>::Codim<codimension>::SubEntityIterator::
      index() const
      {
        return _index;
      }

      template<int dimension, class CoordType>
      template<int codimension>
      typename RefinementImp<dimension, CoordType>::template Codim<codimension>::Geometry
      RefinementImp<dimension, CoordType>::Codim<codimension>::SubEntityIterator::geometry () const
      {
        assert(false && "Not Implemented");
      }

#endif // DOXYGEN

    } // namespace HCube

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
             (GenericGeometry::CubeTopology<dim>::type::id >> 1) ==
             (coerceToId >> 1)
            )>::type
        >
    {
      typedef HCube::RefinementImp<dim, CoordType> Imp;
    };
#endif

  } // namespace RefinementImp

} // namespace Dune

#endif //DUNE_GRID_COMMON_REFINEMENT_HCUBE_CC
