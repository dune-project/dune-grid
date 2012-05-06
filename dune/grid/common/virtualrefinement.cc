// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_COMMON_VIRTUALREFINEMENT_CC
#define DUNE_GRID_COMMON_VIRTUALREFINEMENT_CC

// This file is part of DUNE, a Distributed and Unified Numerics Environment
// This file is copyright (C) 2005 Jorrit Fahlke <jorrit@jorrit.de>
// This file is licensed under version 2 of the GNU General Public License,
// with a special "runtime exception."  See COPYING at the top of the source
// tree for the full licence.

/*! @file

   @brief This file contains the virtual wrapper around refinement.

   @verbatim
   $Id$
   @endverbatim
 */

#include <cassert>
#include <typeinfo>

#include <dune/common/deprecated.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/iteratorfacades.hh>

#include <dune/geometry/type.hh>

#include "refinement.hh"

namespace Dune {

  // //////////////////////////////////////////
  //
  // The virtual base class and its iterators
  //

  //
  // Refinement
  //

  template<int dimension, class CoordType>
  typename VirtualRefinement<dimension, CoordType>::VertexIterator
  VirtualRefinement<dimension, CoordType>::
  vBegin(int level) const
  {
    return VertexIterator(vBeginBack(level));
  }

  template<int dimension, class CoordType>
  typename VirtualRefinement<dimension, CoordType>::VertexIterator
  VirtualRefinement<dimension, CoordType>::
  vEnd(int level) const
  {
    return VertexIterator(vEndBack(level));
  }

  template<int dimension, class CoordType>
  typename VirtualRefinement<dimension, CoordType>::ElementIterator
  VirtualRefinement<dimension, CoordType>::
  eBegin(int level) const
  {
    return ElementIterator(eBeginBack(level));
  }

  template<int dimension, class CoordType>
  typename VirtualRefinement<dimension, CoordType>::ElementIterator
  VirtualRefinement<dimension, CoordType>::
  eEnd(int level) const
  {
    return ElementIterator(eEndBack(level));
  }

  //
  // The iterators
  //

  template<int dimension, class CoordType, int codimension>
  class VirtualRefinementSubEntityIteratorSpecial;

  // The iterator for vertices

  template<int dimension, class CoordType>
  class VirtualRefinementSubEntityIteratorSpecial<dimension, CoordType, dimension>
  {};

  // The iterator for elements

  template<int dimension, class CoordType>
  class VirtualRefinementSubEntityIteratorSpecial<dimension, CoordType, 0>
  {
  public:
    typedef VirtualRefinement<dimension, CoordType> Refinement;
    typedef typename Refinement::template Codim<0>::SubEntityIterator Common;
    typedef typename Refinement::IndexVector IndexVector;

    IndexVector vertexIndices() const;
  };

  template<int dimension, class CoordType>
  typename VirtualRefinementSubEntityIteratorSpecial<dimension, CoordType, 0>::IndexVector
  VirtualRefinementSubEntityIteratorSpecial<dimension, CoordType, 0>::
  vertexIndices() const
  {
    return static_cast<const Common *>(this)->backend->vertexIndices();
  }

  // The iterator common stuff

  template<int dimension, class CoordType>
  template<int codimension>
  class VirtualRefinement<dimension, CoordType>::Codim<codimension>::SubEntityIterator
    : public ForwardIteratorFacade<typename VirtualRefinement<dimension, CoordType>::template Codim<codimension>::SubEntityIterator, int>,
      public VirtualRefinementSubEntityIteratorSpecial<dimension, CoordType, codimension>
  {
  public:
    typedef VirtualRefinement<dimension, CoordType> Refinement;
    typedef typename Refinement::template Codim<codimension>::SubEntityIterator This;
    typedef typename Refinement::template SubEntityIteratorBack<codimension> IteratorBack;
    typedef typename Refinement::CoordVector CoordVector;

    SubEntityIterator(IteratorBack *backend);
    SubEntityIterator(const This &other);
    ~SubEntityIterator();

    This &operator=(const This &other);

    bool equals(const This &other) const;
    void increment();

    int index() const;

    // If you simply use an unqualified CoordVector here g++-4.2 chokes
    typename VirtualRefinement<dimension, CoordType>::template Codim<codimension>::SubEntityIterator::
    CoordVector coords() const;
  private:
    friend class VirtualRefinementSubEntityIteratorSpecial<dimension, CoordType, codimension>;
    IteratorBack *backend;
  };

#ifndef DOXYGEN

  template<int dimension, class CoordType>
  template<int codimension>
  VirtualRefinement<dimension, CoordType>::Codim<codimension>::SubEntityIterator::
  SubEntityIterator(IteratorBack *backend_)
    : backend(backend_)
  {}

  template<int dimension, class CoordType>
  template<int codimension>
  VirtualRefinement<dimension, CoordType>::Codim<codimension>::SubEntityIterator::
  SubEntityIterator(const This &other)
    : backend(other.backend->clone())
  {}

  template<int dimension, class CoordType>
  template<int codimension>
  VirtualRefinement<dimension, CoordType>::Codim<codimension>::SubEntityIterator::
  ~SubEntityIterator()
  {
    delete backend;
  }

  template<int dimension, class CoordType>
  template<int codimension>
  typename VirtualRefinement<dimension, CoordType>::template Codim<codimension>::SubEntityIterator &
  VirtualRefinement<dimension, CoordType>::Codim<codimension>::SubEntityIterator::
  operator=(const This &other)
  {
    delete backend;
    backend = other.backend->clone();
  }

  template<int dimension, class CoordType>
  template<int codimension>
  bool
  VirtualRefinement<dimension, CoordType>::Codim<codimension>::SubEntityIterator::
  equals(const This &other) const
  { return *backend == *(other.backend); }

  template<int dimension, class CoordType>
  template<int codimension>
  void
  VirtualRefinement<dimension, CoordType>::Codim<codimension>::SubEntityIterator::
  increment()
  {
    ++*backend;
  }

  template<int dimension, class CoordType>
  template<int codimension>
  int
  VirtualRefinement<dimension, CoordType>::Codim<codimension>::SubEntityIterator::
  index() const
  { return backend->index(); }

  template<int dimension, class CoordType>
  template<int codimension>
  typename VirtualRefinement<dimension, CoordType>::template Codim<codimension>::SubEntityIterator::CoordVector
  VirtualRefinement<dimension, CoordType>::Codim<codimension>::SubEntityIterator::
  coords() const
  { return backend->coords(); }

#endif // DOXYGEN

  //
  // The iterator backend
  //

  template<int dimension, class CoordType, int codimension>
  class VirtualRefinementSubEntityIteratorBackSpecial;

  // The iterator backend for vertices

  template<int dimension, class CoordType>
  class VirtualRefinementSubEntityIteratorBackSpecial<dimension, CoordType, dimension>
  {
  public:

    virtual ~VirtualRefinementSubEntityIteratorBackSpecial()
    {}

  };

  // The iterator backend for elements

  template<int dimension, class CoordType>
  class VirtualRefinementSubEntityIteratorBackSpecial<dimension, CoordType, 0>
  {
  public:
    typedef VirtualRefinement<dimension, CoordType> Refinement;
    typedef typename Refinement::IndexVector IndexVector;

    virtual IndexVector vertexIndices() const = 0;

    virtual ~VirtualRefinementSubEntityIteratorBackSpecial()
    {}

  };

  // The iterator backend common stuff

  template<int dimension, class CoordType>
  template<int codimension>
  class VirtualRefinement<dimension, CoordType>::SubEntityIteratorBack
    : public VirtualRefinementSubEntityIteratorBackSpecial<dimension, CoordType, codimension>
  {
  public:
    typedef VirtualRefinement<dimension, CoordType> Refinement;
    typedef typename Refinement::template SubEntityIteratorBack<codimension> This;
    typedef typename Refinement::CoordVector CoordVector;

    virtual ~SubEntityIteratorBack() {}

    virtual This *clone() const = 0;

    virtual bool operator==(const This &other) const = 0;
    virtual This &operator++() = 0;

    virtual int index() const = 0;
    virtual CoordVector coords() const = 0;
  };

  // /////////////////////////////////////////////////
  //
  // The derived classes and their iterator backends
  //

  //
  // The refinement implementation
  //

  template<
      unsigned topologyId
      , class CoordType
      , unsigned coerceToId
      , int dimension
      >
  class VirtualRefinementImp
    : public Dune::VirtualRefinement<dimension, CoordType>
  {
  public:
    typedef Dune::StaticRefinement<topologyId, CoordType, coerceToId, dimension> StaticRefinement;
    typedef Dune::VirtualRefinement<dimension, CoordType> VirtualRefinement;

    template<int codimension>
    class SubEntityIteratorBack;

    int nVertices(int level) const;
    int nElements(int level) const;

    static VirtualRefinementImp<topologyId, CoordType, coerceToId, dimension> &instance();
  private:
    VirtualRefinementImp() {}
    static VirtualRefinementImp<topologyId, CoordType, coerceToId, dimension> *instance_;

    typename VirtualRefinement::VertexIteratorBack *vBeginBack(int level) const;
    typename VirtualRefinement::VertexIteratorBack *vEndBack(int level) const;
    typename VirtualRefinement::ElementIteratorBack *eBeginBack(int level) const;
    typename VirtualRefinement::ElementIteratorBack *eEndBack(int level) const;
  };

  template<
      unsigned topologyId
      , class CoordType
      , unsigned coerceToId
      , int dimension
      >
  VirtualRefinementImp<topologyId, CoordType, coerceToId, dimension> &
  VirtualRefinementImp<topologyId, CoordType, coerceToId, dimension>::instance()
  {
    if(instance_ == 0)
      instance_ = new VirtualRefinementImp<topologyId, CoordType, coerceToId, dimension>;
    return *instance_;
  }

  template<
      unsigned topologyId
      , class CoordType
      , unsigned coerceToId
      , int dimension
      >
  VirtualRefinementImp<topologyId, CoordType, coerceToId, dimension> *
  VirtualRefinementImp<topologyId, CoordType, coerceToId, dimension>::instance_ = 0;

  template<
      unsigned topologyId
      , class CoordType
      , unsigned coerceToId
      , int dimension
      >
  int VirtualRefinementImp<topologyId, CoordType, coerceToId, dimension>::
  nVertices(int level) const
  {
    return StaticRefinement::nVertices(level);
  }

  template<
      unsigned topologyId
      , class CoordType
      , unsigned coerceToId
      , int dimension
      >
  typename VirtualRefinementImp<topologyId, CoordType, coerceToId, dimension>::VirtualRefinement::VertexIteratorBack *
  VirtualRefinementImp<topologyId, CoordType, coerceToId, dimension>::
  vBeginBack(int level) const
  { return new SubEntityIteratorBack<dimension>(StaticRefinement::vBegin(level)); }

  template<
      unsigned topologyId
      , class CoordType
      , unsigned coerceToId
      , int dimension
      >
  typename VirtualRefinementImp<topologyId, CoordType, coerceToId, dimension>::VirtualRefinement::VertexIteratorBack *
  VirtualRefinementImp<topologyId, CoordType, coerceToId, dimension>::
  vEndBack(int level) const
  { return new SubEntityIteratorBack<dimension>(StaticRefinement::vEnd(level)); }

  template<
      unsigned topologyId
      , class CoordType
      , unsigned coerceToId
      , int dimension
      >
  int VirtualRefinementImp<topologyId, CoordType, coerceToId, dimension>::nElements(int level) const
  {
    return StaticRefinement::nElements(level);
  }

  template<
      unsigned topologyId
      , class CoordType
      , unsigned coerceToId
      , int dimension
      >
  typename VirtualRefinementImp<topologyId, CoordType, coerceToId, dimension>::VirtualRefinement::ElementIteratorBack *
  VirtualRefinementImp<topologyId, CoordType, coerceToId, dimension>::
  eBeginBack(int level) const
  { return new SubEntityIteratorBack<0>(StaticRefinement::eBegin(level)); }

  template<
      unsigned topologyId
      , class CoordType
      , unsigned coerceToId
      , int dimension
      >
  typename VirtualRefinementImp<topologyId, CoordType, coerceToId, dimension>::VirtualRefinement::ElementIteratorBack *
  VirtualRefinementImp<topologyId, CoordType, coerceToId, dimension>::
  eEndBack(int level) const
  { return new SubEntityIteratorBack<0>(StaticRefinement::eEnd(level)); }

  //
  // The iterator backend implementation
  //

  // The iterator backend implementation specialties

  template<
      unsigned topologyId
      , class CoordType
      , unsigned coerceToId
      , int dimension
      , int codimension
      >
  class VirtualRefinementImpSubEntityIteratorBackSpecial;

  // The iterator backend implementation specialties for vertices

  template<
      unsigned topologyId
      , class CoordType
      , unsigned coerceToId
      , int dimension
      >
  class VirtualRefinementImpSubEntityIteratorBackSpecial<topologyId, CoordType, coerceToId, dimension, dimension>
    : public VirtualRefinement<dimension, CoordType>::template SubEntityIteratorBack<dimension>
  {};

  // The iterator backend implementation specialties for elements

  template<
      unsigned topologyId
      , class CoordType
      , unsigned coerceToId
      , int dimension
      >
  class VirtualRefinementImpSubEntityIteratorBackSpecial<topologyId, CoordType, coerceToId, dimension, 0>
    : public VirtualRefinement<dimension, CoordType>::template SubEntityIteratorBack<0>
  {
  public:
    typedef Dune::VirtualRefinementImp<topologyId, CoordType, coerceToId, dimension> VirtualRefinementImp;
    typedef typename VirtualRefinementImp::template SubEntityIteratorBack<0> Common;
    typedef typename VirtualRefinementImp::StaticRefinement StaticRefinement;
    typedef VirtualRefinement<dimension, CoordType> RefinementBase;
    typedef typename RefinementBase::IndexVector IndexVector;

    IndexVector vertexIndices() const;
  };

  template<
      unsigned topologyId
      , class CoordType
      , unsigned coerceToId
      , int dimension
      >
  typename VirtualRefinementImpSubEntityIteratorBackSpecial<topologyId, CoordType, coerceToId, dimension, 0>::IndexVector
  VirtualRefinementImpSubEntityIteratorBackSpecial<topologyId, CoordType, coerceToId, dimension, 0>::
  vertexIndices() const
  {
    IndexVector vIndices;
    vIndices.reserve(StaticRefinement::IndexVector::dimension);

    typename StaticRefinement::IndexVector sIndices = static_cast<const Common *>(this)->backend.vertexIndices();
    for(int i = 0; i < StaticRefinement::IndexVector::dimension; ++i)
      vIndices.push_back(sIndices[i]);
    return vIndices;
  }

  // The shared iterator backend implementation

  template<
      unsigned topologyId
      , class CoordType
      , unsigned coerceToId
      , int dimension
      >
  template<int codimension>
  class VirtualRefinementImp<topologyId, CoordType, coerceToId, dimension>::SubEntityIteratorBack
    : public VirtualRefinementImpSubEntityIteratorBackSpecial<topologyId, CoordType, coerceToId, dimension, codimension>
  {
  public:
    typedef typename StaticRefinement::template Codim<codimension>::SubEntityIterator BackendIterator;
    typedef typename VirtualRefinementImp<topologyId, CoordType, coerceToId, dimension>::template SubEntityIteratorBack<codimension> This;
    typedef typename VirtualRefinement::template SubEntityIteratorBack<codimension> Base;
    typedef typename VirtualRefinement::CoordVector CoordVector;

    SubEntityIteratorBack(const BackendIterator &backend);
    SubEntityIteratorBack(const This &other);

    Base *clone() const;

    bool operator==(const Base &other) const;
    Base &operator++();

    int index() const;
    CoordVector coords() const;

  private:
    friend class VirtualRefinementImpSubEntityIteratorBackSpecial<topologyId, CoordType, coerceToId, dimension, codimension>;
    BackendIterator backend;
  };

  template<
      unsigned topologyId
      , class CoordType
      , unsigned coerceToId
      , int dimension
      >
  template<int codimension>
  VirtualRefinementImp<topologyId, CoordType, coerceToId, dimension>::SubEntityIteratorBack<codimension>::
  SubEntityIteratorBack(const BackendIterator &backend_)
    : backend(backend_)
  {}

  template<
      unsigned topologyId
      , class CoordType
      , unsigned coerceToId
      , int dimension
      >
  template<int codimension>
  VirtualRefinementImp<topologyId, CoordType, coerceToId, dimension>::SubEntityIteratorBack<codimension>::
  SubEntityIteratorBack(const This &other)
    : VirtualRefinementImpSubEntityIteratorBackSpecial<topologyId, CoordType, coerceToId, dimension, codimension>(other),
      backend(other.backend)
  {}

  template<unsigned topologyId, class CoordType, unsigned coerceToId, int dimension>
  template<int codimension>
  typename VirtualRefinementImp<topologyId, CoordType, coerceToId, dimension>::template SubEntityIteratorBack<codimension>::Base *
  VirtualRefinementImp<topologyId, CoordType, coerceToId, dimension>::SubEntityIteratorBack<codimension>::
  clone() const
  { return new This(*this); }

  template<unsigned topologyId, class CoordType, unsigned coerceToId, int dimension>
  template<int codimension>
  bool
  VirtualRefinementImp<topologyId, CoordType, coerceToId, dimension>::SubEntityIteratorBack<codimension>::
  operator==(const Base &other) const
  {
    try { return backend == dynamic_cast<const This &>(other).backend; }
    catch(std::bad_cast) { return false; }
  }

  template<unsigned topologyId, class CoordType, unsigned coerceToId, int dimension>
  template<int codimension>
  typename VirtualRefinementImp<topologyId, CoordType, coerceToId, dimension>::template SubEntityIteratorBack<codimension>::Base &
  VirtualRefinementImp<topologyId, CoordType, coerceToId, dimension>::SubEntityIteratorBack<codimension>::
  operator++()
  {
    ++backend;
    return *this;
  }

  template<unsigned topologyId, class CoordType, unsigned coerceToId, int dimension>
  template<int codimension>
  int
  VirtualRefinementImp<topologyId, CoordType, coerceToId, dimension>::SubEntityIteratorBack<codimension>::
  index() const
  { return backend.index(); }

  template<unsigned topologyId, class CoordType, unsigned coerceToId, int dimension>
  template<int codimension>
  typename VirtualRefinementImp<topologyId, CoordType, coerceToId, dimension>::template SubEntityIteratorBack<codimension>::CoordVector
  VirtualRefinementImp<topologyId, CoordType, coerceToId, dimension>::SubEntityIteratorBack<codimension>::
  coords() const
  { return backend.coords(); }

  // ////////////////////////
  //
  // The refinement builder
  //

  template<int dimension, class CoordType>
  class RefinementBuilder;

  /*! @brief return a reference to the VirtualRefinement according to
             the parameters

     @tparam dimension Dimension of the element to refine
     @tparam CoordType C++ type of the coordinates

     @throws NotImplemented There is no Refinement implementation for
                           the specified parameters.
   */
  template<int dimension, class CoordType>
  VirtualRefinement<dimension, CoordType> &
  buildRefinement( //! geometry type of the refined element
    GeometryType geometryType,
    //! geometry type of the subelements
    GeometryType coerceTo)
  {
    // Check that the user used valid geometry types
    assert(geometryType.dim() == dimension && coerceTo.dim() == dimension);
    return RefinementBuilder<dimension, CoordType>::build( geometryType.id(), coerceTo.id() );
  }


  // In principle the trick with the class is no longer neccessary,
  // but I'm keeping it in here so it will be easier to specialize
  // buildRefinement when someone implements pyramids and prisms
  template<int dimension, class CoordType>
  class RefinementBuilder
  {
  public:
    static
    VirtualRefinement<dimension, CoordType> &
    build(
      unsigned topologyId
      , unsigned coerceToId
      )
    {
      topologyId &= ~1;
      coerceToId &= ~1;

      const unsigned idSimplex = GenericGeometry::SimplexTopology<dimension>::type::id & ~1;
      const unsigned idCube = GenericGeometry::CubeTopology<dimension>::type::id & ~1;

      switch( topologyId )
      {
      //case GeometryType::simplex:
      case idSimplex :
        //switch( coerceTo )
        switch( coerceToId )
        {
        //case GeometryType::simplex:
        case idSimplex :
          return VirtualRefinementImp< idSimplex, CoordType, idSimplex, dimension>::instance();
        default :
          break;
        }
        break;

      //case GeometryType::cube:
      case idCube :
        switch( coerceToId )
        {
        case idSimplex :
          return VirtualRefinementImp< idCube, CoordType, idSimplex, dimension>::instance();
        case idCube :
          return VirtualRefinementImp< idCube, CoordType, idCube, dimension>::instance();
        default :
          break;
        }
        break;

      default :
        break;
      }
      DUNE_THROW( NotImplemented, "No Refinement<" << topologyId << ", CoordType, "
                                                   << coerceToId << " >.");
    }
  };

  template<class CoordType>
  class RefinementBuilder<1, CoordType>
  {
    static const std::size_t dimension = 1;
  public:
    static
    VirtualRefinement<dimension, CoordType> &
    build(
      unsigned topologyId
      , unsigned coerceToId
      )
    {
      topologyId &= ~1;
      coerceToId &= ~1;

      const unsigned idSimplex = GenericGeometry::SimplexTopology<dimension>::type::id & ~1;

      if(topologyId == 0 && coerceToId == 0)
        return VirtualRefinementImp< idSimplex, CoordType, idSimplex, dimension>::instance();

      DUNE_THROW( NotImplemented, "No Refinement<" << topologyId << ", CoordType, "
                                                   << coerceToId << " >.");
    }
  };

  template<class CoordType>
  class RefinementBuilder<3, CoordType>
  {
    static const std::size_t dimension = 3;
  public:
    static
    VirtualRefinement<dimension, CoordType> &
    build(
      unsigned topologyId
      , unsigned coerceToId
      )
    {
      topologyId &= ~1;
      coerceToId &= ~1;

      const unsigned idSimplex = GenericGeometry::SimplexTopology<dimension>::type::id & ~1;
      const unsigned idCube = GenericGeometry::CubeTopology<dimension>::type::id & ~1;
      const unsigned idPrism = GenericGeometry::PrismTopology<dimension>::type::id & ~1;
      const unsigned idPyramid = GenericGeometry::PyramidTopology<dimension>::type::id & ~1;

      switch( topologyId )
      {
      //case GeometryType::simplex:
      case idSimplex :
        //switch( coerceTo )
        switch( coerceToId )
        {
        //case GeometryType::simplex:
        case idSimplex :
          return VirtualRefinementImp< idSimplex, CoordType, idSimplex, dimension>::instance();
        default :
          break;
        }
        break;

      //case GeometryType::cube:
      case idCube :
        switch( coerceToId )
        {
        case idSimplex :
          return VirtualRefinementImp< idCube, CoordType, idSimplex, dimension>::instance();
        case idCube :
          return VirtualRefinementImp< idCube, CoordType, idCube, dimension>::instance();
        default :
          break;
        }
        break;


      //case GeometryType::prism:
      case idPrism :
        switch( coerceToId )
        {
        case idSimplex :
          return VirtualRefinementImp< idPrism, CoordType, idSimplex, dimension>::instance();
        default :
          break;
        }
        break;

      //case GeometryType::pyramid:
      case idPyramid :
        switch( coerceToId )
        {
        case idSimplex :
          return VirtualRefinementImp< idPyramid, CoordType, idSimplex, dimension>::instance();
        default :
          break;
        }
        break;

      default :
        break;
      }
      DUNE_THROW( NotImplemented, "No Refinement<" << topologyId << ", CoordType, "
                                                   << coerceToId << " >.");
    }
  };

} // namespace Dune

#endif //DUNE_GRID_COMMON_VIRTUALREFINEMENT_CC
