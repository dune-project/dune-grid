// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_TWISTUTILITY_HH
#define DUNE_TWISTUTILITY_HH

// is Alberta was found then also include headers
#ifndef HAVE_ALBERTA
#define HAVE_ALBERTA_FOUND 0
#else
#define HAVE_ALBERTA_FOUND HAVE_ALBERTA
#endif

// is ALU3dGrid was found then also include headers
#ifndef HAVE_ALUGRID
#define HAVE_ALUGRID_FOUND 0
#else
#define HAVE_ALUGRID_FOUND HAVE_ALUGRID
#endif

#if HAVE_ALUGRID_FOUND
#include <dune/grid/alugrid.hh>
#endif

#if HAVE_ALBERTA_FOUND
#include <dune/grid/albertagrid.hh>
#endif

#include <dune/grid/sgrid.hh>
#include <dune/grid/yaspgrid.hh>

namespace Dune {

  template <class GridImp>
  class TwistUtility;
  // for structured grids, the twist is always zero
  // ? is this correct
  template <int dim,int dimw>
  class TwistUtility<SGrid<dim,dimw> >
  {
  public:
    typedef SGrid<dim,dimw> GridType;
  public:
    TwistUtility(const GridType& grid) :
      grid_(grid)
    {}

    // default twist is zero
    template <class IntersectionIterator>
    int twistInSelf(const IntersectionIterator& it) const {
      return 0;
    }

    // default twist is zero
    template <class IntersectionIterator>
    int twistInNeighbor(const IntersectionIterator& it) const {
      return 0;
    }

  private:
    const GridType& grid_;
  };
  template <int dim,int dimw>
  class TwistUtility<YaspGrid<dim,dimw> >
  {
  public:
    typedef YaspGrid<dim,dimw> GridType;
  public:
    TwistUtility(const GridType& grid) :
      grid_(grid)
    {}

    // default twist is zero
    template <class IntersectionIterator>
    int twistInSelf(const IntersectionIterator& it) const {
      return 0;
    }

    // default twist is zero
    template <class IntersectionIterator>
    int twistInNeighbor(const IntersectionIterator& it) const {
      return 0;
    }

  private:
    const GridType& grid_;
  };
#if HAVE_ALBERTA_FOUND
  template <int dim, int dimW>
  class TwistUtility<AlbertaGrid<dim, dimW> >
  {
  public:
    typedef AlbertaGrid<dim, dimW> GridType;
    typedef typename GridType::Traits::LeafIntersectionIterator LeafIntersectionIterator;
    typedef typename GridType::Traits::LevelIntersectionIterator LevelIntersectionIterator;
  public:
    TwistUtility(const GridType& grid) :
      grid_(grid)
    {}

    int twistInSelf(const LeafIntersectionIterator& it) const {
      return grid_.getRealIntersectionIterator(it).twistInSelf();
    }

    //int twistInSelf(const LevelIntersectionIterator& it) const {
    //  return grid_.getRealIntersectionIterator(it).twistInSelf();
    //}

    int twistInNeighbor(const LeafIntersectionIterator& it) const {
      return grid_.getRealIntersectionIterator(it).twistInNeighbor();
    }

    //int twistInNeighbor(const LevelIntersectionIterator& it) const {
    //  return grid_.getRealIntersectionIterator(it).twistInNeighbor();
    //}

  private:
    const GridType& grid_;
  };
#endif

#if HAVE_ALUGRID_FOUND
  template <>
  class TwistUtility<ALUSimplexGrid<3,3>  >
  {
  public:
    typedef ALUSimplexGrid<3,3> GridType;
    typedef GridType::Traits::LeafIntersectionIterator LeafIntersectionIterator;
    typedef GridType::Traits::LevelIntersectionIterator LevelIntersectionIterator;
  public:
    TwistUtility(const GridType& grid) :
      grid_(grid)
    {}

    int twistInSelf(const LeafIntersectionIterator& it) const {
      return grid_.getRealIntersectionIterator(it).twistInSelf();
    }

    int twistInSelf(const LevelIntersectionIterator& it) const {
      return grid_.getRealIntersectionIterator(it).twistInSelf();
    }

    int twistInNeighbor(const LeafIntersectionIterator& it) const {
      return grid_.getRealIntersectionIterator(it).twistInNeighbor();
    }

    int twistInNeighbor(const LevelIntersectionIterator& it) const {
      return grid_.getRealIntersectionIterator(it).twistInNeighbor();
    }

  private:
    TwistUtility(const TwistUtility&);
    TwistUtility& operator=(const TwistUtility&);

  private:
    const GridType& grid_;
  };
  template <>
  class TwistUtility<ALUCubeGrid<3,3>  >
  {
  public:
    typedef ALUCubeGrid<3,3> GridType;
    typedef GridType::Traits::LeafIntersectionIterator LeafIntersectionIterator;
    typedef GridType::Traits::LevelIntersectionIterator LevelIntersectionIterator;
  public:
    TwistUtility(const GridType& grid) :
      grid_(grid)
    {}

    int twistInSelf(const LeafIntersectionIterator& it) const {
      return grid_.getRealIntersectionIterator(it).twistInSelf();
    }

    int twistInSelf(const LevelIntersectionIterator& it) const {
      return grid_.getRealIntersectionIterator(it).twistInSelf();
    }

    int twistInNeighbor(const LeafIntersectionIterator& it) const {
      return grid_.getRealIntersectionIterator(it).twistInNeighbor();
    }

    int twistInNeighbor(const LevelIntersectionIterator& it) const {
      return grid_.getRealIntersectionIterator(it).twistInNeighbor();
    }

  private:
    TwistUtility(const TwistUtility&);
    TwistUtility& operator=(const TwistUtility&);

  private:
    const GridType& grid_;
  };
  template <>
  class TwistUtility<ALUSimplexGrid<2,2>  >
  {
  public:
    typedef ALUSimplexGrid<2, 2> GridType;
    typedef GridType::Traits::LeafIntersectionIterator LeafIntersectionIterator;
    typedef GridType::Traits::LevelIntersectionIterator LevelIntersectionIterator;
  public:
    TwistUtility(const GridType& grid) :
      grid_(grid)
    {}

    int twistInSelf(const LeafIntersectionIterator& it) const {
      return 0;
    }

    int twistInSelf(const LevelIntersectionIterator& it) const {
      return 0;
    }

    int twistInNeighbor(const LeafIntersectionIterator& it) const {
      return 1;
    }

    int twistInNeighbor(const LevelIntersectionIterator& it) const {
      return 1;
    }

  private:
    TwistUtility(const TwistUtility&);
    TwistUtility& operator=(const TwistUtility&);

  private:
    const GridType& grid_;
  };
#endif

#undef HAVE_ALBERTA_FOUND
#undef HAVE_ALUGRID_FOUND
} // end namespace Dune

#endif
