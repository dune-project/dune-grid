// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

// check whether c++20 concept can be used
#if DUNE_HAVE_CXX_UNEVALUATED_CONTEXT_LAMBDA
  #if __has_include(<version>) && __has_include(<concepts>)
    #include <version>
    #if  __cpp_concepts >= 201907L && __cpp_lib_concepts >= 202002L
      #define DUNE_ENABLE_CONCEPTS 1
    #endif
  #endif
#endif

#if DUNE_ENABLE_CONCEPTS

#include <dune/grid/onedgrid.hh>
#include <dune/grid/geometrygrid.hh>
#include <dune/grid/identitygrid.hh>
#include <dune/grid/yaspgrid.hh>

#if HAVE_ALBERTA
#include <dune/grid/albertagrid.hh>
#endif

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

#include <dune/grid/concepts/grid.hh>

int main ( int argc, char **argv )
{
  Dune::MPIHelper::instance( argc, argv );

  static_assert(Dune::Concept::Grid< Dune::OneDGrid >);

  static_assert(Dune::Concept::Grid< Dune::YaspGrid<1> >);
  static_assert(Dune::Concept::Grid< Dune::YaspGrid<2> >);
  static_assert(Dune::Concept::Grid< Dune::YaspGrid<3> >);

#if HAVE_ALBERTA
  static_assert(Dune::Concept::Grid< Dune::AlbertaGrid<1,1> >);
  static_assert(Dune::Concept::Grid< Dune::AlbertaGrid<1,2> >);
  static_assert(Dune::Concept::Grid< Dune::AlbertaGrid<1,3> >);
  static_assert(Dune::Concept::Grid< Dune::AlbertaGrid<2,2> >);
  static_assert(Dune::Concept::Grid< Dune::AlbertaGrid<2,3> >);
  static_assert(Dune::Concept::Grid< Dune::AlbertaGrid<3,3> >);
#endif

#if HAVE_UG
  static_assert(Dune::Concept::Grid< Dune::UGGrid<1> >);
  static_assert(Dune::Concept::Grid< Dune::UGGrid<2> >);
  static_assert(Dune::Concept::Grid< Dune::UGGrid<3> >);
#endif

  // check grid wrappers
  static_assert(Dune::Concept::Grid< Dune::GeometryGrid< Dune::YaspGrid<1> > >);
  static_assert(Dune::Concept::Grid< Dune::IdentityGrid< Dune::YaspGrid<1> > >);

}

#else // DUNE_ENABLE_CONCEPTS

int main()
{
  return 77;
}

#endif // DUNE_ENABLE_CONCEPTS
