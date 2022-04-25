#ifndef DUNE_GRID_CONCEPTS_HH
#define DUNE_GRID_CONCEPTS_HH

/**
 * This file contains a convenience definition that checks if concepts are available
 * If DUNE_GRID_HAVE_CONCEPTS is true, the dune-grid concepts are available and
 * have been included.
 *
 * In order to enable these concepts, you need the following:
 *  (i)   A C++20 compiler with concepts (i.e. __cpp_concepts >= 201907L)
 *  (ii)  The concepts in the standard library (i.e. __cpp_lib_concepts >= 202002L)
 *  (iii) A C++ compiler that supports unevaluated context lambdas (i.e.
 *        can compile `using = decltype([]{})`). This is automatically tested by
 *        CMake and exposed as a macro definition
 *        `DUNE_HAVE_CXX_UNEVALUATED_CONTEXT_LAMBDA`.
 *
 * - `gcc>=9` and `clang>=12` are known to fulfill (i) and (iii).
 * - `libstdc++>10` and `libc++>=13` are known to fulfill (ii).
 */

// check whether c++20 concept can be used
#if __has_include(<version>) && __has_include(<concepts>)
  #include <version>
  #if  __cpp_concepts >= 201907L && __cpp_lib_concepts >= 202002L && DUNE_HAVE_CXX_UNEVALUATED_CONTEXT_LAMBDA
    #ifndef DUNE_GRID_HAVE_CONCEPTS
    #define DUNE_GRID_HAVE_CONCEPTS 1
    #endif
  #endif
#endif

//! Grid concepts are available
#if DUNE_GRID_HAVE_CONCEPTS

// Include all concept headers
#include <dune/grid/concepts/entity.hh>
#include <dune/grid/concepts/entityiterator.hh>
#include <dune/grid/concepts/intersection.hh>
#include <dune/grid/concepts/intersectioniterator.hh>
#include <dune/grid/concepts/geometry.hh>
#include <dune/grid/concepts/indexset.hh>
#include <dune/grid/concepts/gridview.hh>
#include <dune/grid/concepts/grid.hh>

#endif // DUNE_GRID_HAVE_CONCEPTS

#endif // DUNE_GRID_CONCEPTS_HH
