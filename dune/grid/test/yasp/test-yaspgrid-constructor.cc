// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#include "config.h"

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/yaspgrid.hh>

#include <array>
#include <type_traits>

int main(int argc, char** argv)
{
#if __cpp_deduction_guides >= 201611
  using namespace Dune;
  MPIHelper::instance(argc, argv);

  using DVector = FieldVector<double, 2>;
  using FVector = FieldVector<float, 2>;
  std::bitset<2> p{0ULL};

  {
    DVector x{1.0, 1.0};
    std::array N{2, 2};
    YaspGrid grid{x, N, p};
    static_assert(std::is_same_v< decltype(grid), YaspGrid< 2, EquidistantCoordinates<double, 2> > >);
  }
  {
    FVector x{1.0f, 1.0f};
    std::array N{2, 2};
    YaspGrid grid{x, N};
    static_assert(std::is_same_v< decltype(grid), YaspGrid< 2, EquidistantCoordinates<float, 2> > >);
  }

  {
    DVector x1{0.0, 0.0}, x2{1.0, 1.0};
    std::array N{2, 2};
    YaspGrid grid{x1, x2, N};
    static_assert(std::is_same_v< decltype(grid), YaspGrid< 2, EquidistantOffsetCoordinates<double, 2> > >);
  }
  {
    FVector x1{0.0f, 0.0f}, x2{1.0f, 1.0f};
    std::array N{2, 2};
    YaspGrid grid{x1, x2, N};
    static_assert(std::is_same_v< decltype(grid), YaspGrid< 2, EquidistantOffsetCoordinates<float, 2> > >);
  }

  {
    std::vector v{0.0, 1.0, 2.0};
    std::array vs{v, v};
    YaspGrid grid{vs};
    static_assert(std::is_same_v< decltype(grid), YaspGrid< 2, TensorProductCoordinates<double, 2> > >);
  }
  {
    std::vector v{0.0f, 1.0f, 2.0f};
    std::array vs{v, v};
    YaspGrid grid{vs};
    static_assert(std::is_same_v< decltype(grid), YaspGrid< 2, TensorProductCoordinates<float, 2> > >);
  }

  return 0;
#else
  std::cerr << "This test requires full support for C++17's class template argument deduction.\n";
  return 77;
#endif
}
