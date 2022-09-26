// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
/*
 * Test that an implicit copy generator is available for classes derived from
 * - Dune::AnalyticalCoordFunction
 * - Dune::DiscreteCoordFunction
 *
 * See https://gitlab.dune-project.org/flyspray/FS/issues/1463
 */

#include <config.h>

#include <dune/grid/geometrygrid/coordfunction.hh>

class Analytical
  : public Dune::AnalyticalCoordFunction<double, 3, 3, Analytical>
{
  using Base = Dune::AnalyticalCoordFunction<double, 3, 3, Analytical>;
  Base::RangeVector offset_;
public:
  Analytical(Base::RangeVector offset)
    : offset_(offset)
    { /* Nothing. */ }

  void evaluate(const Base::DomainVector& x, Base::RangeVector& y) const
    {
      y = x;
      y += offset_;
    }
};

class Discrete
  : public Dune::DiscreteCoordFunction<double, 3, Discrete>
{
  using Base = Dune::DiscreteCoordFunction<double, 3, Discrete>;
  Base::RangeVector offset_;
public:
  Discrete(Base::RangeVector offset)
    : offset_(offset)
    { /* Nothing. */ }

  template<class HostEntity>
  void evaluate(const HostEntity& hostEntity, unsigned int corner, Base::RangeVector& y) const
    {
      y = hostEntity.geometry().corner(corner);
      y += offset_;
    }
};

int main()
{
  {
    Analytical a({1., 2., 3.});
    Analytical b = a;
    [[maybe_unused]] Analytical c(b);
  }

  {
    Discrete a({1., 2., 3.});
    Discrete b = a;
    [[maybe_unused]] Discrete c(b);
  }

  return 0;
}
