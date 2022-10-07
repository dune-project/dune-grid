// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#include "config.h"

#include <iostream>

#include <dune/common/deprecated.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/grid/yaspgrid.hh>

using Dune::TestSuite;

template<int n>
TestSuite testBinomialTable(std::vector<int> const* reference = nullptr, bool dump = false)
{
  TestSuite t;

  using Table = Dune::Yasp::BinomialTable<n>;

  if (!reference)
    std::cerr << "W: No reference given!\n";

  if (dump)
    std::cout << "    std::vector<int> reference{";

  int i = 0;
  for (int d = 0; d <= n; ++d) {
    for (int c = 0; c <= d; ++c, ++i) {
      const auto value = Table::evaluate(d, c);

      DUNE_NO_DEPRECATED_BEGIN
      t.check(value == Table::binomial(d, c));
      DUNE_NO_DEPRECATED_END
      if (reference)
        t.check(Table::evaluate(d, c) == reference->at(i));

      if (dump)
        std::cout << value << ", ";
    }
  }

  if (dump)
    std::cout << "};\n";

  return t;
}

template<class F, int dim>
TestSuite testEntityShiftTable(std::vector<unsigned long long> const* reference = nullptr, bool dump = false)
{
  TestSuite t;

  using Table = Dune::Yasp::EntityShiftTable<F, dim>;

  if (!reference)
    std::cerr << "W: No reference given!\n";

  if (dump)
    std::cout << "    std::vector<unsigned long long> reference{";

  int i = 0;
  for (int codim = 0; codim <= dim; ++codim) {
    for (int j = 0; j < Dune::Yasp::subEnt<dim>(dim, codim); ++j, ++i) {
      const auto value = Table::evaluate(j, codim);

      t.check(value == F::evaluate(j, codim));
      if (reference)
        t.check(value.to_ullong() == reference->at(i));

      if (dump)
        std::cout << "0b" << value.to_string() << "ull, ";
    }
  }

  if (dump)
    std::cout << "};\n";

  return t;
}

int main(int argc, char** /* argv */)
{
  bool dump = false;
  if (argc > 1)
    dump = true;

  TestSuite t;

  {
    std::vector<int> reference{1, 1, 1};
    t.subTest(testBinomialTable<1>(&reference, dump));
  }
  {
    std::vector<int> reference{1, 1, 1, 1, 2, 1};
    t.subTest(testBinomialTable<2>(&reference, dump));
  }
  {
    std::vector<int> reference{1, 1, 1, 1, 2, 1, 1, 3, 3, 1};
    t.subTest(testBinomialTable<3>(&reference, dump));
  }
  {
    std::vector<int> reference{1, 1, 1, 1, 2, 1, 1, 3, 3, 1, 1, 4, 6, 4, 1};
    t.subTest(testBinomialTable<4>(&reference, dump));
  }

  {
    std::vector<unsigned long long> reference{0b1ull, 0b0ull, 0b0ull};
    t.subTest(testEntityShiftTable<Dune::Yasp::calculate_entity_shift<1>, 1>(&reference, dump));
  }
  {
    std::vector<unsigned long long> reference{
      0b11ull, 0b10ull, 0b10ull, 0b01ull, 0b01ull,
      0b00ull, 0b00ull, 0b00ull, 0b00ull
    };
    t.subTest(testEntityShiftTable<Dune::Yasp::calculate_entity_shift<2>, 2>(&reference, dump));
  }
  {
    std::vector<unsigned long long> reference{
      0b111ull, 0b110ull, 0b110ull, 0b101ull, 0b101ull, 0b011ull, 0b011ull,
      0b100ull, 0b100ull, 0b100ull, 0b100ull, 0b010ull, 0b010ull, 0b001ull,
      0b001ull, 0b010ull, 0b010ull, 0b001ull, 0b001ull, 0b000ull, 0b000ull,
      0b000ull, 0b000ull, 0b000ull, 0b000ull, 0b000ull, 0b000ull
    };
    t.subTest(testEntityShiftTable<Dune::Yasp::calculate_entity_shift<3>, 3>(&reference, dump));
  }
  {
    std::vector<unsigned long long> reference{
      0b1111ull, 0b1110ull, 0b1110ull, 0b1101ull, 0b1101ull, 0b1011ull,
      0b1011ull, 0b0111ull, 0b0111ull, 0b1100ull, 0b1100ull, 0b1100ull,
      0b1100ull, 0b1010ull, 0b1010ull, 0b1001ull, 0b1001ull, 0b1010ull,
      0b1010ull, 0b1001ull, 0b1001ull, 0b0110ull, 0b0110ull, 0b0101ull,
      0b0101ull, 0b0011ull, 0b0011ull, 0b0110ull, 0b0110ull, 0b0101ull,
      0b0101ull, 0b0011ull, 0b0011ull, 0b1000ull, 0b1000ull, 0b1000ull,
      0b1000ull, 0b1000ull, 0b1000ull, 0b1000ull, 0b1000ull, 0b0100ull,
      0b0100ull, 0b0100ull, 0b0100ull, 0b0010ull, 0b0010ull, 0b0001ull,
      0b0001ull, 0b0010ull, 0b0010ull, 0b0001ull, 0b0001ull, 0b0100ull,
      0b0100ull, 0b0100ull, 0b0100ull, 0b0010ull, 0b0010ull, 0b0001ull,
      0b0001ull, 0b0010ull, 0b0010ull, 0b0001ull, 0b0001ull, 0b0000ull,
      0b0000ull, 0b0000ull, 0b0000ull, 0b0000ull, 0b0000ull, 0b0000ull,
      0b0000ull, 0b0000ull, 0b0000ull, 0b0000ull, 0b0000ull, 0b0000ull,
      0b0000ull, 0b0000ull, 0b0000ull
    };
    t.subTest(testEntityShiftTable<Dune::Yasp::calculate_entity_shift<4>, 4>(&reference, dump));
  }

  {
    std::vector<unsigned long long> reference{0b0ull, 0b0ull, 0b1ull};
    t.subTest(testEntityShiftTable<Dune::Yasp::calculate_entity_move<1>, 1>(&reference, dump));
  }
  {
    std::vector<unsigned long long> reference{
      0b00ull, 0b00ull, 0b01ull, 0b00ull, 0b10ull,
      0b00ull, 0b01ull, 0b10ull, 0b11ull
    };
    t.subTest(testEntityShiftTable<Dune::Yasp::calculate_entity_move<2>, 2>(&reference, dump));
  }
  {
    std::vector<unsigned long long> reference{
      0b000ull, 0b000ull, 0b001ull, 0b000ull, 0b010ull, 0b000ull, 0b100ull,
      0b000ull, 0b001ull, 0b010ull, 0b011ull, 0b000ull, 0b001ull, 0b000ull,
      0b010ull, 0b100ull, 0b101ull, 0b100ull, 0b110ull, 0b000ull, 0b001ull,
      0b010ull, 0b011ull, 0b100ull, 0b101ull, 0b110ull, 0b111ull
    };
    t.subTest(testEntityShiftTable<Dune::Yasp::calculate_entity_move<3>, 3>(&reference, dump));
  }
  {
    std::vector<unsigned long long> reference{
      0b0000ull, 0b0000ull, 0b0001ull, 0b0000ull, 0b0010ull, 0b0000ull,
      0b0100ull, 0b0000ull, 0b1000ull, 0b0000ull, 0b0001ull, 0b0010ull,
      0b0011ull, 0b0000ull, 0b0001ull, 0b0000ull, 0b0010ull, 0b0100ull,
      0b0101ull, 0b0100ull, 0b0110ull, 0b0000ull, 0b0001ull, 0b0000ull,
      0b0010ull, 0b0000ull, 0b0100ull, 0b1000ull, 0b1001ull, 0b1000ull,
      0b1010ull, 0b1000ull, 0b1100ull, 0b0000ull, 0b0001ull, 0b0010ull,
      0b0011ull, 0b0100ull, 0b0101ull, 0b0110ull, 0b0111ull, 0b0000ull,
      0b0001ull, 0b0010ull, 0b0011ull, 0b0000ull, 0b0001ull, 0b0000ull,
      0b0010ull, 0b0100ull, 0b0101ull, 0b0100ull, 0b0110ull, 0b1000ull,
      0b1001ull, 0b1010ull, 0b1011ull, 0b1000ull, 0b1001ull, 0b1000ull,
      0b1010ull, 0b1100ull, 0b1101ull, 0b1100ull, 0b1110ull, 0b0000ull,
      0b0001ull, 0b0010ull, 0b0011ull, 0b0100ull, 0b0101ull, 0b0110ull,
      0b0111ull, 0b1000ull, 0b1001ull, 0b1010ull, 0b1011ull, 0b1100ull,
      0b1101ull, 0b1110ull, 0b1111ull
    };
    t.subTest(testEntityShiftTable<Dune::Yasp::calculate_entity_move<4>, 4>(&reference, dump));
  }

  return t.exit();
}
