// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GRID_UTILITY_MULTIINDEX_HH
#define DUNE_GRID_UTILITY_MULTIINDEX_HH

/** \file
 *  \brief Implements a multiindex with arbitrary dimension and fixed index ranges
 *  This is used by various factory classes.
 */

#include<array>

namespace Dune
{
 namespace FactoryUtilities
 {
  template<std::size_t dim>
  class MultiIndex : public std::array<unsigned int,dim>
  {
    // The range of each component
    std::array<unsigned int,dim> limits_;

  public:
    /** \brief Constructor with a given range for each digit */
    MultiIndex(const std::array<unsigned int,dim>& limits) : limits_(limits)
    {
      std::fill(this->begin(), this->end(), 0);
    }

    /** \brief Increment the MultiIndex */
    MultiIndex<dim>& operator++()
    {
      for (std::size_t i=0; i<dim; i++)
      {
        // Augment digit
        (*this)[i]++;

        // If there is no carry-over we can stop here
        if ((*this)[i]<limits_[i])
          break;

        (*this)[i] = 0;
      }
      return *this;
    }

    /** \brief Compute how many times you can call operator++ before getting to (0,...,0) again */
    size_t cycle() const
    {
      size_t result = 1;
      for (std::size_t i=0; i<dim; i++)
        result *= limits_[i];
      return result;
    }
  };
 }
}

#endif
