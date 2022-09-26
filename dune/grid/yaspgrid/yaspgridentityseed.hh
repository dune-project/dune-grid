// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_YASPGRIDENTITYSEED_HH
#define DUNE_GRID_YASPGRIDENTITYSEED_HH

/** \file
 * \brief The YaspEntitySeed class
 */

namespace Dune {

  /** \brief Describes the minimal information necessary to create a fully functional YaspEntity
   */
  template<int codim, class GridImp>
  class YaspEntitySeed
  {
    //! know your own dimension
    constexpr static int dim = GridImp::dimension;

  public:
    //! codimension of entity
    constexpr static int codimension = codim;

    //! default construct an invalid entity seed
    YaspEntitySeed ()
      : _l(-1), _o(0)
    {
      std::fill(_c.begin(),_c.end(),0);
    }

    //! constructor
    YaspEntitySeed (int level, std::array<int, dim> coord, int o = 0)
      : _l(level), _c(coord), _o(o)
    {}

    //! check whether the EntitySeed refers to a valid Entity
    bool isValid() const
    {
      return _l != -1;
    }

    int level () const { return _l; }
    const std::array<int, dim> & coord() const { return _c; }
    int offset () const { return _o; }

  protected:
    int _l;                  // grid level
    std::array<int, dim> _c; // coord in the global grid
    int _o; // the offset: which YGridComponent, does the entity belong to
  };

}  // namespace Dune

#endif   // DUNE_GRID_YASPGRIDENTITYSEED_HH
