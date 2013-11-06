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
    enum { dim=GridImp::dimension };

  public:
    //! codimension of entity pointer
    enum { codimension = codim };

    //! constructor
    YaspEntitySeed (int level, Dune::array<int, dim> coord)
      : _l(level), _c(coord)
    {}

    //! copy constructor
    YaspEntitySeed (const YaspEntitySeed& rhs)
      : _l(rhs._l), _c(rhs._c)
    {}

    int level () const { return _l; }
    const Dune::array<int, dim> & coord() const { return _c; }

  protected:
    int _l;                  // grid level
    Dune::array<int, dim> _c; // coord in the global grid
  };

}  // namespace Dune

#endif   // DUNE_GRID_YASPGRIDENTITYSEED_HH
