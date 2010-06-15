// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_IO_FILE_VTK_SKELETONFUNCTION_HH
#define DUNE_GRID_IO_FILE_VTK_SKELETONFUNCTION_HH

#include <vector>

#include <dune/common/fvector.hh>

namespace Dune {

  //! \addtogroup VTK
  //! \{

  /** @file
      @author Jö Fahlke
      @brief Functions for VTK output on the skeleton
   */

  namespace VTK {

    //////////////////////////////////////////////////////////////////////
    //
    //  Prototype for VTKFunktions on the skeleton
    //

    template<typename GV, typename RF>
    struct SkeletonFunctionTraits {
      typedef GV GridView;
      typedef typename GV::Intersection Cell;

      typedef typename GV::ctype DomainField;
      static const unsigned dimDomain = GV::dimension-1;
      typedef FieldVector<DomainField, dimDomain> Domain;

      typedef RF RangeField;
      typedef std::vector<RangeField> Range;
    };

    //! A prototype for VTKFunctions on the skeleton
    template <typename GV, typename RF>
    class SkeletonFunctionInterface {
    public:
      typedef SkeletonFunctionTraits<GV, RF> Traits;

      //! get dimension of the Range
      unsigned dimRange() const;

      //! evaluate at local point xl in Cell c, store in result
      /**
       * \param c      The cell (intersection) to evaluate in.
       * \param xl     The local coordinate within the cell.
       * \param result Where to store the result.  The vector is resized as
       *               nessecary, and is overwritten.
       */
      void evaluate(const typename Traits::Cell& c,
                    const typename Traits::Domain& xl,
                    typename Traits::Range& result) const;
    };

  } // namespace VTK

  //! \} group VTK

} // namespace Dune

#endif // DUNE_GRID_IO_FILE_VTK_SKELETONFUNCTION_HH
