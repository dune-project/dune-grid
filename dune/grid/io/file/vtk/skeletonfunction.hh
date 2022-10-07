// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_IO_FILE_VTK_SKELETONFUNCTION_HH
#define DUNE_GRID_IO_FILE_VTK_SKELETONFUNCTION_HH

#include <memory>
#include <string>
#include <vector>

#include <dune/common/fvector.hh>

#include <dune/grid/io/file/vtk/functionwriter.hh>
#include <dune/grid/io/file/vtk/pvtuwriter.hh>
#include <dune/grid/io/file/vtk/vtuwriter.hh>

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
       *               necessary, and is overwritten.
       */
      void evaluate(const typename Traits::Cell& c,
                    const typename Traits::Domain& xl,
                    typename Traits::Range& result) const;
    };

    ////////////////////////////////////////////////////////////////////////
    //
    //  Class for writing SkeletonFunctions
    //

    //! function writer for skeleton functions
    /**
     * \tparam Func Function to write.  Must be a model of
     *              SkeletonFunctionInterface.
     */
    template<typename Func>
    class SkeletonFunctionWriter
      : public FunctionWriterBase<typename Func::Traits::Cell>
    {
      typedef typename Func::Traits::RangeField RF;

      std::shared_ptr<const Func> func;
      std::string name_;
      unsigned dimR;
      VTK::Precision precision_;
      std::shared_ptr<DataArrayWriter> arraywriter;

    public:
      SkeletonFunctionWriter(const std::shared_ptr<const Func>& func_,
                             const std::string& name, unsigned dimR_,
                             VTK::Precision prec = VTK::Precision::float32)
        : func(func_), name_(name), dimR(dimR_), precision_(prec)
      { }

      SkeletonFunctionWriter(const std::shared_ptr<const Func>& func_,
                             const std::string& name,
                             VTK::Precision prec = VTK::Precision::float32)
        : func(func_), name_(name), dimR(func->dimRange()), precision_(prec)
      { }

      //! return name
      virtual std::string name() const { return name_; }

      //! return number of components of the vector
      virtual unsigned ncomps() const { return dimR; }

      //! add this field to the given parallel writer
      virtual void addArray(PVTUWriter& writer) {
        writer.addArray(name(), ncomps(), precision_);
      }

      //! start writing with the given writer
      virtual bool beginWrite(VTUWriter& writer, std::size_t nitems) {
        arraywriter.reset(writer.makeArrayWriter(name(), ncomps(),
                                                 nitems, precision_));
        return !arraywriter->writeIsNoop();
      }

      //! write at the given position
      virtual void write(const typename Func::Traits::Cell& cell,
                         const typename Func::Traits::Domain& xl) {
        typename Func::Traits::Range result;
        func->evaluate(cell, xl, result);
        for(unsigned d = 0; d < result.size() && d < dimR; ++d)
          arraywriter->write(result[d]);
        for(unsigned d = result.size(); d < dimR; ++d)
          arraywriter->write(0);
      }

      //! signal end of writing
      virtual void endWrite() {
        arraywriter.reset();
      }
    };

  } // namespace VTK

  //! \} group VTK

} // namespace Dune

#endif // DUNE_GRID_IO_FILE_VTK_SKELETONFUNCTION_HH
