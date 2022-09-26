// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_TEST_CHECKCOMMCORRECTNESS_HH
#define DUNE_GRID_TEST_CHECKCOMMCORRECTNESS_HH

#include <config.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_set>
#include <unordered_map>
#include <unistd.h>

#include <dune/common/hash.hh>
#include <dune/common/float_cmp.hh>
#include <dune/geometry/dimension.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/mcmgmapper.hh>

namespace Dune {

  std::size_t hash_value(const PartitionType& pt)
  {
    return static_cast<std::size_t>(pt);
  }

}

DUNE_DEFINE_HASH(DUNE_HASH_TEMPLATE_ARGS(),DUNE_HASH_TYPE(Dune::PartitionType))

namespace Dune {

  namespace GridCheck {

      struct SymmetryVerifyingDataHandle
        : public Dune::CommDataHandleIF<SymmetryVerifyingDataHandle,std::size_t>
      {

        bool contains(int /* dim */, int codim) const
        {
          return codim == _codim;
        }

        bool fixedSize(int /* dim */, int /* codim */) const
        {
          return _fixed_size;
        }

        template<typename E>
        std::size_t size(const E&) const
        {
          return _writes_per_cell;
        }

        template<typename Buf, typename E>
        void gather(Buf& buf, const E&) const
        {
          for (std::size_t i = 0; i < _writes_per_cell; ++i)
            buf.write(std::size_t(42));
          _writes += _writes_per_cell;
        }

        template<typename Buf, typename E>
        void scatter(Buf& buf, const E&, std::size_t n) const
        {
          assert(_writes_per_cell == n);
          for (std::size_t i = 0; i < _writes_per_cell; ++i)
            {
              std::size_t tmp = 0;
              buf.read(tmp);
              assert(tmp == 42);
            }
          _reads += _writes_per_cell;
        }

        SymmetryVerifyingDataHandle(int codim, bool fixed_size, std::size_t writes_per_cell)
          : _codim(codim)
          , _fixed_size(fixed_size)
          , _writes_per_cell(writes_per_cell)
          , _reads(0)
          , _writes(0)
        {}

        const int _codim;
        const bool _fixed_size;
        const std::size_t _writes_per_cell;
        mutable std::size_t _reads;
        mutable std::size_t _writes;

      };

      struct CodimLayout
      {

        bool operator()(Dune::GeometryType gt, int dim) const
        {
          return (static_cast< int >( gt.dim() ) == dim - _codim);
        }

        int _codim;

      };

      template<typename GV>
      struct CommunicationTestDataHandle
        : public Dune::CommDataHandleIF<CommunicationTestDataHandle<GV>,
                                  typename GV::template Codim<0>::Geometry::GlobalCoordinate>
      {
        using ctype = typename GV::ctype;

        bool contains([[maybe_unused]] int dim, int codim) const
        {
          return codim == _codim;
        }

        bool fixedSize([[maybe_unused]] int dim, [[maybe_unused]] int codim) const
        {
          return true;
        }

        template<typename E>
        std::size_t size([[maybe_unused]] const E& e) const
        {
          return 1;
        }

        template<typename Buf, typename E>
        void gather(Buf& buf, const E& e) const
        {
          assert(_allowed_writes.count(e.partitionType()) > 0);
          auto center = e.geometry().center();
          buf.write(e.geometry().center());
          if (_gv.comm().rank() == 0)
            std::cout << "Gathering from entity " << _mapper.index(e) << " at " << e.geometry().center() << std::endl;
          ++_writes[_mapper.index(e)];
          center -= _coords[_mapper.index(e)];
          assert(Dune::FloatCmp::eq(center.two_norm(),ctype(0)));
        }

        template<typename Buf, typename E>
        void scatter(Buf& buf, const E& e, [[maybe_unused]] std::size_t n) const
        {
          assert(_allowed_reads.count(e.partitionType()) > 0);
          typename E::Geometry::GlobalCoordinate data;
          buf.read(data);
          ++_reads[_mapper.index(e)];
          auto center = e.geometry().center();
          data -= e.geometry().center();
          assert(Dune::FloatCmp::eq(data.two_norm(),ctype(0)));
          center -= _coords[_mapper.index(e)];
          assert(Dune::FloatCmp::eq(center.two_norm(),ctype(0)));
        }

        template<int cd, typename Check>
        void verify(Dune::Codim<cd> codim, Check check)
        {
          assert(codim == _codim);
          for (const auto& e : entities(_gv,codim))
            {
              if (_gv.comm().rank() == 0)
                std::cout << "Entity of codim " << cd
                          << " at " << e.geometry().center()
                          << " with partition " << PartitionName(e.partitionType())
                          << " and " << _writes[_mapper.index(e)] << " / " << _reads[_mapper.index(e)] << " writes / reads"
                          << std::endl;
              check(
                e.partitionType(),
                _allowed_writes.count(e.partitionType()) > 0,
                _allowed_reads.count(e.partitionType()) > 0,
                _writes[_mapper.index(e)],
                _reads[_mapper.index(e)]
                );
            }
        }

        CommunicationTestDataHandle(GV gv, int codim, const std::unordered_set<Dune::PartitionType>& allowed_writes, const std::unordered_set<Dune::PartitionType>& allowed_reads, const std::vector<typename GV::template Codim<0>::Geometry::GlobalCoordinate>& coords)
          : _gv(gv)
          , _codim(codim)
          , _mapper(gv, CodimLayout{_codim})
          , _allowed_writes(allowed_writes)
          , _allowed_reads(allowed_reads)
          , _reads(_mapper.size(),0)
          , _writes(_mapper.size(),0)
          , _coords(coords)
        {}

        GV _gv;
        const int _codim;
        Dune::MultipleCodimMultipleGeomTypeMapper<GV> _mapper;
        const std::unordered_set<Dune::PartitionType> _allowed_writes;
        const std::unordered_set<Dune::PartitionType> _allowed_reads;
        mutable std::vector<std::size_t> _reads;
        mutable std::vector<std::size_t> _writes;
        const std::vector<typename GV::template Codim<0>::Geometry::GlobalCoordinate>& _coords;

      };

      template<typename GV, int cd>
      void check_communication_correctness_do(GV gv, Codim<cd> codim)
      {
        if (gv.grid().comm().rank() == 0)
          {
            std::cout << "Checking codim " << cd << std::endl;
          }

        std::unordered_map<Dune::PartitionType,std::size_t> count;

        Dune::MultipleCodimMultipleGeomTypeMapper<GV> mapper(gv, CodimLayout{codim});

        std::vector<
          typename GV::template Codim<0>::Geometry::GlobalCoordinate
          > coords(mapper.size());

        // start by counting entities by partition type and storing entity positions
        for (const auto& e : entities(gv,codim))
          {
            ++count[e.partitionType()];
            coords[mapper.index(e)] = e.geometry().center();
          }

        {
          SymmetryVerifyingDataHandle dh_forward(cd,false,3);
          gv.communicate(dh_forward,InteriorBorder_InteriorBorder_Interface,ForwardCommunication);

          SymmetryVerifyingDataHandle dh_backward(cd,true,3);
          gv.communicate(dh_backward,InteriorBorder_InteriorBorder_Interface,BackwardCommunication);

          if (count[BorderEntity] > 0)
            {
              assert(dh_forward._writes == dh_forward._reads);
              assert(dh_backward._writes == dh_backward._reads);

              assert(dh_forward._writes == dh_backward._writes);
              assert(dh_backward._writes == dh_backward._writes);
            }

          if (gv.grid().comm().size() == 2)
            {
              if (gv.comm().rank() == 0)
                std::cout << "MPI size == 2, checking writes (" << (dh_forward._writes / 3)
                          << ") against count of border entities (" <<  count[BorderEntity] << ")"
                          << std::endl;
              assert((dh_forward._writes / 3) == count[BorderEntity]);
            }

        }

        {
          SymmetryVerifyingDataHandle dh_forward(cd,true,7);
          gv.communicate(dh_forward,All_All_Interface,ForwardCommunication);

          SymmetryVerifyingDataHandle dh_backward(cd,false,7);
          gv.communicate(dh_backward,All_All_Interface,BackwardCommunication);

          if (count[BorderEntity] > 0)
            {
              assert(dh_forward._writes == dh_forward._reads);
              assert(dh_backward._writes == dh_backward._reads);

              assert(dh_forward._writes == dh_backward._writes);
              assert(dh_backward._writes == dh_backward._writes);
            }

        }

        using PTSet = std::unordered_set<PartitionType>;

        {

          PTSet writers({InteriorEntity,BorderEntity});
          PTSet readers({InteriorEntity,BorderEntity,OverlapEntity,FrontEntity,GhostEntity});

          CommunicationTestDataHandle<GV> dh(gv,codim,writers,readers,coords);
          gv.communicate(dh,InteriorBorder_All_Interface,ForwardCommunication);

        }

        {

          PTSet writers({InteriorEntity,BorderEntity,OverlapEntity,FrontEntity,GhostEntity});
          PTSet readers({InteriorEntity,BorderEntity});

          CommunicationTestDataHandle<GV> dh(gv,codim,writers,readers,coords);
          gv.communicate(dh,InteriorBorder_All_Interface,BackwardCommunication);

        }

        {

          PTSet writers({InteriorEntity,BorderEntity});
          PTSet readers({InteriorEntity,BorderEntity});

          CommunicationTestDataHandle<GV> dh(gv,codim,writers,readers,coords);
          gv.communicate(dh,InteriorBorder_InteriorBorder_Interface,ForwardCommunication);
          dh.verify(
            codim,
            [](PartitionType partition, bool /* write_allowed */, bool /* read_allowed */, [[maybe_unused]] std::size_t writes, [[maybe_unused]] std::size_t reads)
            {
              if (partition == BorderEntity)
                {
                  assert(writes > 0);
                  assert(reads > 0);
                  assert(writes == reads);
                }
              else
                {
                  assert(writes == 0);
                  assert(reads == 0);
                }
            });
        }
      }

      // Need a forward declaration here
      template<typename GV, int cd>
      void check_communication_correctness_iter(GV gv, Codim<cd> codim);

      // Statically iterate over all codimensions
      template<typename GV, int cd>
      void check_communication_correctness_iter(GV, Codim<cd>, std::true_type, std::true_type) {}
      template<typename GV, int cd>
      void check_communication_correctness_iter(GV, Codim<cd>, std::false_type, std::true_type) {}

      template<typename GV, int cd>
      void check_communication_correctness_iter(GV gv, Codim<cd> codim, std::true_type, std::false_type) {
        check_communication_correctness_do(gv, codim);
        check_communication_correctness_iter (gv, Codim<cd+1>());
      }
      template<typename GV, int cd>
      void check_communication_correctness_iter(GV gv, Codim<cd> codim, std::false_type, std::false_type) {
        check_communication_correctness_iter (gv, Codim<cd+1>());
      }

      template<typename GV, int cd>
      void check_communication_correctness_iter(GV gv, Codim<cd> codim) {
        check_communication_correctness_iter (gv, codim,
              std::integral_constant<bool, Dune::Capabilities::hasEntity<typename GV::Grid, cd>::v>(),
              std::integral_constant<bool, cd == GV::dimension + 1>());
      }

      // Start with codim 0
      template<typename GV>
      void check_communication_correctness(GV gv)
      {
        check_communication_correctness_iter (gv, Codim<0>());
      }

  } // namespace GridCheck
} // namespace Dune

#endif
