// SPDX-FileCopyrightText: Copyright (c) DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <config.h>

#include <array>
#include <cassert>
#include <iostream>

#include <dune/common/filledarray.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/common/exceptions.hh>
#include <dune/grid/yaspgrid/partitioning.hh>

template <int d>
void old_optimize_dims (int i, const std::array<int,d>& size, int P, std::array<int,d>& dims,
                        std::array<int,d>& trydims, double& opt)
{
  if (i>0) // test all subdivisions recursively
  {
    for (int k=1; k<=P; k++) {
      if (P%k==0) {
        // P divisible by k
        trydims[i] = k;
        old_optimize_dims<d>(i-1,size,P/k,dims,trydims,opt);
      }
    }
  }
  else
  {
    // found a possible combination
    trydims[0] = P;

    // check for optimality
    double m = -1.0;

    for (int k=0; k<d; k++) {
      double mm=((double)size[k])/((double)trydims[k]);
      if (fmod((double)size[k],(double)trydims[k])>0.0001) mm*=3;
      if ( mm > m ) m = mm;
    }
    if (m<opt) {
      opt = m;
      dims = trydims;
    }
  }
}

template <std::size_t d>
void old_partition (const std::array<int,d>& size, int P, std::array<int,d>& dims, int = 0)
{
  double opt=1E100;
  std::array<int,d> trydims;
  old_optimize_dims<d>(d-1,size,P,dims,trydims,opt);
}


template <std::size_t d>
void run_test (const Dune::Yasp::Partitioning<d>& partitioner, const std::array<int,d>& size)
{
  int maxP = d == 1 ? size[0] :
             d == 2 ? size[0]*size[1] :
             d == 3 ? size[0]*size[1]*size[2] : 0;
  int maxO = d == 1 ? size[0]/2 :
             d == 2 ? std::max({size[0], size[1]})/2 :
             d == 3 ? std::max({size[0], size[1], size[2]})/2 : 0;

  for (int o = 0; o <= maxO; ++o) {
    for (int p = 1; p <= maxP; ++p) {
      std::array<int,d> dims;
      bool failed = false;
      try {
        partitioner.partition(size,p,dims,o);
      } catch(...) {
        failed = true;
      }

      if (!failed) {
        // Check that all successful computed dimensions are positive
        for (std::size_t i = 0; i < d; ++i)
          assert(dims[i] > 0);

        if (o == 0) {
          // Check that we get the same result as in the old implementation for overlap == 0
          std::array<int,d> old;
          old_partition(size,p,old);
          assert(dims == old);
        }
      }
    } // end p
  } // end o
}

template <std::size_t d, class Partitioner>
void test (const Partitioner& partitioner)
{
  int maxS = 4;
  if constexpr (d == 1) {
    for (int s0 = 1; s0 <= maxS; ++s0)
      run_test<1>(partitioner,{s0});
  } else if constexpr (d == 2) {
    for (int s0 = 1; s0 <= maxS; ++s0)
      for (int s1 = 1; s1 <= maxS; ++s1)
        run_test<2>(partitioner,{s0,s1});
  } else if constexpr (d == 3) {
    for (int s0 = 1; s0 <= maxS; ++s0)
      for (int s1 = 1; s1 <= maxS; ++s1)
        for (int s2 = 1; s2 <= maxS; ++s2)
          run_test<3>(partitioner,{s0,s1,s2});
  }
}


int main (int argc , char **argv)
{
  // Initialize MPI, if present
  Dune::MPIHelper::instance(argc, argv);

  Dune::Yasp::DefaultPartitioning<1> p1; test<1>(p1);
  Dune::Yasp::DefaultPartitioning<2> p2; test<2>(p2);
  Dune::Yasp::DefaultPartitioning<3> p3; test<3>(p3);

  // Test construction of partitioners
  {
    Dune::Yasp::PowerDPartitioning<1> ylbp1;
    Dune::Yasp::PowerDPartitioning<2> ylbp2;
    Dune::Yasp::PowerDPartitioning<3> ylbp3;

    Dune::Yasp::FixedSizePartitioning<1> yfsp1(std::array<int,1>{1});
    Dune::Yasp::FixedSizePartitioning<2> yfsp2(std::array<int,2>{1,1});
    Dune::Yasp::FixedSizePartitioning<3> yfsp3(std::array<int,3>{1,1,1});
  }

DUNE_NO_DEPRECATED_BEGIN
  // deprecated partitioners
  {
    Dune::YLoadBalanceDefault<1> ylbd1;
    Dune::YLoadBalanceDefault<2> ylbd2;
    Dune::YLoadBalanceDefault<3> ylbd3;

    Dune::YLoadBalancePowerD<1> ylbp1;
    Dune::YLoadBalancePowerD<2> ylbp2;
    Dune::YLoadBalancePowerD<3> ylbp3;

    Dune::YaspFixedSizePartitioner<1> yfsp1(std::array<int,1>{1});
    Dune::YaspFixedSizePartitioner<2> yfsp2(std::array<int,2>{1,1});
    Dune::YaspFixedSizePartitioner<3> yfsp3(std::array<int,3>{1,1,1});
  }
DUNE_NO_DEPRECATED_END


  return 0;
}
