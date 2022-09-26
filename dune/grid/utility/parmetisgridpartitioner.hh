// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GRID_UTILITY_PARMETISGRIDPARTITIONER_HH
#define DUNE_GRID_UTILITY_PARMETISGRIDPARTITIONER_HH

/** \file
 *  \brief Compute a repartitioning of a Dune grid using ParMetis
 */

#include <algorithm>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/utility/globalindexset.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/rangegenerators.hh>

#if HAVE_PARMETIS

#include <parmetis.h>

// only enable for ParMETIS because the implementation uses functions that
// are not emulated by scotch
#ifdef PARMETIS_MAJOR_VERSION

namespace Dune
{

  /** \brief Compute a repartitioning of a Dune grid using ParMetis
   *
   * \author Benjamin Bykowski
   */
  template<class GridView>
  struct ParMetisGridPartitioner {

    // define index type as provided by ParMETIS
#if PARMETIS_MAJOR_VERSION > 3
    typedef idx_t idx_type;
    typedef ::real_t real_type;
#else
    typedef int idx_type;
    typedef float real_type;
#endif // PARMETIS_MAJOR_VERSION > 3

    constexpr static int dimension = GridView::dimension;


    /** \brief Create an initial partitioning of a Dune grid, i.e., not taking into account communication cost
     *
     * This pipes a Dune grid into the method ParMETIS_V3_PartMeshKway (see the ParMetis documentation for details)
     *
     * \param gv The grid view to be partitioned
     * \param mpihelper The MPIHelper object, needed to get the MPI communicator
     *
     * \return std::vector with one uint per All_Partition element.  For each element, the entry is the
     *    number of the partition the element is assigned to.
     */
    static std::vector<unsigned> partition(const GridView& gv, const Dune::MPIHelper& mpihelper) {
      const unsigned numElements = gv.size(0);

      std::vector<unsigned> part(numElements);

      // Setup parameters for ParMETIS
      idx_type wgtflag = 0;                                  // we don't use weights
      idx_type numflag = 0;                                  // we are using C-style arrays
      idx_type ncon = 1;                                     // number of balance constraints
      idx_type ncommonnodes = 2;                             // number of nodes elements must have in common to be considered adjacent to each other
      idx_type options[4] = {0, 0, 0, 0};                    // use default values for random seed, output and coupling
      idx_type edgecut;                                      // will store number of edges cut by partition
      idx_type nparts = mpihelper.size();                    // number of parts equals number of processes
      std::vector<real_type> tpwgts(ncon*nparts, 1./nparts); // load per subdomain and weight (same load on every process)
      std::vector<real_type> ubvec(ncon, 1.05);              // weight tolerance (same weight tolerance for every weight there is)

      // The difference elmdist[i+1] - elmdist[i] is the number of nodes that are on process i
      std::vector<idx_type> elmdist(nparts+1);
      elmdist[0] = 0;
      std::fill(elmdist.begin()+1, elmdist.end(), gv.size(0)); // all elements are on process zero

      // Create and fill arrays "eptr", where eptr[i] is the number of vertices that belong to the i-th element, and
      // "eind" contains the vertex-numbers of the i-the element in eind[eptr[i]] to eind[eptr[i+1]-1]
      std::vector<idx_type> eptr, eind;
      int numVertices = 0;
      eptr.push_back(numVertices);

      for (const auto& element : elements(gv, Partitions::interior)) {
        const size_t curNumVertices = referenceElement<double, dimension>(element.type()).size(dimension);

        numVertices += curNumVertices;
        eptr.push_back(numVertices);

        for (size_t k = 0; k < curNumVertices; ++k)
          eind.push_back(gv.indexSet().subIndex(element, k, dimension));
      }

      // Partition mesh using ParMETIS
      if (0 == mpihelper.rank()) {
        MPI_Comm comm = Dune::MPIHelper::getLocalCommunicator();

#if PARMETIS_MAJOR_VERSION >= 4
        const int OK =
#endif
        ParMETIS_V3_PartMeshKway(elmdist.data(), eptr.data(), eind.data(), NULL, &wgtflag, &numflag,
                                 &ncon, &ncommonnodes, &nparts, tpwgts.data(), ubvec.data(),
                                 options, &edgecut, reinterpret_cast<idx_type*>(part.data()), &comm);

#if PARMETIS_MAJOR_VERSION >= 4
        if (OK != METIS_OK)
          DUNE_THROW(Dune::Exception, "ParMETIS returned an error code.");
#endif
      }

      return part;
    }

    /** \brief Create a repartitioning of a distributed Dune grid
     *
     * This pipes a Dune grid into the method ParMETIS_V3_AdaptiveRepart (see the ParMetis documentation for details)
     *
     * \param gv The grid view to be partitioned
     * \param mpihelper The MPIHelper object, needed to get the MPI communicator
     * \param itr ParMetis parameter to balance grid transport cost vs. communication cost for the redistributed grid.
     *
     * \return std::vector with one uint per All_Partition element.  For each Interior_Partition element, the entry is the
     *    number of the partition the element is assigned to.
     */
    static std::vector<unsigned> repartition(const GridView& gv, const Dune::MPIHelper& mpihelper, real_type itr = 1000) {

      // Create global index map
      GlobalIndexSet<GridView> globalIndex(gv,0);

      int numElements = std::distance(gv.template begin<0, Interior_Partition>(),
                                                 gv.template end<0, Interior_Partition>());

      std::vector<unsigned> interiorPart(numElements);

      // Setup parameters for ParMETIS
      idx_type wgtflag = 0;                                  // we don't use weights
      idx_type numflag = 0;                                  // we are using C-style arrays
      idx_type ncon = 1;                                     // number of balance constraints
      idx_type options[4] = {0, 0, 0, 0};                    // use default values for random seed, output and coupling
      idx_type edgecut;                                      // will store number of edges cut by partition
      idx_type nparts = mpihelper.size();                    // number of parts equals number of processes
      std::vector<real_type> tpwgts(ncon*nparts, 1./nparts); // load per subdomain and weight (same load on every process)
      std::vector<real_type> ubvec(ncon, 1.05);              // weight tolerance (same weight tolerance for every weight there is)

      MPI_Comm comm = Dune::MPIHelper::getCommunicator();

      // Make the number of interior elements of each processor available to all processors
      std::vector<int> offset(gv.comm().size());
      std::fill(offset.begin(), offset.end(), 0);

      gv.comm().template allgather<int>(&numElements, 1, offset.data());

      // The difference vtxdist[i+1] - vtxdist[i] is the number of elements that are on process i
      std::vector<idx_type> vtxdist(gv.comm().size()+1);
      vtxdist[0] = 0;

      for (unsigned int i=1; i<vtxdist.size(); ++i)
        vtxdist[i] = vtxdist[i-1] + offset[i-1];

      // Set up element adjacency lists
      std::vector<idx_type> xadj, adjncy;
      xadj.push_back(0);

      for (const auto& element : elements(gv, Partitions::interior)) {
        size_t numNeighbors = 0;

        for (const auto& in : intersections(gv, element)) {
          if (in.neighbor()) {
            adjncy.push_back(globalIndex.index(in.outside()));

            ++numNeighbors;
          }
        }

        xadj.push_back(xadj.back() + numNeighbors);
      }

#if PARMETIS_MAJOR_VERSION >= 4
      const int OK =
#endif
        ParMETIS_V3_AdaptiveRepart(vtxdist.data(), xadj.data(), adjncy.data(), NULL, NULL, NULL,
                                   &wgtflag, &numflag, &ncon, &nparts, tpwgts.data(), ubvec.data(),
                                   &itr, options, &edgecut, reinterpret_cast<idx_type*>(interiorPart.data()), &comm);

#if PARMETIS_MAJOR_VERSION >= 4
      if (OK != METIS_OK)
        DUNE_THROW(Dune::Exception, "ParMETIS returned error code " << OK);
#endif

      // At this point, interiorPart contains a target rank for each interior element, and they are sorted
      // by the order in which the grid view traverses them.  Now we need to do two things:
      // a) Add additional dummy entries for the ghost elements
      // b) Use the element index for the actual ordering.  Since there may be different types of elements,
      //    we cannot use the index set directly, but have to go through a Mapper.

      typedef MultipleCodimMultipleGeomTypeMapper<GridView> ElementMapper;
      ElementMapper elementMapper(gv, mcmgElementLayout());

      std::vector<unsigned int> part(gv.size(0));
      std::fill(part.begin(), part.end(), 0);
      unsigned int c = 0;
      for (const auto& element : elements(gv, Partitions::interior))
        part[elementMapper.index(element)] = interiorPart[c++];

      return part;
    }
  };

}  // namespace Dune

#else // PARMETIS_MAJOR_VERSION
#warning "You seem to be using the ParMETIS emulation layer of scotch, which does not work with this file."
#endif

#else // HAVE_PARMETIS
#warning "PARMETIS was not found, please check your configuration"
#endif

#endif // DUNE_GRID_UTILITY_PARMETISGRIDPARTITIONER_HH
