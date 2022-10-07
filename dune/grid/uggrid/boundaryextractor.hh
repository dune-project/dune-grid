// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_BOUNDARY_EXTRACTOR_HH
#define DUNE_BOUNDARY_EXTRACTOR_HH

/** \file
    \brief Contains helper classes for the creation of UGGrid objects
    \author Oliver Sander
 */

#include <iostream>
#include <vector>
#include <set>
#include <array>


namespace Dune {

  /** \brief Boundary segments that can be compared

     This general implementation is empty.  Only specializations for dim==2 and dim==3 exist.
   */
  template <int dim>
  class UGGridBoundarySegment {};

  /** \brief Specialization of the boundary segment class for 2d */
  template <>
  class UGGridBoundarySegment<2> : public std::array<int,2> {

  public:

    /** \brief Always returns 2 */
    int numVertices() const {
      return 2;
    }

    /** \brief Compare the vertex lists modulo permutation */
    bool operator<(const UGGridBoundarySegment<2>& other) const {

      array<int,2> sorted1, sorted2;

      // ////////////////////////////////////////////////////////////////////////////
      // Sort the two arrays to get rid of cyclic permutations in mirror symmetry
      // ////////////////////////////////////////////////////////////////////////////
      if ((*this)[0] < (*this)[1])
        sorted1 = (*this);
      else {
        sorted1[0] = (*this)[1];
        sorted1[1] = (*this)[0];
      }

      if (other[0] < other[1])
        sorted2 = other;
      else {
        sorted2[0] = other[1];
        sorted2[1] = other[0];
      }

      // ////////////////////////////////////////////////////////////////////////////
      //   Compare the two sorted arrays
      // ////////////////////////////////////////////////////////////////////////////

      return sorted1 < sorted2;

    }

  };

  /** \brief Specialization of the boundary segment class for 2d */
  template <>
  class UGGridBoundarySegment<3> : public std::array<int,4> {

  public:

    /** \brief 3 or 4 */
    int numVertices() const {
      return ((*this)[3]==-1) ? 3 : 4;
    }

    /** \brief Compare the vertex lists modulo permutation */
    bool operator<(const UGGridBoundarySegment<3>& other) const {
      UGGridBoundarySegment<3> sorted1 = (*this);
      UGGridBoundarySegment<3> sorted2 = other;

      if (numVertices()<other.numVertices())
        return true;

      if (numVertices()>other.numVertices())
        return false;

      // ////////////////////////////////////////////////////////////////////////////
      // Sort the two arrays to get rid of cyclic permutations and mirror symmetry
      // ////////////////////////////////////////////////////////////////////////////

      // bubble sort: sort and compare together until the first nonmatching digits are found
      for (int i=numVertices()-1; i>=0; i--) {

        for (int j=0; j<i; j++) {

          if (sorted1[j] > sorted1[j+1])
            std::swap(sorted1[j], sorted1[j+1]);

          if (sorted2[j] > sorted2[j+1])
            std::swap(sorted2[j], sorted2[j+1]);

        }

        //
        if (sorted1[i]<sorted2[i])
          return true;
        else if (sorted1[i]>sorted2[i])
          return false;

      }

      // The sorted arrays are identical
      return false;
    }

  };

  //! Output operator for array
  inline std::ostream& operator<< (std::ostream& s, const UGGridBoundarySegment<2>& v)
  {
    return s << "[" << v[0] << ", " << v[1] << "]";
  }

  inline std::ostream& operator<< (std::ostream& s, const UGGridBoundarySegment<3>& v)
  {
    s << "[" << v[0] << ", " << v[1] << ", " << v[2];
    if (v[3]!=-1)     // quadrilateral
      s << ", " << v[3];
    return s << "]";
  }

  /** \brief Extracts the boundary faces and nodes from a set grid given as a set of elements
   */
  class BoundaryExtractor {

    typedef std::set<UGGridBoundarySegment<2> >::iterator SetIterator2d;
    typedef std::set<UGGridBoundarySegment<3> >::iterator SetIterator3d;

  public:

    static void detectBoundarySegments(const std::vector<unsigned char>& elementTypes,
                                       const std::vector<unsigned int>& elementVertices,
                                       std::set<UGGridBoundarySegment<2> >& boundarySegments);

    static void detectBoundarySegments(const std::vector<unsigned char>& elementTypes,
                                       const std::vector<unsigned int>& elementVertices,
                                       std::set<UGGridBoundarySegment<3> >& boundarySegments);

    template <int dim>
    static int detectBoundaryNodes(const std::set<UGGridBoundarySegment<dim> >& boundarySegments,
                                   int noOfNodes,
                                   std::vector<int>& isBoundaryNode);

  };

}

#endif
