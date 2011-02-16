// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef PSURFACE_H
#define PSURFACE_H

/** \file
    \brief This is a wrapper for the standalone psurface library

    \deprecated This is a badly written C-style wrapper to the psurface library.
    It used to be part of libpsurface itself, but it had to be removed from
    there for two reasons.  First of all, C support was to be abandonded
    eventually.  Secondly, there is also an (important) header called PSurface.h
    in libsurface.  The two together created problems on case-insensitive
    file systems such as found, e.g., on certain Macs.  The latter reason
    is the reason why this file had to be removed from libpsurface rather
    quickly.  It didn't have the time to properly adapt the corresponding
    code in dune-grid, therefore I just copied the entire file to here,
    to have working psurface support in dune-grid even with newer versions
    of libpsurface that will not ship psurface.h.  It will be removed
    before the next (2.2) release.
 */

namespace psurface {

  /** \brief Error codes returned by the methods */
  enum {OK, ERROR};

  /** \brief This command loads a boundary parametrization into memory and
   * assigns a label to it.
   */
  int LoadMesh(const char* label, const char* filename);

  /** \brief Removes the parametrization designated by 'label' from memory
   * @return AMIRA_OK on success.
   */
  int RemoveDomain(const char* label);

  /** \brief Removes all parametrizations.
   * @return Number of domains removed
   */
  void RemoveAllDomains();

  /** \brief Sets up a particular parametrization for access.
   *
   * Sets up a particular parametrization for access.  All subsequent
   * access operations will use the P. that has been chosen here.  When you
   * want to use a new boundary, just call this function again with the
   * new label.
   */
  int StartEditingDomain(const char* label);

  /** \brief This is the main access function.
   * @return AMIRA_OK on success
   */
  int CallParametrization(int tri, double* p, int* res, double* coords, int seed);

  /** \brief This is the access function for position.
   * @return AMIRA_OK on success
   */
  int CallPositionParametrization(int tri, double* p, double* res);

  /** \brief This is the access function for position.
   * @return AMIRA_OK on success
   */
  int CallPositionParametrizationForDomain(int domain, int tri, double* p, double* res);

  /** \brief This is the access function for surface normals
   */
  int CallDirectNormalParametrization(int tri, double* p, double* res);

  /** \brief Returns the number of vertices in the currently chosen parametrization.
   */
  int GetNoOfNodes();

  /** \brief Returns the number of base grid triangles in the currently chosen parametrization
   */
  int GetNoOfSegments();

  /** \brief Returns the number of patches in the currently chosen parametrization
   */
  int GetNoOfPatches();

  /** \brief Returns the number of currently loaded domains
   */
  int GetNoOfDomains();

  /** \brief Gets the vertex indices of a given base grid triangle
   */
  void GetNodeNumbersOfSegment(int* points, int tri);

  /** \brief Returns the left and right material of a base grid triangle
   */
  void GetLeftAndRightSideOfSegment(int* left, int* right, int tri);

  /** \brief Returns the boundary Id of a given base grid triangle. */
  int GetBoundaryIdOfSegment(int tri);

  /** \brief Returns the patch number of a given base grid triangle. */
  int GetPatchNoOfSegment(int tri);

}

#endif
