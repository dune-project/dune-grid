// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_IO_FILE_AMIRAMESH_PSURFACE_BOUNDARY_HH
#define DUNE_GRID_IO_FILE_AMIRAMESH_PSURFACE_BOUNDARY_HH

#warning Support for PSurface is deprecated and will be removed after Dune 2.8.

/** \file
 *  \brief A domain boundary implemented by the psurface library
 */
#include <memory>

#include <dune/grid/common/gridfactory.hh>

#if HAVE_PSURFACE
#include <psurface/PSurface.h>
#include "psurface/AmiraMeshIO.h"
#if HAVE_PSURFACE_2_0
#include <psurface/Hdf5IO.h>
#endif

#if HAVE_AMIRAMESH
#include <amiramesh/AmiraMesh.h>
#endif


namespace Dune {

  /** \brief A domain boundary implemented by the psurface library
   *
   * \warning This code is experimental.  It may change in all kinds of unexpected ways
   * without prior notice.  Use it only if you know what you are doing.
   *
   * \tparam dim The dimension of the <b>boundary</b>, not the dimension of the grid.
   *
   * \note The dimension 'dim' must be 1 or 2, the psurface library doesn't implement
   *    anything else.
   */
  template <int dim, class field_type = double>
  class PSurfaceBoundary
  {
    static_assert((dim==1 or dim==2), "PSurfaceBoundaries can only have dimensions 1 or 2!");

  public:

    /** \brief Implementation of a single BoundarySegment using a psurface object */
    class PSurfaceBoundarySegment : public Dune::BoundarySegment<dim+1>
    {
    public:

      /** \brief Construct from a given psurface and triangle number
       *
       * \param psurfaceBoundary The psurface object that implements the segment
       * \param segment The number of the boundary segment in the psurface object
       */
      PSurfaceBoundarySegment(const std::shared_ptr<PSurfaceBoundary<dim, field_type> >& psurfaceBoundary, int segment)
        : psurfaceBoundary_(psurfaceBoundary),
          segment_(segment)
      {}

      /** \brief Evaluate the parametrization function */
      virtual Dune::FieldVector<double, dim+1> operator()(const Dune::FieldVector<double,dim>& local) const {

        Dune::FieldVector<double, dim+1> result;

        // Transform local to barycentric coordinates
        psurface::StaticVector<field_type, dim> barCoords;

        if (dim==2) {
          barCoords[0] = 1 - local[0] - local[1];
          barCoords[1] = local[0];
        } else {          // dim==1
          barCoords[0] = 1 - local[0];
        }

        psurface::StaticVector<field_type,dim+1> r;

        if (!psurfaceBoundary_->getPSurfaceObject()->positionMap(segment_, barCoords, r))
          DUNE_THROW(Dune::GridError, "psurface::positionMap returned error code");

        for (int i=0; i<dim+1; i++)
          result[i] = r[i];

        return result;
      }

      std::shared_ptr<PSurfaceBoundary<dim, field_type> > psurfaceBoundary_;
      int segment_;
    };



    /** \brief Constructor from a given PSurface object */
    PSurfaceBoundary(psurface::PSurface<dim, field_type>* psurface)
      : psurface_(psurface)
    {}

    /** \brief Obtain a pointer to the underlying PSurface object
     *
     * Use this only if you know what you are doing.
     *
     * This class retains control over the memory management.  Do not
     * delete the object you receive.
     */
    psurface::PSurface<dim, field_type>* getPSurfaceObject()
    {
      return psurface_.get();
    }

    /** \brief Read a PSurface boundary description from a file
     *
     * Supported file formats are AmiraMesh, and hdf5 if you have psurface-2.0 or newer.
     * Your psurface library needs to be specially configured to support those file formats.
     * The format is determined by the filename suffix.  If it is .h5, then the file
     * is assumed to be hdf5.  Otherwise it is assumed to be AmiraMesh.
     */
    static std::shared_ptr<PSurfaceBoundary<dim, field_type> > read(const std::string& filename)
    {
      psurface::PSurface<dim, field_type>* newDomain;

#if HAVE_PSURFACE_2_0
      // Try to read the file as an hdf5 file
      if (filename.find(".h5")==filename.length()-3) {
        newDomain = psurface::Hdf5IO<field_type,dim>::read(filename);
        if (newDomain)
          return std::make_shared<PSurfaceBoundary<dim, field_type> >(newDomain);
      }
#endif

#if HAVE_AMIRAMESH
      std::unique_ptr<AmiraMesh> am(AmiraMesh::read(filename.c_str()));

      if (!am.get())
        DUNE_THROW(IOError, "An error has occurred while reading " << filename);

      newDomain
        = (psurface::PSurface<dim, field_type>*) psurface::AmiraMeshIO<field_type>::readAmiraMesh(am.get(), filename.c_str());

      if (!newDomain)
        DUNE_THROW(IOError, "An error has occurred while reading " << filename);

      return std::make_shared<PSurfaceBoundary<dim, field_type> >(newDomain);
#else
      DUNE_THROW(IOError, "The given file is not in a supported format!");
#endif
    }

  private:

    std::unique_ptr<psurface::PSurface<dim, field_type> > psurface_;

  };

}

#endif // #if HAVE_PSURFACE
#endif // #ifndef DUNE_GRID_IO_FILE_AMIRAMESH_PSURFACE_BOUNDARY_HH
