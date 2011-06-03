// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef PSURFACE_BOUNDARY_HH
#define PSURFACE_BOUNDARY_HH

/** \file
 *  \brief A domain boundary implemented by the psurface library
 */

#include <psurface/PSurface.h>

namespace Dune {

  /** \brief A domain boundary implemented by the psurface library
   *
   * \warning This code is experimental.  It may change in all kinds of unexpected ways
   * without prior notice.  Use it only if you know what you are doing.
   *
   * \tparam dim The dimension of the <b>boundary</b>, not the dimension of the grid.
   */
  template <int dim>
  class PSurfaceBoundary
  {
    dune_static_assert((dim==1 or dim==2), "PSurfaceBoundaries can only have dimensions 1 or 2!");

  public:

    /** \brief Implementation of a single BoundarySegment using a psurface object */
    class PSurfaceBoundarySegment : public Dune::BoundarySegment<dim+1>
    {
    public:

      /** \brief Construct from a given psurface and triangle number
       *
       * \param psurface The psurface object that implements the segment
       * \param segment The number of the boundary segment in the psurface object
       */
      PSurfaceBoundarySegment(PSurface<2,float>* psurface, int segment)
        : psurface_(psurface),
          segment_(segment)
      {}

      /** \brief Evaluate the parametrization function */
      virtual Dune::FieldVector<double, dim+1> operator()(const Dune::FieldVector<double,dim>& local) const {

        Dune::FieldVector<double, dim+1> result;

        // Transform local to barycentric coordinates
        StaticVector<float,dim> barCoords;

        if (dim==2) {
          barCoords[0] = 1 - local[0] - local[1];
          barCoords[1] = local[0];
        } else {          // dim==1
          barCoords[0] = 1 - local[0];
        }

        StaticVector<float,dim+1> r;

        if (!psurface_->positionMap(segment_, barCoords, r))
          DUNE_THROW(Dune::GridError, "psurface::positionMap returned error code");

        for (int i=0; i<dim+1; i++)
          result[i] = r[i];

        return result;
      }

      PSurface<dim,float>* psurface_;
      int segment_;
    };



    /** \brief Constructor from a given PSurface object */
    PSurfaceBoundary(PSurface<dim,float>* psurface)
      : psurface_(psurface)
    {}

    /** \brief Obtain a pointer to the underlying PSurface object
     *
     * Use this only if you know what you are doing.
     *
     * This class retains control over the memory management.  Do not
     * delete the object to receive.
     */
    PSurface<dim,float>* getPSurfaceObject()
    {
      return psurface_.get();
    }

  private:

    std::auto_ptr<PSurface<dim,float> > psurface_;

  };

}

#endif
