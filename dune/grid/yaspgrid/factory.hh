// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_YASPGRID_FACTORY_HH
#define DUNE_GRID_YASPGRID_FACTORY_HH

#include <vector>
#include <dune/common/array.hh>
#include <dune/grid/common/exceptions.hh>

#include "coordinates.hh"

/** \file
 *  \brief This file provides a factory class for tensorproduct
 *  YaspGrids. This is a collection of methods to generate monotonous
 *  sequences as needed for a tensorproduct grid. Apart
 *  from easy ones for locally equidistant grids, there are also
 *  more involved methods like splitting a range according to a
 *  geometric series.
 *
 *  \author Dominic Kempf
 */

namespace Dune
{
  /** \brief A factory class for conveniently creating tensorproduct grids
   *
   * \tparam ctype the coordinate type to use
   * \tparam dim the grid dimension
   */
  template<typename ctype, int dim>
  class TensorYaspGridFactory
  {
  public:

    typedef YaspGrid<dim, TensorProductCoordinateContainer<ctype, dim> > Grid;
    typedef typename Grid::Traits::CollectiveCommunication Comm;

    //! initialize the factory with a set of default values
    TensorYaspGridFactory () : _periodic (), _overlap (1)
    {}

    //! finalizes the factory and gives a pointer to the constructed grid
    Grid* createGrid(Comm comm = Comm())
    {
      if (!Dune::Yasp::checkIfMonotonous (_coords))
        DUNE_THROW( Dune::Exception,
          "TensorYaspFactory did not get enough coordinate information to construct a grid!");
      return new grid(_coords, _periodic, _overlap, comm);
    }

    //! set whether the grid is periodic in direction dir
    void setPeriodicity (int dir)
    {
      _periodic[dir] = true;
    }

    //! set the number of overlap cells
    void setOverlap (int overlap)
    {
      _overlap = overlap;
    }

    /** \brief set a starting value in a given direction d
     * \param d the coordinate direction
     * \param value the value to set
     *
     * This resizes the coordinate vector for the given direction to 1.
     * Not using this function will result in 0.0 to be used as a lower
     * bound of the coordinate range.
     */
    void setStart (int d, ctype value)
    {
      _coords[d].resize(1);
      _coords[d][0] = value;
    }

    /** \brief pushs n intervals of length h in direction d
     * \param d the coordinate direction
     * \param n the number of intervals to add
     * \param h the interval length
     *
     * Given a vector with last element \f$x_0\f$, this will add elements
     * \f$x_1,\dots ,x_n\f$ such that \f$x_i=x_0+i*h\f$.
     */
    void fill_intervals (int d, int n, ctype h)
    {
      emptyCheck (d);
      for (int i = 0; i < n; i++)
        _coords[d].push_back (_coords[d].back () + h);
    }

    /** \brief fills the range to end with n intervals in direction d
     * \param d the coordinate direction
     * \param n the number of intervals to add
     * \param end the coordinate on the upper border of the range
     *
     * Given a vector with last element \f$x_0\f$, this will add elements
     * \f$x_1,\dots ,x_n\f$ such that \f$x_i=x_0+i*\frac{end-x_0}{n}\f$.
     */
    void fill_range (int d, int n, ctype end)
    {
      emptyCheck (d);
      const ctype h = (end - _coords[d].back ()) / n;
      for (int i = 0; i < n - 1; i++)
        _coords[d].push_back (_coords[d].back () + h);
      _coords[d].push_back (end);
    }

    /** \brief adds intervals in direction d until a given coordinate is reached
     * \param d the coordinate direction
     * \param h the interval length
     * \param end the coordinate on the upper border of the range
     *
     *  Given a vector with last element \f$x_0\f$, this will add elements
     * \f$x_1,\dots ,x_n\f$ such that \f$x_n < end < x_n + h\f$ and \f$x_{i+1}-x_i = h\f$.
     */
    void fill_until (int d, ctype h, ctype end)
    {
      emptyCheck (d);
      while (_coords[d].back () < end)
        _coords[d].push_back (_coords[d].back () + h);
    }

    /** \brief adds n intervals in direction d with a given length ratio and a given starting interval length.
     *  \param d the coordinate direction
     *  \param n the number of intervals to add
     *  \param ratio the ratio of \f$ h_{i+1}\f$ to \f$h_i\f$
     *  \param h0 the starting interval length (optional)
     *
     *  Given a vector with last element \f$x_0\f$, this will add elements
     *  \f$x_1,\dots ,x_n\f$ such that - with \f$h_i=x_{i+1}-x_i\f$ - \f$h_{i+1}=qh_i$
     *  for a given ratio \f$q\f$. The first interval length can either be explicitly
     *  given or be deduced by multiplying the ratio with the last interval
     *  length in the container
     */
    void geometric_fill_intervals (int d, int n, ctype ratio, ctype h0 =
        static_cast<ctype> (0))
    {
      emptyCheck (d);
      ctype h = h0;
      if (h0 == static_cast<ctype>(0))
        h = lastInterval (d) * ratio;
      for (int i = 0; i < n; i++)
      {
        _coords[d].push_back (_coords[d].back () + h);
        h *= ratio;
      }
    }

    /** \brief adds intervals in direction d according with a given length ratio until a given coordinate is reached
     * \param d the coordinate direction
     * \param ratio the ratio of \f$ h_{i+1}\f$ to \f$h_i\f$
     * \param end the coordinate on the right border of the range
     * \param h0 the starting interval length (optional)
     *
     *  Given a vector with last element \f$x_0\f$, this will add elements
     *  \f$x_1,\dots ,x_n\f$ such that - with \f$h_i=x_{i+1}-x_i\f$ - \f$h_{i+1}=qh_i$
     *  for a given ratio \f$q\f$ and that \f$x_n < end < x_n + h\f$. The first interval
     *  length can either be explicitly given or be deduced by multiplying the ratio with
     *  the last interval length in the container.
     */
    void geometric_fill_until (int d, ctype ratio, ctype end, ctype h0 =
        static_cast<ctype> (0))
    {
      emptyCheck (d);
      ctype h = h0;
      if (h0 == static_cast<ctype>(0))
        h = lastInterval (d) * ratio;
      while (_coords[d].back () < end)
      {
        _coords[d].push_back (_coords[d].back () + h);
        h *= ratio;
      }
    }

    /** \brief fills a coordinate range in direction d with n intervals according to a geometric series
     * \param d the coordinate direction
     * \param n the number of intervals to add
     * \param end the coordinate of the upper border of the range
     * \param h the interval length to start or end with (see below) (optional)
     * \param first true if the given h is to be the first interval, false if last one
     *
     * Given a vector with last element \f$x_0\f$, this will add elements
     * \f$x_1,\dots ,x_n\f$ such that - with \f$h_i=x_{i+1}-x_i\f$ - the ratio
     * \f$h_{i+1} / h_i\f$ is fixed throughout the range and \f$x_n=end\f$.
     * The first interval length can either be explicitly given or be deduced by taking
     * the last interval length in the container. By setting the optional parameter first
     * to false, the given h can instead be used as last interval length in the range.
     */
    void geometric_fill_range (int d, int n, ctype end, ctype h =
        static_cast<ctype> (0),
        bool first = true)
    {
      emptyCheck (d);
      if (h < 1e-8)
        h = lastInterval (d);
      ctype ratio = newton (n, _coords[d].back (), end, h);
      if (!first)
      {
        h = h * pow (ratio, n - 1);
        ratio = 1 / ratio;
      }
      for (int i = 0; i < n - 1; i++)
      {
        _coords[d].push_back (_coords[d].back () + h);
        h *= ratio;
      }
      _coords[d].push_back (end);
    }

    //! print the coordinate information given to the factory so far
    void print()
    {
      for (int i=0; i<dim; i++)
      {
        std::cout << "Container in direction " << i << ":" << std::endl << "Coordinates: ";
        for (auto it = _coords[i].begin(); it != _coords[i].end(); ++it)
          std::cout << *it << "  ";
        std::cout << std::endl << "Interval lengths: ";

        std::vector<ctype> meshsize;
        for (auto it = _coords[i].begin(); it != _coords[i].end()-1;)
        {
          meshsize.push_back(-1.*(*it));
          ++it;
          meshsize.back() += *it;
        }

        for (auto it = meshsize.begin(); it != meshsize.end(); ++it)
          std::cout << *it << "  ";
        std::cout << std::endl << "Ratios between interval lengths: ";

        std::vector<ctype> ratios;
        for (auto it = meshsize.begin(); it != meshsize.end() - 1 ;)
          ratios.push_back((1./(*it)) * *(++it));

        for (auto it = ratios.begin(); it != ratios.end(); ++it)
          std::cout << *it << "  ";
        std::cout << std::endl << std::endl << std::endl;
      }
    }

  private:
    // check whether the ith component is empty and add a 0.0 entry if so
    void emptyCheck (int i)
    {
      if (_coords[i].empty ())
        _coords[i].push_back (static_cast<ctype> (0));
    }

    // returns the last interval length in direction d
    ctype lastInterval (int d)
    {
      if (_coords[d].size () < 2)
        DUNE_THROW(
            GridError,
            "Not enough elements in coordinate container to deduce interval length in TensorYaspFactory");
      else
        return _coords[d].back () - _coords[d][_coords[d].size () - 2];
    }

    /** this implements a simple newton iteration for the function
     *  \f$f(x) = -x^n+\frac{x_e-x_s}{h} (x-1)+1\f$
     */
    ctype newton (int n, ctype x_s, ctype x_e, ctype h)
    {
      ctype m = (x_e - x_s) / h;
      ctype xold = 0.0;
      ctype xnew = x_e - x_s;
      while (std::abs (xnew - xold) > 1E-8)
      {
        xold = xnew;
        xnew = xold
            - (-pow (xold, n) + m * xold - m + 1)
            / (-n * pow (xold, n - 1) + m);
      }
      if (std::abs (xnew - 1) < 1E-6)
      {
        xold = x_e - x_s;
        xnew = 0.0;
        while (std::abs (xnew - xold) > 1E-8)
        {
          xold = xnew;
          xnew = xold
              - (-pow (xold, n) + m * xold - m + 1)
              / (-n * pow (xold, n - 1) + m);
        }
      }
      return xnew;
    }

    Dune::array<std::vector<ctype>, dim> _coords;
    std::bitset<dim> _periodic;
    int _overlap;
  };
}

#endif
