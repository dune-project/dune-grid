// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GRID_UTILITY_TENSORGRIDFACTORY_HH
#define DUNE_GRID_UTILITY_TENSORGRIDFACTORY_HH

/** \file
 *  \brief This file provides a factory class for tensorproduct
 *  grids. This is a collection of methods to generate monotonous
 *  sequences as needed for a tensorproduct grid. Apart
 *  from easy ones for locally equidistant grids, there are also
 *  more involved methods like splitting a range according to a
 *  geometric series.
 *
 *  The grid generation process is implemented for unstructured grids
 *  and for YaspGrid.
 *
 *  \author Dominic Kempf
 */

#include<array>
#include<memory>
#include<vector>

#include <dune/common/fvector.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/yaspgrid.hh>
#include<dune/grid/utility/multiindex.hh>

namespace Dune
{
  // forward declaration of TensorGridFactoryCreator, which is the real factory
  // that should be specialized for each grid.
  template<typename Grid>
  class TensorGridFactoryCreator;

  /** \brief A factory class for conveniently creating tensorproduct grids
   *
   * \tparam Grid the grid type
   */
  template<typename Grid>
  class TensorGridFactory
  {
  public:
    typedef typename Grid::Traits::Communication Comm;
    typedef typename Grid::ctype ctype;
    static const int dim = Grid::dimension;

    std::unique_ptr<Grid> createGrid(Comm comm = Comm())
    {
      TensorGridFactoryCreator<Grid> creator(*this);
      return creator.createGrid(comm);
    }

    std::array<std::vector<ctype> , dim> coords() const
    {
      return _coords;
    }

    //! allow to manually tune the factory by overloading operator[] to export the coordinate vectors in the coordinate directories.
    std::vector<ctype>& operator[](std::size_t d)
    {
      return _coords[d];
    }

    //! allow to manually tune the factory by overloading operator[] to export the coordinate vectors in the coordinate directories.
    const std::vector<ctype>& operator[](std::size_t d) const
    {
      return _coords[d];
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
    void fillIntervals (int d, int n, ctype h)
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
    void fillRange (int d, int n, ctype end)
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
    void fillUntil (int d, ctype h, ctype end)
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
     *  \f$x_1,\dots ,x_n\f$ such that \f$h_{i+1}=qh_i\f$ for a given ratio
     *  \f$q\f$ and interval length \f$h_i=x_{i+1}-x_i\f$. The first interval length
     *  can either be explicitly given or be deduced by multiplying the ratio
     *  with the last interval length in the container.
     */
    void geometricFillIntervals (int d, int n, ctype ratio, ctype h0 =
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
     * \param ratio the ratio of \f$h_{i+1}\f$ to \f$h_i\f$
     * \param end the coordinate on the right border of the range
     * \param h0 the starting interval length (optional)
     *
     *  Given a vector with last element \f$x_0\f$, this will add elements
     *  \f$x_1,\dots ,x_n\f$ such that - with \f$h_i=x_{i+1}-x_i\f$ - \f$h_{i+1}=qh_i\f$
     *  for a given ratio \f$q\f$ and that \f$x_n < end < x_n + h\f$. The first interval
     *  length can either be explicitly given or be deduced by multiplying the ratio with
     *  the last interval length in the container.
     */
    void geometricFillUntil (int d, ctype ratio, ctype end, ctype h0 = static_cast<ctype> (0))
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
     * \f$x_1,\dots ,x_n\f$ such that the ratio \f$h_{i+1} / h_i\f$ is fixed throughout
     * the range and \f$x_n=end\f$, while \f$h_i=x_{i+1}-x_i\f$ is the interval length.
     * The first interval length can either be explicitly given or be deduced by taking
     * the last interval length in the container. By setting the optional parameter first
     * to false, the given h can instead be used as last interval length in the range.
     */
    void geometricFillRange (int d, int n, ctype end, ctype h =
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

    std::array<std::vector<ctype>, dim> _coords;
  };

  // class that implements the actual grid creation process. The default is implementing
  // standard creation for unstructured grids. Provide a specialization for other grids.
  template<typename Grid>
  class TensorGridFactoryCreator
  {
  public:
    typedef typename Grid::Traits::Communication Comm;
    typedef typename Grid::ctype ctype;
    static const int dim = Grid::dimension;

    TensorGridFactoryCreator(const TensorGridFactory<Grid>& factory) : _factory(factory) {}

    std::unique_ptr<Grid> createGrid(Comm comm)
    {
      // The grid factory
      GridFactory<Grid> fac;

      if (comm.rank() == 0)
      {
        // determine the size of the grid
        std::array<unsigned int, dim> vsizes, esizes;
        std::size_t size = 1;
        for (std::size_t i = 0; i<dim; ++i)
        {
          vsizes[i] = _factory[i].size();
          esizes[i] = vsizes[i] - 1;
          size *= vsizes[i];
        }

        // insert all vertices
        FactoryUtilities::MultiIndex<dim> index(vsizes);
        for (std::size_t i=0; i<size; ++i, ++index)
        {
          Dune::FieldVector<ctype, dim> position;
          for (std::size_t j = 0; j<dim; ++j)
            position[j] = _factory[j][index[j]];
          fac.insertVertex(position);
        }

        // compute the offsets
        std::array<std::size_t, dim> offsets;
        offsets[0] = 1;
        for (std::size_t i=1; i<dim; i++)
          offsets[i] = offsets[i-1] * vsizes[i-1];

        // Compute an element template (the cube at (0,...,0).  All
        // other cubes are constructed by moving this template around
        unsigned int nCorners = 1<<dim;

        std::vector<unsigned int> cornersTemplate(nCorners,0);

        for (size_t i=0; i<nCorners; i++)
          for (int j=0; j<dim; j++)
            if ( i & (1<<j) )
              cornersTemplate[i] += offsets[j];

        // Insert elements
        FactoryUtilities::MultiIndex<dim> eindex(esizes);

        // Compute the total number of elements to be created
        int numElements = eindex.cycle();

        for (int i=0; i<numElements; i++, ++eindex)
        {
          // 'base' is the index of the lower left element corner
          unsigned int base = 0;
          for (int j=0; j<dim; j++)
            base += eindex[j] * offsets[j];

          // insert new element
          std::vector<unsigned int> corners = cornersTemplate;
          for (size_t j=0; j<corners.size(); j++)
            corners[j] += base;

          fac.insertElement(GeometryTypes::cube(dim), corners);
        }
      }

      return std::unique_ptr<Grid>(fac.createGrid());
    }

  private:
    const TensorGridFactory<Grid>& _factory;
  };

  template<typename ctype, int dim>
  class TensorGridFactoryCreator<YaspGrid<dim, TensorProductCoordinates<ctype, dim> > >
  {
  public:
    typedef YaspGrid<dim, TensorProductCoordinates<ctype, dim> > Grid;
    typedef typename Grid::Communication Comm;

    TensorGridFactoryCreator(const TensorGridFactory<Grid>& factory) : _factory(factory) {}

    std::unique_ptr<Grid> createGrid(Comm comm)
    {
      return std::make_unique<Grid>(_factory.coords(), std::bitset<dim>(0ULL), 1, comm);
    }
  private:
    const TensorGridFactory<Grid>& _factory;
  };
}

#endif
