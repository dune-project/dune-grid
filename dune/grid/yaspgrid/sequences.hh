// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_YASPGRID_SEQUENCES_HH
#define DUNE_GRID_YASPGRID_SEQUENCES_HH

#include <vector>

#include <dune/grid/common/exceptions.hh>

/** \file
 *  \brief This file provides methods to fill std::vectors with
 *  monotonous sequences as needed for a tensorproduct grid. Apart
 *  from easy ones for locally equidistant grids, there are also
 *  more involved methods like splitting a range according to a
 *  geometric series.
 *
 *  \author Dominic Kempf
 */

namespace Dune
{

  /** \brief pushs n intervals of length n on the vector v
   * \param v the std::vector to operate on
   * \param n the number of intervals to add
   * \param h the interval length
   *
   * Given a vector with last element \f$x_0\f$, this will add elements
   * \f$x_1,\dots ,x_n\f$ such that \f$x_i=x_0+i*h\f$.
   */
  template<class ctype>
  void fill_intervals(std::vector<ctype>& v, int n, ctype h)
  {
    if (v.empty())
      v.push_back(static_cast<ctype>(0));
    for (int i=0; i<n; i++)
      v.push_back(v.back()+h);
  }

  /** \brief fills the range to end with n intervals
   * \param v the std::vector to operate on
   * \param n the number of intervals to add
   * \param end the coordinate on the right border of the range
   *
   * Given a vector with last element \f$x_0\f$, this will add elements
   * \f$x_1,\dots ,x_n\f$ such that \f$x_i=x_0+i*\frac{end-x_0}{n}\f$.
   */
  template<class ctype>
  void fill_range(std::vector<ctype>& v, int n, ctype end)
  {
    if (v.empty())
      v.push_back(static_cast<ctype>(0));
    const ctype h = (end-v.back())/n;
    for (int i=0; i<n-1; i++)
      v.push_back(v.back()+h);
    v.push_back(end);
  }

  /** \brief adds intervals to a vector until a given coordinate is reached
   * \param v the std::vector to operate on
   * \param h the interval length
   * \param end the coordinate on the right border of the range
   *
   *  Given a vector with last element \f$x_0\f$, this will add elements
   * \f$x_1,\dots ,x_n\f$ such that \f$x_n < end < x_n + h\f$ and \f$x_{i+1}-x_i = h\f$.
   */
  template<class ctype>
  void fill_until(std::vector<ctype>& v, ctype h, ctype end)
  {
    if (v.empty())
      v.push_back(static_cast<ctype>(0));
    while (v.back() < end)
      v.push_back(v.back()+h);
  }

  /** \brief adds n intervals with a given length ratio and a given starting interval length.
   *  \param v the std::vector to operate on
   *  \param n the number of intervals to adds
   *  \param ratio the ratio of \f$ h_{i+1}\f$ to \f$h_i\f$
   *  \param h0 the starting interval length (optional)
   *
   *  Given a vector with last element \f$x_0\f$, this will add elements
   *  \f$x_1,\dots ,x_n\f$ such that - with \f$h_i=x_{i+1}-x_i\f$ - \f$h_{i+1}=qh_i$
   *  for a given ratio \f$q\f$. The first interval length can either be explicitly
   *  given or be deduced by multiplying the ratio with the last interval
   *  length in the container
   */
  template<class ctype>
  void geometric_fill_intervals(std::vector<ctype>& v, int n, ctype ratio, ctype h0 = static_cast<ctype>(0))
  {
    if (v.empty())
      v.push_back(static_cast<ctype>(0));
    ctype h = h0;
    if (h0 < 1e-8)
    {
      if (v.size() < 2)
        DUNE_THROW(GridError,"Not enough elements in coordinate container for geometric_fill_intervals");
      else
        h = (v.back() - v[v.size()-2]) * ratio;
    }
    for (int i=0; i<n; i++)
    {
      v.push_back(v.back() + h);
      h *= ratio;
    }
  }

  /** \brief adds intervals according with a given length ratio until a given coordinate is reached
   * \param v the std::vector to operate on
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
  template<class ctype>
  void geometric_fill_until(std::vector<ctype>& v, ctype ratio, ctype end, ctype h0 = static_cast<ctype>(0))
  {
    if (v.empty())
      v.push_back(static_cast<ctype>(0));
    ctype h = h0;
    if (h0 < 1e-8)
    {
      if (v.size() < 2)
        DUNE_THROW(GridError,"Not enough elements in coordinate container for geometric_fill_until");
      else
        h = (v.back() - v[v.size()-2]) * ratio;
    }
    while (v.back() < end)
    {
      v.push_back(v.back()+h);
      h *= ratio;
    }
  }

  /** this implements a simple newton iteration for the function
   *  \f$f(x) = -x^n+\frac{x_e-x_s}{h} (x-1)+1\f$
   */
  template<class ctype>
  ctype newton(int n,ctype x_s,ctype x_e,ctype h)
  {
    ctype m = (x_e - x_s)/h;
    ctype xold = 0.0;
    ctype xnew = x_e-x_s;
    while (fabs(xnew-xold) > 1E-8)
    {
      xold = xnew;
      xnew = xold - (-pow(xold,n)+m*xold-m+1)/(-n*pow(xold,n-1)+m);
    }
    if (fabs(xnew-1) < 1E-6)
    {
      xold = x_e-x_s;
      xnew = 0.0;
      while (fabs(xnew-xold) > 1E-8)
      {
        xold = xnew;
        xnew = xold - (-pow(xold,n)+m*xold-m+1)/(-n*pow(xold,n-1)+m);
      }
    }
    return xnew;
  }

  /** \brief fills a coordinate range with n intervals accoring to a geometric series
   * \param v the std::vector to operate on
   * \param n the number of intervals to add
   * \param end the coordinate of the right border of the range
   * \param h the interval length to start or end with (see below) (optional)
   * \param first true if the given h is to be the first interval, false if last one
   */
  template<class ctype>
  void geometric_fill_range(std::vector<ctype>& v, int n, ctype end, ctype h = static_cast<ctype>(0), bool first=true)
  {

  }
}

#endif
