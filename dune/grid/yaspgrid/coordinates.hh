// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_YASPGRID_COORDINATES_HH
#define DUNE_GRID_YASPGRID_COORDINATES_HH

#include <array>
#include <bitset>
#include <vector>

#include <dune/common/fvector.hh>

/** \file
 *  \brief This provides container classes for the coordinates to be used in YaspGrid
 *  Upon implementation of the tensorproduct feature, the coordinate information has
 *  been encapsulated to keep performance for the equidistant grid. Containers for
 *  equidistant and tensorproduct grids are provided here.
 */

namespace Dune
{
  /** \brief Container for equidistant coordinates in a YaspGrid
   *  @tparam ct the coordinate type
   *  @tparam dim the dimension of the grid
   */
  template<class ct, int dim>
  class EquidistantCoordinates
  {
    public:
    //! export the coordinate type
    typedef ct ctype;
    //! export dimension
    static const int dimension = dim;

    /** \brief default constructor */
    EquidistantCoordinates() {}

    /** \brief construct a container with all necessary information
     *  \param upperRight upper right corner of the domain
     *  \param s the size (in codim 0 elements) of the grid on this processor
     *  the size information is kept with this container, because this is the natural
     *  way to handle this for a tensorproduct grid.
     */
    EquidistantCoordinates(const Dune::FieldVector<ct,dim>& upperRight, const std::array<int,dim>& s)
      : _s(s)
    {
      for (int i=0; i<dim; i++)
        _h[i] = upperRight[i] / _s[i];
    }

    /** \returns the meshsize in given direction at given position
     *  \param d the direction to be used
     *  \param i the global coordinate index where to return the meshsize
     */
    inline ct meshsize(int d, [[maybe_unused]] int i) const
    {
      return _h[d];
    }

    /** \returns a coordinate given a direction and an index
     *  \param d the direction to be used
     *  \param i the global coordinate index
     */
    inline ct coordinate(int d, int i) const
    {
      return i*_h[d];
    }

    /** \returns the size in given direction
     *  \param d the direction to be used
     */
    inline int size(int d) const
    {
      return _s[d];
    }

    /** \returns a container that represents the same grid after one step of uniform refinement
     *  \param ovlp_low whether we have an overlap area at the lower processor boundary
     *  \param ovlp_up whether we have an overlap area at the upper processor boundary
     *  \param overlap the size of the overlap region
     *  \param keep_ovlp the refinement option parameter to be used
     */
    EquidistantCoordinates<ct,dim> refine(std::bitset<dim> ovlp_low, std::bitset<dim> ovlp_up, int overlap, bool keep_ovlp) const
    {
      //determine new size and meshsize
      std::array<int,dim> news;
      Dune::FieldVector<ct,dim> newUpperRight;

      for (int i=0; i<dim; i++)
      {
        news[i] = 2 * _s[i];
        if (!keep_ovlp)
        {
          if (ovlp_low[i])
            news[i] -= overlap;
          if (ovlp_up[i])
            news[i] -= overlap;
        }

        newUpperRight[i] = (_h[i] / ct(2.)) * news[i];
      }
      return EquidistantCoordinates<ct,dim>(newUpperRight,news);
    }

    /** \brief print information on this container */
    void print(std::ostream& s) const
    {
      s << "Printing equidistant coordinate information:" << std::endl;
      s << "Meshsize: " << _h << std::endl << "Size: " << _s << std::endl;
    }

    private:
    Dune::FieldVector<ct,dim> _h;
    std::array<int,dim> _s;
  };

  template<class ct, int dim>
  inline std::ostream& operator<< (std::ostream& s, EquidistantCoordinates<ct,dim>& c)
  {
    c.print(s);
    return s;
  }

  /** \brief Container for equidistant coordinates in a YaspGrid with non-trivial origin
    *  @tparam ct the coordinate type
    *  @tparam dim the dimension of the grid
    */
   template<class ct, int dim>
   class EquidistantOffsetCoordinates
   {
     public:
     //! export the coordinate type
     typedef ct ctype;
     //! export dimension
     static const int dimension = dim;

     /** \brief default constructor */
     EquidistantOffsetCoordinates() {}

     /** \brief construct a container with all necessary information
      *  \param lowerLeft the lower left corner of the grid
      *  \param upperRight the upper right corner of the grid
      *  \param s the size (in codim 0 elements) of the grid on this processor
      *
      *  the size information is kept with this container, because this is the natural
      *  way to handle this for a tensorproduct grid.
      */
     EquidistantOffsetCoordinates(const Dune::FieldVector<ct,dim>& lowerLeft, const Dune::FieldVector<ct,dim>& upperRight, const std::array<int,dim>& s)
       : _origin(lowerLeft), _s(s)
     {
       for (int i=0; i<dim; i++)
         _h[i] = (upperRight[i] - lowerLeft[i]) / s[i];
     }

     /** \returns the meshsize in given direction at given position
      *  \param d the direction to be used
      *  \param i the global coordinate index where to return the meshsize
      */
     inline ct meshsize(int d, [[maybe_unused]] int i) const
     {
       return _h[d];
     }

     /** \returns a coordinate given a direction and an index
      *  \param d the direction to be used
      *  \param i the global coordinate index
      */
     inline ct coordinate(int d, int i) const
     {
       return _origin[d] + i*_h[d];
     }

     /** \returns the size in given direction
      *  \param d the direction to be used
      */
     inline int size(int d) const
     {
       return _s[d];
     }

     /** \returns the dth component of the origin
      *  \param d the direction to be used
      */
     inline ct origin(int d) const
     {
       return _origin[d];
     }

     /** \returns a container that represents the same grid after one step of uniform refinement
      *  \param ovlp_low whether we have an overlap area at the lower processor boundary
      *  \param ovlp_up whether we have an overlap area at the upper processor boundary
      *  \param overlap the size of the overlap region
      *  \param keep_ovlp the refinement option parameter to be used
      */
     EquidistantOffsetCoordinates<ct,dim> refine(std::bitset<dim> ovlp_low, std::bitset<dim> ovlp_up, int overlap, bool keep_ovlp) const
     {
       //determine new size and meshsize
       std::array<int,dim> news;
       Dune::FieldVector<ct,dim> newUpperRight;

       for (int i=0; i<dim; i++)
       {
         news[i] = 2 * _s[i];
         if (!keep_ovlp)
         {
           if (ovlp_low[i])
             news[i] -= overlap;
           if (ovlp_up[i])
             news[i] -= overlap;
         }

         newUpperRight[i] = _origin[i] + (_h[i] / ct(2.)) * news[i];
       }
       return EquidistantOffsetCoordinates<ct,dim>(_origin,newUpperRight,news);
     }

     /** \brief print information on this container */
     void print(std::ostream& s) const
     {
       s << "Printing equidistant coordinate information:" << std::endl;
       s << "Meshsize: " << _h << std::endl << "Size: " << _s << std::endl;
       s << "Offset to origin: " << _origin << std::endl;
     }

     private:
     Dune::FieldVector<ct,dim> _origin;
     Dune::FieldVector<ct,dim> _h;
     std::array<int,dim> _s;
   };

   template<class ct, int dim>
   inline std::ostream& operator<< (std::ostream& s, EquidistantOffsetCoordinates<ct,dim>& c)
   {
     c.print(s);
     return s;
   }

  /** \brief Coordinate container for a tensor product YaspGrid
   *  @tparam ct the coordinate type
   *  @tparam dim the dimension of the grid
   */
  template<class ct, int dim>
  class TensorProductCoordinates
  {
    public:
    //! export the coordinate type
    typedef ct ctype;
    //! export dimension
    static const int dimension = dim;

    /** \brief the default constructor */
    TensorProductCoordinates() {}

    /** \brief construct a container with all necessary information
     *  \param c the array of coordinate vectors
     *  \param offset the offset between global origin and processor origin
     *  the size information is deduced from c. Storing offset allows for use of
     *  global coordinates in the YaspGrid code.
     */
    TensorProductCoordinates(const std::array<std::vector<ct>,dim>& c, const std::array<int,dim>& offset)
      : _c(c),_offset(offset)
    {}

    /** \returns the meshsize in given direction at given position
     *  \param d the direction to be used
     *  \param i the global coordinate index where to return the meshsize
     */
    inline ct meshsize(int d, int i) const
    {
      return _c[d][i+1-_offset[d]] - _c[d][i-_offset[d]];
    }

    /** \returns a coordinate given a direction and an index
     *  \param d the direction to be used
     *  \param i the global coordinate index
     */
    inline ct coordinate(int d, int i) const
    {
      return _c[d][i-_offset[d]];
    }

    /** \returns the size in given direction
     *  \param d the direction to be used
     */
    inline int size(int d) const
    {
      return _c[d].size() - 1;
    }

    /** \returns a container that represents the same grid after one step of uniform refinement
     *  \param ovlp_low whether we have an overlap area at the lower processor boundary
     *  \param ovlp_up whether we have an overlap area at the upper processor boundary
     *  \param overlap the size of the overlap region
     *  \param keep_ovlp the refinement option parameter to be used
     */
    TensorProductCoordinates<ct,dim> refine(std::bitset<dim> ovlp_low, std::bitset<dim> ovlp_up, int overlap, bool keep_ovlp) const
    {
      std::array<std::vector<ct>,dim> newcoords;
      std::array<int,dim> newoffset(_offset);
      for (int i=0; i<dim; i++)
      {
        newoffset[i] *= 2;

        //determine new size
        int newsize = 2 * _c[i].size() - 1;
        if (!keep_ovlp)
        {
          if (ovlp_low[i])
          {
            newoffset[i] += overlap;
            newsize -= overlap;
          }
          if (ovlp_up[i])
            newsize -= overlap;
        }
        newcoords[i].resize(newsize);

        typename std::vector<ct>::const_iterator it = _c[i].begin();
        typename std::vector<ct>::const_iterator end = _c[i].end()-1;
        typename std::vector<ct>::iterator iit = newcoords[i].begin() - 1;
        if (!keep_ovlp)
        {
          if (ovlp_low[i])
          {
            it += overlap/2;
            if (overlap%2)
              *(++iit) = (*it + *(++it)) / ct(2.);
          }
          if (ovlp_up[i])
            end -= overlap/2;
        }

        for (;it!=end;)
        {
          *(++iit) = *it;
          *(++iit) = (*it + *(++it)) / ct(2.);
        }

        if (++iit != newcoords[i].end())
          *iit = *it;
      }
      return TensorProductCoordinates<ct,dim>(newcoords, newoffset);
    }

    /** \brief print information on this container */
    void print(std::ostream& s) const
    {
      s << "Printing TensorProduct Coordinate information:" << std::endl;
      for (int i=0; i<dim; i++)
      {
        s << "Direction " << i << ": " << _c[i].size() << " coordinates" << std::endl;
        for (std::size_t j=0; j<_c[i].size(); j++)
          s << _c[i][j] << std::endl;
      }
    }

    private:
    std::array<std::vector<ct>,dim> _c;
    std::array<int,dim> _offset;
  };

  template<class ct, int dim>
  inline std::ostream& operator<< (std::ostream& s, TensorProductCoordinates<ct,dim>& c)
  {
    c.print(s);
    return s;
  }

 namespace Yasp {
  template<class ctype, std::size_t dim>
  bool checkIfMonotonous(const std::array<std::vector<ctype>, dim>& coords)
  {
    for (std::size_t i=0; i<dim; i++)
    {
      if (coords[i].size() <= 1)
        return false;
      for (std::size_t j=1; j<coords[i].size(); j++)
        if (coords[i][j] < coords[i][j-1])
          return false;
    }
    return true;
  }
 } // namespace Yasp
} // namespace Dune

#endif
