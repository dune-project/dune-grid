// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_YASPGRID_YGRID_HH
#define DUNE_GRID_YASPGRID_YGRID_HH

#include <vector>

#include <dune/common/array.hh>
#include <dune/common/fvector.hh>

/** \file
    \brief This provides a YGrid, the elemental component of the yaspgrid implementation
 */

namespace Dune {

  // forward declarations
  template<class CC, int d, typename ct> class YGrid;

  //! define a tolerance value for coordinate comparisons
  static const double Ytolerance=1E-13;

  /** @returns an array containing the sizes of the grids associated with vectors in given array.
   *  Needed in this form due to the need of such functionality in class initializer lists.
   *  @param v the array of vectors to examine
   *  @param r the shift vector to determine grid type
   */
  template<int d, typename ct>
  Dune::array<int,d> sizeArray(Dune::array<std::vector<ct>,d> v, Dune::FieldVector<ct,d> r)
  {
    Dune::array<int,d> tmp;
    for (int i=0; i<d; ++i)
      if (r[i] < Ytolerance)
        tmp[i] = v[i].size();
      else
        tmp[i] = v[i].size() - 1;
    return tmp;
  }

  /**
     The YGrid considered here describes a finite set \f$d\f$-tupels of the form
     \f[ G = \{ (k_0,\ldots,k_{d-1}) | o_i \leq k_i < o_i+s_i \}  \f]

     together with an affine mapping

     \f[ t : G \to R^d, \ \ \ t(k)_i = k_i h_i + r_i \f].

     Therefore a YGrid is characterized by the following quantities:

     - The origin \f$ o=(o_0,\ldots,o_{d-1}) \in Z^d\f$,
     - the size \f$ s=(s_0,\ldots,s_{d-1}) \in Z^d\f$,
     - The shift \f$ r=(r_0,\ldots,r_{d-1}) \in R^d\f$.
     - a coordinate container (see coordinates.hh)
     The shift can be used to interpret the points of a grid as midpoints of cells, faces, edges, etc.


     The YGrid can be parametrized by the dimension d and the type to be used for the coordinates.

     Here is a graphical illustration of a grid:

     \image html  grid.png "A YGrid."
     \image latex grid.eps "A YGrid." width=\textwidth

     A YGrid allows to iterate over all its cells with an Iterator class.

     A YGrid is always considered as being embedded in a larger grid.
     This embedding is characterized by an offset and an enclosing grid as
     shown in the following picture:

     \image html  subgrid.png "The SubYGrid is shown in red, blue is the enclosing grid."
     \image latex subgrid.eps "The SubYGrid is shown in red, blue is the enclosing grid." width=\textwidth

     The iterator provides also a mapping to the consecutive index in the enclosing grid.

     Note: as of november 2013 there are only YGrid and YGrid::Iterator. These represent
     the functionality of former SubYGrid and SubYGrid::TransformingSubIterator. All other
     classes in the hierarchy have not been used.
   */
  template<class CC, int d, typename ct>
  class YGrid
  {
  public:
    typedef Dune::array<int, d> iTupel;
    typedef FieldVector<ct,d> fTupel;

    //! make uninitialized ygrid
    YGrid () : _shift(0.0)
    {
      std::fill(_origin.begin(), _origin.end(), 0);
      std::fill(_offset.begin(), _offset.end(), 0);
      std::fill(_size.begin(), _size.end(), 0);
    }

    /** @brief make ygrid without coordinate information
     *  @param origin origin of the grid in global coordinates
     *  @param shift shift vector to be used for this grid
     *  @param size size of the grid
     *  Such grid has no coordinate information stored but can be
     *  used to determine an intersection with a grid with coordinate
     *  information. This avoids sending coordinates in the parallel case.
     */
    YGrid(iTupel origin, iTupel size, fTupel shift)
      : _origin(origin), _shift(shift), _size(size)
    {}

    /** @brief make a subgrid by taking coordinates from a larger grid
     *  @param origin origin of the grid to be constructed
     *  @param size size of the grid to be constructed
     *  @param enclosing the grid to take coordinates and shift vector from
     */
    YGrid (iTupel origin, iTupel size, const YGrid<CC,d,ct>& enclosing)
      :  _origin(origin), _shift(enclosing.shift()), _coords(enclosing.getCoords()), _size(size), _supersize(enclosing.supersize())
    {
      for (int i=0; i<d; i++)
        _offset[i] = origin[i] - enclosing.origin(i) + enclosing.offset(i);
    }

    /** @brief Make YGrid by giving all parameters
     *  @param origin the origin of the grid in global coordinates
     *  @param shift the shift vector
     *  @param coords the coordinate vectors to be used
     *  @param size the size vector
     *  @param offset the offset in the enclosing grid
     *  @param supersize size of the enclosing grid
     */
    YGrid (iTupel origin,  fTupel shift, CC* coords, iTupel size, iTupel offset, iTupel supersize)
      : _origin(origin), _shift(shift), _coords(coords), _size(size), _offset(offset), _supersize(supersize)
    {}

    //! Return origin in direction i
    int origin (int i) const
    {
      return _origin[i];
    }

    //! return reference to origin
    const iTupel& origin () const
    {
      return _origin;
    }

    //! Return shift in direction i
    ct shift (int i) const
    {
      return _shift[i];
    }

    //! Return shift tupel
    const fTupel& shift () const
    {
      return _shift;
    }

    CC* getCoords() const
    {
      return _coords;
    }

    //! Return offset to origin of enclosing grid
    int offset (int i) const
    {
      return _offset[i];
    }

    //! Return offset to origin of enclosing grid
    const iTupel & offset () const
    {
      return _offset;
    }

    //! return size of enclosing grid
    int supersize (int i) const
    {
      return _supersize[i];
    }

    //! return size of enclosing grid
    const iTupel & supersize () const
    {
      return _supersize;
    }

    //! return size in direction i
    int size (int i) const
    {
      return _size[i];
    }

    //! retrun size
    iTupel size () const
    {
      return _size;
    }

    //! Return total size of index set which is the product of all size per direction.
    int totalsize () const
    {
      int s=1;
      for (int i=0; i<d; ++i)
        s *= size(i);
      return s;
    }

    //! Return minimum index in direction i
    int min (int i) const
    {
      return _origin[i];
    }

    //! Return maximum index in direction i
    int max (int i) const
    {
      return _origin[i] + size(i) - 1;
    }

    //! Return true if YGrid is empty, i.e. has size 0 in all directions.
    bool empty () const
    {
      for (int i=0; i<d; ++i)
      {
        if (size(i) == 0)
          return true;
      }
      return false;
    }

    //! given a coordinate, return true if it is in the grid
    bool inside (const iTupel& coord) const
    {
      for (int i=0; i<d; i++)
      {
        if ((coord[i]<_origin[i]) || (coord[i]>=_origin[i]+_size[i]))
          return false;
      }
      return true;
    }

    //! given a tupel compute its index in the lexicographic numbering
    int index (const iTupel& coord) const
    {
      int index = (coord[d-1]-_origin[d-1]);

      for (int i=d-2; i>=0; i--)
        index = index*_size[i] + (coord[i]-_origin[i]);

      return index;
    }

    //! return grid moved by the vector v
    YGrid<CC,d,ct> move (iTupel v) const
    {
      for (int i=0; i<d; i++)
        v[i] += _origin[i];
      return YGrid<CC,d,ct>(v,_size,*this);
    }

    //! Return SubYGrid of supergrid of self which is the intersection of self and another YGrid
    YGrid<CC,d,ct> intersection (const YGrid<CC,d,ct>& r) const
    {
      for (int i=0; i<d; i++)
      {
        //empty coordinate vectors result in empty intersections
        if (empty() || r.empty())
          return YGrid<CC,d,ct>();

        //intersectable grids must have the same shift
        if (std::abs(shift(i)-r.shift(i)) > Ytolerance)
          return YGrid<CC,d,ct>();
      }

      iTupel neworigin;
      iTupel newsize;
      for (int i=0; i<d; ++i)
      {
        neworigin[i] = std::max(origin(i),r.origin(i));
        newsize[i] = std::min(max(i),r.max(i)) - neworigin[i] + 1;
      }

      return YGrid<CC,d,ct>(neworigin,newsize,*this);
    }


   /*! Iterator class allows one to run over all cells of a grid.
       The cells of the grid to iterate over are numbered consecutively starting
       with zero. Via the index() method the iterator provides a mapping of the
       cells of the grid to a one-dimensional array. The number of entries
       in this array must be the size of the grid.
     */
    class Iterator {
    public:
      //! Make iterator pointing to first cell in a grid.
      Iterator (const YGrid<CC,d,ct>& r) : _grid(&r)
      {
        iTupel coord(r.origin());
        reinit(r,coord);
      }

      //! Make iterator pointing to given cell in a grid.
      Iterator (const YGrid<CC,d,ct>& r, const iTupel& coord)
      {
        reinit(r,coord);
      }

      //! reinitialize iterator to given position
      void reinit (const YGrid<CC,d,ct>& r, const iTupel& coord)
      {
        // compute increments;
        int inc = 1;
        for (int i=0; i<d; ++i)
        {
          _increment[i] = inc;
          inc *= r.size(i);
        }

        // initialize to given position in index set
        for (int i=0; i<d; ++i)
          _coord[i] = coord[i];
        _index = r.index(coord);

        // compute superincrements
        inc = 1;
        for (int i=0; i<d; ++i)
        {
          _superincrement[i] = inc;
          inc *= r.supersize(i);
        }

        // move superindex to first cell in subgrid
        _superindex = 0;
        for (int i=0; i<d; ++i)
          _superindex += (r.offset(i)+coord[i]-r.origin(i))*_superincrement[i];

        _grid = &r;
        if (_grid->inside(coord))
        {
          for (int i=0; i<d; ++i)
          {
            if (!_grid->empty())
              _begin[i] = _grid->getCoords()->coordinate(i,_grid->origin(i));
            if ((_grid->getCoords()->size(i) > 0) && (_grid->shift(i) > Ytolerance))
              _begin[i] += _grid->shift(i) * _grid->getCoords()->meshsize(i,_grid->origin(i));
            _position[i] = _grid->getCoords()->coordinate(i,_coord[i]);
            if ((_grid->getCoords()->size(i) > 0) && (_grid->shift(i) > Ytolerance))
              _position[i] += _grid->shift(i)* _grid->getCoords()->meshsize(i,_coord[i]);
          }
        }
      }

      //! Return true when two iterators over the same grid are equal (!).
      bool operator== (const Iterator& i) const
      {
        return _superindex == i._superindex;
      }

      //! Return true when two iterators over the same grid are not equal (!).
      bool operator!= (const Iterator& i) const
      {
        return _superindex != i._superindex;
      }

      //! Return index of the current cell in the consecutive numbering.
      int index () const
      {
        return _index;
      }

      //! Return consecutive index in enclosing grid
      int superindex () const
      {
        return _superindex;
      }

      //! Return coordinate of the cell in direction i.
      int coord (int i) const
      {
        return _coord[i];
      }

      //! Return coordinate of the cell as reference (do not modify).
      const iTupel& coord () const
      {
        return _coord;
      }

      //! move this iterator dist cells in direction i
      void move (int i, int dist)
      {
        _coord[i] += dist;
        _index += dist*_increment[i];
        _superindex += dist*_superincrement[i];
        if (_grid->inside(_coord))
        {
          _position[i] = _grid->getCoords()->coordinate(i,_coord[i]);
          if (_grid->shift(i) > Ytolerance)
            _position[i] += _grid->shift(i) * meshsize(i);
        }
      }

      //! Increment iterator to next cell with position.
      Iterator& operator++ ()
      {
        ++_index;               // update consecutive index in subgrid
        for (int i=0; i<d; i++)         // check for wrap around
        {
          _superindex += _superincrement[i];   // move on cell in direction i
          if (++_coord[i] <= _grid->max(i))
          {
            _position[i] = _grid->getCoords()->coordinate(i,_coord[i]);
            if (_grid->shift(i) > Ytolerance)
              _position[i] += _grid->shift(i) * _grid->getCoords()->meshsize(i,_coord[i]);
            return *this;
          }
          else
          {
            _coord[i] = _grid->origin(i);         // move back to origin in direction i
            _superindex -= _grid->size(i) * _superincrement[i];
            _position[i] = _begin[i];
          }
        }
        // if we wrapped around, back to to begin(), we must put the iterator to end()
        if (_coord == _grid->origin())
        {
          for (int i=0; i<d; i++)
            _superindex += (_grid->size(i)-1) * _superincrement[i];
          _superindex += _superincrement[0];
        }
        return *this;
      }

      //! Return position of current cell in direction i.
      ct position (int i) const
      {
        return _position[i];
      }

      //! Return position of current cell as reference.
      const fTupel& position () const
      {
        return _position;
      }

      //! Return meshsize in direction i
      ct meshsize (int i) const
      {
        return _grid->getCoords()->meshsize(i,_coord[i]);
      }

      //! Return meshsize of current cell as reference.
      fTupel meshsize () const
      {
        fTupel h;
        for (int i=0; i<d; i++)
          h[i] = meshsize(i);
        return h;
      }

    protected:
      int _index;          //!< current lexicographic position in index set
      iTupel _coord;       //!< current position in index set
      iTupel _increment;   //!< increment for next neighbor in direction i
      int _superindex;        //!< consecutive index in enclosing grid
      iTupel _superincrement; //!< moves consecutive index by one in this direction in supergrid
      const YGrid<CC,d,ct>* _grid;
      fTupel _begin;    //!< position of origin of grid
      fTupel _position; //!< current position
    };

    //! return iterator to first element of index set
    Iterator begin () const
    {
      return Iterator(*this);
    }

    //! return iterator to given element of index set
    Iterator begin (iTupel& co) const
    {
      return Iterator(*this,co);
    }

    //! return subiterator to last element of index set
    Iterator end () const
    {
      iTupel last;
      for (int i=0; i<d; i++)
        last[i] = max(i);
      last[0] += 1;
      return Iterator(*this,last);
    }

  private:
    iTupel _origin;
    fTupel _shift;
    CC* _coords;
    iTupel _size;
    iTupel _offset;    //!< offset to origin of the enclosing grid
    iTupel _supersize; //!< size of the enclosing grid
  };


  //! Output operator for ygrids
  template <class CC, int d, typename ct>
  inline std::ostream& operator<< (std::ostream& s, YGrid<CC,d,ct> e)
  {
    s << "Printing YGrid structure:" << std::endl;
    s << "Origin: " << e.origin() << std::endl;
    s << "Shift: " << e.shift() << std::endl;
    s << "Size: " << e.size() << std::endl;
    s << "Offset: " << e.offset() << std::endl;
    s << "Supersize: " << e.supersize() << std::endl;
    return s;
  }

  //! Output operator for ygrids
  template <class CC, int d, typename ct>
  inline std::ostream& operator<< (std::ostream& s, typename YGrid<CC,d,ct>::Iterator& e)
  {
    s << "please reimplement this" << std::endl;
    return s;
  }


} // namespace Dune

#endif
