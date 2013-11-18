// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_YGRIDS_HH
#define DUNE_YGRIDS_HH

// C++ includes
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <deque>
#include <bitset>

// C includes
#if HAVE_MPI
#include <mpi.h>
#endif
#include <string.h>

// local includes
#include <dune/common/fvector.hh>
#include <dune/common/stdstreams.hh>
#include <dune/common/power.hh>
#include <dune/grid/common/grid.hh>

/** \file
    \brief This is the basis for the yaspgrid implementation of the Dune grid interface.
 */

namespace Dune {

  // forward declarations
  template<int d, typename ct> class YGrid;
  template<int d, typename ct> class SubYGrid;

  static const double Ytolerance=1E-13;

  /**
     This is the basis of a parallel implementation of the dune grid interface
     supporting codim 0 and dim.

     You can also use the structured interface and write really fast code.

     The YGrid considered here describes a finite set \f$d\f$-tupels of the form
     \f[ G = \{ (k_0,\ldots,k_{d-1}) | o_i \leq k_i < o_i+s_i \}  \f]

     together with an affine mapping

     \f[ t : G \to R^d, \ \ \ t(k)_i = k_i h_i + r_i \f].

     Therefore a YGrid is characterized by the following four quantities:

     - The origin \f$ o=(o_0,\ldots,o_{d-1}) \in Z^d\f$,
     - the size \f$ s=(s_0,\ldots,s_{d-1}) \in Z^d\f$,
     - the mesh width \f$ h=(h_0,\ldots,h_{d-1}) \in R^d\f$,
     - The shift \f$ r=(r_0,\ldots,r_{d-1}) \in R^d\f$. The shift can be used to interpret the
     points of a grid as midpoints of cells, faces, edges, etc.

     The YGrid can be parametrized by the dimension d and the type to be used for the coordinates.

     Here is a graphical illustration of a grid:

     \image html  grid.png "A YGrid."
     \image latex grid.eps "A YGrid." width=\textwidth

     A grid can be manipulated either in the origin/size representation or in the
     min index / max index representation.

     A YGrid allows to iterate over all its cells with an Iterator class.
   */
  template<int d, typename ct>
  class YGrid {
  public:
    //! define types used for arguments
    typedef FieldVector<int, d>  iTupel;
    typedef FieldVector<ct, d> fTupel;

    //! Destructor
    virtual ~YGrid()
    {}

    //! Make an empty YGrid with origin 0
    YGrid () :
      _origin(0), _size(0), _h(0.0), _r(0.0)
    {}

    //! Make YGrid from origin and size arrays
    YGrid (iTupel o, iTupel s, fTupel h, fTupel r) :
      _origin(o), _size(s), _h(h), _r(r)
    {
#ifndef NDEBUG
      for (int i=0; i<d; ++i)
        assert (_size[i] >= 0);
#endif
    }

    //! Return origin in direction i
    int origin (int i) const
    {
      return _origin[i];
    }

    //! Set origin in direction i
    void origin (int i, int oi) const
    {
      _origin[i] = oi;
    }

    //! return reference to origin
    const iTupel& origin () const
    {
      return _origin;
    }

    //! Return size in direction i
    int size (int i) const
    {
      return _size[i];
    }

    //! Set size in direction i
    void size (int i, int si) const
    {
      _size[i] = si;
      if (_size[i]<0) _size[i] = 0;
    }

    //! Return reference to size tupel
    const iTupel& size () const
    {
      return _size;
    }

    //! Return total size of index set which is the product of all size per direction.
    int totalsize () const
    {
      int s=1;
      for (int i=0; i<d; ++i) s=s*_size[i];
      return s;
    }

    //! Return minimum index in direction i
    int min (int i) const
    {
      return _origin[i];
    }

    //! Set minimum index in direction i
    void min (int i, int mi) const
    {
      _size[i] = max(i)-mi+1;
      _origin[i] = mi;
      if (_size[i]<0) _size[i] = 0;
    }

    //! Return maximum index in direction i
    int max (int i) const
    {
      return _origin[i]+_size[i]-1;
    }

    //! Set maximum index in direction i
    void max (int i, int mi) const
    {
      _size[i] = mi-min(i)+1;
      if (_size[i]<0) _size[i] = 0;
    }

    //! Return reference to mesh size tupel for read write access
    const fTupel& meshsize () const
    {
      return _h;
    }

    //! Return mesh size in direction i
    ct meshsize (int i) const
    {
      return _h[i];
    }

    //! Set mesh size in direction i
    void meshsize (int i, int hi) const
    {
      _h[i] = hi;
    }

    //! Return shift tupel
    const fTupel& shift () const
    {
      return _r;
    }

    //! Return shift in direction i
    ct shift (int i) const
    {
      return _r[i];
    }

    //! Set shift in direction i
    void shift (int i, int ri) const
    {
      _r[i] = ri;
    }

    //! Return true if YGrid is empty, i.e. has size 0 in all directions.
    bool empty () const
    {
      for (int i=0; i<d; ++i) if (_size[i]<=0) return true;
      return false;
    }

    //! given a tupel compute its index in the lexicographic numbering
    int index (const iTupel& coord) const
    {
      int index = (coord[d-1]-_origin[d-1]);

      for (int i=d-2; i>=0; i--)
        index = index*_size[i] + (coord[i]-_origin[i]);

      return index;
    }

    //! given a tupel compute its index in the lexicographic numbering
    int index (const array<int,d>& coord) const
    {
      int index = (coord[d-1]-_origin[d-1]);

      for (int i=d-2; i>=0; i--)
        index = index*_size[i] + (coord[i]-_origin[i]);

      return index;
    }

    //! given a coordinate, return true if it is in the grid
    bool inside (const iTupel& coord) const
    {
      for (int i=0; i<d; i++)
      {
        if (coord[i]<_origin[i] || coord[i]>=_origin[i]+_size[i]) return false;
      }
      return true;
    }

    //! return grid moved by the vector v
    YGrid<d,ct> move (iTupel v) const
    {
      for (int i=0; i<d; i++) v[i] += _origin[i];
      return YGrid<d,ct>(v,_size,_h,_r);
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
      Iterator (const YGrid<d,ct>& r)
      {
        // copy data coming from grid to iterate over
        for (int i=0; i<d; ++i) _origin[i] = r.origin(i);
        for (int i=0; i<d; ++i) _end[i] = r.origin(i)+r.size(i)-1;

        // initialize to first position in index set
        for (int i=0; i<d; ++i) _coord[i] = _origin[i];
        _index = 0;

        // compute increments;
        int inc = 1;
        for (int i=0; i<d; ++i)
        {
          _increment[i] = inc;
          inc *= r.size(i);
        }
      }

      //! Make iterator pointing to given cell in a grid.
      Iterator (const YGrid<d,ct>& r, const iTupel& coord)
      {
        // copy data coming from grid to iterate over
        for (int i=0; i<d; ++i) _origin[i] = r.origin(i);
        for (int i=0; i<d; ++i) _end[i] = r.origin(i)+r.size(i)-1;

        // compute increments;
        int inc = 1;
        for (int i=0; i<d; ++i)
        {
          _increment[i] = inc;
          inc *= r.size(i);
        }

        // initialize to given position in index set
        for (int i=0; i<d; ++i) _coord[i] = coord[i];
        _index = r.index(coord);
      }

      //! reinitialize iterator to given position
      void reinit (const YGrid<d,ct>& r, const iTupel& coord)
      {
        // copy data coming from grid to iterate over
        for (int i=0; i<d; ++i) _origin[i] = r.origin(i);
        for (int i=0; i<d; ++i) _end[i] = r.origin(i)+r.size(i)-1;

        // compute increments;
        int inc = 1;
        for (int i=0; i<d; ++i)
        {
          _increment[i] = inc;
          inc *= r.size(i);
        }

        // initialize to given position in index set
        for (int i=0; i<d; ++i) _coord[i] = coord[i];
        _index = r.index(coord);
      }

      //! reinitialize iterator to given position
      void reinit (const YGrid<d,ct>& r, const array<int,d>& coord)
      {
        // copy data coming from grid to iterate over
        for (int i=0; i<d; ++i) _origin[i] = r.origin(i);
        for (int i=0; i<d; ++i) _end[i] = r.origin(i)+r.size(i)-1;

        // compute increments;
        int inc = 1;
        for (int i=0; i<d; ++i)
        {
          _increment[i] = inc;
          inc *= r.size(i);
        }

        // initialize to given position in index set
        for (int i=0; i<d; ++i) _coord[i] = coord[i];
        _index = r.index(coord);
      }

      //! Return true when two iterators over the same grid are equal (!).
      bool operator== (const Iterator& i) const
      {
        return _index == i._index;
      }

      //! Return true when two iterators over the same grid are not equal (!).
      bool operator!= (const Iterator& i) const
      {
        return _index != i._index;
      }

      //! Return index of the current cell in the consecutive numbering.
      int index () const
      {
        return _index;
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

      //! Get index of cell which is dist cells away in direction i.
      int neighbor (int i, int dist) const
      {
        return _index+dist*_increment[i];
      }

      //! Get index of neighboring cell which is -1 away in direction i.
      int down (int i) const
      {
        return _index-_increment[i];
      }

      //! Get index of neighboring cell which is +1 away in direction i.
      int up (int i) const
      {
        return _index+_increment[i];
      }

      //! move this iterator dist cells in direction i
      void move (int i, int dist)
      {
        _coord[i] += dist;
        _index += dist*_increment[i];
      }

      //! Increment iterator to next cell.
      Iterator& operator++ ()
      {
        ++_index;
        for (int i=0; i<d; i++)
          if (++(_coord[i])<=_end[i])
            return *this;
          else
            _coord[i]=_origin[i];
        return *this;
      }

      //! Print position of iterator
      void print (std::ostream& s) const
      {
        s << index() << " : [";
        for (int i=0; i<d-1; i++) s << coord(i) << ",";
        s << coord(d-1) << "]";
      }

    protected:
      int _index;          //< current lexicographic position in index set
      iTupel _coord;       //< current position in index set
      iTupel _increment;   //< increment for next neighbor in direction i
      iTupel _origin;      //< origin and
      iTupel _end;         //< last index in direction i
    };

    //! return iterator to first element of index set
    Iterator begin () const
    {
      return Iterator(*this);
    }

    //! return iterator to one past the last element of index set
    Iterator end () const {
      iTupel last;
      for (int i=0; i<d; i++)
        last[i] = max(i);
      last[0] += 1;
      return Iterator(*this,last);
    }

    /*! TransformingIterator is an Iterator providing in addition a linear transformation
       of the coordinates of the grid in the form \f$ y_i = x_i h_i + s_i \f$.
       This can be used to interpret the grid cells as vertices, edges, faces, etc.
     */
    class TransformingIterator : public Iterator {
    public:
      //! Make iterator pointing to first cell in a grid.
      TransformingIterator (const YGrid<d,ct>& r) : Iterator(r)
      {
        for (int i=0; i<d; ++i) _h[i] = r.meshsize(i);
        for (int i=0; i<d; ++i) _begin[i] = r.origin(i)*r.meshsize(i)+r.shift(i);
        for (int i=0; i<d; ++i) _position[i] = _begin[i];
      }

      //! Make iterator pointing to given cell in a grid.
      TransformingIterator (const YGrid<d,ct>& r, iTupel& coord) : Iterator(r,coord)
      {
        for (int i=0; i<d; ++i) _h[i] = r.meshsize(i);
        for (int i=0; i<d; ++i) _begin[i] = r.origin(i)*r.meshsize(i)+r.shift(i);
        for (int i=0; i<d; ++i) _position[i] = coord[i]*r.meshsize(i)+r.shift(i);
      }

      //! Make transforming iterator from iterator (used for automatic conversion of end)
      TransformingIterator (Iterator i) : Iterator(i)
      {}

      //! Increment iterator to next cell with position.
      TransformingIterator& operator++ ()
      {
        ++(this->_index);
        for (int i=0; i<d; i++)
          if (++(this->_coord[i])<=this->_end[i])
          {
            _position[i] += _h[i];
            return *this;
          }
          else
          {
            this->_coord[i]=this->_origin[i];
            _position[i] = _begin[i];
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
        return _h[i];
      }

      //! Return meshsize of current cell as reference.
      const fTupel& meshsize () const
      {
        return _h;
      }

      //! Move cell position by dist cells in direction i.
      void move (int i, int dist)
      {
        Iterator::move(i,dist);
        _position[i] += dist*_h[i];
      }

      //! Print contents of iterator
      void print (std::ostream& s) const
      {
        Iterator::print(s);
        s << " " << _position;
      }

    private:
      fTupel _h;        //!< mesh size per direction
      fTupel _begin;    //!< position of origin of grid
      fTupel _position; //!< current position
    };

    //! return iterator to first element of index set
    TransformingIterator tbegin () const
    {
      return TransformingIterator(*this);
    }

    //! return iterator to one past the last element of the grid
    TransformingIterator tend () const
    {
      iTupel last;
      for (int i=0; i<d; i++) last = max(i);
      last[0] += 1;
      return TransformingIterator(*this,last);
    }

  protected:
    //! internal representation uses origin/size
    iTupel _origin;
    iTupel _size;
    fTupel _h;        //!< mesh size per direction
    fTupel _r;        //!< shift per direction
  };

  //! Output operator for grids
  template <int d, typename ct>
  inline std::ostream& operator<< (std::ostream& s, YGrid<d,ct> e)
  {
    s << "{";
    for (int i=0; i<d-1; i++)
      s << "[" << e.min(i) << "," << e.max(i) << "]x";
    s << "[" << e.min(d-1) << "," << e.max(d-1) << "]";
    s << " = [";
    for (int i=0; i<d-1; i++) s << e.origin(i) << ",";
    s << e.origin(d-1) << "]x[";
    for (int i=0; i<d-1; i++) s << e.size(i) << ",";
    s << e.size(d-1) << "]";
    s << " h=[";
    for (int i=0; i<d-1; i++) s << e.meshsize(i) << ",";
    s << e.meshsize(d-1) << "]";
    s << " r=[";
    for (int i=0; i<d-1; i++) s << e.shift(i) << ",";
    s << e.shift(d-1) << "]";
    s << "}";
    return s;
  }

  //! Output operator for Iterators
  template <int d, typename ct>
  inline std::ostream& operator<< (std::ostream& s, typename YGrid<d,ct>::Iterator& e)
  {
    e.print(s);
    return s;
  }


  /*! A SubYGrid is a grid that is embedded in a larger grid
     It is characterized by an offset and an enclosing grid as
     shown in the following picture:

     \image html  subgrid.png "The SubYGrid is shown in red, blue is the enclosing grid."
     \image latex subgrid.eps "The SubYGrid is shown in red, blue is the enclosing grid." width=\textwidth

     SubYGrid has additional iterators that provide a mapping to
     the consecutive index in the enclosing grid.
   */
  template<int d, typename ct>
  class SubYGrid : public YGrid<d,ct> {
  public:
    typedef typename YGrid<d,ct>::iTupel iTupel;
    typedef typename YGrid<d,ct>::fTupel fTupel;

    //! Destructor
    virtual ~SubYGrid()
    {}

    //! make uninitialized subgrid
    SubYGrid () {}

    //! Make SubYGrid from origin, size, offset and supersize
    SubYGrid (iTupel origin, iTupel size, iTupel offset, iTupel supersize, fTupel h, fTupel r)
      : YGrid<d,ct>::YGrid(origin,size,h,r), _offset(offset), _supersize(supersize)
    {
      for (int i=0; i<d; ++i)
      {
        if (offset[i]<0)
          std::cout << "warning: offset["
          << i <<"] negative in SubYGrid"
          << std::endl;
        if (-offset[i]+supersize[i]<size[i])
          std::cout << "warning: subgrid larger than enclosing grid in direction "
          << i <<" in SubYGrid"
          << std::endl;
      }
    }

    //! Make SubYGrid from YGrid
    SubYGrid (YGrid<d,ct> base) : YGrid<d,ct>(base)
    {
      for (int i=0; i<d; ++i)
      {
        _offset[i] = 0;
        _supersize[i] = this->size(i);
      }
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

    //! Return SubYGrid of supergrid of self which is the intersection of self and another YGrid
    virtual SubYGrid<d,ct> intersection (const YGrid<d,ct>& r) const
    {
      // check if the two grids can be intersected, must have same mesh size and shift
      for (int i=0; i<d; i++)
        if (fabs(this->meshsize(i)-r.meshsize(i))>Ytolerance) return SubYGrid<d,ct>();
      for (int i=0; i<d; i++)
        if (fabs(this->shift(i)-r.shift(i))>Ytolerance) return SubYGrid<d,ct>();

      iTupel neworigin;
      iTupel newsize;
      iTupel offset;

      for (int i=0; i<d; ++i)
      {
        // intersect
        neworigin[i] = std::max(this->min(i),r.min(i));
        newsize[i] = std::min(this->max(i),r.max(i))-neworigin[i]+1;
        if (newsize[i]<0) {
          newsize[i] = 0;
          neworigin[i] = this->min(i);
        }

        // offset to my supergrid
        offset[i] = _offset[i]+neworigin[i]-this->origin(i);
      }
      return SubYGrid<d,ct>(neworigin,newsize,offset,_supersize,this->meshsize(),this->shift());
    }

    /*! SubIterator is an Iterator that provides in addition the consecutive
       index in the enclosing grid.
     */
    class SubIterator : public YGrid<d,ct>::Iterator {
    public:
      //! Make iterator pointing to first cell in subgrid.
      SubIterator (const SubYGrid<d,ct>& r) : YGrid<d,ct>::Iterator::Iterator (r)
      {
        //! store some grid information
        for (int i=0; i<d; ++i) _size[i] = r.size(i);

        // compute superincrements
        int inc = 1;
        for (int i=0; i<d; ++i)
        {
          _superincrement[i] = inc;
          inc *= r.supersize(i);
        }

        // move superindex to first cell in subgrid
        _superindex = 0;
        for (int i=0; i<d; ++i)
          _superindex += r.offset(i)*_superincrement[i];
      }

      //! Make iterator pointing to given cell in subgrid.
      SubIterator (const SubYGrid<d,ct>& r, const iTupel& coord) : YGrid<d,ct>::Iterator::Iterator (r,coord)
      {
        //! store some grid information
        for (int i=0; i<d; ++i) _size[i] = r.size(i);

        // compute superincrements
        int inc = 1;
        for (int i=0; i<d; ++i)
        {
          _superincrement[i] = inc;
          inc *= r.supersize(i);
        }

        // move superindex to first cell in subgrid
        _superindex = 0;
        for (int i=0; i<d; ++i)
          _superindex += (r.offset(i)+coord[i]-r.origin(i))*_superincrement[i];
      }

      //! Make transforming iterator from iterator (used for automatic conversion of end)
      SubIterator (const typename YGrid<d,ct>::Iterator& i) : YGrid<d,ct>::Iterator::Iterator(i)
      {}

      //! Make iterator pointing to given cell in subgrid.
      void reinit (const SubYGrid<d,ct>& r, const iTupel& coord)
      {
        YGrid<d,ct>::Iterator::reinit(r,coord);

        //! store some grid information
        for (int i=0; i<d; ++i) _size[i] = r.size(i);

        // compute superincrements
        int inc = 1;
        for (int i=0; i<d; ++i)
        {
          _superincrement[i] = inc;
          inc *= r.supersize(i);
        }

        // move superindex to first cell in subgrid
        _superindex = 0;
        for (int i=0; i<d; ++i)
          _superindex += (r.offset(i)+coord[i]-r.origin(i))*_superincrement[i];
      }

      //! Make iterator pointing to given cell in subgrid.
      void reinit (const SubYGrid<d,ct>& r, const array<int,d>& coord)
      {
        YGrid<d,ct>::Iterator::reinit(r,coord);

        //! store some grid information
        for (int i=0; i<d; ++i) _size[i] = r.size(i);

        // compute superincrements
        int inc = 1;
        for (int i=0; i<d; ++i)
        {
          _superincrement[i] = inc;
          inc *= r.supersize(i);
        }

        // move superindex to first cell in subgrid
        _superindex = 0;
        for (int i=0; i<d; ++i)
          _superindex += (r.offset(i)+coord[i]-r.origin(i))*_superincrement[i];
      }

      //! Return true when two iterators over the same grid are equal (!).
      bool operator== (const SubIterator& i) const
      {
        return _superindex == i._superindex;
      }

      //! Return true when two iterators over the same grid are not equal (!).
      bool operator!= (const SubIterator& i) const
      {
        return _superindex != i._superindex;
      }

      //! Return consecutive index in enclosing grid
      int superindex () const
      {
        return _superindex;
      }

      //! Get index of cell which is dist cells away in direction i in enclosing grid.
      int superneighbor (int i, int dist) const
      {
        return _superindex+dist*_superincrement[i];
      }

      //! Get index of neighboring cell which is -1 away in direction i in enclosing grid.
      int superdown (int i) const
      {
        return _superindex-_superincrement[i];
      }

      //! Get index of neighboring cell which is +1 away in direction i in enclosing grid.
      int superup (int i) const
      {
        return _superindex+_superincrement[i];
      }

      //! move this iterator dist cells in direction i
      void move (int i, int dist)
      {
        YGrid<d,ct>::Iterator::move(i,dist);    // move base iterator
        _superindex += dist*_superincrement[i]; // move superindex
      }

      //! Increment iterator to next cell in subgrid
      SubIterator& operator++ ()
      {
        ++(this->_index);               // update consecutive index in grid
        for (int i=0; i<d; i++)         // check for wrap around
        {
          _superindex += _superincrement[i];   // move on cell in direction i
          if (++(this->_coord[i])<=this->_end[i])
            return *this;
          else
          {
            this->_coord[i]=this->_origin[i];         // move back to origin in direction i
            _superindex -= _size[i]*_superincrement[i];
          }
        }
        // if we wrapped around, back to to begin(), we must put the iterator to end()
        if (this->_coord == this->_origin)
        {
          for (int i=0; i<d; i++)
            this->_superindex += (this->_size[i]-1)*this->_superincrement[i];
          this->_superindex += this->_superincrement[0];
        }
        return *this;
      }

      //! Print position of iterator
      void print (std::ostream& s) const
      {
        YGrid<d,ct>::Iterator::print(s);
        s << " super=" << superindex();
      }

    protected:
      int _superindex;        //!< consecutive index in enclosing grid
      iTupel _superincrement; //!< moves consecutive index by one in this direction in supergrid
      iTupel _size;           //!< size of subgrid
    };

    //! return subiterator to first element of index set
    SubIterator subbegin () const { return SubIterator(*this); }

    //! return subiterator to last element of index set
    SubIterator subend () const
    {
      iTupel last;
      for (int i=0; i<d; i++) last[i] = this->max(i);
      last[0] += 1;
      return SubIterator(*this,last);
    }

    /*! TransformingSubIterator is a SubIterator providing in addition a linear transformation
       of the coordinates of the grid in the form \f$ y_i = x_i h_i + s_i \f$.
       This can be used to interpret the grid cells as vertices, edges, faces, etc.
     */
    class TransformingSubIterator : public SubIterator {
    public:
      //! Make iterator pointing to first cell in a grid.
      TransformingSubIterator (const SubYGrid<d,ct>& r) : SubIterator(r)
      {
        for (int i=0; i<d; ++i) _h[i] = r.meshsize(i);
        for (int i=0; i<d; ++i) _begin[i] = r.origin(i)*r.meshsize(i)+r.shift(i);
        for (int i=0; i<d; ++i) _position[i] = _begin[i];
      }

      //! Make iterator pointing to given cell in a grid.
      TransformingSubIterator (const SubYGrid<d,ct>& r, const iTupel& coord) : SubIterator(r,coord)
      {
        for (int i=0; i<d; ++i) _h[i] = r.meshsize(i);
        for (int i=0; i<d; ++i) _begin[i] = r.origin(i)*r.meshsize(i)+r.shift(i);
        for (int i=0; i<d; ++i) _position[i] = coord[i]*r.meshsize(i)+r.shift(i);
      }

      //! Make transforming iterator from iterator (used for automatic conversion of end)
      TransformingSubIterator (const SubIterator& i) :
        SubIterator(i)
      {}

      TransformingSubIterator (const TransformingSubIterator & t) :
        SubIterator(t), _h(t._h), _begin(t._begin), _position(t._position)
      {}

      //! Make iterator pointing to given cell in a grid.
      void reinit (const SubYGrid<d,ct>& r, const iTupel& coord)
      {
        SubIterator::reinit(r,coord);
        for (int i=0; i<d; ++i) _h[i] = r.meshsize(i);
        for (int i=0; i<d; ++i) _begin[i] = r.origin(i)*r.meshsize(i)+r.shift(i);
        for (int i=0; i<d; ++i) _position[i] = coord[i]*r.meshsize(i)+r.shift(i);
      }

      //! Make iterator pointing to given cell in a grid.
      void reinit (const SubYGrid<d,ct>& r, const array<int,d>& coord)
      {
        SubIterator::reinit(r,coord);
        for (int i=0; i<d; ++i) _h[i] = r.meshsize(i);
        for (int i=0; i<d; ++i) _begin[i] = r.origin(i)*r.meshsize(i)+r.shift(i);
        for (int i=0; i<d; ++i) _position[i] = coord[i]*r.meshsize(i)+r.shift(i);
      }

      //! Increment iterator to next cell with position.
      TransformingSubIterator& operator++ ()
      {
        ++(this->_index);               // update consecutive index in subgrid
        for (int i=0; i<d; i++)         // check for wrap around
        {
          this->_superindex += this->_superincrement[i];   // move on cell in direction i
          if (++(this->_coord[i])<=this->_end[i])
          {
            _position[i] += _h[i];
            return *this;
          }
          else
          {
            this->_coord[i]=this->_origin[i];         // move back to origin in direction i
            this->_superindex -= this->_size[i]*this->_superincrement[i];
            _position[i] = _begin[i];
          }
        }
        // if we wrapped around, back to to begin(), we must put the iterator to end()
        if (this->_coord == this->_origin)
        {
          for (int i=0; i<d; i++)
            this->_superindex += (this->_size[i]-1)*this->_superincrement[i];
          this->_superindex += this->_superincrement[0];
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
        return _h[i];
      }

      //! Return meshsize of current cell as reference.
      const fTupel& meshsize () const
      {
        return _h;
      }

      //! Move cell position by dist cells in direction i.
      void move (int i, int dist)
      {
        SubIterator::move(i,dist);
        _position[i] += dist*_h[i];
      }

      //! Print contents of iterator
      void print (std::ostream& s) const
      {
        SubIterator::print(s);
        s << " [";
        for (int i=0; i<d-1; i++) s << position(i) << ",";
        s << position(d-1) << "]";
      }

    private:
      fTupel _h;        //!< mesh size per direction
      fTupel _begin;    //!< position of origin of grid
      fTupel _position; //!< current position
    };

    //! return iterator to first element of index set
    TransformingSubIterator tsubbegin () const
    {
      return TransformingSubIterator(*this);
    }

    //! return iterator to given element of index set
    TransformingSubIterator tsubbegin (iTupel& co) const
    {
      return TransformingSubIterator(*this,co);
    }

    //! return subiterator to last element of index set
    TransformingSubIterator tsubend () const
    {
      SubIterator endit = subend();
      return TransformingSubIterator(endit);
    }

  private:
    iTupel _offset;    //!< offset to origin of the enclosing grid
    iTupel _supersize; //!< size of the enclosing grid
  };


  //! Output operator for subgrids
  template <int d, typename ct>
  inline std::ostream& operator<< (std::ostream& s, SubYGrid<d,ct> e)
  {
    YGrid<d,ct> x = e;
    s << x << " ofs=" << e.offset() << " ss=" << e.supersize();
    return s;
  }

  //! Output operator for subgrids
  template <int d, typename ct>
  inline std::ostream& operator<< (std::ostream& s, typename SubYGrid<d,ct>::TransformingSubIterator& e)
  {
    e.print(s);
    return s;
  }

  /** \brief Implement the default load balance strategy of yaspgrid
   */
  template<int d>
  class YLoadBalance
  {
  public:
    typedef FieldVector<int, d>  iTupel;
    virtual ~YLoadBalance() {}
    virtual void loadbalance (const iTupel& size, int P, iTupel& dims) const
    {
      double opt=1E100;
      iTupel trydims;

      optimize_dims(d-1,size,P,dims,trydims,opt);
    }
  private:
    void optimize_dims (int i, const iTupel& size, int P, iTupel& dims, iTupel& trydims, double &opt ) const
    {
      if (i>0) // test all subdivisions recursively
      {
        for (int k=1; k<=P; k++)
          if (P%k==0)
          {
            // P divisible by k
            trydims[i] = k;
            optimize_dims(i-1,size,P/k,dims,trydims,opt);
          }
      }
      else
      {
        // found a possible combination
        trydims[0] = P;

        // check for optimality
        double m = -1.0;

        for (int k=0; k<d; k++)
        {
          double mm=((double)size[k])/((double)trydims[k]);
          if (fmod((double)size[k],(double)trydims[k])>0.0001) mm*=3;
          if ( mm > m ) m = mm;
        }
        //if (_rank==0) std::cout << "optimize_dims: " << size << " | " << trydims << " norm=" << m << std::endl;
        if (m<opt)
        {
          opt = m;
          dims = trydims;
        }
      }
    }
  };

  /** \brief Implement yaspgrid load balance strategy for P=x^{dim} processors
   */
  template<int d>
  class YLoadBalancePowerD : public YLoadBalance<d>
  {
  public:
    typedef FieldVector<int, d>  iTupel;
    virtual void loadbalance (const iTupel& size, int P, iTupel& dims) const
    {
      bool found=false;
      for(int i=1; i<P; ++i)
        if(Power<d>::eval(i)==P) {
          for(int j=0; j<d; ++j)
            dims[j]=i;
          found=true;
        }
      if(!found)
        DUNE_THROW(GridError, "Loadbalancing failed\n");
    }
  };

  /*! Torus provides all the functionality to handle a toroidal communication structure:

     - Map a set of processes (given by an MPI communicator) to a torus of dimension d. The "optimal"
     torus dimensions are determined by a coarse mesh. The maximum side length is minimized.

     - Provide lists of neighboring processes and a method for nearest neighbor exchange
     using asynchronous communication with MPI. The periodic case is handled where one process
     might have to exchange several messages with the same process. (Logically, a process has always
     \f$3^d-1\f$ neighbors, but several of these logical neighbors might be identical)

     - Provide means to partition a grid to the torus.

   */
  template<int d>
  class Torus {
  public:
    //! type used to pass tupels in and out
    typedef FieldVector<int, d> iTupel;


  private:
    struct CommPartner {
      int rank;
      iTupel delta;
      int index;
    };

    struct CommTask {
      int rank;      // process to send to / receive from
      void *buffer;  // buffer to send / receive
      int size;      // size of buffer
#if HAVE_MPI
      MPI_Request request; // used by MPI to handle request
#else
      int request;
#endif
      int flag;      // used by MPI
    };

  public:
    //! constructor making uninitialized object
    Torus ()
    {}

    //! make partitioner from communicator and coarse mesh size
#if HAVE_MPI
    Torus (MPI_Comm comm, int tag, iTupel size, const YLoadBalance<d>* lb)
#else
    Torus (int tag, iTupel size, const YLoadBalance<d>* lb)
#endif
    {
      // MPI stuff
#if HAVE_MPI
      _comm = comm;
      MPI_Comm_size(comm,&_procs);
      MPI_Comm_rank(comm,&_rank);
#else
      _procs=1; _rank=0;
#endif
      _tag = tag;

      // determine dimensions
      lb->loadbalance(size, _procs, _dims);
      // if (_rank==0) std::cout << "Torus<" << d
      //                         << ">: mapping " << _procs << " processes onto "
      //                         << _dims << " torus." << std::endl;

      // compute increments for lexicographic ordering
      int inc = 1;
      for (int i=0; i<d; i++)
      {
        _increment[i] = inc;
        inc *= _dims[i];
      }

      // make full schedule
      proclists();
    }

    //! make partitioner from communicator and coarse mesh size
#if HAVE_MPI
    Torus (MPI_Comm comm, int tag, Dune::array<int,d> size, const YLoadBalance<d>* lb)
#else
    Torus (int tag, Dune::array<int,d> size, const YLoadBalance<d>* lb)
#endif
    {
      // MPI stuff
#if HAVE_MPI
      _comm = comm;
      MPI_Comm_size(comm,&_procs);
      MPI_Comm_rank(comm,&_rank);
#else
      _procs=1; _rank=0;
#endif
      _tag = tag;

      // determine dimensions
      iTupel sizeITupel;
      std::copy(size.begin(), size.end(), sizeITupel.begin());
      lb->loadbalance(sizeITupel, _procs, _dims);

      // compute increments for lexicographic ordering
      int inc = 1;
      for (int i=0; i<d; i++)
      {
        _increment[i] = inc;
        inc *= _dims[i];
      }

      // make full schedule
      proclists();
    }


    //! return own rank
    int rank () const
    {
      return _rank;
    }

    //! return own coordinates
    iTupel coord () const
    {
      return rank_to_coord(_rank);
    }

    //! return number of processes
    int procs () const
    {
      return _procs;
    }

    //! return dimensions of torus
    const iTupel & dims () const
    {
      return _dims;
    }

    //! return dimensions of torus in direction i
    int dims (int i) const
    {
      return _dims[i];
    }

    //! return MPI communicator
#if HAVE_MPI
    MPI_Comm comm () const
    {
      return _comm;
    }
#endif

    //! return tag used by torus
    int tag () const
    {
      return _tag;
    }

    //! return true if coordinate is inside torus
    bool inside (iTupel c) const
    {
      for (int i=d-1; i>=0; i--)
        if (c[i]<0 || c[i]>=_dims[i]) return false;
      return true;
    }

    //! map rank to coordinate in torus using lexicographic ordering
    iTupel rank_to_coord (int rank) const
    {
      iTupel coord;
      rank = rank%_procs;
      for (int i=d-1; i>=0; i--)
      {
        coord[i] = rank/_increment[i];
        rank = rank%_increment[i];
      }
      return coord;
    }

    //! map coordinate in torus to rank using lexicographic ordering
    int coord_to_rank (iTupel coord) const
    {
      for (int i=0; i<d; i++) coord[i] = coord[i]%_dims[i];
      int rank = 0;
      for (int i=0; i<d; i++) rank += coord[i]*_increment[i];
      return rank;
    }

    //! return rank of process where its coordinate in direction dir has offset cnt (handles periodic case)
    int rank_relative (int rank, int dir, int cnt) const
    {
      iTupel coord = rank_to_coord(rank);
      coord[dir] = (coord[dir]+_dims[dir]+cnt)%_dims[dir];
      return coord_to_rank(coord);
    }

    //! assign color to given coordinate
    int color (const iTupel & coord) const
    {
      int c = 0;
      int power = 1;

      // interior coloring
      for (int i=0; i<d; i++)
      {
        if (coord[i]%2==1) c += power;
        power *= 2;
      }

      // extra colors for boundary processes
      for (int i=0; i<d; i++)
      {
        if (_dims[i]>1 && coord[i]==_dims[i]-1) c += power;
        power *= 2;
      }

      return c;
    }

    //! assign color to given rank
    int color (int rank) const
    {
      return color(rank_to_coord(rank));
    }

    //! return the number of neighbors, which is \f$3^d-1\f$
    int neighbors () const
    {
      int n=1;
      for (int i=0; i<d; ++i)
        n *= 3;
      return n-1;
    }

    //! return true if neighbor with given delta is a neighbor under the given periodicity
    bool is_neighbor (iTupel delta, std::bitset<d> periodic) const
    {
      iTupel coord = rank_to_coord(_rank); // my own coordinate with 0 <= c_i < dims_i


      for (int i=0; i<d; i++)
      {
        if (delta[i]<0)
        {
          // if I am on the boundary and domain is not periodic => no neighbor
          if (coord[i]==0 && periodic[i]==false) return false;
        }
        if (delta[i]>0)
        {
          // if I am on the boundary and domain is not periodic => no neighbor
          if (coord[i]==_dims[i]-1 && periodic[i]==false) return false;
        }
      }
      return true;
    }

    //! partition the given grid onto the torus and return the piece of the process with given rank; returns load imbalance
    double partition (int rank, iTupel origin_in, Dune::array<int,d> size_in, iTupel& origin_out, iTupel& size_out) const
    {
      iTupel coord = rank_to_coord(rank);
      double maxsize = 1;
      double sz = 1;

      // make a tensor product partition
      for (int i=0; i<d; i++)
      {
        // determine
        int m = size_in[i]/_dims[i];
        int r = size_in[i]%_dims[i];

        sz *= size_in[i];

        if (coord[i]<_dims[i]-r)
        {
          origin_out[i] = origin_in[i] + coord[i]*m;
          size_out[i] = m;
          maxsize *= m;
        }
        else
        {
          origin_out[i] = origin_in[i] + (_dims[i]-r)*m + (coord[i]-(_dims[i]-r))*(m+1);
          size_out[i] = m+1;
          maxsize *= m+1;
        }
      }
      return maxsize/(sz/_procs);
    }

    /*!
       ProcListIterator provides access to a list of neighboring processes. There are always
       \f$ 3^d-1 \f$ entries in such a list. Two lists are maintained, one for sending and one for
       receiving. The lists are sorted in such a way that in sequence message delivery ensures that
       e.g. a message send to the left neighbor is received as a message from the right neighbor.
     */
    class ProcListIterator {
    public:
      //! make an iterator
      ProcListIterator (typename std::deque<CommPartner>::const_iterator iter)
      {
        i = iter;
      }

      //! return rank of neighboring process
      int rank () const
      {
        return i->rank;
      }

      //! return distance vector
      iTupel delta () const
      {
        return i->delta;
      }

      //! return index in proclist
      int index () const
      {
        return i->index;
      }

      //! return 1-norm of distance vector
      int distance () const
      {
        int dist = 0;
        iTupel delta=i->delta;
        for (int j=0; j<d; ++j)
          dist += std::abs(delta[j]);
        return dist;
      }

      //! Return true when two iterators point to same member
      bool operator== (const ProcListIterator& iter)
      {
        return i == iter.i;
      }


      //! Return true when two iterators do not point to same member
      bool operator!= (const ProcListIterator& iter)
      {
        return i != iter.i;
      }

      //! Increment iterator to next cell.
      ProcListIterator& operator++ ()
      {
        ++i;
        return *this;
      }

    private:
      typename std::deque<CommPartner>::const_iterator i;
    };

    //! first process in send list
    ProcListIterator sendbegin () const
    {
      return ProcListIterator(_sendlist.begin());
    }

    //! end of send list
    ProcListIterator sendend () const
    {
      return ProcListIterator(_sendlist.end());
    }

    //! first process in receive list
    ProcListIterator recvbegin () const
    {
      return ProcListIterator(_recvlist.begin());
    }

    //! last process in receive list
    ProcListIterator recvend () const
    {
      return ProcListIterator(_recvlist.end());
    }

    //! store a send request; buffers are sent in order; handles also local requests with memcpy
    void send (int rank, void* buffer, int size) const
    {
      CommTask task;
      task.rank = rank;
      task.buffer = buffer;
      task.size = size;
      if (rank!=_rank)
        _sendrequests.push_back(task);
      else
        _localsendrequests.push_back(task);
    }

    //! store a receive request; buffers are received in order; handles also local requests with memcpy
    void recv (int rank, void* buffer, int size) const
    {
      CommTask task;
      task.rank = rank;
      task.buffer = buffer;
      task.size = size;
      if (rank!=_rank)
        _recvrequests.push_back(task);
      else
        _localrecvrequests.push_back(task);
    }

    //! exchange messages stored in request buffers; clear request buffers afterwards
    void exchange () const
    {
      // handle local requests first
      if (_localsendrequests.size()!=_localrecvrequests.size())
      {
        std::cout << "[" << rank() << "]: ERROR: local sends/receives do not match in exchange!" << std::endl;
        return;
      }
      for (unsigned int i=0; i<_localsendrequests.size(); i++)
      {
        if (_localsendrequests[i].size!=_localrecvrequests[i].size)
        {
          std::cout << "[" << rank() << "]: ERROR: size in local sends/receive does not match in exchange!" << std::endl;
          return;
        }
        memcpy(_localrecvrequests[i].buffer,_localsendrequests[i].buffer,_localsendrequests[i].size);
      }
      _localsendrequests.clear();
      _localrecvrequests.clear();

#if HAVE_MPI
      // handle foreign requests
      int sends=0;
      int recvs=0;

      // issue sends to foreign processes
      for (unsigned int i=0; i<_sendrequests.size(); i++)
        if (_sendrequests[i].rank!=rank())
        {
          //          std::cout << "[" << rank() << "]" << " send " << _sendrequests[i].size << " bytes "
          //                    << "to " << _sendrequests[i].rank << " p=" << _sendrequests[i].buffer << std::endl;
          MPI_Isend(_sendrequests[i].buffer, _sendrequests[i].size, MPI_BYTE,
                    _sendrequests[i].rank, _tag, _comm, &(_sendrequests[i].request));
          _sendrequests[i].flag = false;
          sends++;
        }

      // issue receives from foreign processes
      for (unsigned int i=0; i<_recvrequests.size(); i++)
        if (_recvrequests[i].rank!=rank())
        {
          //          std::cout << "[" << rank() << "]"  << " recv " << _recvrequests[i].size << " bytes "
          //                    << "fm " << _recvrequests[i].rank << " p=" << _recvrequests[i].buffer << std::endl;
          MPI_Irecv(_recvrequests[i].buffer, _recvrequests[i].size, MPI_BYTE,
                    _recvrequests[i].rank, _tag, _comm, &(_recvrequests[i].request));
          _recvrequests[i].flag = false;
          recvs++;
        }

      // poll sends
      while (sends>0)
      {
        for (unsigned int i=0; i<_sendrequests.size(); i++)
          if (!_sendrequests[i].flag)
          {
            MPI_Status status;
            MPI_Test( &(_sendrequests[i].request), &(_sendrequests[i].flag), &status);
            if (_sendrequests[i].flag)
            {
              sends--;
              //                  std::cout << "[" << rank() << "]"  << " send to " << _sendrequests[i].rank << " OK" << std::endl;
            }
          }
      }

      // poll receives
      while (recvs>0)
      {
        for (unsigned int i=0; i<_recvrequests.size(); i++)
          if (!_recvrequests[i].flag)
          {
            MPI_Status status;
            MPI_Test( &(_recvrequests[i].request), &(_recvrequests[i].flag), &status);
            if (_recvrequests[i].flag)
            {
              recvs--;
              //                  std::cout << "[" << rank() << "]"  << " recv fm " << _recvrequests[i].rank << " OK" << std::endl;
            }

          }
      }

      // clear request buffers
      _sendrequests.clear();
      _recvrequests.clear();
#endif
    }

    //! global max
    double global_max (double x) const
    {
      double res = 0.0;

      if (_procs==1) return x;
#if HAVE_MPI
      MPI_Allreduce(&x,&res,1,MPI_DOUBLE,MPI_MAX,_comm);
#endif
      return res;
    }

    //! print contents of torus object
    void print (std::ostream& s) const
    {
      s << "[" << rank() <<  "]: Torus " << procs() << " processor(s) arranged as " << dims() << std::endl;
      for (ProcListIterator i=sendbegin(); i!=sendend(); ++i)
      {
        s << "[" << rank() <<  "]: send to   "
          << "rank=" << i.rank()
          << " index=" << i.index()
          << " delta=" << i.delta() << " dist=" << i.distance() << std::endl;
      }
      for (ProcListIterator i=recvbegin(); i!=recvend(); ++i)
      {
        s << "[" << rank() <<  "]: recv from "
          << "rank=" << i.rank()
          << " index=" << i.index()
          << " delta=" << i.delta() << " dist=" << i.distance() << std::endl;
      }
    }

  private:

    void proclists ()
    {
      // compile the full neighbor list
      CommPartner cp;
      iTupel delta;

      delta = -1;
      bool ready = false;
      iTupel me, nb;
      me = rank_to_coord(_rank);
      int index = 0;
      int last = neighbors()-1;
      while (!ready)
      {
        // find neighbors coordinates
        for (int i=0; i<d; i++)
          nb[i] = ( me[i]+_dims[i]+delta[i] ) % _dims[i];

        // find neighbors rank
        int nbrank = coord_to_rank(nb);

        // check if delta is not zero
        for (int i=0; i<d; i++)
          if (delta[i]!=0)
          {
            cp.rank = nbrank;
            cp.delta = delta;
            cp.index = index;
            _recvlist.push_back(cp);
            cp.index = last-index;
            _sendlist.push_front(cp);
            index++;
            break;
          }

        // next neighbor
        ready = true;
        for (int i=0; i<d; i++)
          if (delta[i]<1)
          {
            (delta[i])++;
            ready=false;
            break;
          }
          else
          {
            delta[i] = -1;
          }
      }

    }

#if HAVE_MPI
    MPI_Comm _comm;
#endif
    int _rank;
    int _procs;
    iTupel _dims;
    iTupel _increment;
    int _tag;
    std::deque<CommPartner> _sendlist;
    std::deque<CommPartner> _recvlist;

    mutable std::vector<CommTask> _sendrequests;
    mutable std::vector<CommTask> _recvrequests;
    mutable std::vector<CommTask> _localsendrequests;
    mutable std::vector<CommTask> _localrecvrequests;

  };

  //! Output operator for Torus
  template <int d>
  inline std::ostream& operator<< (std::ostream& s, const Torus<d> & t)
  {
    t.print(s);
    return s;
  }

} // namespace Dune

#endif
