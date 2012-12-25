// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_YASPGRIDGEOMETRY_HH
#define DUNE_GRID_YASPGRIDGEOMETRY_HH

/** \file
 * \brief The YaspGeometry class and its specializations

   YaspGeometry realizes the concept of the geometric part of a mesh entity.

   We have specializations for dim == dimworld (elements) and dim == 0
   (vertices).  The general version implements dim == dimworld-1 (faces)
   and otherwise throws a GridError.
 */

namespace Dune {

  //! The general version can do any dimension, but constructors currently exist only for dim==dimworld-1
  template<int mydim,int cdim, class GridImp>
  class YaspGeometry : public AxisAlignedCubeGeometry<typename GridImp::ctype,mydim,cdim>
  {
  public:
    //! define type used for coordinates in grid module
    typedef typename GridImp::ctype ctype;

    //! default constructor
    YaspGeometry ()
      : AxisAlignedCubeGeometry<ctype,mydim,cdim>(FieldVector<ctype,cdim>(0),FieldVector<ctype,cdim>(0)) // anything
    {}

    //! constructor from midpoint and extension and missing direction number
    YaspGeometry (const FieldVector<ctype, cdim>& p, const FieldVector<ctype, cdim>& h, uint8_t& m)
      : AxisAlignedCubeGeometry<ctype,mydim,cdim>(FieldVector<ctype,cdim>(0),FieldVector<ctype,cdim>(0)) // anything
    {
      if (cdim!=mydim+1)
        DUNE_THROW(GridError, "This YaspGeometry constructor assumes cdim=mydim+1");

      FieldVector<ctype, cdim> lower = p;
      FieldVector<ctype, cdim> upper = p;
      lower.axpy(-0.5,h);
      upper.axpy( 0.5,h);

      lower[m] = upper[m] = p[m];

      std::bitset<cdim> axes((1<<cdim)-1);    // all bits set
      axes[m] = false; // except the one at 'missing'

      // set up base class
      static_cast< AxisAlignedCubeGeometry<ctype,mydim,cdim> & >( *this ) = AxisAlignedCubeGeometry<ctype,mydim,cdim>(lower, upper, axes);
    }

    //! copy constructor
    YaspGeometry (const YaspGeometry& other)
      : AxisAlignedCubeGeometry<ctype,mydim,cdim>(other)
    {}

    //! print function
    void print (std::ostream& s) const
    {
      s << "YaspGeometry<"<<mydim<<","<<cdim<< "> ";
      s << "midpoint";
      for (int i=0; i<cdim; i++)
        s << " " << 0.5*(this->lower_[i] + this->upper_[i]);
      s << " extension";
      for (int i=0; i<cdim; i++)
        s << " " << (this->upper_[i] - this->lower_[i]);
      s << " coordinates: " << this->axes_;
    }

  };



  //! specialize for dim=dimworld, i.e. a volume element
  template<int mydim, class GridImp>
  class YaspGeometry<mydim,mydim,GridImp> : public AxisAlignedCubeGeometry<typename GridImp::ctype,mydim,mydim>
  {
  public:
    typedef typename GridImp::ctype ctype;

    //! default constructor
    YaspGeometry ()
      : AxisAlignedCubeGeometry<ctype,mydim,mydim>(FieldVector<ctype,mydim>(0),FieldVector<ctype,mydim>(0)) // anything
    {}

    //! constructor from midpoint and extension
    YaspGeometry (const FieldVector<ctype, mydim>& p, const FieldVector<ctype, mydim>& h)
      : AxisAlignedCubeGeometry<ctype,mydim,mydim>(FieldVector<ctype,mydim>(0),FieldVector<ctype,mydim>(0)) // anything
    {
      FieldVector<ctype, mydim> lower = p;
      FieldVector<ctype, mydim> upper = p;
      lower.axpy(-0.5,h);
      upper.axpy( 0.5,h);
      // set up base class
      static_cast< AxisAlignedCubeGeometry<ctype,mydim,mydim> & >( *this ) = AxisAlignedCubeGeometry<ctype,mydim,mydim>(lower, upper);
    }

    //! copy constructor (skipping temporary variables)
    YaspGeometry (const YaspGeometry& other)
      : AxisAlignedCubeGeometry<ctype,mydim,mydim>(other)
    {}

    //! print function
    void print (std::ostream& s) const
    {
      s << "YaspGeometry<"<<mydim<<","<<mydim<< "> ";
      s << "midpoint";
      for (int i=0; i<mydim; i++)
        s << " " << 0.5 * (this->lower_[i] + this->upper_[i]);
      s << " extension";
      for (int i=0; i<mydim; i++)
        s << " " << (this->upper_[i] + this->lower_[i]);
    }

  };

  //! specialization for dim=0, this is a vertex
  template<int cdim, class GridImp>
  class YaspGeometry<0,cdim,GridImp> : public AxisAlignedCubeGeometry<typename GridImp::ctype,0,cdim>
  {
  public:
    typedef typename GridImp::ctype ctype;

    //! default constructor
    YaspGeometry ()
      : AxisAlignedCubeGeometry<typename GridImp::ctype,0,cdim>(FieldVector<ctype,cdim>(0)) // anything
    {}

    //! constructor
    explicit YaspGeometry ( const FieldVector< ctype, cdim > &p )
      : AxisAlignedCubeGeometry<typename GridImp::ctype,0,cdim>( p )
    {}

    YaspGeometry ( const FieldVector< ctype, cdim > &p, const FieldVector< ctype, cdim > &, uint8_t &)
      : AxisAlignedCubeGeometry<typename GridImp::ctype,0,cdim>( p )
    {}

    //! print function
    void print (std::ostream& s) const
    {
      s << "YaspGeometry<"<<0<<","<<cdim<< "> ";
      s << "position " << this->lower_;
    }
  };

  // operator<< for all YaspGeometrys
  template <int mydim, int cdim, class GridImp>
  inline
  std::ostream& operator<< (std::ostream& s, YaspGeometry<mydim,cdim,GridImp>& e)
  {
    e.print(s);
    return s;
  }

}  // namespace Dune

#endif // DUNE_GRID_YASPGRIDGEOMETRY_HH
