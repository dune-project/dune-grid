// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_IO_FILE_VTK_FUNCTION_HH
#define DUNE_GRID_IO_FILE_VTK_FUNCTION_HH

#include <string>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/common.hh>

/** @file
    @author Peter Bastian, Christian Engwer
    @brief Functions for VTK output
 */

namespace Dune
{
  //! \addtogroup VTK
  //! \{

  //////////////////////////////////////////////////////////////////////
  //
  //  Base VTKFunction
  //

  /** \brief A base class for grid functions with any return type and dimension

      Trick : use double as return type
   */
  template< class GridView >
  class VTKFunction
  {
  public:
    typedef typename GridView::ctype ctype;
    constexpr static int dim = GridView::dimension;
    typedef typename GridView::template Codim< 0 >::Entity Entity;

    //! return number of components (1 for scalar valued functions, 3 for
    //! vector valued function in 3D etc.)
    virtual int ncomps () const = 0;

    //! evaluate single component comp in the entity e at local coordinates xi
    /*! Evaluate the function in an entity at local coordinates.
       @param[in]  comp   number of component to be evaluated
       @param[in]  e      reference to grid entity of codimension 0
       @param[in]  xi     point in local coordinates of the reference element
                         of e
       \return            value of the component
     */
    virtual double evaluate (int comp, const Entity& e,
                             const Dune::FieldVector<ctype,dim>& xi) const = 0;

    //! get name
    virtual std::string name () const = 0;

    //! get output precision for the field
    virtual VTK::Precision precision() const
    { return VTK::Precision::float32; }

    //! virtual destructor
    virtual ~VTKFunction () {}
  };

  //////////////////////////////////////////////////////////////////////
  //
  //  P0VTKFunction
  //

  //! Take a vector and interpret it as cell data for the VTKWriter
  /**
   * This class turns a generic vector containing cell data into a
   * VTKFunction.  The vector must allow read access to the data via
   * operator[]() and store the data in the order given by
   * MultipleCodimMultipleGeomTypeMapper with a layout class that allows only
   * elements.  Also, it must support the method size().
   *
   * While the number of components of the function is always 1, the vector
   * may represent a field with multiple components of which one may be
   * selected.
   *
   * \tparam GV Type of GridView the vector applies to.
   * \tparam V  Type of vector.
   */
  template<typename GV, typename V>
  class P0VTKFunction
    : public VTKFunction< GV >
  {
    //! Base class
    typedef VTKFunction< GV > Base;
    //! Mapper for elements
    typedef MultipleCodimMultipleGeomTypeMapper<GV> Mapper;

    //! store a reference to the vector
    const V& v;
    //! name of this function
    std::string s;
    //! number of components of the field stored in the vector
    int ncomps_;
    //! index of the component of the field in the vector this function is
    //! responsible for
    int mycomp_;
    //! precision with which to output the field
    VTK::Precision prec_;
    //! mapper used to map elements to indices
    Mapper mapper;

  public:
    typedef typename Base::Entity Entity;
    typedef typename Base::ctype ctype;
    using Base::dim;

    //! return number of components
    int ncomps () const override
    {
      return 1;
    }

    //! evaluate
    double evaluate (int, const Entity& e,
                     const Dune::FieldVector<ctype,dim>&) const override
    {
      return v[mapper.index(e)*ncomps_+mycomp_];
    }

    //! get name
    std::string name () const override
    {
      return s;
    }

    //! get output precision for the field
    VTK::Precision precision() const override
    {
      return prec_;
    }

    //! construct from a vector and a name
    /**
     * \param gv     GridView to operate on (used to instantiate a
     *               MultipleCodimMultipleGeomeTypeMapper, otherwise no
     *               reference or copy is stored).  Note that this must be the
     *               GridView the vector applies to as well as the GridView
     *               later used by the VTKWriter -- i.e. we do not implicitly
     *               restrict or prolongate the data.
     * \param v_     Reference to the vector holding the data.  The reference
     *               is stored internally and must be valid for as long as
     *               this functions evaluate method is used.
     * \param s_     Name of this function in the VTK file.
     * \param ncomps Number of components of the field represented by the
     *               vector.
     * \param mycomp Number of the field component this function is
     *               responsible for.
     * \param prec   the precision with which to output the field
     */
    P0VTKFunction(const GV &gv, const V &v_, const std::string &s_,
                  int ncomps=1, int mycomp=0, VTK::Precision prec = VTK::Precision::float32)
      : v( v_ ),
        s( s_ ),
        ncomps_(ncomps),
        mycomp_(mycomp),
        prec_(prec),
        mapper( gv, mcmgElementLayout() )
    {
      if (v.size()!=(unsigned int)(mapper.size()*ncomps_))
        DUNE_THROW(IOError, "P0VTKFunction: size mismatch");
    }

    //! destructor
    virtual ~P0VTKFunction() {}
  };

  //////////////////////////////////////////////////////////////////////
  //
  //  P1VTKFunction
  //

  //! Take a vector and interpret it as point data for the VTKWriter
  /**
   * This class turns a generic vector containing point data into a
   * VTKFunction.  The vector must allow read access to the data via
   * operator[]() and store the data in the order given by
   * MultipleCodimMultipleGeomTypeMapper with a layout class that allows only
   * vertices.  Also, it must support the method size().
   *
   * While the number of components of the function is always 1, the vector
   * may represent a field with multiple components of which one may be
   * selected.
   *
   * \tparam GV Type of GridView the vector applies to.
   * \tparam V  Type of vector.
   */
  template<typename GV, typename V>
  class P1VTKFunction
    : public VTKFunction< GV >
  {
    //! Base class
    typedef VTKFunction< GV > Base;
    //! Mapper for vertices
    typedef MultipleCodimMultipleGeomTypeMapper<GV> Mapper;

    //! store a reference to the vector
    const V& v;
    //! name of this function
    std::string s;
    //! number of components of the field stored in the vector
    int ncomps_;
    //! index of the component of the field in the vector this function is
    //! responsible for
    int mycomp_;
    //! precision with which to output the field
    VTK::Precision prec_;
    //! mapper used to map elements to indices
    Mapper mapper;

  public:
    typedef typename Base::Entity Entity;
    typedef typename Base::ctype ctype;
    using Base::dim;

    //! return number of components
    int ncomps () const override
    {
      return 1;
    }

    //! evaluate
    double evaluate ([[maybe_unused]] int comp, const Entity& e,
                     const Dune::FieldVector<ctype,dim>& xi) const override
    {
      const unsigned int myDim = Entity::mydimension;
      const unsigned int nVertices = e.subEntities(dim);

      std::vector<FieldVector<ctype,1> > cornerValues(nVertices);
      for (unsigned i=0; i<nVertices; ++i)
        cornerValues[i] = v[mapper.subIndex(e,i,myDim)*ncomps_+mycomp_];

      // (Ab)use the MultiLinearGeometry class to do multi-linear interpolation between scalars
      const MultiLinearGeometry<ctype,dim,1> interpolation(e.type(), cornerValues);
      return interpolation.global(xi);
    }

    //! get name
    std::string name () const override
    {
      return s;
    }

    //! get output precision for the field
    VTK::Precision precision() const override
    {
      return prec_;
    }

    //! construct from a vector and a name
    /**
     * \param gv     GridView to operate on (used to instantiate a
     *               MultipleCodimMultipleGeomTypeMapper, otherwise no
     *               reference or copy is stored).  Note that this must be the
     *               GridView the vector applies to as well as the GridView
     *               later used by the VTKWriter -- i.e. we do not implicitly
     *               restrict or prolongate the data.
     * \param v_     Reference to the vector holding the data.  The reference
     *               is stored internally and must be valid for as long as
     *               this functions evaluate method is used.
     * \param s_     Name of this function in the VTK file.
     * \param ncomps Number of components of the field represented by the
     *               vector.
     * \param mycomp Number of the field component this function is
     *               responsible for.
     * \param prec   the precision with which to output the field
     */
    P1VTKFunction(const GV& gv, const V &v_, const std::string &s_,
                  int ncomps=1, int mycomp=0, VTK::Precision prec = VTK::Precision::float32)
      : v( v_ ),
        s( s_ ),
        ncomps_(ncomps),
        mycomp_(mycomp),
        prec_(prec),
        mapper( gv, mcmgVertexLayout() )
    {
      if (v.size()!=(unsigned int)(mapper.size()*ncomps_))
        DUNE_THROW(IOError,"P1VTKFunction: size mismatch");
    }

    //! destructor
    virtual ~P1VTKFunction() {}
  };

  //! \} group VTK

} // namespace Dune

#endif // DUNE_GRID_IO_FILE_VTK_FUNCTION_HH
