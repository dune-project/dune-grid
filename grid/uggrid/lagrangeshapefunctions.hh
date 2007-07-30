// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id: lagrangeshapefunctions.hh 411 2006-11-20 16:46:35Z oliver $

#ifndef DUNE_UGGRID_LAGRANGESHAPEFUNCTIONS_HH
#define DUNE_UGGRID_LAGRANGESHAPEFUNCTIONS_HH

#include <iostream>
#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/misc.hh>
#include <dune/common/geometrytype.hh>

#include "shapefunctions.hh"
#include "cubeshapefunctions.hh"
#include "prismshapefunctions.hh"
#include "pyramidshapefunctions.hh"
#include "simplexshapefunctions.hh"

/**
 * @file
 * @brief  define Lagrange type shape functions
 * @author Peter Bastian
 */
namespace Dune {

  namespace UGShapeFunctions {

    /***********************************************************
    * The interface for Lagrange shape functions of arbitrary
    * order and element type.
    ***********************************************************/


    /** \brief A scalar (N=1) ShapeFunction extended by a method providing a position
        \param C Type used for coordinates in the reference element
        \param T Type used for the shape function values
        \param d Dimension of the element
     */
    template<typename T, int d>
    class LagrangeShapeFunction : public ShapeFunction<T,d,1>
    {
    public:
      // compile time sizes
      enum { dim=d };
      enum { comps=1 };
      typedef T CoordType;
      typedef T ResultType;

      //! interpolation point associated with shape function
      virtual const FieldVector<CoordType,dim>& position () const = 0;
    };



    /** \brief A scalar (N=1) ShapeFunctionSet that returns a LagrangeShapeFunction
       \param T Type used for the shape function values and coordinates
       \param d Dimension of the element
     */
    template<typename T, int d>
    class LagrangeShapeFunctionSet : public ShapeFunctionSet<T,d,1>
    {
    public:
      // compile time sizes
      enum { dim=d };
      enum { comps=1 };

      // exported types
      typedef T CoordType;
      typedef T ResultType;
      typedef LagrangeShapeFunction<ResultType,dim> value_type;

      //! random access to i'th LagrangeShapeFunction
      virtual const value_type& operator[] (int i) const = 0;
    };


    /***********************************************************
    * Wrappers
    *
    ***********************************************************/

    /*! wrap inlinable implementation into class that is derived
     * from abstract base class and implement the functions with
     * the given implementation
     */
    template<typename Imp>
    class LagrangeShapeFunctionWrapper :
      public LagrangeShapeFunction<typename Imp::ResultType,Imp::dim>,
      private Imp
    {
    public:

      // compile time sizes
      enum { dim=Imp::dim };
      enum { comps=1 };

      // exported types
      typedef typename Imp::CoordType CoordType;
      typedef typename Imp::ResultType ResultType;
      typedef Imp ImplementationType;

      //! assignment from implementation type (this class has no data)
      LagrangeShapeFunctionWrapper& operator= (const Imp& imp)
      {
        Imp::operator=(imp);
        return *this;
      }

      //! evaluate component comp at point x
      virtual ResultType evaluateFunction (int comp, const FieldVector<CoordType,dim>& x) const
      {
        return Imp::evaluateFunction(comp,x);
      }

      //! evaluate derivative of component comp in direction dir at point x
      virtual ResultType evaluateDerivative (int comp, int dir, const FieldVector<CoordType,dim>& x) const
      {
        return Imp::evaluateDerivative(comp,dir,x);
      }

      //! consecutive number of associated dof within element
      virtual int localindex (int comp) const
      {
        return Imp::localindex(comp);
      }

      //! codim of associated dof
      virtual int codim () const
      {
        return Imp::codim();
      }

      //! entity (of codim) of associated dof
      virtual int entity () const
      {
        return Imp::entity();
      }

      //! consecutive number of dof within entity
      virtual int entityindex () const
      {
        return Imp::entityindex();
      }

      //! interpolation point associated with shape function
      virtual const FieldVector<CoordType,dim>& position () const
      {
        return Imp::position();
      }
    };




    /*! wrap inlinable implementation into class that is derived
     * from abstract base class and implement the functions with
     * the given implementation
     */
    template<typename Imp>
    class LagrangeShapeFunctionSetWrapper :
      public LagrangeShapeFunctionSet<typename Imp::ResultType,Imp::dim>,
      private Imp
    {
    public:

      // compile time sizes
      enum { dim=Imp::dim };
      enum { comps=1 };       // must be available at compile time

      // exported types
      typedef typename Imp::CoordType CoordType;
      typedef typename Imp::ResultType ResultType;
      typedef Imp ImplementationType;
      typedef LagrangeShapeFunction<ResultType,dim> value_type;   // Note: Imp::value_type references
      // must be convertible to references to
      // this type !
      //! return total number of shape functions
      virtual int size () const
      {
        return Imp::size();
      }

      //! total number of shape functions associated with entity in codim
      virtual int size (int entity, int codim) const
      {
        return Imp::size(entity,codim);
      }

      //! random access to i'th ShapeFunction
      virtual const value_type& operator[] (int i) const
      {
        return Imp::operator[](i);
      }

      //! return order
      virtual int order () const
      {
        return Imp::order();
      }

      //! return type of element
      virtual GeometryType type () const
      {
        return Imp::type();
      }
    };



    /***********************************************************
    * The general container for Lagrange shape functions of
    * any order and element type (if someone implements them
    * in finite time). All containers are accessible in a singleton.
    ***********************************************************/

    /** \brief This are Lagrange shape functions for any element type and order
        \param C Type used for coordinates in the reference element
        \param T Type used for the shape function values
        \param d Dimension of the element
     */
    template<typename T, int d>
    class LagrangeShapeFunctionSetContainer : public ShapeFunctionSetContainer<T,d,1,Power_m_p<3,d>::power >
    {
    public:
      // compile time sizes
      enum { dim=d };
      enum { comps=1 };
      enum { maxsize=Power_m_p<3,dim>::power };

      // exported types
      typedef T CoordType;
      typedef T ResultType;

      //! type of objects in the container
      typedef LagrangeShapeFunctionSet<T,d> value_type;

      const value_type& operator() (GeometryType type, int order) const
      {
        if ( type.isCube() )
          return wrappedp2cube;

        if ( type.isSimplex() )
          return wrappedp2simplex;

        if ( type.isPyramid() )
          DUNE_THROW(RangeError, "No pyramid for this dimension");

        if ( type.isPrism() )
          DUNE_THROW(RangeError, "No prism for this dimension ");

        DUNE_THROW(RangeError, "type not available");
      }

    private:
      // the cubes
      typedef LagrangeShapeFunctionWrapper<P2CubeShapeFunction<T,d> > WrappedP2CubeShapeFunction;
      typedef P2CubeShapeFunctionSet<T,d,WrappedP2CubeShapeFunction> P2CubeWrappedShapeFunctionSet;
      typedef LagrangeShapeFunctionSetWrapper<P2CubeWrappedShapeFunctionSet> WrappedP2CubeShapeFunctionSet;
      WrappedP2CubeShapeFunctionSet wrappedp2cube;

      // the simplices
      typedef LagrangeShapeFunctionWrapper<P2SimplexShapeFunction<T,d> > WrappedP2SimplexShapeFunction;
      typedef P2SimplexShapeFunctionSet<T,d,WrappedP2SimplexShapeFunction> P2SimplexWrappedShapeFunctionSet;
      typedef LagrangeShapeFunctionSetWrapper<P2SimplexWrappedShapeFunctionSet> WrappedP2SimplexShapeFunctionSet;
      WrappedP2SimplexShapeFunctionSet wrappedp2simplex;
    };

    //! This are Lagrange shape functions for any element type and order (in the future ... )
    template<typename T>
    class LagrangeShapeFunctionSetContainer<T,3> : public ShapeFunctionSetContainer<T,3,1,Power_m_p<3,3>::power >
    {
    public:
      // compile time sizes
      enum { dim=3 };
      enum { comps=1 };
      enum { maxsize=Power_m_p<3,dim>::power };

      // exported types
      typedef T CoordType;
      typedef T ResultType;

      //! type of objects in the container
      typedef LagrangeShapeFunctionSet<T,3> value_type;

      const value_type& operator() (GeometryType type, int order) const
      {
        if ( type.isCube() )
          return wrappedp2cube;

        if ( type.isSimplex() )
          return wrappedp2simplex;

        if (type.isPyramid() )
        {
          //return wrappedp2pyramid;
          DUNE_THROW(RangeError, "order not available for pyramid");
        }

        if (type.isPrism())
          return wrappedp2prism;

        DUNE_THROW(RangeError, "type or order not available");
      }

    private:
      // the cubes
      typedef LagrangeShapeFunctionWrapper<P2CubeShapeFunction<T,dim> > WrappedP2CubeShapeFunction;
      typedef P2CubeShapeFunctionSet<T,dim,WrappedP2CubeShapeFunction> P2CubeWrappedShapeFunctionSet;
      typedef LagrangeShapeFunctionSetWrapper<P2CubeWrappedShapeFunctionSet> WrappedP2CubeShapeFunctionSet;
      WrappedP2CubeShapeFunctionSet wrappedp2cube;

      // the simplices
      typedef LagrangeShapeFunctionWrapper<P2SimplexShapeFunction<T,dim> > WrappedP2SimplexShapeFunction;
      typedef P2SimplexShapeFunctionSet<T,dim,WrappedP2SimplexShapeFunction> P2SimplexWrappedShapeFunctionSet;
      typedef LagrangeShapeFunctionSetWrapper<P2SimplexWrappedShapeFunctionSet> WrappedP2SimplexShapeFunctionSet;
      WrappedP2SimplexShapeFunctionSet wrappedp2simplex;

      // Pyramid
      //     typedef LagrangeShapeFunctionWrapper<P2PyramidShapeFunction<T> > WrappedP2PyramidShapeFunction;
      //     typedef P2PyramidShapeFunctionSet<C,T,WrappedP2PyramidShapeFunction> P2PyramidWrappedShapeFunctionSet;
      //     typedef LagrangeShapeFunctionSetWrapper<P2PyramidWrappedShapeFunctionSet> WrappedP2PyramidShapeFunctionSet;
      //     WrappedP2PyramidShapeFunctionSet wrappedp2pyramid;

      // Prism
      typedef LagrangeShapeFunctionWrapper<P2PrismShapeFunction<T> > WrappedP2PrismShapeFunction;
      typedef P2PrismShapeFunctionSet<T,WrappedP2PrismShapeFunction> P2PrismWrappedShapeFunctionSet;
      typedef LagrangeShapeFunctionSetWrapper<P2PrismWrappedShapeFunctionSet> WrappedP2PrismShapeFunctionSet;
      WrappedP2PrismShapeFunctionSet wrappedp2prism;
    };


    // singleton holding the reference element container
    template<typename T, int d>
    struct LagrangeShapeFunctions {
      static LagrangeShapeFunctionSetContainer<T,d> general;
    };

    /** @} */

  }
}
#endif
