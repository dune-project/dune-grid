// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id: cubeshapefunctions.hh 336 2006-05-03 13:09:05Z oliver $

#ifndef DUNE_UGGRID_CUBESHAPEFUNCTIONS_HH
#define DUNE_UGGRID_CUBESHAPEFUNCTIONS_HH

#include <iostream>
#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/misc.hh>
#include <dune/common/geometrytype.hh>
#include <dune/grid/common/referenceelements.hh>

/**
 * @file
 * @brief  define Lagrange type shape functions
 * @author Peter Bastian
 */
namespace Dune {

  namespace UGShapeFunctions {

    /***********************************************************
    * P2 shape functions for the cube of any dimension
    ***********************************************************/

    /*!
     * A class for piecewise quadratic shape functions in the cube.
     *
     * Shape function no. i is given by
     *
     *    phi_i = \prod_{j=0}^{dim-1} { a[j_i] + b[j_i]x_{j_i} + c[j_i]x_{j_i}^2 }
     *
     * where (j_{d-1},...,j_0) is the integer coordinate of the i'th node
     * with j_i \in {0,1,2}.
     */
    template<typename T, int d>
    class P2CubeShapeFunction
    {
    public:

      // compile time sizes
      enum { dim=d };      // maps from R^d
      enum { comps=1 };        // to R^1

      enum { m=Power_m_p<3,d>::power };   // total number of basis functions

      // export types
      typedef T CoordType;
      typedef T ResultType;
      typedef P2CubeShapeFunction ImplementationType;

      //! make a shape function object
      P2CubeShapeFunction (int no, int en, int co, const FieldVector<int,dim>& ipos)   // make it the i'th shape function
      {
        num = no;
        ent = en;
        cod = co;
        for (int j=0; j<dim; j++)
        {
          if (ipos[j]==0)
          {
            a[j] = 1.0;
            b[j] = -3.0;
            c[j] = 2.0;
            pos[j] = 0;
          }
          if (ipos[j]==1)
          {
            a[j] = 0.0;
            b[j] = 4.0;
            c[j] = -4.0;
            pos[j] = 0.5;
          }
          if (ipos[j]==2)
          {
            a[j] = 0.0;
            b[j] = -1.0;
            c[j] = 2.0;
            pos[j] = 1.0;
          }
        }
      }

      //! must be defaultconstructible
      P2CubeShapeFunction ()
      {       }

      //! evaluate shape function in local coordinates
      ResultType evaluateFunction (int comp, const FieldVector<CoordType,d>& x) const
      {
        ResultType phi = a[0]+x[0]*b[0] +x[0]*x[0]*c[0];
        for (int j=1; j<dim; j++) phi *= a[j]+x[j]*b[j]+x[j]*x[j]*c[j];
        return phi;
      }

      //! evaluate gradient in local coordinates
      ResultType evaluateDerivative (int comp, int dir, const FieldVector<CoordType,d>& x) const
      {
        ResultType deriv = b[dir]+2*c[dir]*x[dir];
        for (int j=0; j<dim; j++)
          if (j!=dir) deriv *= a[j]+x[j]*b[j]+x[j]*x[j]*c[j];
        return deriv;
      }

      //! consecutive number of associated dof within element
      int localindex (int comp) const
      {
        return num;
      }

      //! codim of associated dof
      int codim () const
      {
        return cod;
      }

      //! entity (of codim) of associated dof
      int entity () const
      {
        return ent;
      }

      //! consecutive number of dof within entity
      int entityindex () const
      {
        return 0;
      }

      //! interpolation point associated with shape function
      const FieldVector<CoordType,dim>& position () const
      {
        return pos;
      }

    private:
      int num;
      int ent;
      int cod;
      ResultType a[dim], b[dim], c[dim];   // store coefficients for this shape function
      FieldVector<CoordType,d> pos;
    };


    /*! P2CubeShapeFunctionSet implements the interface of
     * LagrangeShapeFunctionSet but is NOT derived from it.
     * S is either P1CubeShapeFunction<..> or LagrangeShapeFunctionWrapper<P1CubeShapeFunction<..> >
     */
    template<typename T, int d, typename S>
    class P2CubeShapeFunctionSet
    {
    public:

      // compile time sizes
      enum { dim=d };      // maps from R^d
      enum { comps=1 };        // to R^1

      enum { m=Power_m_p<3,d>::power };   // total number of basis functions

      // export types
      typedef T CoordType;
      typedef T ResultType;
      typedef S value_type;
      typedef typename S::ImplementationType Imp;   // Imp is either S or derived from S

      //! make a shape function object
      P2CubeShapeFunctionSet ()
      {
        ReferenceCube<T,d> cube;
        int i=0;
        for (int c=0; c<=dim; c++)
          for (int e=0; e<cube.size(c); e++)
          {
            sf[i] = Imp(i,e,c,cube.iposition(e,c));
            i++;
          }
      }

      //! return total number of shape functions
      int size () const
      {
        return m;
      }

      //! total number of shape functions associated with entity in codim
      int size (int entity, int codim) const
      {
        return 1;
      }

      //! random access to shape functions
      const value_type& operator[] (int i) const
      {
        return sf[i];   // ok derived class reference goes for base class reference
      }

      //! return order
      int order () const
      {
        return 2;
      }

      //! return type of element
      GeometryType type () const
      {
        static GeometryType cube(GeometryType::cube, dim);
        return cube;
      }

    private:
      S sf[m];
    };


    //! This are P1 shape functions in the cube without virtual functions
    template<typename T, int d>
    class P2CubeShapeFunctionSetContainer
    {
    public:
      // compile time sizes
      enum { dim=d };
      enum { comps=1 };
      enum { maxsize=Power_m_p<3,d>::power };

      // exported types
      typedef T CoordType;
      typedef T ResultType;
      typedef P2CubeShapeFunctionSet<T,d,P2CubeShapeFunction<T,d> > value_type;

      const value_type& operator() (GeometryType type, int order) const
      {
        if (type.isCube()) return p2cube;
        DUNE_THROW(NotImplemented, "type not implemented yet");
      }
    private:
      value_type p2cube;
    };

  }
}
#endif
