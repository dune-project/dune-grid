// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_UGGRID_PRISMSHAPEFUNCTIONS_HH
#define DUNE_UGGRID_PRISMSHAPEFUNCTIONS_HH

#include <iostream>
#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/geometrytype.hh>

/**
 * @file
 * @brief  define Lagrange type shape functions for prism elements
 * @author Sreejith Pulloor Kuttanikkad
 */

namespace Dune {

  namespace UGShapeFunctions {

    /***********************************************************
    * P2 shape functions for prism
    ***********************************************************/

    /*!
     * A class for piecewise linear shape functions in a prism
     * Evaluation is done at the quadrature points of the prism
     */


    /*
          polynomial of the form: ( a[0] + b[0]*x + b[1]*y + b[2]z) * (a[1] + c[0]*x + c[1]*y + c[2]*z)
          derrivative of the form:( aa[0][dir] + bb[0][dir]*x + bb[1][dir]*y + bb[2][dir]*z)
          dir=0 for x, dir=1 for y, dir=2 for z
     */
    template<typename T>
    class P2PrismShapeFunction
    {
    public:

      enum {dim=3};
      enum {comps=1};
      enum {m=18};   // number of basis functions

      typedef T CoordType;
      typedef T ResultType;
      typedef P2PrismShapeFunction ImplementationType;

      // i'th shape function
      P2PrismShapeFunction(int i)
      {
        number=i;

        // A renumbering of the shape functions might simplify this constructor

        if (i<3 || (i>=12 && i<15))
          codim_=3;
        else if (i>=9 && i<12)
          codim_=1;
        else
          codim_=2;

        static const int entities[] = {0,1,2,        // vertices of the lower triangle
                                       1,2,0,        // edges of the lower triangle
                                       3, 4, 5,      // vertical edges
                                       2,3,1,        // quadrilateral sides
                                       3,4,5,        // vertices of the upper triangle
                                       7,8,6};       // edges of the upper triangle

        entity_ = entities[i];

        if (i<6) {

          a1d = 1.0;
          b1d = -3.0;
          c1d = 2.0;
          pos[2] = 0;

        } else if (i<12) {

          a1d = 0.0;
          b1d = 4.0;
          c1d = -4.0;
          pos[2] = 0.5;

        } else {
          a1d = 0.0;
          b1d = -1.0;
          c1d = 2.0;
          pos[2] = 1.0;
        }

        switch(i) {

        case 0 :
        case 6 :
        case 12 :

          //--interpolation point associated with shape fn
          pos[0]=0.0;
          pos[1]=0.0;
          //--
          coeff=2;
          a[0]=1.0;
          a[1]=0.5;
          b[0]=-1.0;
          b[1]=-1.0;
          c[0]=-1.0;
          c[1]=-1.0;
          //x derivative
          aa[0][0]=-3;
          bb[0][0]=4;
          bb[1][0]=4;
          //y derivative
          aa[0][1]=-3;
          bb[0][1]=4;
          bb[1][1]=4;
          break;
        case 1 :
        case 7 :
        case 13 :

          //--interpolation point associated with shape fn
          pos[0]=1.0;
          pos[1]=0.0;
          //--
          coeff=2;
          a[0]=0.0;
          a[1]=-0.5;
          b[0]=1.0;
          b[1]=0.0;
          c[0]=1.0;
          c[1]=0.0;
          //x derivative
          aa[0][0]=-1;
          bb[0][0]=4;
          bb[1][0]=0;
          //y derivative
          aa[0][1]=0;
          bb[0][1]=0;
          bb[1][1]=0;
          break;
        case 2 :
        case 8 :
        case 14 :
          //--interpolation point associated with shape fn
          pos[0]=0.0;
          pos[1]=1.0;
          //--
          coeff=2;
          a[0]=0.0;
          a[1]=-0.5;
          b[0]=0.0;
          b[1]=1.0;
          c[0]=0.0;
          c[1]=1.0;
          //x derivative
          aa[0][0]=0;
          bb[0][0]=0;
          bb[1][0]=0;
          //y derivative
          aa[0][1]=-1;
          bb[0][1]=0;
          bb[1][1]=4;
          break;
        case 3 :
        case 9 :
        case 15 :
          //--interpolation point associated with shape fn
          pos[0]=0.5;
          pos[1]=0.5;
          //--
          coeff=4;
          a[0]=0.0;
          a[1]=0.0;
          b[0]=1.0;
          b[1]=0.0;
          c[0]=0.0;
          c[1]=1.0;
          //x derivative
          aa[0][0]=0;
          bb[0][0]=0;
          bb[1][0]=4;
          //y derivative
          aa[0][1]=0;
          bb[0][1]=4;
          bb[1][1]=0;
          break;
        case 4 :
        case 10 :
        case 16 :
          //--interpolation point associated with shape fn
          pos[0]=0.0;
          pos[1]=0.5;
          //--
          coeff=4;
          a[0]=0.0;
          a[1]=1.0;
          b[0]=0.0;
          b[1]=1.0;
          c[0]=-1.0;
          c[1]=-1.0;
          //x derivative
          aa[0][0]=0;
          bb[0][0]=0;
          bb[1][0]=-4;
          //y derivative
          aa[0][1]=4;
          bb[0][1]=-4;
          bb[1][1]=-8;
          break;
        case 5 :
        case 11 :
        case 17 :
          //--interpolation point associated with shape fn
          pos[0]=0.5;
          pos[1]=0.0;
          //--
          coeff=4;
          a[0]=0.0;
          a[1]=1.0;
          b[0]=1.0;
          b[1]=0.0;
          c[0]=-1.0;
          c[1]=-1.0;
          //x derivative
          aa[0][0]=4;
          bb[0][0]=-8;
          bb[1][0]=-4;
          //y derivative
          aa[0][1]=0;
          bb[0][1]=-4;
          bb[1][1]=0;
          break;
        default :
          DUNE_THROW(RangeError, "wrong no of shape fns in Prism?");
          break;
        }
      }
      //! must be defaultconstructible
      P2PrismShapeFunction ()
      {}

      //! evaluate shape function in local coordinates
      ResultType evaluateFunction (int comp, const FieldVector<CoordType,3>& x) const
      {
        // Part I: Evaluate a P2 function on the triangle
        ResultType phi1=a[0];
        ResultType phi2=a[1];

        for (int j=0; j<2; ++j)
        {
          phi1 +=b[j]*x[j];
          phi2 +=c[j]*x[j];
        }

        // Part II: Evaluate a P2 function on the segment, multiply the results
        ResultType phi3 = a1d + x[2]*b1d + x[2]*x[2]*c1d;

        return coeff * phi1*phi2*phi3;


      }


      ResultType evaluateDerivative (int comp, int dir, const FieldVector<CoordType,3>& x) const

      {
        if (dir==2) {
          ResultType phi1=a[0];
          ResultType phi2=a[1];

          for (int j=0; j<2; ++j) {
            phi1 +=b[j]*x[j];
            phi2 +=c[j]*x[j];
          }

          ResultType phi = coeff * phi1 * phi2;
          return phi * (b1d + 2*c1d*x[2]);

        }

        ResultType deriv=aa[0][dir];

        for (int j=0; j<2; ++j)
          deriv +=bb[j][dir]*x[j];

        return deriv * (a1d + x[2]*b1d + x[2]*x[2]*c1d);

      }



      //! consecutive number of associated dof within element
      int localindex (int comp) const
      {
        return number;
      }

      //! codim of associated dof
      int codim () const
      {
        return codim_;
      }

      //! entity (of codim) of associated dof
      int entity () const
      {
        return entity_;
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
      int number,coeff,entity_,codim_;
      ResultType a[2], b[2], c[2],aa[2][2], bb[2][2];

      // Coefficients for the 1d-P2 shape functions
      ResultType a1d, b1d, c1d;
      FieldVector<CoordType,3> pos;


    };


    template<typename T, typename S>
    class P2PrismShapeFunctionSet
    {
    public:

      // compile time sizes
      enum { dim=3 };      // maps from R^d
      enum { comps=1 };        // to R^1

      enum { m=18 };   // total number of basis functions

      // export types
      typedef T CoordType;
      typedef T ResultType;
      typedef S value_type;
      typedef typename S::ImplementationType Imp;   // Imp is either S or derived from S

      //! make a shape function object
      P2PrismShapeFunctionSet ()
      {
        for (int i=0; i<m; i++)
          sf[i] = Imp(i);       // assignment of derived class objects defined in wrapper

      }

      //! return total number of shape functions
      int size () const
      {
        return m;
      }

      //! total number of shape functions associated with entity in codim
      int size (int entity, int codim) const
      {
        if (codim==dim)
          return 6;
        else if (codim==2)
          return 9;
        else if (codim==1)
          return 3;

        return 0;
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
        static GeometryType prism(GeometryType::prism, dim);
        return prism;
      }

    private:
      S sf[m];
    };


    //! P2 shape functions in the prism without virtual functions
    template<typename T>
    class P2PrismShapeFunctionSetContainer
    {
    public:
      // compile time sizes
      enum { dim=3 };
      enum { comps=1 };
      enum { maxsize=18 };

      // exported types
      typedef T CoordType;
      typedef T ResultType;
      typedef P2PrismShapeFunctionSet<T,P2PrismShapeFunction<T> > value_type;

      const value_type& operator() (GeometryType type, int order) const
      {
        if(type.isPrism()) return p2prism;
        DUNE_THROW(NotImplemented, "type not implemented yet");
      }

    private:
      value_type p2prism;
    };

  }
}
#endif
