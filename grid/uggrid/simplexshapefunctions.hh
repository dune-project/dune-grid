// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_UGGRID_SIMPLEXSHAPEFUNCTIONS_HH
#define DUNE_UGGRID_SIMPLEXSHAPEFUNCTIONS_HH

#include <iostream>
#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/geometrytype.hh>

/**
 * @file
 * @brief  define Lagrange type shape functions for simplicial elements
 * @author Sreejith Pulloor Kuttanikkad
 */

namespace Dune {

  namespace UGShapeFunctions {

    /***********************************************************
    * P2 shape functions for the simplex
    ***********************************************************/
    /*!
     * A class for piecewise quadratic shape functions in a simplex geometry(lines/triangle/tetrehedron).
     * Evaluation is done at the quadrature points
     */

    /*  specialization for 1D
        polynomials of the form coeff[0]*x^2 + coeff[1]*x + coeff[2]
     */


    template<typename T, int dim>
    class P2SimplexShapeFunction;

    template<typename T>
    class P2SimplexShapeFunction<T,1>
      : public P2CubeShapeFunction<T,1>
    {
    public:
      enum {dim=1};
      enum {comps=1};
      enum {m=3};
      typedef T CoordType;
      typedef T ResultType;
      typedef P2SimplexShapeFunction ImplementationType;
      P2SimplexShapeFunction(int i,int en,int co)
        : P2CubeShapeFunction<T,1>(i,en,co,FieldVector<int,1>(i))
      {}

      P2SimplexShapeFunction()
      {}
    };

    //specialization 1D
    //-------------------------------------------------
    template<typename T, int d, typename S>
    class P2SimplexShapeFunctionSet;

    template<typename T,typename S>
    class P2SimplexShapeFunctionSet<T,1,S>

    {

    public:
      enum {dim=1};
      enum {comps=1};
      enum {m=3};

      typedef T CoordType;
      typedef T ResultType;
      typedef S value_type;
      typedef typename S::ImplementationType Imp;
      //! make a shape fn object
      P2SimplexShapeFunctionSet()
      {
        ReferenceSimplex<T,1> simpline;
        int j=0;
        for (int c=2; c>=0; --c)
          for(int e=0; e<simpline.size(c); ++e)
          {
            sf[j] = Imp(j,e,c);
            j++;
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
        return sf[i]; // ok derived class reference goes for base class reference
      }

      //! return order
      int order () const
      {
        return 2;
      }

      //! return type of element
      GeometryType type () const
      {
        static GeometryType simplex(GeometryType::simplex, dim);
        return simplex;
      }
    private:
      S sf[m];
    };

    //! These are P2 shape functions in the simplex without virtual functions
    // specialization for line
    template<typename T, int d>
    class P2SimplexShapeFunctionSetContainer;

    template<typename T>
    class P2SimplexShapeFunctionSetContainer<T,1>
    {

    public:

      enum {dim=1};
      enum {comps=1};
      enum {maxsize=3};

      // exported types
      typedef T CoordType;
      typedef T ResultType;
      typedef P2SimplexShapeFunctionSet<T,dim,P2SimplexShapeFunction<T,dim> > value_type;

      const value_type& operator() (GeometryType type, int order) const
      {
        if (type.isSimplex()) return p2simplex;
        DUNE_THROW(NotImplemented,"type not yet implemented");
      }
    private:
      value_type p2simplex;
    };

    /*  specialization for 2D
        A class for piecewise quadratic shape functions in a simplex geometry for triangles
        Evaluation is done at the quadrature points
     */
    /*
       polynomial of the form: coeff * ( a[0] + b[0]*x + b[1]*y )(a[1] + c[0]*x + c[1]*y )
       derrivative of the form:( aa[0][dir] + bb[0][dir]*x + bb[1][dir]*y)
       dir=0 for x, dir=1 for y
     */

    template<typename T, int dim>
    class P2SimplexShapeFunction;

    template<typename T>
    class P2SimplexShapeFunction<T,2>
    {

    public:

      enum {dim=2};
      enum {comps=1};
      enum {m=6}; // no of basis functions=6

      typedef T CoordType;
      typedef T ResultType;
      typedef P2SimplexShapeFunction ImplementationType;

      P2SimplexShapeFunction(int i,int en,int co)
      {
        number=i;
        ent = en;
        cod = co;

        switch(i)
        {
        case 0 :
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
          DUNE_THROW(RangeError, "wrong no of shape fns? check simplexshapefunctins.hh");
          break;
        }
      }

      P2SimplexShapeFunction()
      {}

      //! evaluate shape function in local coordinates
      ResultType evaluateFunction (int comp, const FieldVector<CoordType,2>& x) const
      {

        ResultType phi1=a[0];
        ResultType phi2=a[1];

        for (int j=0; j<2; ++j)
        {
          phi1 +=b[j]*x[j];
          phi2 +=c[j]*x[j];
        }
        ResultType phi=coeff*phi1*phi2;
        return phi;
      }


      ResultType evaluateDerivative (int comp, int dir, const FieldVector<CoordType,2>& x) const

      {

        ResultType deriv=aa[0][dir];

        for (int j=0; j<2; ++j)
        {
          deriv +=bb[j][dir]*x[j];
        }

        return deriv;
      }



      //! consecutive number of associated dof within element
      int localindex (int comp) const
      {
        return number;
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
      const FieldVector<CoordType,2>& position () const
      {
        return pos;
      }

    private:
      int number,coeff,ent,cod;
      ResultType a[2], b[2], c[2],aa[2][2], bb[2][2];
      FieldVector<CoordType,2> pos;

    };


    //specialization 2D
    //-------------------------------------------------
    template<typename T, int d, typename S>
    class P2SimplexShapeFunctionSet;

    template<typename T,typename S>
    class P2SimplexShapeFunctionSet<T,2,S>

    {

    public:
      enum {dim=2};
      enum {comps=1};
      enum {m=6};

      typedef T CoordType;
      typedef T ResultType;
      typedef S value_type;
      typedef typename S::ImplementationType Imp; // Imp is either S or derived from S
      //! make a shape fn object
      P2SimplexShapeFunctionSet()
      {
        ReferenceSimplex<T,2> triang;
        int j=0;
        for (int c=2; c>=1; --c)
          for(int e=0; e<triang.size(c); ++e)
          {
            sf[j] = Imp(j,e,c);
            j++;
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
        return sf[i]; // ok derived class reference goes for base class reference
      }

      //! return order
      int order () const
      {
        return 2;
      }

      //! return type of element
      GeometryType type () const
      {
        static GeometryType simplex(GeometryType::simplex, dim);
        return simplex;
      }
    private:
      S sf[m];
    };

    //! These are P2 shape functions in the simplex without virtual functions
    // specialization for triangle
    template<typename T, int d>
    class P2SimplexShapeFunctionSetContainer;

    template<typename T>
    class P2SimplexShapeFunctionSetContainer<T,2>
    {

    public:

      enum {dim=2};
      enum {comps=1};
      enum {maxsize=6};

      // exported types
      typedef T CoordType;
      typedef T ResultType;
      typedef P2SimplexShapeFunctionSet<T,dim,P2SimplexShapeFunction<T,dim> > value_type;

      const value_type& operator() (GeometryType type, int order) const
      {
        if (type.isSimplex()) return p2simplex;
        DUNE_THROW(NotImplemented,"type not yet implemented");
      }
    private:
      value_type p2simplex;
    };

    /*  specialization for 3D
        A class for piecewise quadratic shape functions in a simplex geometry for tetrahedron
        Evaluation is done at the quadrature points
     */
    template<typename T, int dim>
    class P2SimplexShapeFunction;

    template<typename T>
    class P2SimplexShapeFunction<T,3>
    {
    public:
      enum {dim=3};
      enum {comps=1};
      enum {m=10}; // no of basis functions=10

      typedef T CoordType;
      typedef T ResultType;
      typedef P2SimplexShapeFunction ImplementationType;


      P2SimplexShapeFunction(int i,int en,int co)
      {
        number=i;
        ent = en;
        cod = co;

        switch(i)
        {
        case 0 :
          //--interpolation point associated with shape fn
          pos[0]=0.0;
          pos[1]=0.0;
          pos[2]=0.0;
          //--
          coeff=2;
          a[0]=1.0;
          a[1]=0.5;
          b[0]=-1.0;
          b[1]=-1.0;
          b[2]=-1.0;
          c[0]=-1.0;
          c[1]=-1.0;
          c[2]=-1.0;
          //x derivative
          aa[0][0]=-3.0;
          bb[0][0]=4.0;
          bb[1][0]=4.0;
          bb[2][0]=4.0;
          //y derivative
          aa[0][1]=-3.0;
          bb[0][1]=4.0;
          bb[1][1]=4.0;
          bb[2][1]=4.0;
          // z derivative
          aa[0][2]=-3.0;
          bb[0][2]=4.0;
          bb[1][2]=4.0;
          bb[2][2]=4.0;
          break;
        case 1 :
          //--interpolation point associated with shape fn
          pos[0]=1.0;
          pos[1]=0.0;
          pos[2]=0.0;
          //--
          coeff=2;
          a[0]=0.0;
          a[1]=-0.5;
          b[0]=1.0;
          b[1]=0.0;
          b[2]=0.0;
          c[0]=1.0;
          c[1]=0.0;
          c[2]=0.0;
          //x derivative
          aa[0][0]=-1.0;
          bb[0][0]=4.0;
          bb[1][0]=0.0;
          bb[2][0]=0.0;
          //y derivative
          aa[0][1]=0.0;
          bb[0][1]=0.0;
          bb[1][1]=0.0;
          bb[2][1]=0.0;
          // z derivative
          aa[0][2]=0.0;
          bb[0][2]=0.0;
          bb[1][2]=0.0;
          bb[2][2]=0.0;
          break;
        case 2 :
          //--interpolation point associated with shape fn
          pos[0]=0.0;
          pos[1]=1.0;
          pos[2]=0.0;
          //--
          coeff=2;
          a[0]=0.0;
          a[1]=-0.5;
          b[0]=0.0;
          b[1]=1.0;
          b[2]=0.0;
          c[0]=0.0;
          c[1]=1.0;
          c[2]=0.0;
          //x derivative
          aa[0][0]=0.0;
          bb[0][0]=0.0;
          bb[1][0]=0.0;
          bb[2][0]=0.0;
          //y derivative
          aa[0][1]=-1.0;
          bb[0][1]=0.0;
          bb[1][1]=4.0;
          bb[2][1]=0.0;
          // z derivative
          aa[0][2]=0.0;
          bb[0][2]=0.0;
          bb[1][2]=0.0;
          bb[2][2]=0.0;
          break;
        case 3 :
          //--interpolation point associated with shape fn
          pos[0]=0.0;
          pos[1]=0.0;
          pos[2]=1.0;
          //--
          coeff=2;
          a[0]=0.0;
          a[1]=-0.5;
          b[0]=0.0;
          b[1]=0.0;
          b[2]=1.0;
          c[0]=0.0;
          c[1]=0.0;
          c[2]=1.0;
          //x derivative
          aa[0][0]=0.0;
          bb[0][0]=0.0;
          bb[1][0]=0.0;
          bb[2][0]=0.0;
          //y derivative
          aa[0][1]=0.0;
          bb[0][1]=0.0;
          bb[1][1]=0.0;
          bb[2][1]=0.0;
          // z derivative
          aa[0][2]=-1.0;
          bb[0][2]=0.0;
          bb[1][2]=0.0;
          bb[2][2]=4.0;
          break;
        case 4 :
          //--interpolation point associated with shape fn
          pos[0]=0.5;
          pos[1]=0.0;
          pos[2]=0.0;
          //--
          coeff=4;
          a[0]=0.0;
          a[1]=1.0;
          b[0]=1.0;
          b[1]=0.0;
          b[2]=0.0;
          c[0]=-1.0;
          c[1]=-1.0;
          c[2]=-1.0;
          //x derivative
          aa[0][0]=4.0;
          bb[0][0]=-8.0;
          bb[1][0]=-4.0;
          bb[2][0]=-4.0;
          //y derivative
          aa[0][1]=0.0;
          bb[0][1]=-4.0;
          bb[1][1]=0.0;
          bb[2][1]=0.0;
          // z derivative
          aa[0][2]=0.0;
          bb[0][2]=-4.0;
          bb[1][2]=0.0;
          bb[2][2]=0.0;
          break;
        case 5 :
          //--interpolation point associated with shape fn
          pos[0]=0.5;
          pos[1]=0.5;
          pos[2]=0.0;
          //--
          coeff=4;
          a[0]=0.0;
          a[1]=0.0;
          b[0]=1.0;
          b[1]=0.0;
          b[2]=0.0;
          c[0]=0.0;
          c[1]=1.0;
          c[2]=0.0;
          //x derivative
          aa[0][0]=0.0;
          bb[0][0]=0.0;
          bb[1][0]=4.0;
          bb[2][0]=0.0;
          //y derivative
          aa[0][1]=0.0;
          bb[0][1]=4.0;
          bb[1][1]=0.0;
          bb[2][1]=0.0;
          // z derivative
          aa[0][2]=0.0;
          bb[0][2]=0.0;
          bb[1][2]=0.0;
          bb[2][2]=0.0;
          break;
        case 6 :
          //--interpolation point associated with shape fn
          pos[0]=0.0;
          pos[1]=0.5;
          pos[2]=0.0;
          //--
          coeff=4;
          a[0]=0.0;
          a[1]=1.0;
          b[0]=0.0;
          b[1]=1.0;
          b[2]=0.0;
          c[0]=-1.0;
          c[1]=-1.0;
          c[2]=-1.0;
          //x derivative
          aa[0][0]=0.0;
          bb[0][0]=0.0;
          bb[1][0]=-4.0;
          bb[2][0]=0.0;
          //y derivative
          aa[0][1]=4.0;
          bb[0][1]=-4.0;
          bb[1][1]=-8.0;
          bb[2][1]=-4.0;
          // z derivative
          aa[0][2]=0.0;
          bb[0][2]=0.0;
          bb[1][2]=-4.0;
          bb[2][2]=0.0;
          break;
        case 7 :
          //--interpolation point associated with shape fn
          pos[0]=0.0;
          pos[1]=0.0;
          pos[2]=0.5;
          //--
          coeff=4;
          a[0]=0.0;
          a[1]=1.0;
          b[0]=0.0;
          b[1]=0.0;
          b[2]=1.0;
          c[0]=-1.0;
          c[1]=-1.0;
          c[2]=-1.0;
          //x derivative
          aa[0][0]=0.0;
          bb[0][0]=0.0;
          bb[1][0]=0.0;
          bb[2][0]=-4.0;
          //y derivative
          aa[0][1]=0.0;
          bb[0][1]=0.0;
          bb[1][1]=0.0;
          bb[2][1]=-4.0;
          // z derivative
          aa[0][2]=4.0;
          bb[0][2]=-4.0;
          bb[1][2]=-4.0;
          bb[2][2]=-8.0;
          break;
        case 8 :
          //--interpolation point associated with shape fn
          pos[0]=0.5;
          pos[1]=0.0;
          pos[2]=0.5;
          //--
          coeff=4;
          a[0]=0.0;
          a[1]=0.0;
          b[0]=1.0;
          b[1]=0.0;
          b[2]=0.0;
          c[0]=0.0;
          c[1]=0.0;
          c[2]=1.0;
          //x derivative
          aa[0][0]=0.0;
          bb[0][0]=0.0;
          bb[1][0]=0.0;
          bb[2][0]=4.0;
          //y derivative
          aa[0][1]=0.0;
          bb[0][1]=0.0;
          bb[1][1]=0.0;
          bb[2][1]=0.0;
          // z derivative
          aa[0][2]=0.0;
          bb[0][2]=4.0;
          bb[1][2]=0.0;
          bb[2][2]=0.0;
          break;
        case 9 :
          //--interpolation point associated with shape fn
          pos[0]=0.0;
          pos[1]=0.5;
          pos[2]=0.5;
          //--
          coeff=4;
          a[0]=0.0;
          a[1]=0.0;
          b[0]=0.0;
          b[1]=1.0;
          b[2]=0.0;
          c[0]=0.0;
          c[1]=0.0;
          c[2]=1.0;
          //x derivative
          aa[0][0]=0.0;
          bb[0][0]=0.0;
          bb[1][0]=0.0;
          bb[2][0]=0.0;
          //y derivative
          aa[0][1]=0.0;
          bb[0][1]=0.0;
          bb[1][1]=0.0;
          bb[2][1]=4.0;
          // z derivative
          aa[0][2]=0.0;
          bb[0][2]=0.0;
          bb[1][2]=4.0;
          bb[2][2]=0.0;
          break;
        default :
          DUNE_THROW(RangeError, "wrong no of shape fns? check out please");
          break;
        }
      }


      // default constructor
      P2SimplexShapeFunction()
      {}

      //! evaluate shape function in local coordinates
      ResultType evaluateFunction (int comp, const FieldVector<CoordType,3>& x) const
      {

        ResultType phi1=a[0];
        ResultType phi2=a[1];

        for (int j=0; j<3; ++j)
        {
          phi1 +=b[j]*x[j];
          phi2 +=c[j]*x[j];
        }
        ResultType phi=coeff*phi1*phi2;
        return phi;
      }


      ResultType evaluateDerivative (int comp, int dir, const FieldVector<CoordType,3>& x) const

      {

        ResultType deriv=aa[0][dir];

        for (int j=0; j<3; ++j)
        {
          deriv +=bb[j][dir]*x[j];
        }

        return deriv;

      }


      //! consecutive number of associated dof within element
      int localindex (int comp) const
      {
        return number;
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

      int number,coeff,ent,cod;
      ResultType a[2], b[3], c[3],aa[3][3], bb[3][3];
      FieldVector<CoordType,3> pos;


    };


    //specialization 3D
    //-----------------------
    template<typename T, int d, typename S>
    class P2SimplexShapeFunctionSet;

    template<typename T,typename S>
    class P2SimplexShapeFunctionSet<T,3,S>

    {

    public:
      enum {dim=3};
      enum {comps=1};
      enum {m=10};

      typedef T CoordType;
      typedef T ResultType;
      typedef S value_type;
      typedef typename S::ImplementationType Imp; // Imp is either S or derived from S
      //! make a shape fn object
      P2SimplexShapeFunctionSet()
      {
        ReferenceSimplex<T,3> tetrahed;
        int j=0;
        for (int c=3; c>=2; --c)
          for(int e=0; e<tetrahed.size(c); ++e)
          {
            sf[j] = Imp(j,e,c);
            j++;
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
        return sf[i]; // ok derived class reference goes for base class reference
      }

      //! return order
      int order () const
      {
        return 2;
      }

      //! return type of element
      GeometryType type () const
      {
        static GeometryType simplex(GeometryType::simplex, dim);
        return simplex;
      }
    private:
      S sf[m];
    };

  }
}
#endif
