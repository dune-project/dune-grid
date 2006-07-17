// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_QUADRATURERULES_HH
#define DUNE_QUADRATURERULES_HH

#include <iostream>
#include <vector>
#include <map>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/stdstreams.hh>
#include <dune/common/geometrytype.hh>

/**
   \file
   Interface for quadrature points and rules
 */

namespace Dune {

  /**
     Exception thrown if an Aquired QuadratureRule is not available,
     because the order is to high
   */
  class QuadratureOrderOutOfRange : public NotImplemented {};

  /**
     quadrature rules for cubes of any dimension based on Gauss quadrature
   */
  template<typename ct, int dim>
  class QuadraturePoint {
  public:
    // compile time parameters
    enum { d=dim };
    typedef ct CoordType;

    //! set up quadrature of given order in d dimensions
    QuadraturePoint (const FieldVector<ct, dim>& x, double w) : local(x)
    {
      wght = w;
    }

    //! return local coordinates of integration point i
    const FieldVector<ct, dim>& position () const
    {
      return local;
    }

    //! return weight associated with integration point i
    double weight () const
    {
      return wght;
    }
    virtual ~QuadraturePoint(){}

  protected:
    FieldVector<ct, dim> local;
    double wght;
  };

  /**
     Defines an \p enum for currently available quadrature rules.
   */
  namespace QuadratureType {
    enum Enum {
      Gauss      = 0,

      Jacobian_1_0 = 1,
      Jacobian_2_0 = 2,

      Simpson    = 3,
      Trap       = 4,
      Grid       = 5,

      Clough     = 21,

      Invalid_Rule = 127
    };
  };

  /**
     abstract base class for quadrature rules
   */
  template<typename ct, int dim>
  class QuadratureRule : public std::vector<QuadraturePoint<ct,dim> >
  {
  public:
    QuadratureRule() : delivered_order(-1) {}
    QuadratureRule(GeometryType t) : geometry_type(t), delivered_order(-1) {}
    QuadratureRule(GeometryType t, int order) : geometry_type(t), delivered_order(order) {}
    // compile time parameters
    enum { d=dim };
    typedef ct CoordType;

    //! return order
    virtual int order () const { return delivered_order; }

    //! return type of element
    virtual GeometryType type () const { return geometry_type; }
    virtual ~QuadratureRule(){}

    //! this container is always a const container,
    //! therefor iterator is the same as const_iterator
    typedef typename std::vector<QuadraturePoint<ct,dim> >::const_iterator iterator;

  protected:
    GeometryType geometry_type;
    int delivered_order;

    void tensor_product (const QuadratureRule<ct,1> & q1d)
    {
      // fill in all the gauss points
      int m = q1d.size();
      int n = power(m,dim);
      for (int i=0; i<n; i++)
      {
        // compute multi index for Gauss point
        int x[dim];
        int z = i;
        for (int k=0; k<dim; ++k)
        {
          x[k] = z%m;
          z = z/m;
        }

        // compute coordinates and weight
        double weight = 1.0;
        FieldVector<ct, dim> local;
        for (int k=0; k<dim; k++)
        {
          local[k] = q1d[x[k]].position()[0];
          weight *= q1d[x[k]].weight();
        }

        // put in container
        push_back(QuadraturePoint<ct,dim>(local,weight));
      }
    }

    int power (int y, int d)
    {
      int m=1;
      for (int i=0; i<d; i++) m *= y;
      return m;
    }
  };

  // Forward declaration of the factory class,
  // needed internally by the QuadratureRules container class.
  template<typename ctype, int dim> struct QuadratureRuleFactory;

  /**
     A container for all quadrature rules of dimension dim
   */
  template<typename ctype, int dim>
  struct QuadratureRules {
    typedef std::pair<GeometryType,int> QuadratureRuleKey;
    typedef QuadratureRule<ctype, dim> QuadratureRule;
    //! real rule creator
    const QuadratureRule& _rule(const GeometryType& t, int p)
    {
      static std::map<QuadratureRuleKey, QuadratureRule> _quadratureMap;
      QuadratureRuleKey key(t,p);
      if (_quadratureMap.find(key) == _quadratureMap.end()) {
        _quadratureMap[key] = QuadratureRuleFactory<ctype,dim>::rule(t,p,QuadratureType::Gauss);
      }
      return _quadratureMap[key];
    }
    //! singleton provider
    static QuadratureRules& instance()
    {
      static QuadratureRules instance;
      return instance;
    }
    //! private constructor
    QuadratureRules () {};
  public:
    //! select the appropriate QuadratureRule for GeometryType t and order p
    static const QuadratureRule& rule(const GeometryType& t, int p)
    {
      return instance()._rule(t,p);
    }
    //! @copydoc rule
    static const QuadratureRule& rule(const GeometryType::BasicType t, int p)
    {
      GeometryType gt(t,dim);
      return instance()._rule(gt,p);
    }
  };

  /***********************************************************/

  //! Gauss quadrature for the n-dimensional cube
  template<typename ct, int dim>
  class CubeQuadratureRule :
    public QuadratureRule<ct,dim>
  {
  public:
    // compile time parameters
    enum { d=dim };
    enum { highest_order=CubeQuadratureRule<ct,1>::highest_order };
    typedef ct CoordType;
    typedef CubeQuadratureRule value_type;

    ~CubeQuadratureRule(){}
  private:
    friend class QuadratureRuleFactory<ct,dim>;
    //! set up quadrature of given order in d dimensions
    CubeQuadratureRule (int p) : QuadratureRule<ct,dim>(GeometryType(GeometryType::cube, d))
    {
      const QuadratureRule<ct,1> & q1D = QuadratureRules<ct,1>::rule(GeometryType::cube, p);
      tensor_product( q1D );
      this->delivered_order = q1D.order();
    }

  };

  //! @copydoc CubeQuadratureRule
  //! Specialization for 1D.
  template<typename ct>
  class CubeQuadratureRule<ct,1> :
    public QuadratureRule<ct,1>
  {
  public:
    // compile time parameters
    enum { d=1 };
    enum { dim=1 };
    enum { highest_order=43 };
    typedef ct CoordType;
    typedef CubeQuadratureRule value_type;

    ~CubeQuadratureRule(){}
  private:
    friend class QuadratureRuleFactory<ct,dim>;
    CubeQuadratureRule (int p);
  };

  /************************************************
   * Quadraturerule for Simplices/Triangle
   *************************************************/

  template<typename ct, int dim>
  class SimplexQuadratureRule;

  template<typename ct>
  class SimplexQuadratureRule<ct,2> : public QuadratureRule<ct,2>
  {
  public:
    enum {d=2};
    enum { highest_order=12 };
    typedef ct CoordType;
    typedef SimplexQuadratureRule value_type;
    ~SimplexQuadratureRule(){}
  private:
    friend class QuadratureRuleFactory<ct,d>;
    SimplexQuadratureRule (int p);
  };

  template<typename ct>
  class SimplexQuadratureRule<ct,3> : public QuadratureRule<ct,3>
  {

  public:
    enum {d=3};
    enum { highest_order=5 };
    typedef ct CoordType;
    typedef SimplexQuadratureRule<ct,3> value_type;
    ~SimplexQuadratureRule(){}
  private:
    friend class QuadratureRuleFactory<ct,d>;
    SimplexQuadratureRule (int p);
  };

  /***********************************
   * quadrature for Prism
   **********************************/

  //!
  template<int dim>
  class PrismQuadraturePoints;

  template<>
  class PrismQuadraturePoints<3>
  {
  public:
    enum { MAXP=6};
    enum { highest_order=2 };

    //! initialize quadrature points on the interval for all orders
    PrismQuadraturePoints ()
    {
      int m = 0;
      O[m] = 0;

      // polynom degree 0  ???
      m = 6;
      G[m][0][0] = 0.0;
      G[m][0][1] = 0.0;
      G[m][0][2] = 0.0;

      G[m][1][0] = 1.0;
      G[m][1][1] = 0.0;
      G[m][1][2] = 0.0;

      G[m][2][0] = 0.0;
      G[m][2][1] = 1.0;
      G[m][2][2] = 0.0;

      G[m][3][0] = 0.0;
      G[m][3][1] = 0.0;
      G[m][3][2] = 1.0;

      G[m][4][0] = 1.0;
      G[m][4][1] = 0.0;
      G[m][4][2] = 1.0;

      G[m][5][0] = 0.0;
      G[m][5][1] = 0.1;
      G[m][5][2] = 1.0;

      W[m][0] = 0.16666666666666666 / 2.0;
      W[m][1] = 0.16666666666666666 / 2.0;
      W[m][2] = 0.16666666666666666 / 2.0;
      W[m][3] = 0.16666666666666666 / 2.0;
      W[m][4] = 0.16666666666666666 / 2.0;
      W[m][5] = 0.16666666666666666 / 2.0;

      O[m] = 0;  // verify ????????


      // polynom degree 2  ???
      m = 6;
      G[m][0][0] =0.66666666666666666 ;
      G[m][0][1] =0.16666666666666666 ;
      G[m][0][2] =0.211324865405187 ;

      G[m][1][0] = 0.16666666666666666;
      G[m][1][1] =0.66666666666666666 ;
      G[m][1][2] = 0.211324865405187;

      G[m][2][0] = 0.16666666666666666;
      G[m][2][1] = 0.16666666666666666;
      G[m][2][2] = 0.211324865405187;

      G[m][3][0] = 0.66666666666666666;
      G[m][3][1] = 0.16666666666666666;
      G[m][3][2] = 0.788675134594813;

      G[m][4][0] = 0.16666666666666666;
      G[m][4][1] = 0.66666666666666666;
      G[m][4][2] = 0.788675134594813;

      G[m][5][0] = 0.16666666666666666;
      G[m][5][1] = 0.16666666666666666;
      G[m][5][2] = 0.788675134594813;

      W[m][0] = 0.16666666666666666 / 2.0;
      W[m][1] = 0.16666666666666666 / 2.0;
      W[m][2] = 0.16666666666666666 / 2.0;
      W[m][3] = 0.16666666666666666 / 2.0;
      W[m][4] = 0.16666666666666666 / 2.0;
      W[m][5] = 0.16666666666666666 / 2.0;

      O[m] = 2;  // verify ????????

    }

    FieldVector<double, 3> point(int m, int i)
    {
      return G[m][i];
    }

    double weight (int m, int i)
    {
      return W[m][i];
    }

    int order (int m)
    {
      return O[m];
    }

  private:
    FieldVector<double, 3> G[MAXP+1][MAXP]; //positions

    double W[MAXP+1][MAXP];     // weights associated with points
    int O[MAXP+1];              // order of the rule
  };


  /** Singleton holding the Prism Quadrature points  */
  template<int dim>
  struct PrismQuadraturePointsSingleton {
    static PrismQuadraturePoints<3> prqp;
  };
  template<>
  struct PrismQuadraturePointsSingleton<3> {
    static PrismQuadraturePoints<3> prqp;
  };

  template<typename ct, int dim>
  class PrismQuadratureRule;

  template<typename ct>
  class PrismQuadratureRule<ct,3> : public QuadratureRule<ct,3>
  {
  public:
    enum { d=3 };
    enum {
      /* min(Line::order, Triangle::order) */
      highest_order =
        (int)CubeQuadratureRule<ct,1>::highest_order
        < (int)SimplexQuadratureRule<ct,2>::highest_order
        ? (int)CubeQuadratureRule<ct,1>::highest_order
        : (int)SimplexQuadratureRule<ct,2>::highest_order
    };
    typedef ct CoordType;
    typedef PrismQuadratureRule<ct,3> value_type;

    ~PrismQuadratureRule(){}
  private:
    friend class QuadratureRuleFactory<ct,d>;
    PrismQuadratureRule(int p) : QuadratureRule<ct,3>(GeometryType(GeometryType::prism, d))
    {
      if (p>highest_order)
        DUNE_THROW(QuadratureOrderOutOfRange,
                   "QuadratureRule for order " << p << " and GeometryType "
                                               << this->type() << " not available");

      if (p<=2) {
        int m=6;
        this->delivered_order = PrismQuadraturePointsSingleton<3>::prqp.order(m);
        for(int i=0; i<m; ++i)
        {
          FieldVector<ct,3> local;
          for (int k=0; k<d; k++)
            local[k] = PrismQuadraturePointsSingleton<3>::prqp.point(m,i)[k];
          double weight =
            PrismQuadraturePointsSingleton<3>::prqp.weight(m,i);
          // put in container
          push_back(QuadraturePoint<ct,d>(local,weight));
        }
      }
      else {
        const QuadratureRule<ct,2> & triangle = QuadratureRules<ct,2>::rule(GeometryType::simplex, p);
        const QuadratureRule<ct,1> & line = QuadratureRules<ct,1>::rule(GeometryType::cube, p);

        this->delivered_order = std::min(triangle.order(),line.order());

        for (typename QuadratureRule<ct,1>::const_iterator
             lit = line.begin(); lit != line.end(); ++lit)
        {
          for (typename QuadratureRule<ct,2>::const_iterator
               tit = triangle.begin(); tit != triangle.end(); ++tit)
          {
            FieldVector<ct, d> local;
            local[0] = tit->position()[0];
            local[1] = tit->position()[1];
            local[2] = lit->position()[0];

            double weight = tit->weight() * lit->weight();

            // put in container
            push_back(QuadraturePoint<ct,d>(local,weight));
          }
        }
      }
    }
  };

  /** Singleton holding the Pyramid Quadrature points  */

  class PyramidQuadraturePoints
  {
  public:
    enum { MAXP=8};
    enum { highest_order=2 };

    //! initialize quadrature points on the interval for all orders
    PyramidQuadraturePoints()
    {
      int m = 0;
      O[m] = 0;


      // polynom degree 2  ???
      m = 8;
      G[m][0][0] =0.58541020;
      G[m][0][1] =0.72819660;
      G[m][0][2] =0.13819660;

      G[m][1][0] =0.13819660;
      G[m][1][1] =0.72819660;
      G[m][1][2] =0.13819660;

      G[m][2][0] =0.13819660;
      G[m][2][1] =0.27630920;
      G[m][2][2] =0.58541020;

      G[m][3][0] =0.13819660;
      G[m][3][1] =0.27630920;
      G[m][3][2] =0.13819660;

      G[m][4][0] =0.72819660;
      G[m][4][1] =0.13819660;
      G[m][4][2] =0.13819660;

      G[m][5][0] =0.72819660;
      G[m][5][1] =0.58541020;
      G[m][5][2] =0.13819660;

      G[m][6][0] =0.27630920;
      G[m][6][1] =0.13819660;
      G[m][6][2] =0.58541020;

      G[m][7][0] =0.27630920;
      G[m][7][1] =0.13819660;
      G[m][7][2] =0.13819660;

      W[m][0] = 0.125 / 3.0;
      W[m][1] = 0.125 / 3.0;
      W[m][2] = 0.125 / 3.0;
      W[m][3] = 0.125 / 3.0;
      W[m][4] = 0.125 / 3.0;
      W[m][5] = 0.125 / 3.0;
      W[m][6] = 0.125 / 3.0;
      W[m][7] = 0.125 / 3.0;

      O[m] = 2;  // verify ????????

    }

    FieldVector<double, 3> point(int m, int i)
    {
      return G[m][i];
    }

    double weight (int m, int i)
    {
      return W[m][i];
    }

    int order (int m)
    {
      return O[m];
    }

  private:
    FieldVector<double, 3> G[MAXP+1][MAXP];
    double W[MAXP+1][MAXP];     // weights associated with points
    int O[MAXP+1];              // order of the rule
  };

  template<int dim>
  struct PyramidQuadraturePointsSingleton {};
  template<>
  struct PyramidQuadraturePointsSingleton<3> {
    static PyramidQuadraturePoints pyqp;
  };

  template<typename ct, int dim>
  class PyramidQuadratureRule;

  template<typename ct>
  class PyramidQuadratureRule<ct,3> : public QuadratureRule<ct,3>
  {
  public:
    enum {d=3};
    enum {highest_order=PyramidQuadraturePoints::highest_order};
    typedef ct CoordType;
    typedef PyramidQuadratureRule<ct,3> value_type;

    ~PyramidQuadratureRule(){}
  private:
    friend class QuadratureRuleFactory<ct,d>;
    PyramidQuadratureRule(int p) : QuadratureRule<ct,3>(GeometryType(GeometryType::pyramid, d))
    {
      int m;

      if (p>highest_order)
        DUNE_THROW(QuadratureOrderOutOfRange,
                   "QuadratureRule for order " << p << " and GeometryType "
                                               << this->type() << " not available");

      switch(p)
      {
      default : m=8;
      }

      this->delivered_order = PyramidQuadraturePointsSingleton<3>::pyqp.order(m);
      FieldVector<ct, d> local;
      double weight;
      for(int i=0; i<m; ++i)
      {
        for(int k=0; k<d; ++k)
          local[k]=PyramidQuadraturePointsSingleton<3>::pyqp.point(m,i)[k];
        weight=PyramidQuadraturePointsSingleton<3>::pyqp.weight(m,i);
        // put in container
        push_back(QuadraturePoint<ct,d>(local,weight));
      }
    }
  };

  /**
     Factory class for creation of quadrature rules,
     depending on GeometryType, order and QuadratureType.

     The whole class is private and can only be accessed
     by the singleton container class QuadratureRules.
   */

  template<typename ctype, int dim>
  class QuadratureRuleFactory {
  private:
    friend class QuadratureRules<ctype, dim>;
    static QuadratureRule<ctype, dim> rule(const GeometryType& t, int p, QuadratureType::Enum qt)
    {
      if (t.isCube())
      {
        return CubeQuadratureRule<ctype,dim>(p);
      }
      if (t.isSimplex())
      {
        return SimplexQuadratureRule<ctype,dim>(p);
      }
      DUNE_THROW(Exception, "Unknown GeometryType");
    }
  };

  template<typename ctype>
  class QuadratureRuleFactory<ctype, 1> {
  private:
    enum { dim = 1 };
    friend class QuadratureRules<ctype, dim>;
    static QuadratureRule<ctype, dim> rule(const GeometryType& t, int p, QuadratureType::Enum qt)
    {
      if (t.isLine())
      {
        switch (qt) {
        case QuadratureType::Gauss :
          return CubeQuadratureRule<ctype,dim>(p);
        default :
          DUNE_THROW(Exception, "Unknown QuadratureType");
        }
      }
      DUNE_THROW(Exception, "Unknown GeometryType");
    }
  };

  template<typename ctype>
  class QuadratureRuleFactory<ctype, 3> {
  private:
    enum { dim = 3 };
    friend class QuadratureRules<ctype, dim>;
    static QuadratureRule<ctype, dim> rule(const GeometryType& t, int p, QuadratureType::Enum qt)
    {
      if (t.isCube())
      {
        return CubeQuadratureRule<ctype,dim>(p);
      }
      if (t.isSimplex())
      {
        return SimplexQuadratureRule<ctype,dim>(p);
      }
      if (t.isPrism())
      {
        return PrismQuadratureRule<ctype,dim>(p);
      }
      if (t.isPyramid())
      {
        return PyramidQuadratureRule<ctype,dim>(p);
      }
      DUNE_THROW(Exception, "Unknown GeometryType");
    }
  };

} // end namespace

#endif
