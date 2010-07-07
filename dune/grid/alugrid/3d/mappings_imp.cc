// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALUGRID_MAPPINGS_IMP_CC
#define DUNE_ALUGRID_MAPPINGS_IMP_CC

#include "mappings.hh"

namespace Dune {

  //- Trilinear mapping (from alu3dmappings.hh)
  alu_inline TrilinearMapping ::
  TrilinearMapping (const coord_t& p0, const coord_t& p1,
                    const coord_t& p2, const coord_t& p3,
                    const coord_t& p4, const coord_t& p5,
                    const coord_t& p6, const coord_t& p7)
  {
    buildMapping(p0,p1,p2,p3,p4,p5,p6,p7);
    return ;
  }

  template <class vector_t>
  alu_inline void TrilinearMapping ::
  buildMapping(const vector_t& p0, const vector_t& p1,
               const vector_t& p2, const vector_t& p3,
               const vector_t& p4, const vector_t& p5,
               const vector_t& p6, const vector_t& p7)
  {
    // build mapping
    a [0][0] = p0 [0] ;
    a [0][1] = p0 [1] ;
    a [0][2] = p0 [2] ;
    a [1][0] = p1 [0] - p0 [0] ;
    a [1][1] = p1 [1] - p0 [1] ;
    a [1][2] = p1 [2] - p0 [2] ;
    a [2][0] = p2 [0] - p0 [0] ;
    a [2][1] = p2 [1] - p0 [1] ;
    a [2][2] = p2 [2] - p0 [2] ;
    a [3][0] = p4 [0] - p0 [0] ;
    a [3][1] = p4 [1] - p0 [1] ;
    a [3][2] = p4 [2] - p0 [2] ;
    a [4][0] = p3 [0] - p2 [0] - a [1][0] ;
    a [4][1] = p3 [1] - p2 [1] - a [1][1] ;
    a [4][2] = p3 [2] - p2 [2] - a [1][2] ;
    a [5][0] = p6 [0] - p4 [0] - a [2][0] ;
    a [5][1] = p6 [1] - p4 [1] - a [2][1] ;
    a [5][2] = p6 [2] - p4 [2] - a [2][2] ;
    a [6][0] = p5 [0] - p1 [0] - a [3][0] ;
    a [6][1] = p5 [1] - p1 [1] - a [3][1] ;
    a [6][2] = p5 [2] - p1 [2] - a [3][2] ;
    a [7][0] = p7 [0] - p5 [0] + p4 [0] - p6 [0] - p3 [0] + p1 [0] + a [2][0] ;
    a [7][1] = p7 [1] - p5 [1] + p4 [1] - p6 [1] - p3 [1] + p1 [1] + a [2][1] ;
    a [7][2] = p7 [2] - p5 [2] + p4 [2] - p6 [2] - p3 [2] + p1 [2] + a [2][2] ;

    {
      alu3d_ctype sum = 0.0;
      // sum all factor from non-linear terms
      for(int i=4; i<8; ++i)
      {
        for(int j=0; j<3; ++j)
        {
          sum += fabs(a[i][j]);
        }
      }

      // mapping is affine when all higher terms are zero
      affine_ = (sum < _epsilon);
    }

    // initialize flags
    calcedDet_ = calcedLinear_ = calcedInv_ = false;

    return ;
  }

  alu_inline TrilinearMapping :: TrilinearMapping (const TrilinearMapping & map)
  {
    // copy mapping
    for (int i = 0 ; i < 8 ; ++i)
      for (int j = 0 ; j < 3 ; ++j)
        a [i][j] = map.a [i][j] ;
    // copy flags
    affine_ = map.affine_;
    calcedDet_ = calcedLinear_ = calcedInv_ = false;
    return ;
  }

  alu_inline const FieldMatrix<alu3d_ctype, 3, 3>&
  TrilinearMapping::jacobianTransposed(const coord_t& p)
  {
    linear( p );
    return Df;
  }

  alu_inline const FieldMatrix<alu3d_ctype, 3, 3>&
  TrilinearMapping::jacobianInverseTransposed(const coord_t& p)
  {
    // calculate inverse if not calculated or not affine
    inverse (p);

    // return reference when already calculated
    if( calcedInv_ )
    {
      return Dfi;
    }
    else
    {
      // make a copy since Dfi could change during world2map
      invTransposed_ = Dfi;
      return invTransposed_;
    }
  }

  alu_inline void TrilinearMapping ::
  map2world(const coord_t& p, coord_t& world) const
  {
    map2world(p[0], p[1], p[2], world);
    return ;
  }

  alu_inline void TrilinearMapping ::
  map2world(const alu3d_ctype x, const alu3d_ctype y,
            const alu3d_ctype z, coord_t& world ) const
  {
    const alu3d_ctype yz  = y * z ;
    const alu3d_ctype xz  = x * z ;
    const alu3d_ctype xy  = x * y ;
    const alu3d_ctype xyz = x * yz ;
    world [0] = a [0][0] + a [1][0] * x + a [2][0] * y + a [3][0] * z + a [4][0] * xy + a [5][0] * yz + a [6][0] * xz + a [7][0] * xyz ;
    world [1] = a [0][1] + a [1][1] * x + a [2][1] * y + a [3][1] * z + a [4][1] * xy + a [5][1] * yz + a [6][1] * xz + a [7][1] * xyz ;
    world [2] = a [0][2] + a [1][2] * x + a [2][2] * y + a [3][2] * z + a [4][2] * xy + a [5][2] * yz + a [6][2] * xz + a [7][2] * xyz ;
    return ;
  }

  alu_inline void TrilinearMapping :: linear(const coord_t& p )
  {
    linear(p[0], p[1], p[2]);
  }

  alu_inline void TrilinearMapping :: linear(const alu3d_ctype x,
                                             const alu3d_ctype y,
                                             const alu3d_ctype z)
  {
    if( ! calcedLinear_ )
    {
      const alu3d_ctype yz = y * z ;
      const alu3d_ctype xz = x * z ;
      const alu3d_ctype xy = x * y ;

      // derivatives with respect to x
      Df[0][0] = a[1][0] + y * a[4][0] + z * a[6][0] + yz * a[7][0] ;
      Df[0][1] = a[1][1] + y * a[4][1] + z * a[6][1] + yz * a[7][1] ;
      Df[0][2] = a[1][2] + y * a[4][2] + z * a[6][2] + yz * a[7][2] ;

      // derivatives with respect to y
      Df[1][0] = a[2][0] + x * a[4][0] + z * a[5][0] + xz * a[7][0] ;
      Df[1][1] = a[2][1] + x * a[4][1] + z * a[5][1] + xz * a[7][1] ;
      Df[1][2] = a[2][2] + x * a[4][2] + z * a[5][2] + xz * a[7][2] ;

      // derivatives with respect to z
      Df[2][0] = a[3][0] + y * a[5][0] + x * a[6][0] + xy * a[7][0] ;
      Df[2][1] = a[3][1] + y * a[5][1] + x * a[6][1] + xy * a[7][1] ;
      Df[2][2] = a[3][2] + y * a[5][2] + x * a[6][2] + xy * a[7][2] ;

      // set calced det to affine (true if affine false otherwise)
      calcedLinear_ = affine_ ;
    }
  }

  alu_inline alu3d_ctype TrilinearMapping :: det(const coord_t& point )
  {
    // use cached value of determinant
    if( calcedDet_ ) return DetDf;

    //  Determinante der Abbildung f:[-1,1]^3 -> Hexaeder im Punkt point.
    linear (point) ;

    // code generated by maple
    const alu3d_ctype t4  = Df[0][0] * Df[1][1];
    const alu3d_ctype t6  = Df[0][0] * Df[1][2];
    const alu3d_ctype t8  = Df[0][1] * Df[1][0];
    const alu3d_ctype t10 = Df[0][2] * Df[1][0];
    const alu3d_ctype t12 = Df[0][1] * Df[2][0];
    const alu3d_ctype t14 = Df[0][2] * Df[2][0];

    // determinant
    DetDf = (t4*Df[2][2]-t6*Df[2][1]-t8*Df[2][2]+
             t10*Df[2][1]+t12*Df[1][2]-t14*Df[1][1]);

    assert( DetDf > 0 );

    // set calced det to affine (true if affine false otherwise)
    calcedDet_ = affine_ ;
    return DetDf;
  }

  alu_inline void TrilinearMapping :: inverse(const coord_t& point)
  {
    // return when inverse already calculated
    if( calcedInv_ ) return ;

    //  Kramer - Regel, det() rechnet Df und DetDf neu aus.
    const alu3d_ctype val = 1.0 / det(point) ;

    // calculate inverse^T
    Dfi[0][0] = ( Df[1][1] * Df[2][2] - Df[2][1] * Df[1][2] ) * val ;
    Dfi[1][0] = ( Df[2][0] * Df[1][2] - Df[1][0] * Df[2][2] ) * val ;
    Dfi[2][0] = ( Df[1][0] * Df[2][1] - Df[2][0] * Df[1][1] ) * val ;
    Dfi[0][1] = ( Df[2][1] * Df[0][2] - Df[0][1] * Df[2][2] ) * val ;
    Dfi[1][1] = ( Df[0][0] * Df[2][2] - Df[2][0] * Df[0][2] ) * val ;
    Dfi[2][1] = ( Df[2][0] * Df[0][1] - Df[0][0] * Df[2][1] ) * val ;
    Dfi[0][2] = ( Df[0][1] * Df[1][2] - Df[1][1] * Df[0][2] ) * val ;
    Dfi[1][2] = ( Df[1][0] * Df[0][2] - Df[0][0] * Df[1][2] ) * val ;
    Dfi[2][2] = ( Df[0][0] * Df[1][1] - Df[1][0] * Df[0][1] ) * val ;

    // set calcedInv_ to affine (true if affine false otherwise)
    calcedInv_ = affine_;
    return ;
  }

  alu_inline void TrilinearMapping::world2map (const coord_t& wld , coord_t& map )
  {
    //  Newton - Iteration zum Invertieren der Abbildung f.
    double err = 10.0 * _epsilon ;
#ifndef NDEBUG
    int count = 0 ;
#endif
    map [0] = map [1] = map [2] = .0 ;
    coord_t upd ;
    do {
      // do mapping
      map2world (map, upd) ;
      // get inverse
      inverse ( map ) ;
      const alu3d_ctype u0 = upd [0] - wld [0] ;
      const alu3d_ctype u1 = upd [1] - wld [1] ;
      const alu3d_ctype u2 = upd [2] - wld [2] ;

      // jacobian is stored as transposed
      const alu3d_ctype c0 = Dfi [0][0] * u0 + Dfi [1][0] * u1 + Dfi [2][0] * u2 ;
      const alu3d_ctype c1 = Dfi [0][1] * u0 + Dfi [1][1] * u1 + Dfi [2][1] * u2 ;
      const alu3d_ctype c2 = Dfi [0][2] * u0 + Dfi [1][2] * u1 + Dfi [2][2] * u2 ;
      map [0] -= c0 ;
      map [1] -= c1 ;
      map [2] -= c2 ;
      err = std::fabs (c0) + std::fabs (c1) + std::fabs (c2) ;
      assert (count ++ < 1000) ;
    } while (err > _epsilon) ;
    return ;
  }

  //- Bilinear surface mapping
  // Constructor for FieldVectors
  alu_inline SurfaceNormalCalculator :: SurfaceNormalCalculator()
  {
    alu3d_ctype p[3] = {0.0,0.0,0.0};
    //initialize with zero
    buildMapping(p,p,p,p);
  }

  // the real constructor, this can be called for FieldVectors
  // and double[3], we dont have to convert one type
  template <class vector_t>
  alu_inline void SurfaceNormalCalculator ::
  buildMapping  (const vector_t & _p0, const vector_t & _p1,
                 const vector_t & _p2, const vector_t & _p3)
  {
    alu3d_ctype b[4][3];
    buildMapping( _p0, _p1, _p2, _p3, b );
  }

  // the real constructor, this can be called for FieldVectors
  // and double[3], we dont have to convert one type
  template <class vector_t>
  alu_inline void SurfaceNormalCalculator ::
  buildMapping  (const vector_t & _p0, const vector_t & _p1,
                 const vector_t & _p2, const vector_t & _p3,
                 alu3d_ctype (&_b)[4][3])
  {

    _b [0][0] = _p0 [0] ;
    _b [0][1] = _p0 [1] ;
    _b [0][2] = _p0 [2] ;
    _b [1][0] = _p1 [0] - _p0 [0] ;
    _b [1][1] = _p1 [1] - _p0 [1] ;
    _b [1][2] = _p1 [2] - _p0 [2] ;
    _b [2][0] = _p2 [0] - _p0 [0] ;
    _b [2][1] = _p2 [1] - _p0 [1] ;
    _b [2][2] = _p2 [2] - _p0 [2] ;
    _b [3][0] = _p3 [0] - _p2 [0] - _b [1][0] ;
    _b [3][1] = _p3 [1] - _p2 [1] - _b [1][1] ;
    _b [3][2] = _p3 [2] - _p2 [2] - _b [1][2] ;

    _n [0][0] = _b [1][1] * _b [2][2] - _b [1][2] * _b [2][1] ;
    _n [0][1] = _b [1][2] * _b [2][0] - _b [1][0] * _b [2][2] ;
    _n [0][2] = _b [1][0] * _b [2][1] - _b [1][1] * _b [2][0] ;
    _n [1][0] = _b [1][1] * _b [3][2] - _b [1][2] * _b [3][1] ;
    _n [1][1] = _b [1][2] * _b [3][0] - _b [1][0] * _b [3][2] ;
    _n [1][2] = _b [1][0] * _b [3][1] - _b [1][1] * _b [3][0] ;
    _n [2][0] = _b [3][1] * _b [2][2] - _b [3][2] * _b [2][1] ;
    _n [2][1] = _b [3][2] * _b [2][0] - _b [3][0] * _b [2][2] ;
    _n [2][2] = _b [3][0] * _b [2][1] - _b [3][1] * _b [2][0] ;


    {
      alu3d_ctype sum = 0.0;
      // sum all factor from non-linear terms
      for(int j=0; j<3; ++j)
      {
        sum += fabs(_b[3][j]);
      }

      // mapping is affine when all higher terms are zero
      _affine = (sum < _epsilon);
    }

    return ;
  }

  alu_inline SurfaceNormalCalculator ::
  SurfaceNormalCalculator(const SurfaceNormalCalculator & m)
  {
    // copy n
    {
      for (int i = 0 ; i < 3 ; ++i)
        for (int j = 0 ; j < 3 ; ++j)
          _n [i][j] = m._n [i][j] ;
    }
    _affine = m._affine;
    return ;
  }

  alu_inline void SurfaceNormalCalculator::
  normal (const coord2_t& map, coord3_t& norm) const
  {
    normal(map[0],map[1],norm);
    return ;
  }

  alu_inline void SurfaceNormalCalculator ::
  normal (const alu3d_ctype x, const alu3d_ctype y, coord3_t& norm) const {
    norm [0] = -(_n [0][0] + _n [1][0] * x + _n [2][0] * y);
    norm [1] = -(_n [0][1] + _n [1][1] * x + _n [2][1] * y);
    norm [2] = -(_n [0][2] + _n [1][2] * x + _n [2][2] * y);
    return ;
  }

  alu_inline void SurfaceNormalCalculator ::
  negativeNormal (const coord2_t& map, coord3_t& norm) const
  {
    negativeNormal(map[0],map[1],norm);
    return ;
  }

  alu_inline void SurfaceNormalCalculator ::
  negativeNormal(const alu3d_ctype x, const alu3d_ctype y, coord3_t& norm) const {
    norm [0] = (_n [0][0] + _n [1][0] * x + _n [2][0] * y);
    norm [1] = (_n [0][1] + _n [1][1] * x + _n [2][1] * y);
    norm [2] = (_n [0][2] + _n [1][2] * x + _n [2][2] * y);
    return ;
  }



  // BilinearSurfaceMapping
  // ----------------------

  // Constructor for FieldVectors
  alu_inline BilinearSurfaceMapping ::
  BilinearSurfaceMapping ()
  {
    alu3d_ctype p[3] = {0.0,0.0,0.0};
    //initialize with zero
    buildMapping(p,p,p,p);
  }

  //- Bilinear surface mapping
  // Constructor for FieldVectors
  alu_inline BilinearSurfaceMapping ::
  BilinearSurfaceMapping (const coord3_t& x0, const coord3_t& x1,
                          const coord3_t& x2, const coord3_t& x3)
  {
    buildMapping(x0,x1,x2,x3);
  }

  // Constructor for double[3]
  alu_inline BilinearSurfaceMapping ::
  BilinearSurfaceMapping (const double3_t & x0, const double3_t & x1,
                          const double3_t & x2, const double3_t & x3)
  {
    buildMapping(x0,x1,x2,x3);
  }

  // the real constructor, this can be called for FieldVectors
  // and double[3], we dont have to convert one type
  template <class vector_t>
  alu_inline void BilinearSurfaceMapping ::
  buildMapping  (const vector_t & _p0, const vector_t & _p1,
                 const vector_t & _p2, const vector_t & _p3)
  {
    BaseType :: buildMapping( _p0, _p1, _p2, _p3, _b );
    // initialize flags
    _calcedInv = _calcedTransposed = _calcedMatrix = false ;
    return ;
  }

  alu_inline BilinearSurfaceMapping ::
  BilinearSurfaceMapping (const BilinearSurfaceMapping & m)
    : BaseType(m)
  {
    // copy _b
    {
      for (int i = 0 ; i < 4 ; ++i)
        for (int j = 0 ; j < 3 ; ++j)
          _b [i][j] = m._b [i][j] ;
    }

    // initialize flags
    _calcedInv = _calcedTransposed = _calcedMatrix = false ;
    return ;
  }

  alu_inline void BilinearSurfaceMapping ::
  map2world (const coord2_t& map, coord3_t& wld) const
  {
    map2world(map[0],map[1],wld);
  }

  alu_inline void BilinearSurfaceMapping ::
  map2world (const alu3d_ctype x, const alu3d_ctype y, coord3_t& w) const
  {
    const alu3d_ctype xy = x * y ;
    w[0] = _b [0][0] + x * _b [1][0] + y * _b [2][0] + xy * _b [3][0] ;
    w[1] = _b [0][1] + x * _b [1][1] + y * _b [2][1] + xy * _b [3][1] ;
    w[2] = _b [0][2] + x * _b [1][2] + y * _b [2][2] + xy * _b [3][2] ;
    return ;
  }


  alu_inline void BilinearSurfaceMapping ::
  map2worldnormal (const alu3d_ctype x,
                   const alu3d_ctype y,
                   const alu3d_ctype z,
                   coord3_t& w) const
  {
    normal(x,y,normal_);

    const alu3d_ctype xy = x * y ;
    w[0] = _b [0][0] + x * _b [1][0] + y * _b [2][0] + xy * _b [3][0] + z*normal_[0];
    w[1] = _b [0][1] + x * _b [1][1] + y * _b [2][1] + xy * _b [3][1] + z*normal_[1];
    w[2] = _b [0][2] + x * _b [1][2] + y * _b [2][2] + xy * _b [3][2] + z*normal_[2];
    return ;
  }

  alu_inline void BilinearSurfaceMapping ::
  map2worldlinear(const alu3d_ctype x, const alu3d_ctype y, const alu3d_ctype z) const
  {
    normal(x,y,normal_);

    Df[0][0] = _b [1][0] + y * _b [3][0] + z * _n[1][0] ;
    Df[1][0] = _b [1][1] + y * _b [3][1] + z * _n[1][1] ;
    Df[2][0] = _b [1][2] + y * _b [3][2] + z * _n[1][2] ;

    Df[0][1] = _b [2][0] + x * _b [3][0] + z * _n[2][0] ;
    Df[1][1] = _b [2][1] + x * _b [3][1] + z * _n[2][1] ;
    Df[2][1] = _b [2][2] + x * _b [3][2] + z * _n[2][2] ;

    Df[0][2] = normal_[0];
    Df[1][2] = normal_[1];
    Df[2][2] = normal_[2];

    return ;
  }

  alu_inline const BilinearSurfaceMapping:: matrix_t&
  BilinearSurfaceMapping::jacobianTransposed(const coord2_t & local) const
  {
    if( ! _calcedMatrix )
    {
      const alu3d_ctype x = local[0];
      const alu3d_ctype y = local[1];

      matrix_[0][0] = _b [1][0] + y * _b [3][0] ;
      matrix_[0][1] = _b [1][1] + y * _b [3][1] ;
      matrix_[0][2] = _b [1][2] + y * _b [3][2] ;

      matrix_[1][0] = _b [2][0] + x * _b [3][0] ;
      matrix_[1][1] = _b [2][1] + x * _b [3][1] ;
      matrix_[1][2] = _b [2][2] + x * _b [3][2] ;

      // only true for affine mappings
      _calcedMatrix = _affine ;
    }

    return matrix_;
  }

  // calculates determinant of face mapping
  alu_inline alu3d_ctype BilinearSurfaceMapping :: det(const coord2_t& point ) const
  {
    // calculate normal
    normal(point[0], point[1], normal_);

    // return length
    return normal_.two_norm();
  }

  alu_inline void BilinearSurfaceMapping :: inverse(const coord3_t& point ) const
  {
    if( _calcedInv ) return ;

    //  Determinante der Abbildung f:[-1,1]^3 -> Hexaeder im Punkt point.
    map2worldlinear (point[0],point[1],point[2]) ;

    //  Kramer - Regel, det() rechnet Df und DetDf neu aus.
    const alu3d_ctype val = 1.0 / Df.determinant();

    Dfi[0][0] = ( Df[1][1] * Df[2][2] - Df[1][2] * Df[2][1] ) * val ;
    Dfi[0][1] = ( Df[0][2] * Df[2][1] - Df[0][1] * Df[2][2] ) * val ;
    Dfi[0][2] = ( Df[0][1] * Df[1][2] - Df[0][2] * Df[1][1] ) * val ;
    Dfi[1][0] = ( Df[1][2] * Df[2][0] - Df[1][0] * Df[2][2] ) * val ;
    Dfi[1][1] = ( Df[0][0] * Df[2][2] - Df[0][2] * Df[2][0] ) * val ;
    Dfi[1][2] = ( Df[0][2] * Df[1][0] - Df[0][0] * Df[1][2] ) * val ;
    Dfi[2][0] = ( Df[1][0] * Df[2][1] - Df[1][1] * Df[2][0] ) * val ;
    Dfi[2][1] = ( Df[0][1] * Df[2][0] - Df[0][0] * Df[2][1] ) * val ;
    Dfi[2][2] = ( Df[0][0] * Df[1][1] - Df[0][1] * Df[1][0] ) * val ;

    // only true for affine mappings
    _calcedInv = _affine ;
    return ;
  }

  alu_inline const BilinearSurfaceMapping:: inv_t&
  BilinearSurfaceMapping::jacobianInverseTransposed(const coord2_t & local) const
  {
    // if calculated return
    if( _calcedTransposed) return invTransposed_;

    tmp_[0] = local[0];
    tmp_[1] = local[1];
    tmp_[2] = 0.0;

    inverse (tmp_) ;

    // calculate transposed inverse
    invTransposed_[0][0] = Dfi[0][0];
    invTransposed_[0][1] = Dfi[1][0];

    invTransposed_[1][0] = Dfi[0][1];
    invTransposed_[1][1] = Dfi[1][1];

    invTransposed_[2][0] = Dfi[0][2];
    invTransposed_[2][1] = Dfi[1][2];

    // only true for affine mappings
    _calcedTransposed = _affine ;

    return invTransposed_;
  }

  alu_inline void BilinearSurfaceMapping::world2map (const coord3_t& wld , coord2_t& map ) const
  {
    //  Newton - Iteration zum Invertieren der Abbildung f.
    double err = 10.0 * _epsilon ;
    coord3_t map_ (0);
#ifndef NDEBUG
    int count = 0 ;
#endif
    coord3_t upd ;
    do {
      // apply mapping
      map2worldnormal (map_[0],map_[1],map_[2], upd) ;
      // calculate inverse
      inverse (map_) ;
      const alu3d_ctype u0 = upd [0] - wld [0] ;
      const alu3d_ctype u1 = upd [1] - wld [1] ;
      const alu3d_ctype u2 = upd [2] - wld [2] ;
      const alu3d_ctype c0 = Dfi [0][0] * u0 + Dfi [0][1] * u1 + Dfi [0][2] * u2 ;
      const alu3d_ctype c1 = Dfi [1][0] * u0 + Dfi [1][1] * u1 + Dfi [1][2] * u2 ;
      const alu3d_ctype c2 = Dfi [2][0] * u0 + Dfi [2][1] * u1 + Dfi [2][2] * u2 ;
      map_ [0] -= c0 ;
      map_ [1] -= c1 ;
      map_ [2] -= c2 ;
      err = fabs (c0) + fabs (c1) + fabs (c2) ;
      assert (count ++ < 3000);
    }
    while (err > _epsilon) ;

    // get local coordinates
    map[0] = map_[0];
    map[1] = map_[1];
    return ;
  }



  // BilinearMapping
  // ---------------

  template< int cdim >
  alu_inline BilinearMapping< cdim >::BilinearMapping ()
    : calcedMatrix_( false ),
      calcedDet_( false ),
      calcedInv_( false )
  {
    for( int i = 0; i < 4; ++i )
      for( int j = 0; j < cdim; ++j )
        _b[ i ][ j ] = ctype( 0 );
    affine_ = true;
  }


  template< int cdim >
  alu_inline BilinearMapping< cdim >
  ::BilinearMapping ( const world_t &p0, const world_t &p1,
                      const world_t &p2, const world_t &p3 )
  {
    buildMapping( p0, p1, p2, p3 );
  }


  template< int cdim >
  alu_inline BilinearMapping< cdim >
  ::BilinearMapping ( const ctype (&p0)[ cdim ], const ctype (&p1)[ cdim ],
                      const ctype (&p2)[ cdim ], const ctype (&p3)[ cdim ] )
  {
    buildMapping( p0, p1, p2, p3 );
  }


  template< int cdim >
  alu_inline bool BilinearMapping< cdim >::affine () const
  {
    return affine_;
  }


  template< int cdim >
  alu_inline void BilinearMapping< cdim >::map2world ( const map_t &m, world_t &w ) const
  {
    map2world( m[0], m[1], w );
  }


  template< int cdim >
  alu_inline void BilinearMapping< cdim >::map2world ( const ctype x, const ctype y, world_t &w ) const
  {
    const alu3d_ctype xy = x * y;
    for( int i = 0; i < cdim; ++i )
      w[i] = _b [0][i] + x * _b [1][i] + y * _b [2][i] + xy * _b [3][i];
  }


  template< int cdim >
  alu_inline void BilinearMapping< cdim >::world2map ( const world_t &w, map_t &m ) const
  {
    world_t dw;
    m[0] = m[1] = 0.5;
    map2world( m, dw );
    dw -= w;

    alu3d_ctype step;
    do {
      step = 0;
      for( int j = 0; j < 2; ++j )
      {
        world_t d;
        for( int i = 0; i < cdim; ++i )
          d[i] = _b [j+1][i] + m[1-j] * _b [3][i];
        const ctype dd = d*d;
        if( dd < ALUnumericEpsilon*ALUnumericEpsilon )
          continue;

        const ctype alpha = -(dw*d) / dd;
        dw.axpy( alpha, d );
        m[j] += alpha;
        step += std::abs( alpha );
      }
    }
    while( step > ALUnumericEpsilon );
    //while( dw.two_norm2() > ALUnumericEpsilon*ALUnumericEpsilon );
  }


  template< int cdim >
  alu_inline typename BilinearMapping< cdim >::ctype
  BilinearMapping< cdim >::det ( const map_t &m ) const
  {
    inverse( m );
    return det_;
  }

  template<>
  alu_inline BilinearMapping< 2 >::ctype
  BilinearMapping< 2 >::det ( const map_t &m ) const
  {
    if( !calcedDet_ )
    {
      map2worldlinear( m[0], m[1] );

      const world_t &u = matrix_[0];
      const world_t &v = matrix_[1];

      det_ = fabs( u[0] * v[1] - u[1] * v[0] );
      calcedDet_ = affine();
    }
    return det_;
  }

  template<>
  alu_inline BilinearMapping< 3 >::ctype
  BilinearMapping< 3 >::det ( const map_t &m ) const
  {
    if( !calcedDet_ )
    {
      map2worldlinear( m[0], m[1] );

      const world_t &u = matrix_[0];
      const world_t &v = matrix_[1];

      world_t n;
      for ( int i = 0; i < 3; ++i )
        n[i] = v[(i+1)%3] * u[(i+2)%3] - v[(i+2)%3] * u[(i+1)%3];

      det_ = n.two_norm();
      calcedDet_ = affine();
    }
    return det_;
  }


  template< int cdim >
  alu_inline const typename BilinearMapping< cdim >::matrix_t &
  BilinearMapping< cdim >::jacobianTransposed ( const map_t &m ) const
  {
    map2worldlinear( m[0], m[1] );
    return matrix_;
  }


  template< int cdim >
  alu_inline const typename BilinearMapping< cdim >::inv_t &
  BilinearMapping< cdim >::jacobianInverseTransposed ( const map_t &m ) const
  {
    inverse( m );
    return invTransposed_;
  }


  template< int cdim >
  template< class vector_t >
  alu_inline void BilinearMapping< cdim >
  ::buildMapping ( const vector_t &p0, const vector_t &p1,
                   const vector_t &p2, const vector_t &p3 )
  {
    for( int i = 0; i < cdim; ++i )
    {
      _b [0][i] = p0 [i];
      _b [1][i] = p1 [i] - p0 [i];
      _b [2][i] = p2 [i] - p0 [i];
      _b [3][i] = p3 [i] - p2 [i] - _b [1][i];
    }

    ctype sum = fabs( _b [3][0] );
    for( int i = 1; i < cdim; ++i )
      sum += fabs( _b [3][i] );
    affine_ = (sum < ALUnumericEpsilon );

    calcedMatrix_ = calcedDet_ = calcedInv_ = false;
  }


  template< int cdim >
  alu_inline void BilinearMapping< cdim >
  ::multTransposedMatrix ( const matrix_t &A, FieldMatrix< ctype, 2, 2 > &C )
  {
    for( int i = 0; i < 2; ++i )
    {
      for( int j = 0; j < 2; ++j )
      {
        C[i][j] = A[i][0] * A[j][0];
        for( int k=1; k < cdim; ++k )
          C[i][j] += A[i][k] * A[j][k];
      }
    }
  }


  template< int cdim >
  alu_inline void BilinearMapping< cdim >
  ::multMatrix ( const matrix_t &A, const FieldMatrix< ctype, 2, 2 > &B, inv_t &C )
  {
    for( int i = 0; i < cdim; ++i )
      for( int j = 0; j < 2; ++j )
        C[i][j] = A[0][i] * B[0][j] + A[1][i] * B[1][j];
  }


  template< int cdim >
  alu_inline void BilinearMapping< cdim >
  ::map2worldlinear ( const alu3d_ctype x, const alu3d_ctype y ) const
  {
    if( !calcedMatrix_ )
    {
      for( int i = 0; i < cdim; ++i )
      {
        matrix_[0][i] = _b [1][i] + y * _b [3][i];
        matrix_[1][i] = _b [2][i] + x * _b [3][i];
      }
      calcedMatrix_ = affine();
    }
  }


  template<>
  alu_inline void BilinearMapping< 2 >::inverse ( const map_t &m ) const
  {
    if( !calcedInv_ )
    {
      map2worldlinear ( m[0], m[1] );
      det_ = std::abs( FMatrixHelp::invertMatrix( matrix_, invTransposed_ ) );
      calcedDet_ = calcedInv_ = affine();
    }
  }

  template< int cdim >
  alu_inline void BilinearMapping< cdim >::inverse ( const map_t &m ) const
  {
    // use least squares approach
    if( !calcedInv_ )
    {
      FieldMatrix< ctype, 2, 2 > AT_A, inv_AT_A;
      map2worldlinear ( m[0], m[1] );
      multTransposedMatrix( matrix_, AT_A );
      FMatrixHelp::invertMatrix( AT_A, inv_AT_A );
      multMatrix( matrix_, inv_AT_A, invTransposed_ );
      calcedInv_ = affine();
    }
  }



  ////////////////////////////////////////////////////////////////////////
  //
  //  Tetra specializations
  //  -- tetra spec
  //
  ////////////////////////////////////////////////////////////////////////

  //- Bilinear surface mapping
  // Constructor for FieldVectors
  template <int cdim, int mydim>
  alu_inline LinearMapping<cdim, mydim> ::
  LinearMapping ()
    : _matrix ( 0.0 )
      , _invTransposed( 0.0 )
      , _p0( 0.0 )
      , _calcedInv( false )
      , _calcedDet( false )
  {}

  // copy constructor
  template <int cdim, int mydim>
  alu_inline LinearMapping<cdim, mydim> ::
  LinearMapping (const LinearMapping & m)
    : _matrix ( m._matrix )
      , _invTransposed( m._invTransposed )
      , _p0( m._p0 )
      , _calcedInv( m._calcedInv )
      , _calcedDet( m._calcedDet )
  {}

  // the real constructor, this can be called for FieldVectors
  // and double[3], we dont have to convert one type
  template <>
  template <class vector_t>
  alu_inline void LinearMapping<3, 3> ::
  buildMapping  (const vector_t & p0, const vector_t & p1,
                 const vector_t & p2, const vector_t & p3 )
  {
    _matrix [0][0] = p1[0] - p0 [0] ;
    _matrix [0][1] = p1[1] - p0 [1] ;
    _matrix [0][2] = p1[2] - p0 [2] ;

    _matrix [1][0] = p2[0] - p0 [0] ;
    _matrix [1][1] = p2[1] - p0 [1] ;
    _matrix [1][2] = p2[2] - p0 [2] ;

    _matrix [2][0] = p3[0] - p0 [0] ;
    _matrix [2][1] = p3[1] - p0 [1] ;
    _matrix [2][2] = p3[2] - p0 [2] ;

    _p0[0] = p0[0];
    _p0[1] = p0[1];
    _p0[2] = p0[2];

    // initialize flags
    _calcedDet = _calcedInv = false ;
    return ;
  }

  // the real constructor, this can be called for FieldVectors
  // and double[3], we dont have to convert one type
  template <>
  template <class vector_t>
  alu_inline void LinearMapping<3, 2> ::
  buildMapping  (const vector_t & p0, const vector_t & p1,
                 const vector_t & p2)
  {
    _matrix [0][0] = p1[0] - p0 [0] ;
    _matrix [0][1] = p1[1] - p0 [1] ;
    _matrix [0][2] = p1[2] - p0 [2] ;

    _matrix [1][0] = p2[0] - p0 [0] ;
    _matrix [1][1] = p2[1] - p0 [1] ;
    _matrix [1][2] = p2[2] - p0 [2] ;

    _p0[0] = p0[0];
    _p0[1] = p0[1];
    _p0[2] = p0[2];

    // initialize flags
    _calcedDet = _calcedInv = false ;
    return ;
  }

  // the real constructor,
  // this can be called for FieldVectors
  // and double[3], we dont have to convert one type
  template <int cdim, int mydim>
  template <class vector_t>
  alu_inline void LinearMapping<cdim, mydim> ::
  buildMapping  (const vector_t & p0, const vector_t & p1)
  {
    assert( mydim == 1 );

    _matrix [0][0] = p1[0] - p0 [0] ;
    _matrix [0][1] = p1[1] - p0 [1] ;
    _matrix [0][2] = p1[2] - p0 [2] ;

    _p0[0] = p0[0];
    _p0[1] = p0[1];
    _p0[2] = p0[2];

    // initialize flags
    _calcedDet = _calcedInv = false ;
    return ;
  }


  // the real constructor, this can be called for FieldVectors
  // and double[3], we dont have to convert one type
  template <>
  template <class vector_t>
  alu_inline void LinearMapping<2, 2> ::
  buildMapping  (const vector_t & p0,
                 const vector_t & p1,
                 const vector_t & p2)
  {
    _matrix [0][0] = p1[0] - p0 [0] ;
    _matrix [0][1] = p1[1] - p0 [1] ;

    _matrix [1][0] = p2[0] - p0 [0] ;
    _matrix [1][1] = p2[1] - p0 [1] ;

    _p0[0] = p0[0];
    _p0[1] = p0[1];

    // initialize flags
    _calcedDet = _calcedInv = false ;
    return ;
  }

  // the real constructor, this can be called for FieldVectors
  // and double[3], we dont have to convert one type
  template <>
  template <class vector_t>
  alu_inline void LinearMapping<2, 1> ::
  buildMapping  (const vector_t & p0, const vector_t & p1)
  {
    _matrix [0][0] = p1[0] - p0 [0] ;
    _matrix [0][1] = p1[1] - p0 [1] ;

    _p0[0] = p0[0];
    _p0[1] = p0[1];

    _det = std::sqrt( (_matrix [0][0] * _matrix [0][0]) +
                      (_matrix [0][1] * _matrix [0][1]) );

    // initialize flags
    _calcedDet = true;
    _calcedInv = false ;
    return ;
  }

  template<>
  template< class vector_t >
  alu_inline void LinearMapping< 2, 0 >::buildMapping ( const vector_t &p0 )
  {
    _p0[0] = p0[0];
    _p0[1] = p0[1];

    _det = 1.0;

    // initialize flags
    _calcedDet = _calcedInv = true;
  }

  template<>
  template< class vector_t >
  alu_inline void LinearMapping< 3, 0 >::buildMapping ( const vector_t &p0 )
  {
    _p0[0] = p0[0];
    _p0[1] = p0[1];
    _p0[2] = p0[2];

    _det = 1.0;

    // initialize flags
    _calcedDet = _calcedInv = true;
  }

  // local --> global
  template <int cdim, int mydim>
  alu_inline void LinearMapping<cdim, mydim> ::
  map2world (const map_t& local, world_t& global) const
  {
    // initialize
    global = _p0;

    // multiply with (transposed)
    _matrix.umtv(local, global);
  }

  // global --> local
  template <int cdim, int mydim>
  alu_inline void LinearMapping<cdim, mydim> ::
  world2map (const world_t& global, map_t& local) const
  {
    // initialize
    world_t globalCoord( global );
    // substract p0
    globalCoord -= _p0;

    // multiply with jacobian inverse transposed
    jacobianInverseTransposed( local ).mtv(globalCoord, local);
  }


  // tetra mapping
  template <int cdim, int mydim>
  alu_inline void LinearMapping<cdim, mydim> ::
  inverse(const map_t& local) const
  {
    // invert transposed matrix and return determinant
    _det = std::abs( FMatrixHelp::invertMatrix(_matrix , _invTransposed ) );
    // set flag
    _calcedDet = _calcedInv = true ;
  }

  template <>
  alu_inline void LinearMapping<3, 0> ::
  inverse(const map_t& local) const
  {
    // invert transposed matrix and return determinant
    _det = 1.;
    // set flag
    _calcedDet = _calcedInv = true ;
  }
  template <>
  alu_inline void LinearMapping<2, 0> ::
  inverse(const map_t& local) const
  {
    // invert transposed matrix and return determinant
    _det = 1.;
    // set flag
    _calcedDet = _calcedInv = true ;
  }

  // tetra mapping
  template <int cdim, int mydim>
  alu_inline void LinearMapping<cdim, mydim> ::
  calculateDeterminant(const map_t& local) const
  {
    inverse( local );
  }

  // triangle mapping
  template <>
  alu_inline void LinearMapping<3, 2> ::
  inverse(const map_t& local) const
  {
    inverseCodimOne( local );
  }

  // edge mapping
  template <>
  alu_inline void LinearMapping<2, 1> ::
  inverse(const map_t& local) const
  {
    inverseCodimOne( local );
  }

  template <>
  alu_inline void LinearMapping<3, 0> ::
  inverseCodimOne(const map_t& local) const
  {}

  template <>
  alu_inline void LinearMapping<2, 0> ::
  inverseCodimOne(const map_t& local) const
  {}

  template <int cdim, int mydim>
  alu_inline void LinearMapping<cdim, mydim> ::
  inverseCodimOne(const map_t& local) const
  {
    // use least squares approach
    FieldMatrix<ctype, mydim, mydim> AT_A;

    /*
       inv_t matrix;
       for( int i=0; i<cdim; ++i)
       for( int j=0; j<mydim; ++j)
       {
        matrix[i][j] = _matrix[j][i];
       }
     */

    // calc ret = A^T*A
    //FMatrixHelp::multTransposedMatrix( matrix, AT_A);
    multTransposedMatrix(_matrix, AT_A );

    // calc Jinv_ = A (A^T*A)^-1
    FieldMatrix< ctype, mydim, mydim> inv_AT_A;

    FMatrixHelp :: invertMatrix( AT_A, inv_AT_A );
    //FMatrixHelp :: multMatrix( matrix, inv_AT_A, _invTransposed );
    multMatrix( _matrix, inv_AT_A, _invTransposed );

    // set flag
    _calcedInv = true ;
  }

  // triangle mapping
  template <>
  alu_inline void LinearMapping<3, 2> ::
  calculateDeterminant(const map_t& local) const
  {
    enum { cdim  = 3 };
    world_t tmpV; //! temporary memory
    world_t tmpU; //! temporary memory

    for(int i=0; i<cdim; ++i)
    {
      // p1 - p0 (see buildMapping method)
      tmpV[i] = _matrix[0][i];

      // p2 - p1 = (p2 - p0) - (p1 - p0)
      tmpU[i] = _matrix[1][i] - _matrix[0][i];
    }

    world_t globalCoord;

    // calculate scaled outer normal
    for(int i=0; i<cdim; ++i)
    {
      globalCoord[i] = (  tmpV[(i+1)%cdim] * tmpU[(i+2)%cdim]
                          - tmpV[(i+2)%cdim] * tmpU[(i+1)%cdim] );
    }

    // calculate determinant
    _det = globalCoord.two_norm();

    // set flag
    _calcedDet = true ;
  }

  template <int cdim, int mydim>
  alu_inline void LinearMapping<cdim, mydim> ::
  multTransposedMatrix(const matrix_t& matrix,
                       FieldMatrix<ctype, mydim, mydim>& result) const
  {
    typedef typename matrix_t::size_type size_type;
    for(size_type i=0; i<mydim; ++i)
    {
      for(size_type j=0; j<mydim; ++j)
      {
        result[i][j] = 0.0;
        for(size_type k=0; k<cdim; ++k)
        {
          result[i][j] += matrix[i][k] * matrix[j][k];
        }
      }
    }
  }

  template <int cdim, int mydim>
  alu_inline void LinearMapping<cdim, mydim> ::
  multMatrix ( const matrix_t &A,
               const FieldMatrix< ctype, mydim, mydim > &B,
               inv_t& ret ) const
  {
    //! calculates ret = A * B
    typedef typename matrix_t :: size_type size_type;

    for( size_type i = 0; i < cdim; ++i )
    {
      for( size_type j = 0; j < mydim; ++j )
      {
        ret[ i ][ j ] = 0 ;
        for( size_type k = 0; k < mydim; ++k )
          ret[ i ][ j ] += A[ k ][ i ] * B[ k ][ j ];
      }
    }
  }

  // edge mapping
  template <>
  alu_inline void LinearMapping<3, 1> ::
  inverse(const map_t& local) const
  {
    FieldMatrix<ctype, 1, 1> AT_A_;

    // calc ret = A^T*A
    multTransposedMatrix(_matrix, AT_A_ );

    // calc Jinv_ = A (A^T*A)^-1
    FieldMatrix< ctype, 1, 1 > inv_AT_A;
    FMatrixHelp :: invertMatrix( AT_A_, inv_AT_A );
    multMatrix( _matrix, inv_AT_A, _invTransposed );

    // set flag
    _calcedInv = true ;
  }

  // triangle mapping
  template <>
  alu_inline void LinearMapping<3, 1> ::
  calculateDeterminant(const map_t& local) const
  {
    // calculate length
    _det = std::sqrt( (_matrix[0][0] * _matrix[0][0]) +
                      (_matrix[0][1] * _matrix[0][1]) +
                      (_matrix[0][2] * _matrix[0][2]) );

    // set flag
    _calcedDet = true ;
  }
  // triangle mapping
  template <>
  alu_inline void LinearMapping<2, 1> ::
  calculateDeterminant(const map_t& local) const
  {
    // calculate length
    _det = std::sqrt( (_matrix[0][0] * _matrix[0][0]) +
                      (_matrix[0][1] * _matrix[0][1]) );

    // set flag
    _calcedDet = true ;
  }

  template <int cdim, int mydim>
  alu_inline typename LinearMapping< cdim, mydim >::ctype
  LinearMapping<cdim, mydim >::det( const map_t& local ) const
  {
    // return det if already calculated
    if( _calcedDet ) return _det;

    // calculate inverse
    calculateDeterminant( local );

    return _det;
  }

  template <int cdim, int mydim>
  alu_inline const typename LinearMapping<cdim, mydim> :: matrix_t&
  LinearMapping<cdim, mydim> ::
  jacobianTransposed(const map_t & local) const
  {
    return _matrix;
  }

  template <int cdim, int mydim>
  alu_inline const typename LinearMapping<cdim, mydim> :: inv_t&
  LinearMapping<cdim, mydim> ::
  jacobianInverseTransposed(const map_t & local) const
  {
    // if calculated return
    if( _calcedInv ) return _invTransposed;

    // calculate
    inverse ( local ) ;

    return _invTransposed;
  }

  template <>
  alu_inline const LinearMapping<3, 0> :: inv_t&
  LinearMapping<3, 0> ::
  jacobianInverseTransposed(const map_t & local) const
  {
    return _invTransposed;
  }

  template <>
  alu_inline const LinearMapping<2, 0> :: inv_t&
  LinearMapping<2, 0> ::
  jacobianInverseTransposed(const map_t & local) const
  {
    return _invTransposed;
  }

#if COMPILE_ALUGRID_LIB
  // Instantiation
  class TrilinearMapping ;
  template void TrilinearMapping::buildMapping< TrilinearMapping::coord_t >
    (const TrilinearMapping::coord_t& ,
    const TrilinearMapping::coord_t& ,
    const TrilinearMapping::coord_t& ,
    const TrilinearMapping::coord_t& ,
    const TrilinearMapping::coord_t& ,
    const TrilinearMapping::coord_t& ,
    const TrilinearMapping::coord_t& ,
    const TrilinearMapping::coord_t& );
  template void TrilinearMapping::buildMapping< TrilinearMapping::double_t >
    (const TrilinearMapping::double_t& ,
    const TrilinearMapping::double_t& ,
    const TrilinearMapping::double_t& ,
    const TrilinearMapping::double_t& ,
    const TrilinearMapping::double_t& ,
    const TrilinearMapping::double_t& ,
    const TrilinearMapping::double_t& ,
    const TrilinearMapping::double_t& );

  class SurfaceNormalCalculator ;

  template void SurfaceNormalCalculator::buildMapping< SurfaceNormalCalculator::coord3_t >
    (const SurfaceNormalCalculator::coord3_t & ,
    const SurfaceNormalCalculator::coord3_t & ,
    const SurfaceNormalCalculator::coord3_t & ,
    const SurfaceNormalCalculator::coord3_t & );

  template void SurfaceNormalCalculator::buildMapping< SurfaceNormalCalculator::double3_t >
    (const SurfaceNormalCalculator::double3_t & ,
    const SurfaceNormalCalculator::double3_t & ,
    const SurfaceNormalCalculator::double3_t & ,
    const SurfaceNormalCalculator::double3_t & );

  class BilinearSurfaceMapping ;

  template void BilinearSurfaceMapping::buildMapping< BilinearSurfaceMapping::coord3_t >
    (const BilinearSurfaceMapping::coord3_t & ,
    const BilinearSurfaceMapping::coord3_t & ,
    const BilinearSurfaceMapping::coord3_t & ,
    const BilinearSurfaceMapping::coord3_t & );

  template void BilinearSurfaceMapping::buildMapping< BilinearSurfaceMapping::double3_t >
    (const BilinearSurfaceMapping::double3_t & ,
    const BilinearSurfaceMapping::double3_t & ,
    const BilinearSurfaceMapping::double3_t & ,
    const BilinearSurfaceMapping::double3_t & );

  template class LinearMapping<3, 3> ;
  template void LinearMapping<3, 3>::buildMapping< LinearMapping<3, 3>::world_t >
    ( const LinearMapping<3, 3>::world_t&,
    const LinearMapping<3, 3>::world_t&,
    const LinearMapping<3, 3>::world_t&,
    const LinearMapping<3, 3>::world_t& );
  template void LinearMapping<3, 3>::buildMapping< LinearMapping<3, 3>::double_t >
    ( const LinearMapping<3, 3>::double_t&,
    const LinearMapping<3, 3>::double_t&,
    const LinearMapping<3, 3>::double_t&,
    const LinearMapping<3, 3>::double_t& );

  template class LinearMapping<3, 2> ;
  template void LinearMapping<3, 2>::buildMapping< LinearMapping<3, 2>::world_t >
    ( const LinearMapping<3, 2>::world_t&,
    const LinearMapping<3, 2>::world_t&,
    const LinearMapping<3, 2>::world_t&);
  template void LinearMapping<3, 2>::buildMapping< LinearMapping<3, 2>::double_t >
    ( const LinearMapping<3, 2>::double_t&,
    const LinearMapping<3, 2>::double_t&,
    const LinearMapping<3, 2>::double_t&);

  template class LinearMapping<2, 2> ;
  template void LinearMapping<2, 2>::buildMapping< LinearMapping<2, 2>::world_t >
    ( const LinearMapping<2, 2>::world_t&,
    const LinearMapping<2, 2>::world_t&,
    const LinearMapping<2, 2>::world_t&);
  template void LinearMapping<2, 2>::buildMapping< LinearMapping<2, 2>::double_t >
    ( const LinearMapping<2, 2>::double_t&,
    const LinearMapping<2, 2>::double_t&,
    const LinearMapping<2, 2>::double_t&);

  template class LinearMapping<3, 1> ;
  template void LinearMapping<3, 1>::buildMapping< LinearMapping<3, 1>::world_t >
    ( const LinearMapping<3, 1>::world_t&,
    const LinearMapping<3, 1>::world_t& );
  template void LinearMapping<3, 1>::buildMapping< LinearMapping<3, 1>::double_t >
    ( const LinearMapping<3, 1>::double_t&,
    const LinearMapping<3, 1>::double_t& );

  template class LinearMapping<2, 1> ;
  template void LinearMapping<2, 1>::buildMapping< LinearMapping<2, 1>::world_t >
    ( const LinearMapping<2, 1>::world_t&,
    const LinearMapping<2, 1>::world_t& );
  template void LinearMapping<2, 1>::buildMapping< LinearMapping<2, 1>::double_t >
    ( const LinearMapping<2, 1>::double_t&,
    const LinearMapping<2, 1>::double_t& );
  /// wtf?
  template void LinearMapping<2, 1>::buildMapping< LinearMapping<3, 1>::world_t >
    ( const LinearMapping<3, 1>::world_t&,
    const LinearMapping<3, 1>::world_t& );

  template class LinearMapping<3, 0> ;
  template void LinearMapping<3, 0>::buildMapping< LinearMapping<3, 0>::double_t >
    ( const LinearMapping<3, 0>::double_t& );
  template void LinearMapping<3, 0>::buildMapping< LinearMapping<3, 0>::world_t >
    ( const LinearMapping<3, 0>::world_t& );

  template class LinearMapping<2, 0> ;
  template void LinearMapping<2, 0>::buildMapping< LinearMapping<2, 0>::world_t >
    ( const LinearMapping<2, 0>::world_t& );
  template void LinearMapping<2, 0>::buildMapping< LinearMapping<2, 0>::double_t >
    ( const LinearMapping<2, 0>::double_t& );

  template class BilinearMapping< 2 > ;
  template void BilinearMapping< 2 >::buildMapping< BilinearMapping< 2 >::world_t >
    ( const BilinearMapping< 2 >::world_t&,
    const BilinearMapping< 2 >::world_t&,
    const BilinearMapping< 2 >::world_t&,
    const BilinearMapping< 2 >::world_t& );

  template class BilinearMapping< 3 > ;
  template void BilinearMapping< 3 >::buildMapping< BilinearMapping< 3 >::world_t >
    ( const BilinearMapping< 3 >::world_t&,
    const BilinearMapping< 3 >::world_t&,
    const BilinearMapping< 3 >::world_t&,
    const BilinearMapping< 3 >::world_t& );
#endif

} // end namespace Dune
#endif // end DUNE_ALUGRID_MAPPINGS_IMP_CC
