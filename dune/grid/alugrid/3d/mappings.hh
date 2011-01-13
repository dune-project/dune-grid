// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU3DGRIDMAPPINGS_HH
#define DUNE_ALU3DGRIDMAPPINGS_HH

// System includes
#include <limits>
#include <cmath>

// Dune includes
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/exceptions.hh>

// Local includes
#include "alu3dinclude.hh"

namespace Dune {

  static const alu3d_ctype ALUnumericEpsilon = 10.0 * std::numeric_limits< alu3d_ctype >::epsilon();

  template<int mydim, int coorddim, class GridImp>
  class ALU3dGridGeometry;

  template< ALU3dGridElementType, class >
  class ALU3dGrid;

  //! A trilinear mapping from the Dune reference hexahedron into the physical
  //! space (same as in mapp_cube_3d.h, but for a different reference hexahedron)
  class TrilinearMapping
  {
  public:
    typedef alu3d_ctype double_t[3];
    typedef FieldVector<alu3d_ctype, 3> coord_t;
    typedef FieldMatrix<alu3d_ctype, 3, 3> mat_t;
  private:
    static const double _epsilon ;

    // the internal mapping
    alu3d_ctype a [8][3] ;
    mat_t Df;
    mat_t Dfi;
    mat_t invTransposed_;
    alu3d_ctype DetDf ;

    bool calcedDet_;
    bool calcedLinear_;
    bool calcedInv_;
    bool affine_;

    void linear (const alu3d_ctype, const alu3d_ctype, const alu3d_ctype) ;
    void linear (const coord_t&) ;
    void inverse (const coord_t&, const bool checkDet = true ) ;
  public:
    TrilinearMapping (const coord_t&, const coord_t&,
                      const coord_t&, const coord_t&,
                      const coord_t&, const coord_t&,
                      const coord_t&, const coord_t&);

    // only to call from geometry class
    TrilinearMapping () {}

    TrilinearMapping (const TrilinearMapping &) ;

    ~TrilinearMapping () {}
    alu3d_ctype det (const coord_t&, const bool checkDet = true ) ;
    const mat_t& jacobianInverseTransposed(const coord_t&);
    const mat_t& jacobianTransposed(const coord_t&);
    void map2world (const coord_t&, coord_t&) const ;
    void map2world (const alu3d_ctype , const alu3d_ctype , const alu3d_ctype ,
                    coord_t&) const ;
    void world2map (const coord_t&, coord_t&) ;

    template <class vector_t>
    void buildMapping(const vector_t&, const vector_t&,
                      const vector_t&, const vector_t&,
                      const vector_t&, const vector_t&,
                      const vector_t&, const vector_t&);

    // returns true if mapping is affine
    inline bool affine () const { return affine_; }
  };

  //! A bilinear surface mapping
  // NOTE: this class is different to the BilinearSurfaceMapping in
  // ALUGrid, for example the reference elements differ
  // here we have [0,1]^2 and in ALUGrid its [-1,1]^2
  // also the point numbering is different
  class SurfaceNormalCalculator
  {
  public:
    // our coordinate types
    typedef FieldVector<alu3d_ctype, 3> coord3_t;
    typedef FieldVector<alu3d_ctype, 2> coord2_t;

    // type of coordinate vectors from elements
    typedef alu3d_ctype double3_t[3];
  protected:

    alu3d_ctype _n [3][3] ;

    static const double _epsilon ;

    bool _affine;

  public:
    //! Constructor creating empty mapping with double , i.e. zero
    SurfaceNormalCalculator();

    SurfaceNormalCalculator (const SurfaceNormalCalculator &) ;
    ~SurfaceNormalCalculator () {}

    // returns true if mapping is affine
    inline bool affine () const { return _affine ; }

    // calcuates normal
    void normal(const coord2_t&, coord3_t&) const ;
    void normal(const alu3d_ctype, const alu3d_ctype, coord3_t&) const;

    void negativeNormal(const coord2_t&, coord3_t&) const ;
    void negativeNormal(const alu3d_ctype, const alu3d_ctype, coord3_t&) const;

  public:
    // builds _b and _n, called from the constructors
    // public because also used in faceutility
    template <class vector_t>
    void buildMapping (const vector_t & , const vector_t & ,
                       const vector_t & , const vector_t & );
  protected:
    // builds _b and _n, called from the constructors
    // public because also used in faceutility
    template <class vector_t>
    void buildMapping (const vector_t & , const vector_t & ,
                       const vector_t & , const vector_t & ,
                       alu3d_ctype (&_b)[4][3] );
  } ;


  //! A bilinear surface mapping
  // NOTE: this class is different to the BilinearSurfaceMapping in
  // ALUGrid, for example the reference elements differ
  // here we have [0,1]^2 and in ALUGrid its [-1,1]^2
  // also the point numbering is different
  class BilinearSurfaceMapping : public SurfaceNormalCalculator
  {
  protected:
    typedef SurfaceNormalCalculator BaseType;

    using BaseType :: _n;
    static const double _epsilon;

    // our coordinate types
    typedef FieldVector<alu3d_ctype, 3> coord3_t;
    typedef FieldVector<alu3d_ctype, 2> coord2_t;

    // type of coordinate vectors from elements
    typedef alu3d_ctype double3_t[3];

    // type for helper matrices
    typedef FieldMatrix<alu3d_ctype,3,3> mat3_t;

    // type for inverse matrices
    typedef FieldMatrix<alu3d_ctype,2,3> matrix_t;

    // type for inverse matrices
    typedef FieldMatrix<alu3d_ctype,3,2> inv_t;

    alu3d_ctype _b [4][3] ;

    mutable mat3_t Df,Dfi;
    mutable inv_t invTransposed_;
    mutable matrix_t matrix_;
    mutable alu3d_ctype DetDf;

    mutable coord3_t normal_;
    mutable coord3_t tmp_;

    mutable bool _calcedInv;
    mutable bool _calcedTransposed;
    mutable bool _calcedMatrix;

  public:
    //! Constructor creating empty mapping with double , i.e. zero
    BilinearSurfaceMapping ();

    //! Constructor getting FieldVectors
    BilinearSurfaceMapping (const coord3_t&, const coord3_t&,
                            const coord3_t&, const coord3_t&) ;
    //! Constructor for double[3]
    BilinearSurfaceMapping (const double3_t &, const double3_t &,
                            const double3_t &, const double3_t &) ;
    BilinearSurfaceMapping (const BilinearSurfaceMapping &) ;
    ~BilinearSurfaceMapping () {}

    void inverse (const coord3_t&) const;
    const inv_t& jacobianInverseTransposed(const coord2_t&) const ;

    const matrix_t& jacobianTransposed(const coord2_t&) const ;

    // calculates determinant of face mapping using the normal
    alu3d_ctype det(const coord2_t&) const;

    // maps from local coordinates to global coordinates
    void world2map(const coord3_t &, coord2_t & ) const;

    // maps form global coordinates to local (within reference element)
    // coordinates
    void map2world(const coord2_t&, coord3_t&) const ;
    void map2world(const alu3d_ctype ,const alu3d_ctype , coord3_t&) const ;

  private:
    void map2worldnormal(const alu3d_ctype, const alu3d_ctype, const alu3d_ctype , coord3_t&) const;
    void map2worldlinear(const alu3d_ctype, const alu3d_ctype, const alu3d_ctype ) const;

  public:
    // builds _b and _n, called from the constructors
    // public because also used in faceutility
    template <class vector_t>
    void buildMapping (const vector_t & , const vector_t & ,
                       const vector_t & , const vector_t & );
  } ;



  //! A bilinear mapping
  template< int cdim >
  class BilinearMapping
  {
  public:
    typedef alu3d_ctype ctype;

    typedef FieldVector< ctype, cdim > world_t;
    typedef FieldVector< ctype, 2 > map_t;

    typedef FieldMatrix< ctype, 2, cdim > matrix_t;
    typedef FieldMatrix< ctype, cdim, 2 > inv_t;

  protected:
    ctype _b [4][cdim];

    mutable ctype det_;
    mutable matrix_t matrix_;
    mutable inv_t invTransposed_;

    mutable bool affine_;
    mutable bool calcedMatrix_;
    mutable bool calcedDet_;
    mutable bool calcedInv_;

  public:
    BilinearMapping ();
    BilinearMapping ( const world_t &p0, const world_t &p1,
                      const world_t &p2, const world_t &p3 );
    BilinearMapping ( const ctype (&p0)[ cdim ], const ctype (&p1)[ cdim ],
                      const ctype (&p2)[ cdim ], const ctype (&p3)[ cdim ] );

    bool affine () const;

    void world2map ( const world_t &, map_t & ) const;
    void map2world ( const ctype x, const ctype y, world_t &w ) const;
    void map2world ( const map_t &, world_t & ) const;

    ctype det ( const map_t & ) const;

    const matrix_t &jacobianTransposed ( const map_t & ) const;
    const inv_t &jacobianInverseTransposed ( const map_t & ) const;

    template< class vector_t >
    void buildMapping ( const vector_t &, const vector_t &,
                        const vector_t &, const vector_t & );

  protected:
    static void multTransposedMatrix ( const matrix_t &, FieldMatrix< ctype, 2, 2 > & );
    static void multMatrix ( const matrix_t &, const FieldMatrix< ctype, 2, 2 > &, inv_t & );

    void map2worldlinear ( const ctype, const ctype ) const;
    void inverse ( const map_t & ) const;
  };



  //! A linear mapping
  template< int cdim, int mydim >
  class LinearMapping
  {
  public:
    typedef alu3d_ctype ctype;

    typedef ctype double_t[ cdim ];

    typedef FieldVector< ctype, cdim > world_t;
    typedef FieldVector< ctype, mydim > map_t;

    typedef FieldMatrix< ctype, mydim, cdim > matrix_t;
    typedef FieldMatrix< ctype, cdim, mydim > inv_t;

  protected:
    matrix_t _matrix;             //!< transformation matrix (transposed)
    mutable inv_t _invTransposed; //!< storage for inverse of jacobian (transposed)
    world_t _p0;                  //! P[0]

    //! stores the determinant of the inverse
    mutable ctype _det;

    //! true if inverse has been calculated
    mutable bool _calcedInv;

    //! true if determinant has been calculated
    mutable bool _calcedDet;

  public:
    //! Constructor creating empty mapping with double , i.e. zero
    LinearMapping ();

    //! copy constructor
    LinearMapping (const LinearMapping &) ;

    // returns true if mapping is affine (which is always true)
    inline bool affine () const { return true ; }

    // return reference to transposed jacobian
    const matrix_t& jacobianTransposed(const map_t &) const ;

    // return reference to transposed jacobian inverse
    const inv_t& jacobianInverseTransposed(const map_t &) const ;

    // calculates determinant of mapping
    ctype det(const map_t&) const;

    // maps from local coordinates to global coordinates
    void world2map(const world_t &, map_t &) const;

    // maps form global coordinates to local (within reference element)
    // coordinates
    void map2world(const map_t &, world_t &) const ;

  protected:
    // calculate inverse
    void inverse (const map_t&) const;

    // calculate inverse one codim one entity
    void inverseCodimOne (const map_t&) const;

    // calculate determinant
    void calculateDeterminant (const map_t&) const;

    void multTransposedMatrix(const matrix_t& matrix,
                              FieldMatrix<ctype, mydim, mydim>& result) const;

    void multMatrix ( const matrix_t& A,
                      const FieldMatrix< ctype, mydim, mydim> &B,
                      inv_t& ret ) const ;

  public:
    // builds _b and _n, called from the constructors
    // public because also used in faceutility
    template <class vector_t>
    void buildMapping (const vector_t & , const vector_t & ,
                       const vector_t & , const vector_t & );

    // builds _b and _n, called from the constructors
    // public because also used in faceutility
    template <class vector_t>
    void buildMapping (const vector_t & , const vector_t & ,
                       const vector_t & );

    // builds _b and _n, called from the constructors
    // public because also used in faceutility
    template <class vector_t>
    void buildMapping (const vector_t & , const vector_t & );

    template <class vector_t>
    void buildMapping (const vector_t & );
  } ;


  ///////////////////////////////////////////////////////////////////
  //
  // NonConforming Mappings
  //
  ///////////////////////////////////////////////////////////////////


  //! General form of non-conforming face mapping
  //! This class is empty and needs to be specialised
  template< ALU3dGridElementType type, class Comm >
  class NonConformingFaceMapping;

  //! Non-conforming face mappings for tetrahedra
  template< class Comm >
  struct NonConformingFaceMapping< tetra, Comm >
  {
    typedef FieldVector< alu3d_ctype, 3 > CoordinateType;
    typedef typename ALU3dImplTraits< tetra, Comm >::HfaceRuleType RefinementRuleType;

    NonConformingFaceMapping ( RefinementRuleType rule, int nChild )
      : rule_( rule ), nChild_( nChild )
    {}

    void child2parent ( const CoordinateType &childCoordinates,
                        CoordinateType &parentCoordinates) const;

    CoordinateType child2parent ( const FieldVector< alu3d_ctype, 2 > &childCoordinates ) const;

  private:
    void child2parentNosplit(const CoordinateType& childCoordinates,
                             CoordinateType& parentCoordinates) const;
    void child2parentE01(const CoordinateType& childCoordinates,
                         CoordinateType& parentCoordinates) const;
    void child2parentE12(const CoordinateType& childCoordinates,
                         CoordinateType& parentCoordinates) const;
    void child2parentE20(const CoordinateType& childCoordinates,
                         CoordinateType& parentCoordinates) const;
    void child2parentIso4(const CoordinateType& childCoordinates,
                          CoordinateType& parentCoordinates) const;

    RefinementRuleType rule_;
    int nChild_;
  };

  //! Non-conforming face mappings for hexahedra
  template< class Comm >
  struct NonConformingFaceMapping< hexa, Comm >
  {
    typedef FieldVector< alu3d_ctype, 2 > CoordinateType;
    typedef typename ALU3dImplTraits< hexa, Comm >::HfaceRuleType RefinementRuleType;

    NonConformingFaceMapping ( RefinementRuleType rule, int nChild )
      : rule_( rule ), nChild_( nChild )
    {}

    void child2parent ( const CoordinateType &childCoordinates,
                        CoordinateType &parentCoordinates) const;

    CoordinateType child2parent ( const FieldVector< alu3d_ctype, 2 > &childCoordinates ) const;

  private:
    void child2parentNosplit(const CoordinateType& childCoordinates,
                             CoordinateType& parentCoordinates) const;
    void child2parentIso4(const CoordinateType& childCoordinates,
                          CoordinateType& parentCoordinates) const;

    RefinementRuleType rule_;
    int nChild_;
  };

} // end namespace Dune

#if COMPILE_ALUGRID_INLINE
  #include "mappings_imp.cc"
#endif
#endif
