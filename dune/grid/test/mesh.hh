#ifndef DUNE_GRID_TEST_MESH_HH
#define DUNE_GRID_TEST_MESH_HH

#include <array>
#include <initializer_list>
#include <utility>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>


namespace Dune
{

  template< int dim >
  struct BasicSimplexElement
  {
    BasicSimplexElement () {}
    BasicSimplexElement ( const BasicSimplexElement & ) = default;

    BasicSimplexElement ( std::initializer_list< unsigned int > const &l ) : vertices_( l ) {}


    GeometryType type () const { return GeometryType( Dune::GeometryType::simplex, dim ); }
    const std::vector< unsigned int > &vertices () const { return vertices_; }

  private:
    std::vector< unsigned int > vertices_;
  };



  struct Kuhn2dSimplexMesh
  {
    typedef Dune::FieldVector< double, 2 > Vertex;
    typedef std::vector< Vertex > Vertices;
    typedef BasicSimplexElement< 2 > Element;
    typedef std::vector< Element > Elements;

    Kuhn2dSimplexMesh ()
    {
      vertices_.push_back( {{ 0.0, 0.0 }} );
      vertices_.push_back( {{ 1.0, 0.0 }} );
      vertices_.push_back( {{ 0.0, 1.0 }} );
      vertices_.push_back( {{ 1.0, 1.0 }} );

      elements_.push_back( {{ 0, 1, 2 }} );
      elements_.push_back( {{ 3, 2, 0 }} );
    }

    const Vertices &vertices () const { return vertices_; }
    const Elements &elements () const { return elements_; }

  protected:
    Vertices vertices_;
    Elements elements_;
  };


  struct Kuhn3dSimplexMesh
  {
    typedef Dune::FieldVector< double, 3 > Vertex;
    typedef std::vector< Vertex > Vertices;
    typedef BasicSimplexElement< 3 > Element;
    typedef std::vector< Element > Elements;

    Kuhn3dSimplexMesh ()
    {
      vertices_.push_back( {{0.0, 0.0, 0.0}} );
      vertices_.push_back( {{0.0, 0.0, 1.0}} );
      vertices_.push_back( {{0.0, 1.0, 0.0}} );
      vertices_.push_back( {{0.0, 1.0, 1.0}} );
      vertices_.push_back( {{1.0, 0.0, 0.0}} );
      vertices_.push_back( {{1.0, 0.0, 1.0}} );
      vertices_.push_back( {{1.0, 1.0, 0.0}} );
      vertices_.push_back( {{1.0, 1.0, 1.0}} );

      elements_.push_back( {{ 0, 1, 3, 7 }} );
      elements_.push_back( {{ 0, 1, 5, 7 }} );
      elements_.push_back( {{ 0, 2, 3, 7 }} );
      elements_.push_back( {{ 0, 2, 6, 7 }} );
      elements_.push_back( {{ 0, 4, 6, 7 }} );
      elements_.push_back( {{ 0, 4, 5, 7 }} );
    }

    const Vertices &vertices () const { return vertices_; }
    const Elements &elements () const { return elements_; }

  protected:
    Vertices vertices_;
    Elements elements_;
  };

} // namespace Dune

#endif // #ifndef DUNE_GRID_TEST_MESH_HH
