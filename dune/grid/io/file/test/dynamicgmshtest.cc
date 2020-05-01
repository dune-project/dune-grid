// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include "config.h"

#include <exception>
#include <iostream>
#include <string>
#include <utility>

#include <dune/common/classname.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/io/file/gmshreader.hh>

#if GMSH_ALBERTAGRID
#include <dune/grid/albertagrid.hh>
#ifndef ALBERTA_DIM
#define ALBERTA_DIM 2
#endif
#ifndef GRIDDIM
#define GRIDDIM ALBERTA_DIM
#endif
#endif

#if GMSH_ONEDGRID
#include <dune/grid/onedgrid.hh>
#endif

#if GMSH_UGGRID
#include <dune/grid/uggrid.hh>
#endif


template<class T>
T &discarded(T &&v) { return v; }

void check_set(Dune::TestSuite &suite, bool result, const char *test,
               const char *file, unsigned line)
{
  suite.check(result, test) << file << ':' << line << ": Here";
}

void check_exception(Dune::TestSuite &suite, std::exception_ptr exception,
                     const char *test, const char *file, unsigned line)
{
  try {
    std::rethrow_exception(exception);
  }
  catch(const std::exception &e) {
    suite.check(false, test) << file << ':' << line << ": Exception: "
                             << Dune::className(e) << ": " << e.what();
  }
  catch(...) {
    suite.check(false, test) << file << ':' << line
                             << ": Exception: (unknown type)";
  }
}

#define CHECK(suite, ...)                                               \
  [&]{                                                                  \
    bool check_result;                                                  \
    try {                                                               \
      check_result = (__VA_ARGS__);                                     \
    }                                                                   \
    catch(...) {                                                        \
      check_exception((suite), std::current_exception(),                \
                      #__VA_ARGS__, __FILE__, __LINE__);                \
      return;                                                           \
    }                                                                   \
    check_set((suite), check_result, #__VA_ARGS__, __FILE__, __LINE__); \
  }()

template<class E, class F, class... Args>
bool throws(F &&f, Args &&... args)
{
  try {
    std::forward<F>(f)(std::forward<Args>(args)...);
    return false;
  }
  catch(const E &) {
    return true;
  }
  catch(...) {
    return false;
  }
}

template<class Factory, class Grid>
bool hasInsertedBoundarySegments(const Factory &factory, const Grid &grid)
{
  auto gv = leafGridView(grid);
  for(auto &&e : elements(gv))
    for(auto &&is : intersections(gv, e))
      if(is.boundary() && factory.wasInserted(is))
        return true;
  return false;
}

// test DynamicGmshReader Methods
template<class Factory>
struct DynamicReaderTest {
  std::string inputName;
  bool expectNoBoundary = false;

private:
  template<class Setter, class CheckSet, class CheckUnset>
  static void testSetFlagHelper(Setter setter, CheckSet checkSet,
                                CheckUnset checkUnset)
  {
    Dune::DynamicGmshReader reader{};
    checkUnset(reader);

    setter(reader);
    checkSet(reader);

    setter(reader, true);
    checkSet(reader);

    setter(reader, false);
    checkUnset(reader);
  }

  template<class Setter, class Getter>
  void testDataHelper(Dune::TestSuite &suite, Setter setter, Getter getter,
                      bool expectMissingData)
  {
    Dune::DynamicGmshReader reader{};
    const auto &creader = reader;
    CHECK(suite, throws<Dune::DynamicGmshReader::StateError>(getter, reader));
    CHECK(suite, throws<Dune::DynamicGmshReader::StateError>(getter, creader));

    reader.read(discarded(Factory{}), inputName);
    CHECK(suite, throws<Dune::DynamicGmshReader::StateError>(getter, reader));
    CHECK(suite, throws<Dune::DynamicGmshReader::StateError>(getter, creader));

    setter(reader);
    reader.read(discarded(Factory{}), inputName);
    CHECK(suite, getter(reader).size() > 0 || expectMissingData);
    CHECK(suite, getter(creader).size() > 0 || expectMissingData);

    setter(reader);
    CHECK(suite, throws<Dune::DynamicGmshReader::StateError>(getter, reader));
    CHECK(suite, throws<Dune::DynamicGmshReader::StateError>(getter, creader));
  }

  template<class Setter, class Getter>
  void testDataOrEmptyHelper(Dune::TestSuite &suite, Setter setter,
                             Getter getter, bool expectMissingData)
  {
    Dune::DynamicGmshReader reader{};
    const auto &creader = reader;
    CHECK(suite, throws<Dune::DynamicGmshReader::StateError>(getter, reader));
    CHECK(suite, throws<Dune::DynamicGmshReader::StateError>(getter, creader));

    reader.read(discarded(Factory{}), inputName);
    CHECK(suite, getter(reader).size() == 0);
    CHECK(suite, getter(creader).size() == 0);

    setter(reader);
    reader.read(discarded(Factory{}), inputName);
    CHECK(suite, getter(reader).size() > 0 || expectMissingData);
    CHECK(suite, getter(creader).size() > 0 || expectMissingData);

    setter(reader);
    CHECK(suite, throws<Dune::DynamicGmshReader::StateError>(getter, reader));
    CHECK(suite, throws<Dune::DynamicGmshReader::StateError>(getter, creader));
  }

public:
  // test read()
  void testRead(Dune::TestSuite &suite) const
  {
    Dune::TestSuite test("read()");
    Dune::DynamicGmshReader reader{};
    Factory factory;
    reader.read(factory, inputName);
    auto gridp = factory.createGrid();
    CHECK(test, leafGridView(*gridp).size(0) > 0);
    suite.subTest(test);
  }

  // test setBoundaryData()
  void testSetBoundaryData(Dune::TestSuite &suite)
  {
    Dune::TestSuite test("setBoundaryData()");
    auto setter = [](auto &reader, auto... args) {
      reader.setReadBoundaryData(args...);
    };
    auto checkSet = [&](auto reader) {
      Factory factory;
      reader.read(factory, inputName);
      auto gridp = factory.createGrid();
      CHECK(test,
            expectNoBoundary || hasInsertedBoundarySegments(factory, *gridp));
      CHECK(test,
            reader.boundaryDataOrEmpty().size() > 0 || expectNoBoundary);
    };
    auto checkUnset = [&](auto reader) {
      Factory factory;
      reader.read(factory, inputName);
      auto gridp = factory.createGrid();
      CHECK(test, !hasInsertedBoundarySegments(factory, *gridp));
      CHECK(test, reader.boundaryDataOrEmpty().size() == 0);
    };

    testSetFlagHelper(setter, checkSet, checkUnset);

    suite.subTest(test);
  }

  // test setElementData()
  void testSetElementData(Dune::TestSuite &suite)
  {
    Dune::TestSuite test("setElementData()");
    auto setter = [](auto &reader, auto... args) {
      reader.setReadElementData(args...);
    };
    auto checkSet = [&](auto reader) {
      reader.read(discarded(Factory{}), inputName);
      CHECK(test, reader.elementDataOrEmpty().size() > 0);
    };
    auto checkUnset = [&](auto reader) {
      reader.read(discarded(Factory{}), inputName);
      CHECK(test, reader.elementDataOrEmpty().size() == 0);
    };

    testSetFlagHelper(setter, checkSet, checkUnset);

    suite.subTest(test);
  }

  // test setForceInsertBoundarySegments()
  void testSetForceInsertBoundarySegments(Dune::TestSuite &suite)
  {
    Dune::TestSuite test("setForceInsertBoundarySegments()");
    auto setter = [](auto &reader, auto... args) {
      reader.setForceInsertBoundarySegments(args...);
    };
    auto checkSet = [&](auto reader) {
      Factory factory;
      reader.read(factory, inputName);
      auto gridp = factory.createGrid();
      CHECK(test,
            expectNoBoundary || hasInsertedBoundarySegments(factory, *gridp));
      CHECK(test, reader.boundaryDataOrEmpty().size() == 0);
    };
    auto checkUnset = [&](auto reader) {
      Factory factory;
      reader.read(factory, inputName);
      auto gridp = factory.createGrid();
      CHECK(test, !hasInsertedBoundarySegments(factory, *gridp));
      CHECK(test, reader.boundaryDataOrEmpty().size() == 0);
    };

    testSetFlagHelper(setter, checkSet, checkUnset);

    suite.subTest(test);
  }

  // test setVerbose()
  void testSetVerbose()
  {
    // not much possibility to test beyond that we can call it
    Dune::DynamicGmshReader reader{};
    reader.setVerbose();
    reader.setVerbose(true);
    reader.setVerbose(false);
  }

  // test boundaryData()
  void testBoundaryData(Dune::TestSuite &suite)
  {
    auto setter = [](auto &reader) { reader.setReadBoundaryData(); };
    auto getter = [](auto &reader) -> decltype(auto) {
      return reader.boundaryData();
    };
    auto defaultGetter = [](auto &reader) -> decltype(auto) {
      return reader.boundaryDataOrEmpty();
    };

    {
      Dune::TestSuite test("boundaryData()");
      testDataHelper(test, setter, getter, expectNoBoundary);
      suite.subTest(test);
    }
    {
      Dune::TestSuite test("boundaryDataOrEmpty()");
      testDataOrEmptyHelper(test, setter, defaultGetter, expectNoBoundary);
      suite.subTest(test);
    }
  }

  // test elementData()
  void testElementData(Dune::TestSuite &suite)
  {
    auto setter = [](auto &reader) { reader.setReadElementData(); };
    auto getter = [](auto &reader) -> decltype(auto) {
      return reader.elementData();
    };
    auto defaultGetter = [](auto &reader) -> decltype(auto) {
      return reader.elementDataOrEmpty();
    };

    {
      Dune::TestSuite test("boundaryData()");
      testDataHelper(test, setter, getter, true);
      suite.subTest(test);
    }
    {
      Dune::TestSuite test("boundaryDataOrEmpty()");
      testDataOrEmptyHelper(test, setter, defaultGetter, true);
      suite.subTest(test);
    }
  }

  // test reading with dynamic interface
  void test(Dune::TestSuite &suite)
  {
    Dune::TestSuite test(inputName);

    std::cout << "testing input " << inputName << std::endl;
    testRead(test);
    testSetBoundaryData(test);
    testSetElementData(test);
    testSetForceInsertBoundarySegments(test);
    testSetVerbose();
    testBoundaryData(test);
    testElementData(test);

    suite.subTest(test);
  }
};

template <typename Grid>
void testReader(Dune::TestSuite &suite, const std::string& path,
                const std::string& gridName, bool expectNoBoundary = false)
{
  DynamicReaderTest<Dune::GridFactory<Grid>>
    test{ path+gridName+".msh", expectNoBoundary };
  test.test(suite);
}


int main( int argc, char** argv )
{
  Dune::MPIHelper::instance( argc, argv );
  std::string path(std::string(DUNE_GRID_EXAMPLE_GRIDS_PATH)+"gmsh/");
  Dune::TestSuite suite;

#if GMSH_UGGRID
  {
    Dune::TestSuite test("UGGRID-2D");
    testReader<Dune::UGGrid<2> >( test, path, "curved2d" );
    testReader<Dune::UGGrid<2> >( test, path, "circle2ndorder" );
    testReader<Dune::UGGrid<2> >( test, path, "unitsquare_quads_2x2" );
    testReader<Dune::UGGrid<2> >( test, path, "hybrid-testgrid-2d", true );
    suite.subTest(test);
  }
  {
    Dune::TestSuite test("UGGrid-3D");
    testReader<Dune::UGGrid<3> >( test, path, "pyramid" );
    testReader<Dune::UGGrid<3> >( test, path, "pyramid2ndorder" );
    testReader<Dune::UGGrid<3> >( test, path, "hybrid-testgrid-3d", true );
    suite.subTest(test);
  }
#endif

#if GMSH_ALBERTAGRID
#if ALBERTA_DIM==2
  {
    Dune::TestSuite test("AlbertaGrid-2D");
    testReader<Dune::AlbertaGrid<2> >( test, path, "curved2d" );
    suite.subTest(test);
  }
#endif
#if ALBERTA_DIM==3
  {
    Dune::TestSuite test("AlbertaGrid-2D");
    testReader<Dune::AlbertaGrid<2> >( test, path, "sphere", true );
    suite.subTest(test);
  }
  {
    Dune::TestSuite test("AlbertaGrid-3D");
    testReader<Dune::AlbertaGrid<3> >( test, path, "pyramid" );
    suite.subTest(test);
  }
#endif
#endif

#if GMSH_ONEDGRID
  {
    Dune::TestSuite test("OneDGrid");
    testReader<Dune::OneDGrid>( test, path, "oned-testgrid" );
    suite.subTest(test);
  }
#endif

  return suite.exit();
}
