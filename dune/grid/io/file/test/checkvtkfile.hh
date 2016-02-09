#ifndef DUNE_GRID_IO_FILE_TEST_CHECKVTKFILE_HH
#define DUNE_GRID_IO_FILE_TEST_CHECKVTKFILE_HH

#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <ostream>
#include <sstream>
#include <string>

#include <sys/wait.h>

#include <dune/common/exceptions.hh>

// quote so the result can be used inside '...' in python
// quotes not included in the result
inline std::string pyq(const std::string &s)
{
  std::ostringstream result;
  for(std::size_t i = 0; i < s.size(); ++i)
  {
    char c = s[i];
    switch(c) {
    case '\'': result << "\\'";  break;
    case '\\': result << "\\\\"; break;
    case '\n': result << "\\n";  break;
    default:
      if(c < 32 || c >= 127)
        result << "\\x" << std::hex << std::setfill('0') << std::setw(2)
               << static_cast<int>(c);
      else
        result << c;
    }
  }
  return result.str();
}

// quote so the result can be used inside '...' in the bourne shell
// quotes not included in the result
inline std::string shq(const std::string &s)
{
  std::ostringstream result;
  bool pend = false;
  for(std::size_t i = 0; i < s.size(); ++i)
  {
    char c = s[i];
    switch(c) {
    case '\0': DUNE_THROW(Dune::NotImplemented,
                          "Can't pass \\0 through the shell");
    case '\'': result << (pend ? "" : "'") << "\\'"; pend = true;  break;
    default:   result << (pend ? "'" : "") << c;     pend = false; break;
    }
  }
  if(pend) result << "'";
  return result.str();
}

inline int runShell(const std::string &code)
{
  int result = std::system(code.c_str());
  // Avoid returning anything that is a multiple of 256, unless the return
  // value was 0.  This way the return value can be directly used as an
  // argument to exit(), which usually interprets its argument modulo 256.
  if(WIFEXITED(result))
    return WEXITSTATUS(result);
  if(WIFSIGNALED(result))
    return WTERMSIG(result) + 256;
  else
    return 513;
}

inline int runPython(const std::string &code)
{
  return runShell("python -c '"+shq(code)+"'");
}

inline bool is_suffix(const std::string &s, const std::string &suffix)
{
  return s.size() >= suffix.size() &&
    s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

inline int checkVTKFile(const std::string &name)
{
  static const bool havePythonVTK = [] {
    // This check is invoked only once, even in a multithreading environment,
    // since it is invoked in the initializer of a static variable.
    if(runPython("from vtk import *") == 0)
      return true;
    std::cerr << "warning: python or python vtk module not available.  This "
              << "will result in skipped tests, since we cannot check that "
              << "vtk can read the files we wrote." << std::endl;
    return false;
  } ();
  if(!havePythonVTK)
  {
    std::cerr << "skip: " << name << std::endl;
    return 77;
  }

  std::string reader;
  if     (is_suffix(name, ".vtu"))  reader = "vtkXMLUnstructuredGridReader";
  else if(is_suffix(name, ".pvtu")) reader = "vtkXMLPUnstructuredGridReader";
  else if(is_suffix(name, ".vtp"))  reader = "vtkXMLPolyDataReader";
  else if(is_suffix(name, ".pvtp")) reader = "vtkXMLPPolyDataReader";
  else DUNE_THROW(Dune::NotImplemented,
                  "Unknown vtk file extension: " << name);

  std::cout << "Loading " << name << " using python vtk" << std::endl;
  std::string pycode =
    "from vtk import *;"
    "import sys;"
    "reader = "+reader+"();"
    "reader.SetFileName('"+pyq(name)+"');"
    "reader.Update();"
    // check that the number of of cells is > 0
    "sys.exit(not (reader.GetOutput().GetNumberOfCells() > 0));";
  int result = runPython(pycode);
  std::cout << (result == 0 ? "pass: " : "fail: ") << name << std::endl;
  return result;
}

#endif // DUNE_GRID_IO_FILE_TEST_CHECKVTKFILE_HH
