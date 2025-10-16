// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GRID_IO_FILE_TEST_CHECKVTKFILE_HH
#define DUNE_GRID_IO_FILE_TEST_CHECKVTKFILE_HH

#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <ostream>
#include <set>
#include <sstream>
#include <string>

#include <dune/common/exceptions.hh>
#include <dune/common/test/testsuite.hh>

#ifdef PYTHON_INTERPRETER
#include <Python.h>
#endif

namespace Dune {

namespace Impl {

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

inline int runPython(const std::string &code)
{
#ifdef PYTHON_INTERPRETER
#if (PY_MAJOR_VERSION >= 3) && (PY_MINOR_VERSION >= 8)
  PyConfig config;
  PyConfig_InitPythonConfig(&config);
  PyConfig_SetString(&config, &config.program_name, PYTHON_INTERPRETER);
  PyStatus status = Py_InitializeFromConfig(&config);
  PyConfig_Clear(&config);
  if (PyStatus_Exception(status)) {
    Py_ExitStatusException(status);
    return PyStatus_IsExit(status) ? status.exitcode : EXIT_FAILURE;
  }
#else
  Py_SetProgramName(PYTHON_INTERPRETER);
  Py_Initialize();
#endif

  int exitcode = PyRun_SimpleString(code.c_str());
  int finalizecode = Py_FinalizeEx();
  return (exitcode + finalizecode) ? EXIT_FAILURE : EXIT_SUCCESS;
#else
  return EXIT_FAILURE;
#endif // PYTHON_INTERPRETER
}

inline bool is_suffix(const std::string &s, const std::string &suffix)
{
  return s.size() >= suffix.size() &&
    s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

inline bool havePythonVTK()
{
  static const bool result = [] {
    // This check is invoked only once, even in a multithreading environment,
    // since it is invoked in the initializer of a static variable.
    if(runPython("from vtk import *") == EXIT_SUCCESS)
      return true;
    std::cerr << "warning: python or python vtk module not available.  This "
              << "will result in skipped tests, since we cannot check that "
              << "vtk can read the files we wrote." << std::endl;
    return false;
  } ();

  return result;
}

inline std::string pythonVTKReader(const std::string& filename)
{
  if     (is_suffix(filename, ".vtu"))  return "vtkXMLUnstructuredGridReader";
  else if(is_suffix(filename, ".pvtu")) return "vtkXMLPUnstructuredGridReader";
  else if(is_suffix(filename, ".vtp"))  return "vtkXMLPolyDataReader";
  else if(is_suffix(filename, ".pvtp")) return "vtkXMLPPolyDataReader";
  else DUNE_THROW(Dune::NotImplemented,
                  "Unknown vtk file extension: " << filename);
}

} /* namespace Impl */

class VTKChecker
{
public:
  void push(const std::string& file)
    {
      auto res = files_.insert(file);
      if (not res.second) {
        testSuite_.check(false, "VTKChecker")
          << "'" << file << "' was added multiple times";
      }
    }

  int check()
    {
      if (not Impl::havePythonVTK()) {
        return 77;
      }
      else if (not files_.empty()) {
        const int result = Impl::runPython(generatePythonCode());
        testSuite_.check(result == EXIT_SUCCESS);
      }
      return testSuite_.exit();
    }

  const TestSuite& testSuite() const
    {
      return testSuite_;
    }

private:
  std::string generatePythonCode() const
    {
      std::stringstream code;

      code << "from vtk import *\n"
           << "import sys\n"
           << "passed = True\n";

      for (const auto& file : files_) {
        code << "reader = " << Impl::pythonVTKReader(file) << "()\n"
             << "reader.SetFileName('" << Impl::pyq(file) << "')\n"
             << "reader.Update()\n"
             << "if (not (reader.GetOutput().GetNumberOfCells() > 0)):\n"
             << "    print('ERROR in {}'.format('" << Impl::pyq(file) << "'))\n"
             << "    passed = False\n";
      }

      code << "sys.exit(0 if passed else 1)\n";

      return code.str();
    }

  std::set< std::string > files_;
  TestSuite testSuite_;
};

} /* namespace Dune */

#endif // DUNE_GRID_IO_FILE_TEST_CHECKVTKFILE_HH
