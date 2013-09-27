#Do not follow symlinks during FILE GLOB_RECURSE
cmake_policy(SET CMP0009 NEW)
file(GLOB_RECURSE makefiles RELATIVE ${RELPATH} "CMakeLists.txt")
foreach(_file ${makefiles})
  string(REGEX MATCH ".*/test$" _testdir ${_file})
  if(_testdir)
    list(APPEND _makefiles ${_testdir})
  endif(_testdir)
endforeach(_file ${makefiles})
if(_makefiles)
  message("${_makefiles}")
endif(_makefiles)
