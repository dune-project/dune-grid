if(HAVE_STD_THREAD AND ALUGRID_FOUND)
  set(TESTS ${TESTS} finitevolume_alusimplex)
  add_executable(finitevolume_alusimplex finitevolume_alusimplex.cc)
  add_dune_alugrid_flags(finitevolume_alusimplex)
endif(HAVE_STD_THREAD AND ALUGRID_FOUND)
set(CLEANFILES ${CLEANFILES}
  alusimplex_*_concentration-*.vtu
  alusimplex_*_concentration.series)

if(HAVE_STD_THREAD)
  set(TESTS ${TESTS} finitevolume_yasp)
  add_executable(finitevolume_yasp finitevolume_yasp.cc)
  add_dune_mpi_flags(finitevolume_yasp)
endif(HAVE_STD_THREAD)
set(CLEANFILES ${CLEANFILES}
  yasp_*_concentration-*.vtu
  yasp_*_concentration.series)

foreach(_T ${TESTS})
  target_link_libraries(${_T} dunegrid ${DUNE_LIBS})
  add_test(${_T} ${CMAKE_CURRENT_SOURCE_DIR}/${_T}.sh)
endforeach(_T ${TESTS})

# We do not want want to build the tests during make all,
# but just build them on demand and only if we've
# actually created tests here
if(${TESTS})
  add_directory_test_target(_test_target)
  add_dependencies(${_test_target} ${TESTS})
endif(${TESTS})

file(GLOB CLEANFILES ${CLEANFILES})
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${CLEANFILES}")
