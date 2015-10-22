dune_add_test(NAME test-dgf-yasp
              SOURCES test-dgf.cc
              COMPILE_DEFINITIONS GRIDDIM=3
                                  YASPGRID
                                  HAVE_DUNE_GRID=1
                                  DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\")

dune_add_test(NAME test-dgf-yasp-offset
              SOURCES test-dgf.cc
              COMPILE_DEFINITIONS GRIDDIM=3
                                  YASPGRID
                                  HAVE_DUNE_GRID=1
                                  DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\"
                                  COMPLETE_GRID_TYPE=YaspGrid<2,EquidistantOffsetCoordinates<double,2>>
              CMD_ARGS ${PROJECT_SOURCE_DIR}/doc/grids/dgf/test2d_offset.dgf)

dune_add_test(NAME test-dgf-oned
              SOURCES test-dgf.cc
              COMPILE_DEFINITIONS GRIDDIM=1
                                  ONEDGRID
                                  HAVE_DUNE_GRID=1
                                  DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\")

if(ALUGRID_FOUND)
  dune_add_test(NAME test-dgf-alu
                SOURCES test-dgf.cc
                COMPILE_DEFINITIONS ALUGRID_CUBE
                                    GRIDDIM=3
                                    HAVE_DUNE_GRID=1)
endif(ALUGRID_FOUND)

if(ALBERTA_FOUND)
  add_executable(test-dgf-alberta test-dgf.cc)
  add_dune_alberta_flags(test-dgf-alberta GRIDDIM 2)
  target_compile_definitions(test-dgf-alberta PUBLIC
    ALBERTAGRID
    GRIDDIM=2
    HAVE_DUNE_GRID=1
    DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\")
  dune_add_test(TARGET test-dgf-alberta)
endif(ALBERTA_FOUND)

if(UG_FOUND)
  dune_add_test(NAME test-dgf-ug
                SOURCES test-dgf.cc
                COMPILE_DEFINITIONS UGGRID
                                    GRIDDIM=3
                                    HAVE_DUNE_GRID=1
                                    DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\")
endif(UG_FOUND)
