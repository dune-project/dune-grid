dune_python_add_test(NAME pytest_gf1
                     COMMAND ${PYTHON_EXECUTABLE} test_gf1.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     LABELS quick)
dune_python_add_test(NAME pytest_gf2
                     COMMAND ${PYTHON_EXECUTABLE} test_gf2.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     LABELS quick)
dune_python_add_test(NAME pytest_indexset
                     COMMAND ${PYTHON_EXECUTABLE} test_indexset.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     LABELS quick)
dune_python_add_test(NAME pyinterpolate
                     COMMAND ${PYTHON_EXECUTABLE} interpolate.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     LABELS quick)
