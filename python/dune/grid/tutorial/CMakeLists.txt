# make sure the example works and run in during tests
dune_python_add_test(NAME pyexample
                     SCRIPT example.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     LABELS quick)
