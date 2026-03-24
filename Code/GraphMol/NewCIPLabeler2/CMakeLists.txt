rdkit_library(NewCIPLabeler2
              NewCIPLabeler.cpp
              LINK_LIBRARIES GraphMol)
target_compile_definitions(NewCIPLabeler2 PRIVATE RDKIT_NEWCIPLABELER2_BUILD)

rdkit_headers(NewCIPLabeler.h
              DEST GraphMol/NewCIPLabeler2)

rdkit_catch_test(testNewCIPLabeler2 catch_tests.cpp
                 NewCIPLabeler.cpp
                 LINK_LIBRARIES SmilesParse GraphMol FileParsers)

if(RDK_BUILD_CPP_TESTS)
    target_compile_definitions(testNewCIPLabeler2 PRIVATE RDKIT_NEWCIPLABELER2_BUILD)
endif()
