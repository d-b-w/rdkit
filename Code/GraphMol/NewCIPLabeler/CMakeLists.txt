rdkit_library(NewCIPLabeler
              NewCIPLabeler.cpp
              Descriptor.cpp
              Substituent.cpp
              Rules.cpp
              Ranker.cpp
              Tetrahedral.cpp
              DoubleBond.cpp
              LINK_LIBRARIES GraphMol)
target_compile_definitions(NewCIPLabeler PRIVATE RDKIT_NEWCIPLABELER_BUILD)

rdkit_headers(NewCIPLabeler.h
              Descriptor.h
              DEST GraphMol/NewCIPLabeler)

rdkit_catch_test(testNewCIPLabeler catch_tests.cpp
                 NewCIPLabeler.cpp
                 Descriptor.cpp
                 Substituent.cpp
                 Rules.cpp
                 Ranker.cpp
                 Tetrahedral.cpp
                 DoubleBond.cpp
                 LINK_LIBRARIES SmilesParse GraphMol FileParsers)

# Set the BUILD macro for DLL linkage consistency
if(RDK_BUILD_CPP_TESTS)
    target_compile_definitions(testNewCIPLabeler PRIVATE RDKIT_NEWCIPLABELER_BUILD)
endif()
