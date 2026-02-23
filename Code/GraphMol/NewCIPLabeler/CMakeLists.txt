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

# PDB test executable
if(RDK_BUILD_CPP_TESTS)
    add_executable(testNewCIPLabeler_pdb test_pdb.cpp)
    target_link_libraries(testNewCIPLabeler_pdb PRIVATE
                          NewCIPLabeler
                          SmilesParse
                          GraphMol
                          FileParsers)
    target_compile_definitions(testNewCIPLabeler_pdb PRIVATE RDKIT_NEWCIPLABELER_BUILD)
endif()

# Benchmark executable
if(RDK_BUILD_CPP_TESTS)
    add_executable(benchmarkNewCIPLabeler quick_benchmark.cpp)
    target_link_libraries(benchmarkNewCIPLabeler PRIVATE
                          NewCIPLabeler
                          SmilesParse
                          GraphMol
                          FileParsers)
    target_compile_definitions(benchmarkNewCIPLabeler PRIVATE RDKIT_NEWCIPLABELER_BUILD)
endif()
