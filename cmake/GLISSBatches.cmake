# Shared pointwise energy batch products; ordinary builds do not need a CAS.
target_sources(gliss PRIVATE src/generated/ideal_mhd_batches.f90)
if(BUILD_TESTING)
    add_executable(test_ideal_mhd_batches test/test_ideal_mhd_batches.f90)
    target_link_libraries(test_ideal_mhd_batches PRIVATE gliss)
    add_test(NAME ideal_mhd_batch_oracles COMMAND test_ideal_mhd_batches)
endif()
if(GLISS_ENABLE_CODEGEN)
    add_executable(gliss_batch_gen gen/ideal_mhd/gen_ideal_mhd_batches.f90
        gen/ideal_mhd/ideal_mhd_symbolic.f90)
    set_target_properties(gliss_batch_gen PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/batch-codegen-modules)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/batch-codegen-modules)
    target_link_libraries(gliss_batch_gen PRIVATE fortsym::fortsym)
    file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/gen/fortsym-revision.txt" batch_revision)
    add_test(NAME generated_batch_reproduction COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gliss_batch_gen>
        -DREVISION=${batch_revision}
        -DWORK=${CMAKE_CURRENT_BINARY_DIR}/batch-regeneration
        -DSOURCE=${CMAKE_CURRENT_SOURCE_DIR}/src/generated
        -P ${CMAKE_CURRENT_SOURCE_DIR}/gen/ideal_mhd/check_batches.cmake)
    set_tests_properties(generated_batch_reproduction PROPERTIES LABELS codegen)
endif()
