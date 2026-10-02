# The oracle is DCON's independently qualified Solov'ev sign, not a GLISS
# eigenvalue reference. Both controls use the current conforming FEEC space.
execute_process(
    COMMAND "${EXECUTABLE}" "${EXPORT_FILE}" 3 64 8 1e-8 1
        --count-shifts=0 0,1 1,-1 1,1 2,-1 2,1 3,-1 3,1
        4,-1 4,1 5,-1 5,1 6,-1 6,1
    RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "count diagnostic failed: ${error}")
endif()
string(REGEX MATCH "[\n\r] *0\\.0+[eE]\\+0+,([0-9]+)" row "${output}")
if(NOT row)
    message(FATAL_ERROR "zero-shift count row was not emitted")
endif()
if(NOT CMAKE_MATCH_1 EQUAL EXPECTED_COUNT)
    message(FATAL_ERROR "DCON sign differs: expected ${EXPECTED_COUNT}, got ${CMAKE_MATCH_1}")
endif()
