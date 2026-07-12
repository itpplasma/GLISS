execute_process(
    COMMAND "${PROGRAM}" ${ARGUMENTS}
    RESULT_VARIABLE status
    OUTPUT_VARIABLE standard_output
    ERROR_VARIABLE standard_error)

if(status EQUAL 0)
    message(FATAL_ERROR "invalid input was accepted")
endif()

set(combined_output "${standard_output}${standard_error}")
string(FIND "${combined_output}" "${EXPECTED}" expected_position)
if(expected_position EQUAL -1)
    message(FATAL_ERROR
        "missing diagnostic '${EXPECTED}' in output:\n${combined_output}")
endif()

string(FIND "${combined_output}" "Backtrace" backtrace_position)
if(NOT backtrace_position EQUAL -1)
    message(FATAL_ERROR "input error emitted a backtrace:\n${combined_output}")
endif()
