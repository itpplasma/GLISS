if(NOT DEFINED COMMAND_PATH OR NOT DEFINED EXPECTED)
    message(FATAL_ERROR "COMMAND_PATH and EXPECTED are required")
endif()

set(arguments)
foreach(index RANGE 1 3)
    if(DEFINED ARG${index})
        list(APPEND arguments "${ARG${index}}")
    endif()
endforeach()

execute_process(
    COMMAND "${COMMAND_PATH}" ${arguments}
    RESULT_VARIABLE result
    OUTPUT_VARIABLE output
    ERROR_VARIABLE error)

if(result EQUAL 0)
    message(FATAL_ERROR "command unexpectedly succeeded")
endif()

string(FIND "${output}${error}" "${EXPECTED}" match)
if(match EQUAL -1)
    message(FATAL_ERROR "expected diagnostic was not emitted")
endif()
