# The CLI contract gives eta the final normal-power table unless explicitly
# overridden. A non-default m=0 power must propagate through parsing.
execute_process(
    COMMAND "${EXECUTABLE}" "${EXPORT_FILE}" 3 64 8 1e-8 1
        --stored-powers=0.5,0.5,0 --count-shifts=0 0,1 1,1 2,1
    RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "power-table diagnostic failed: ${error}")
endif()
string(REGEX MATCH
    "eta_stored_power_index,m,n,power[\n\r]+1,0,1, +5\\.0+[eE]-01"
    eta_default "${output}")
if(NOT eta_default)
    message(FATAL_ERROR "eta did not inherit the requested normal power 0.5")
endif()
