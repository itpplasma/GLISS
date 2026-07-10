include_guard(GLOBAL)

if(NOT CMAKE_Fortran_COMPILER_ID MATCHES "LLVMFlang|Flang")
    message(FATAL_ERROR
        "GLISS_ENABLE_ENZYME requires flang-new, got ${CMAKE_Fortran_COMPILER_ID}")
endif()

string(REGEX MATCH "^[0-9]+" _enzyme_compiler_major
    "${CMAKE_Fortran_COMPILER_VERSION}")
if(NOT _enzyme_compiler_major)
    message(FATAL_ERROR "Cannot determine the Flang LLVM major version")
endif()

find_program(ENZYME_OPT_EXECUTABLE
    NAMES "opt-${_enzyme_compiler_major}" opt REQUIRED)
find_program(ENZYME_LLVM_LINK_EXECUTABLE
    NAMES "llvm-link-${_enzyme_compiler_major}" llvm-link REQUIRED)

execute_process(
    COMMAND "${ENZYME_OPT_EXECUTABLE}" --version
    OUTPUT_VARIABLE _enzyme_opt_version
    OUTPUT_STRIP_TRAILING_WHITESPACE)
string(REGEX MATCH "LLVM version ([0-9]+)" _enzyme_opt_match
    "${_enzyme_opt_version}")
if(NOT CMAKE_MATCH_1 STREQUAL _enzyme_compiler_major)
    message(FATAL_ERROR
        "Flang LLVM ${_enzyme_compiler_major} and opt LLVM ${CMAKE_MATCH_1} differ")
endif()

set(ENZYME_PLUGIN "" CACHE FILEPATH "Path to LLVMEnzyme plugin")
if(NOT ENZYME_PLUGIN AND DEFINED ENV{ENZYME_PLUGIN})
    set(ENZYME_PLUGIN "$ENV{ENZYME_PLUGIN}" CACHE FILEPATH
        "Path to LLVMEnzyme plugin" FORCE)
endif()
if(NOT ENZYME_PLUGIN)
    find_library(_enzyme_plugin
        NAMES "LLVMEnzyme-${_enzyme_compiler_major}" LLVMEnzyme
        PATHS ENV ENZYME_PLUGIN_DIR ENV LD_LIBRARY_PATH
        PATH_SUFFIXES lib enzyme Enzyme)
    if(_enzyme_plugin)
        set(ENZYME_PLUGIN "${_enzyme_plugin}" CACHE FILEPATH
            "Path to LLVMEnzyme plugin" FORCE)
    endif()
endif()
if(NOT EXISTS "${ENZYME_PLUGIN}")
    message(FATAL_ERROR
        "Set ENZYME_PLUGIN to LLVMEnzyme-${_enzyme_compiler_major}.so")
endif()

function(add_enzyme_fortran_test target)
    cmake_parse_arguments(ARG "" "" "SOURCES" ${ARGN})
    if(NOT ARG_SOURCES)
        message(FATAL_ERROR "add_enzyme_fortran_test requires SOURCES")
    endif()

    set(_absolute_sources "")
    foreach(_source IN LISTS ARG_SOURCES)
        if(IS_ABSOLUTE "${_source}")
            list(APPEND _absolute_sources "${_source}")
        else()
            list(APPEND _absolute_sources
                "${CMAKE_CURRENT_SOURCE_DIR}/${_source}")
        endif()
    endforeach()

    set(_work "${CMAKE_CURRENT_BINARY_DIR}/${target}.enzyme")
    set(_modules "${_work}/modules")
    set(_executable "${_work}/${target}")
    set(_driver "${_work}/pipeline.cmake")
    file(MAKE_DIRECTORY "${_work}" "${_modules}")
    string(REPLACE ";" "|" _packed_sources "${_absolute_sources}")

    file(WRITE "${_driver}" "
set(FLANG \"${CMAKE_Fortran_COMPILER}\")
set(OPT \"${ENZYME_OPT_EXECUTABLE}\")
set(LLVM_LINK \"${ENZYME_LLVM_LINK_EXECUTABLE}\")
set(PLUGIN \"${ENZYME_PLUGIN}\")
set(WORK \"${_work}\")
set(MODULES \"${_modules}\")
set(EXECUTABLE \"${_executable}\")
string(REPLACE \"|\" \";\" SOURCES \"${_packed_sources}\")

function(run_checked)
    execute_process(COMMAND \${ARGN} RESULT_VARIABLE result)
    if(NOT result EQUAL 0)
        message(FATAL_ERROR \"command failed: \${ARGN}\")
    endif()
endfunction()
set(IR_FILES \"\")
set(index 0)
foreach(source IN LISTS SOURCES)
    set(ir \"\${WORK}/stage_\${index}.ll\")
    run_checked(\${FLANG} -O2 -fPIC -module-dir \${MODULES} -I\${MODULES}
        -S -emit-llvm \${source} -o \${ir})
    list(APPEND IR_FILES \${ir})
    math(EXPR index \"\${index} + 1\")
endforeach()

set(linked \"\${WORK}/linked.ll\")
run_checked(\${LLVM_LINK} -S \${IR_FILES} -o \${linked})
set(differentiated \"\${WORK}/enzyme.ll\")
run_checked(\${OPT} -load-pass-plugin=\${PLUGIN} -passes=enzyme
    \${linked} -S -o \${differentiated})
run_checked(\${FLANG} -O2 \${differentiated} -o \${EXECUTABLE})
")

    add_custom_command(
        OUTPUT "${_executable}"
        COMMAND "${CMAKE_COMMAND}" -P "${_driver}"
        DEPENDS ${_absolute_sources} "${_driver}"
        VERBATIM)
    add_custom_target(${target}_build ALL DEPENDS "${_executable}")
    add_test(NAME ${target} COMMAND "${_executable}")
    set_tests_properties(${target} PROPERTIES LABELS enzyme)
endfunction()
