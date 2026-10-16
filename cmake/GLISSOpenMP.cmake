# Keep BLAS and GLISS on one OpenMP runtime. This bounded compatibility path
# covers macOS GNU Fortran with a dynamically linked LLVM-OpenMP BLAS; ordinary
# pthread/GNU-OpenMP BLAS builds retain FindOpenMP's compiler runtime.
# LLVM implements the GNU GOMP ABI, but compatibility is established below by
# compiling and running the actual compiler/runtime/BLAS combination, not by
# assuming every GNU/LLVM version is interchangeable.
# https://cmake.org/cmake/help/latest/module/FindOpenMP.html
# https://gcc.gnu.org/onlinedocs/libgomp/Enabling-OpenMP.html
# https://github.com/llvm/llvm-project/blob/main/openmp/runtime/src/kmp_gsupport.cpp
include_guard(GLOBAL)

set(GLISS_OPENMP_RUNTIME AUTO CACHE STRING
    "OpenMP runtime policy: AUTO or COMPILER (reject incompatible BLAS runtimes)")
set_property(CACHE GLISS_OPENMP_RUNTIME PROPERTY STRINGS AUTO COMPILER)
if(NOT GLISS_OPENMP_RUNTIME MATCHES "^(AUTO|COMPILER)$")
    message(FATAL_ERROR "GLISS_OPENMP_RUNTIME must be AUTO or COMPILER")
endif()

# Recognize a command-line/native-language override before FindOpenMP fills its
# normal cache. Do not silently replace an explicitly requested GNU runtime.
set(_gliss_explicit_runtime FALSE)
if(DEFINED OpenMP_Fortran_LIB_NAMES)
    get_property(_gliss_openmp_help CACHE OpenMP_Fortran_LIB_NAMES PROPERTY HELPSTRING)
    if(NOT _gliss_openmp_help STREQUAL
            "Fortran compiler libraries for OpenMP parallelization")
        set(_gliss_explicit_runtime TRUE)
    endif()
endif()
if(DEFINED OpenMP_gomp_LIBRARY)
    get_property(_gliss_openmp_help CACHE OpenMP_gomp_LIBRARY PROPERTY HELPSTRING)
    if(NOT _gliss_openmp_help STREQUAL "Path to the gomp library for OpenMP")
        set(_gliss_explicit_runtime TRUE)
    endif()
endif()
find_package(OpenMP REQUIRED COMPONENTS Fortran)

function(_gliss_link_library_paths output)
    set(queue ${ARGN})
    set(visited)
    set(paths)
    while(queue)
        list(POP_FRONT queue item)
        if(item IN_LIST visited)
            continue()
        endif()
        list(APPEND visited "${item}")
        if(TARGET "${item}")
            foreach(property INTERFACE_LINK_LIBRARIES IMPORTED_LOCATION
                    IMPORTED_LOCATION_RELEASE IMPORTED_LOCATION_DEBUG)
                get_target_property(values "${item}" "${property}")
                if(values)
                    list(APPEND queue ${values})
                endif()
            endforeach()
        elseif(IS_ABSOLUTE "${item}" AND EXISTS "${item}" AND NOT IS_DIRECTORY "${item}")
            list(APPEND paths "${item}")
        endif()
    endwhile()
    list(REMOVE_DUPLICATES paths)
    set(${output} "${paths}" PARENT_SCOPE)
endfunction()

function(_gliss_macos_resolve_dependency output dependency library)
    get_filename_component(directory "${library}" DIRECTORY)
    if(dependency MATCHES "^@loader_path/")
        string(REPLACE "@loader_path" "${directory}" dependency "${dependency}")
    elseif(dependency MATCHES "^@rpath/(.+)$")
        set(suffix "${CMAKE_MATCH_1}")
        execute_process(COMMAND "${CMAKE_OTOOL}" -l "${library}"
            RESULT_VARIABLE status OUTPUT_VARIABLE commands)
        if(NOT status EQUAL 0)
            message(FATAL_ERROR "Cannot inspect LC_RPATH in ${library}")
        endif()
        string(REPLACE "\n" ";" commands "${commands}")
        set(in_rpath FALSE)
        set(dependency "")
        foreach(line IN LISTS commands)
            if(line MATCHES "^[ \t]*cmd LC_RPATH$")
                set(in_rpath TRUE)
            elseif(in_rpath AND line MATCHES "^[ \t]*path (.+) \\(offset [0-9]+\\)$")
                set(rpath "${CMAKE_MATCH_1}")
                string(REPLACE "@loader_path" "${directory}" rpath "${rpath}")
                if(IS_ABSOLUTE "${rpath}" AND EXISTS "${rpath}/${suffix}")
                    set(dependency "${rpath}/${suffix}")
                    break()
                endif()
                set(in_rpath FALSE)
            endif()
        endforeach()
    endif()
    if(IS_ABSOLUTE "${dependency}" AND EXISTS "${dependency}")
        get_filename_component(dependency "${dependency}" REALPATH)
        set(${output} "${dependency}" PARENT_SCOPE)
    else()
        set(${output} "" PARENT_SCOPE)
    endif()
endfunction()

function(_gliss_macos_openmp_dependencies output)
    _gliss_link_library_paths(queue ${ARGN})
    set(visited)
    set(runtimes)
    while(queue)
        list(POP_FRONT queue library)
        get_filename_component(library "${library}" REALPATH)
        if(library IN_LIST visited OR library MATCHES "\\.(a|tbd)$")
            continue()
        endif()
        list(APPEND visited "${library}")
        execute_process(COMMAND "${CMAKE_OTOOL}" -L "${library}"
            RESULT_VARIABLE status OUTPUT_VARIABLE dependencies ERROR_VARIABLE error)
        if(NOT status EQUAL 0)
            message(FATAL_ERROR "Cannot inspect BLAS dependency ${library}: ${error}")
        endif()
        string(REPLACE "\n" ";" dependencies "${dependencies}")
        foreach(line IN LISTS dependencies)
            if(NOT line MATCHES "^[ \t]+(.+)[ \t]+\\(compatibility version ")
                continue()
            endif()
            string(STRIP "${CMAKE_MATCH_1}" dependency)
            _gliss_macos_resolve_dependency(resolved "${dependency}" "${library}")
            if(dependency MATCHES "/lib(omp|gomp|iomp5)(\\.[0-9]+)*\\.dylib$")
                if(NOT resolved)
                    message(FATAL_ERROR
                        "Cannot resolve BLAS OpenMP dependency ${dependency} from ${library}; "
                        "its absolute/@loader_path identity or declared LC_RPATH must resolve")
                endif()
                list(APPEND runtimes "${resolved}")
            endif()
            if(resolved)
                list(APPEND queue "${resolved}")
            endif()
        endforeach()
    endwhile()
    list(REMOVE_DUPLICATES runtimes)
    set(${output} "${runtimes}" PARENT_SCOPE)
endfunction()

function(_gliss_check_gnu_link_flags output)
    set(configurations DEBUG RELEASE RELWITHDEBINFO MINSIZEREL
        ${CMAKE_CONFIGURATION_TYPES} ${CMAKE_BUILD_TYPE} ${CMAKE_TRY_COMPILE_CONFIGURATION})
    list(TRANSFORM configurations TOUPPER)
    list(REMOVE_DUPLICATES configurations)
    set(values)
    foreach(prefix CMAKE_Fortran_FLAGS CMAKE_EXE_LINKER_FLAGS
            CMAKE_SHARED_LINKER_FLAGS CMAKE_MODULE_LINKER_FLAGS)
        set(variables "${prefix}")
        foreach(configuration IN LISTS configurations)
            list(APPEND variables "${prefix}_${configuration}")
        endforeach()
        foreach(variable IN LISTS variables)
            list(APPEND values "${variable}=${${variable}}")
            if("${${variable}}" MATCHES "(^|[ ;,])(-fopenmp|-lgomp|-l[ \t]+gomp)([ ;,=]|$)"
                    OR "${${variable}}" MATCHES "libgomp[.](dylib|a)")
                message(FATAL_ERROR
                    "${variable} contains a GNU runtime link flag; remove it and "
                    "let GLISS::OpenMP_Fortran supply compile-only OpenMP flags")
            endif()
        endforeach()
    endforeach()
    set(${output} "${values}" PARENT_SCOPE)
endfunction()

add_library(GLISS::OpenMP_Fortran INTERFACE IMPORTED)
set_property(TARGET GLISS::OpenMP_Fortran PROPERTY
    INTERFACE_LINK_LIBRARIES OpenMP::OpenMP_Fortran)
set(GLISS_OPENMP_SELECTED_LIBRARY "${OpenMP_Fortran_LIBRARIES}" CACHE INTERNAL
    "Actual runtime linked by GLISS::OpenMP_Fortran" FORCE)

if(APPLE)
    if(NOT CMAKE_OTOOL)
        find_program(CMAKE_OTOOL otool REQUIRED)
    endif()
    _gliss_macos_openmp_dependencies(_gliss_blas_runtimes
        LAPACK::LAPACK BLAS::BLAS ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
    list(LENGTH _gliss_blas_runtimes _gliss_runtime_count)
    if(_gliss_runtime_count GREATER 1)
        message(FATAL_ERROR "BLAS loads multiple OpenMP runtimes: ${_gliss_blas_runtimes}")
    endif()
    if(_gliss_runtime_count EQUAL 1)
        list(GET _gliss_blas_runtimes 0 _gliss_blas_runtime)
        get_filename_component(_gliss_blas_runtime_name "${_gliss_blas_runtime}" NAME)
        set(_gliss_llvm_gomp FALSE)
        if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU"
                AND _gliss_blas_runtime_name MATCHES "^libomp(\\.[0-9]+)*\\.dylib$")
            # At link time GNU -fopenmp adds libgomp, even when FindOpenMP
            # itself has already been explicitly configured to use libomp.
            _gliss_check_gnu_link_flags(_gliss_runtime_flags)
            set(_gliss_llvm_gomp TRUE)
        endif()
        set(_gliss_native_runtimes)
        foreach(library IN LISTS OpenMP_Fortran_LIBRARIES)
            if(EXISTS "${library}")
                get_filename_component(library "${library}" REALPATH)
            endif()
            list(APPEND _gliss_native_runtimes "${library}")
        endforeach()
        if(NOT _gliss_blas_runtime IN_LIST _gliss_native_runtimes OR _gliss_llvm_gomp)
            if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "GNU"
                    OR NOT _gliss_blas_runtime_name MATCHES "^libomp(\\.[0-9]+)*\\.dylib$")
                message(FATAL_ERROR
                    "BLAS runtime ${_gliss_blas_runtime} differs from Fortran OpenMP "
                    "${OpenMP_Fortran_LIBRARIES}; configure one compatible runtime")
            endif()
            if(NOT _gliss_blas_runtime IN_LIST _gliss_native_runtimes
                    AND (GLISS_OPENMP_RUNTIME STREQUAL "COMPILER" OR _gliss_explicit_runtime))
                message(FATAL_ERROR
                    "The requested Fortran OpenMP runtime conflicts with BLAS's LLVM libomp. "
                    "Use GLISS_OPENMP_RUNTIME=AUTO and remove conflicting OpenMP overrides")
            endif()
            if(CMAKE_CROSSCOMPILING AND NOT CMAKE_CROSSCOMPILING_EMULATOR)
                message(FATAL_ERROR
                    "GNU/LLVM OpenMP compatibility requires execution; configure natively "
                    "or provide CMAKE_CROSSCOMPILING_EMULATOR")
            endif()
            set_target_properties(GLISS::OpenMP_Fortran PROPERTIES
                INTERFACE_COMPILE_OPTIONS
                    "$<$<COMPILE_LANGUAGE:Fortran>:SHELL:${OpenMP_Fortran_FLAGS}>"
                INTERFACE_LINK_LIBRARIES "${_gliss_blas_runtime}")
            if(OpenMP_Fortran_INCLUDE_DIRS)
                set_property(TARGET GLISS::OpenMP_Fortran PROPERTY
                    INTERFACE_INCLUDE_DIRECTORIES "${OpenMP_Fortran_INCLUDE_DIRS}")
            endif()
            # Probe the actual GOMP feature set and concurrent disjoint BLAS
            # consumer. No BLAS thread-setting API or global lock is used.
            set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
                "${CMAKE_CURRENT_LIST_DIR}/../test/test_parallel_blas.f90")
            file(READ "${CMAKE_CURRENT_LIST_DIR}/../test/test_parallel_blas.f90"
                _gliss_openmp_probe)
            _gliss_link_library_paths(_gliss_blas_paths
                LAPACK::LAPACK BLAS::BLAS ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
            set(_gliss_blas_digests)
            foreach(library IN LISTS _gliss_blas_paths)
                file(SHA256 "${library}" digest)
                list(APPEND _gliss_blas_digests "${library}:${digest}")
            endforeach()
            set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
                "${_gliss_blas_runtime}" ${_gliss_blas_paths})
            file(SHA256 "${_gliss_blas_runtime}" _gliss_runtime_digest)
            set(_gliss_probe_inputs "${CMAKE_Fortran_COMPILER}"
                "${CMAKE_Fortran_COMPILER_VERSION}" "${CMAKE_BUILD_TYPE}"
                "${CMAKE_TRY_COMPILE_CONFIGURATION}" "${_gliss_runtime_flags}"
                "${OpenMP_Fortran_FLAGS}" "${_gliss_runtime_digest}"
                "${_gliss_blas_digests}" "${_gliss_openmp_probe}")
            string(SHA256 _gliss_probe_digest
                "${_gliss_probe_inputs}")
            include(CMakePushCheckState)
            include(CheckFortranSourceRuns)
            cmake_push_check_state(RESET)
            set(CMAKE_REQUIRED_LIBRARIES GLISS::OpenMP_Fortran LAPACK::LAPACK)
            # A separate incorrectly loaded BLAS runtime would retain the
            # environment's setting, even when the probe sets its own team.
            # One inherited thread can therefore mask the incompatibility.
            set(_gliss_had_omp_threads FALSE)
            if(DEFINED ENV{OMP_NUM_THREADS})
                set(_gliss_had_omp_threads TRUE)
                set(_gliss_previous_omp_threads "$ENV{OMP_NUM_THREADS}")
            endif()
            set(ENV{OMP_NUM_THREADS} 2)
            check_fortran_source_runs("${_gliss_openmp_probe}"
                "GLISS_OPENMP_COMPATIBLE_${_gliss_probe_digest}")
            if(_gliss_had_omp_threads)
                set(ENV{OMP_NUM_THREADS} "${_gliss_previous_omp_threads}")
            else()
                unset(ENV{OMP_NUM_THREADS})
            endif()
            cmake_pop_check_state()
            if(NOT "${GLISS_OPENMP_COMPATIBLE_${_gliss_probe_digest}}")
                message(FATAL_ERROR
                    "GNU Fortran/LLVM libomp failed the two-thread feature and concurrent "
                    "BLAS consumer probe; choose a compatible compiler/BLAS configuration")
            endif()
            set(GLISS_OPENMP_SELECTED_LIBRARY "${_gliss_blas_runtime}" CACHE INTERNAL
                "Actual runtime linked by GLISS::OpenMP_Fortran" FORCE)
            message(STATUS "GLISS OpenMP runtime: ${_gliss_blas_runtime} (GOMP/BLAS gates passed)")
        else()
            message(STATUS "GLISS OpenMP runtime matches BLAS: ${_gliss_blas_runtime}")
        endif()
    endif()
endif()
