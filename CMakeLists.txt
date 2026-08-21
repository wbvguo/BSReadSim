cmake_minimum_required(VERSION 3.20)

project(HTSIM LANGUAGES C CXX)

include(CTest)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_SKIP_RPATH ON)

find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)

# HTSlib owns the generic genomics I/O boundary (VCF, compressed textual
# FASTA/BED-family inputs, and BAM serialization). The default is the
# version-pinned Git submodule:
# unlike a configure-time download, it is auditable and works offline after
# `git submodule update --init --recursive`.  A current system installation is
# an explicit opt-in for distro/institutional builds.
set(
    BSREADSIM_HTSLIB_PROVIDER "SUBMODULE"
    CACHE STRING "HTSlib provider: SUBMODULE or SYSTEM"
)
set_property(
    CACHE BSREADSIM_HTSLIB_PROVIDER PROPERTY STRINGS SUBMODULE SYSTEM
)
string(TOUPPER "${BSREADSIM_HTSLIB_PROVIDER}" BSREADSIM_HTSLIB_PROVIDER)
if(NOT BSREADSIM_HTSLIB_PROVIDER MATCHES "^(SUBMODULE|SYSTEM)$")
    message(
        FATAL_ERROR
        "BSREADSIM_HTSLIB_PROVIDER must be SUBMODULE or SYSTEM"
    )
endif()

set(BSREADSIM_HTSLIB_MINIMUM_VERSION "1.24")
set(BSREADSIM_HTSLIB_TARGET "")
if(BSREADSIM_HTSLIB_PROVIDER STREQUAL "SYSTEM")
    find_package(PkgConfig QUIET)
    if(PkgConfig_FOUND)
        pkg_check_modules(
            BSREADSIM_SYSTEM_HTSLIB
            QUIET
            IMPORTED_TARGET
            "htslib>=${BSREADSIM_HTSLIB_MINIMUM_VERSION}"
        )
        if(TARGET PkgConfig::BSREADSIM_SYSTEM_HTSLIB)
            set(BSREADSIM_HTSLIB_TARGET PkgConfig::BSREADSIM_SYSTEM_HTSLIB)
            message(
                STATUS
                "Using system HTSlib ${BSREADSIM_SYSTEM_HTSLIB_VERSION}"
            )
        endif()
    endif()
    if(NOT BSREADSIM_HTSLIB_TARGET)
        message(
            FATAL_ERROR
            "HTSlib ${BSREADSIM_HTSLIB_MINIMUM_VERSION}+ was requested from "
            "the system but was not found through pkg-config"
        )
    endif()
endif()

if(BSREADSIM_HTSLIB_PROVIDER STREQUAL "SUBMODULE")
    include(ExternalProject)
    find_program(BSREADSIM_MAKE_EXECUTABLE NAMES gmake make REQUIRED)
    find_program(BSREADSIM_AUTORECONF_EXECUTABLE NAMES autoreconf REQUIRED)

    set(
        BSREADSIM_HTSLIB_SUBMODULE_DIR
        "${CMAKE_CURRENT_SOURCE_DIR}/HTSLIB"
    )
    set(
        BSREADSIM_HTSLIB_BUILD_SOURCE_DIR
        "${CMAKE_CURRENT_BINARY_DIR}/_deps/htslib-src"
    )
    set(
        BSREADSIM_HTSLIB_INSTALL_DIR
        "${CMAKE_CURRENT_BINARY_DIR}/_deps/htslib-install"
    )
    set(
        BSREADSIM_HTSLIB_ARCHIVE
        "${BSREADSIM_HTSLIB_INSTALL_DIR}/lib/libhts.a"
    )
    if(NOT EXISTS "${BSREADSIM_HTSLIB_SUBMODULE_DIR}/htslib/hts.h"
       OR NOT EXISTS
          "${BSREADSIM_HTSLIB_SUBMODULE_DIR}/htscodecs/htscodecs/htscodecs.h")
        message(
            FATAL_ERROR
            "The pinned HTSlib submodule (including htscodecs) is missing. "
            "Run: git submodule update --init --recursive"
        )
    endif()
    file(
        STRINGS "${BSREADSIM_HTSLIB_SUBMODULE_DIR}/version.sh"
        BSREADSIM_HTSLIB_VERSION_DECLARATION
        REGEX "^VERSION=1\\.24$"
    )
    if(NOT BSREADSIM_HTSLIB_VERSION_DECLARATION)
        message(
            FATAL_ERROR
            "HTSLIB is not the expected HTSlib 1.24 source tree"
        )
    endif()
    file(
        STRINGS
        "${BSREADSIM_HTSLIB_SUBMODULE_DIR}/htscodecs/configure.ac"
        BSREADSIM_HTSCODECS_VERSION_DECLARATION
        REGEX "^AC_INIT\\(htscodecs, 1\\.6\\.7\\)$"
    )
    if(NOT BSREADSIM_HTSCODECS_VERSION_DECLARATION)
        message(
            FATAL_ERROR
            "HTSLIB does not contain the expected htscodecs 1.6.7"
        )
    endif()
    # HTSlib's Autotools bootstrap writes generated files.  Stage a clean copy
    # below the build tree so the pinned submodule remains byte-for-byte clean.
    # Excluding Git administrative files also makes version discovery behave
    # exactly like an official source release rather than following stale
    # relative gitdir pointers from the copied submodule.
    file(MAKE_DIRECTORY "${BSREADSIM_HTSLIB_INSTALL_DIR}/include")
    file(MAKE_DIRECTORY "${BSREADSIM_HTSLIB_INSTALL_DIR}/lib")
    file(MAKE_DIRECTORY "${BSREADSIM_HTSLIB_BUILD_SOURCE_DIR}")
    file(
        COPY "${BSREADSIM_HTSLIB_SUBMODULE_DIR}/"
        DESTINATION "${BSREADSIM_HTSLIB_BUILD_SOURCE_DIR}"
        PATTERN ".git" EXCLUDE
    )
    file(
        WRITE
        "${BSREADSIM_HTSLIB_BUILD_SOURCE_DIR}/htscodecs/htscodecs/version.h"
        "#define HTSCODECS_VERSION_TEXT \"1.6.7\"\n"
    )

    ExternalProject_Add(
        bsreadsim_htslib_build
        DOWNLOAD_COMMAND ""
        SOURCE_DIR "${BSREADSIM_HTSLIB_BUILD_SOURCE_DIR}"
        BUILD_IN_SOURCE TRUE
        CONFIGURE_COMMAND
            "${BSREADSIM_AUTORECONF_EXECUTABLE}" -i
        COMMAND
            "${CMAKE_COMMAND}" -E env
            "CC=${CMAKE_C_COMPILER}"
            "AR=${CMAKE_AR}"
            "RANLIB=${CMAKE_RANLIB}"
            <SOURCE_DIR>/configure
            --disable-bz2
            --disable-gcs
            --disable-libcurl
            --disable-lzma
            --disable-plugins
            --disable-ref-cache
            --disable-s3
            --without-libdeflate
        BUILD_COMMAND "${BSREADSIM_MAKE_EXECUTABLE}" lib-static
        INSTALL_COMMAND
            "${CMAKE_COMMAND}" -E copy
            <SOURCE_DIR>/libhts.a
            "${BSREADSIM_HTSLIB_ARCHIVE}"
        COMMAND
            "${CMAKE_COMMAND}" -E copy_directory
            <SOURCE_DIR>/htslib
            "${BSREADSIM_HTSLIB_INSTALL_DIR}/include/htslib"
        UPDATE_COMMAND ""
        BUILD_BYPRODUCTS "${BSREADSIM_HTSLIB_ARCHIVE}"
    )

    add_library(bsreadsim_htslib STATIC IMPORTED GLOBAL)
    set_target_properties(
        bsreadsim_htslib
        PROPERTIES
            IMPORTED_LOCATION "${BSREADSIM_HTSLIB_ARCHIVE}"
            INTERFACE_INCLUDE_DIRECTORIES
                "${BSREADSIM_HTSLIB_INSTALL_DIR}/include"
            INTERFACE_LINK_LIBRARIES
                "ZLIB::ZLIB;Threads::Threads;${CMAKE_DL_LIBS};m"
    )
    add_dependencies(bsreadsim_htslib bsreadsim_htslib_build)
    set(BSREADSIM_HTSLIB_TARGET bsreadsim_htslib)
    message(STATUS "Using pinned HTSlib 1.24 Git submodule")
endif()

add_library(htsim_compile_options INTERFACE)
target_compile_options(
    htsim_compile_options
    INTERFACE
        $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall;-Wextra;-Wpedantic>
)

# One private library is the complete native simulator. The flat source pairs
# are architectural boundaries; consumers do not need dozens of link targets.
add_library(
    htsim_native STATIC
    HTSIM/src/bam.cpp
    HTSIM/src/core.cpp
    HTSIM/src/fragment.cpp
    HTSIM/src/methdb.cpp
    HTSIM/src/protocol.cpp
    HTSIM/src/reference.cpp
    HTSIM/src/rrbs.cpp
    HTSIM/src/tbs.cpp
    HTSIM/src/utilities.cpp
    HTSIM/src/variant.cpp
    HTSIM/src/wgbs.cpp
)
target_include_directories(
    htsim_native PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/HTSIM/src"
)
target_link_libraries(
    htsim_native
    PUBLIC ZLIB::ZLIB Threads::Threads
    PRIVATE "${BSREADSIM_HTSLIB_TARGET}"
    PRIVATE htsim_compile_options
)
if(TARGET bsreadsim_htslib_build)
    add_dependencies(htsim_native bsreadsim_htslib_build)
endif()

add_executable(htsim_core HTSIM/src/htsim.cpp)
target_link_libraries(
    htsim_core PRIVATE htsim_compile_options htsim_native
)
set_target_properties(
    htsim_core
    PROPERTIES
        OUTPUT_NAME "htsim-core"
        BUILD_RPATH ""
        INSTALL_RPATH ""
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
)
install(TARGETS htsim_core RUNTIME DESTINATION bin)

if(BUILD_TESTING)
    file(
        GLOB_RECURSE HTSIM_CPP_TEST_SOURCES CONFIGURE_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/HTSIM/tests/*_test.cpp"
    )
    foreach(test_source IN LISTS HTSIM_CPP_TEST_SOURCES)
        get_filename_component(test_target "${test_source}" NAME_WE)
        string(REGEX REPLACE "_test$" "" test_name "${test_target}")
        if(test_target STREQUAL "rng_test")
            set(test_name "rng_contract")
        endif()
        add_executable("${test_target}" "${test_source}")
        target_link_libraries(
            "${test_target}"
            PRIVATE htsim_compile_options htsim_native
        )
        set_target_properties(
            "${test_target}"
            PROPERTIES
                RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/tests"
        )
        add_test(NAME "${test_name}" COMMAND "${test_target}")
    endforeach()

    find_package(Python3 COMPONENTS Interpreter)
    if(Python3_Interpreter_FOUND)
        set(
            BSREADSIM_PYTHON_TEST_ROOT
            "${CMAKE_CURRENT_BINARY_DIR}/python-test"
        )
        set(
            BSREADSIM_PYTHON_TEST_SOURCE
            "${BSREADSIM_PYTHON_TEST_ROOT}/src"
        )
        add_custom_target(
            bsreadsim_python_test_package ALL
            COMMAND
                "${CMAKE_COMMAND}" -E rm -rf
                "${BSREADSIM_PYTHON_TEST_ROOT}"
            COMMAND
                "${CMAKE_COMMAND}" -E make_directory
                "${BSREADSIM_PYTHON_TEST_SOURCE}"
            COMMAND
                "${CMAKE_COMMAND}" -E copy_directory
                "${CMAKE_CURRENT_SOURCE_DIR}/src/bsreadsim"
                "${BSREADSIM_PYTHON_TEST_SOURCE}/bsreadsim"
            COMMAND
                "${Python3_EXECUTABLE}" setup.py build_ext
                --build-lib "${BSREADSIM_PYTHON_TEST_SOURCE}"
                --build-temp "${BSREADSIM_PYTHON_TEST_ROOT}/build"
            WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
            COMMENT "Staging the Python package and required native extension"
            VERBATIM
        )

        function(htsim_add_python_test test_name script_name)
            add_test(
                NAME "${test_name}"
                COMMAND
                    "${Python3_EXECUTABLE}"
                    "${CMAKE_CURRENT_SOURCE_DIR}/tests/${script_name}"
                    "$<TARGET_FILE:htsim_core>"
            )
            set_tests_properties(
                "${test_name}"
                PROPERTIES
                    ENVIRONMENT
                        "PYTHONPATH=${BSREADSIM_PYTHON_TEST_SOURCE}"
            )
        endfunction()

        add_test(
            NAME protocol_cross_language
            COMMAND
                "${Python3_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/tests/check_cpp_protocol.py"
                "$<TARGET_FILE:protocol_writer_test>"
                "${CMAKE_CURRENT_BINARY_DIR}/tests/protocol-v2-golden"
                "${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/protocol-v2"
        )
        set_tests_properties(
            protocol_cross_language
            PROPERTIES
                ENVIRONMENT "PYTHONPATH=${BSREADSIM_PYTHON_TEST_SOURCE}"
        )

        htsim_add_python_test(core_cross_language check_core_stream.py)
        htsim_add_python_test(
            pipeline_cross_language check_pipeline_e2e.py
        )
        htsim_add_python_test(
            bam_cross_language check_bam.py
        )
        htsim_add_python_test(
            methdb_cross_language check_methdb.py
        )
        htsim_add_python_test(
            technology_outputs_cross_language check_technology_outputs.py
        )
    endif()
endif()
