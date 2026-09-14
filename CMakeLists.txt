cmake_minimum_required(VERSION 3.20)

project(BSReadSim LANGUAGES C CXX)

file(
    STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/src/bsreadsim/_version.py"
    BSREADSIM_VERSION_ASSIGNMENT
    REGEX
        "^__version__ = \"[0-9]+\\.[0-9]+\\.[0-9]+([-.+][0-9A-Za-z.-]+)?\"$"
)
list(LENGTH BSREADSIM_VERSION_ASSIGNMENT BSREADSIM_VERSION_ASSIGNMENT_COUNT)
if(NOT BSREADSIM_VERSION_ASSIGNMENT_COUNT EQUAL 1)
    message(FATAL_ERROR "src/bsreadsim/_version.py must define one semantic version")
endif()
string(
    REGEX REPLACE
    "^__version__ = \"([^\"]+)\"$"
    "\\1"
    BSREADSIM_VERSION
    "${BSREADSIM_VERSION_ASSIGNMENT}"
)

include(CTest)
include(ExternalProject)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(
    $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang>:-Wall>
    $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang>:-Wextra>
    $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang>:-Wpedantic>
)

find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)
find_program(BSREADSIM_MAKE NAMES gmake make REQUIRED)
find_program(BSREADSIM_AUTORECONF NAMES autoreconf REQUIRED)

# Build the pinned HTSlib source in the build tree. Autotools writes generated
# files, so the submodule itself must remain untouched.
set(HTSLIB_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/htslib")
set(HTSLIB_BUILD "${CMAKE_CURRENT_BINARY_DIR}/_deps/htslib-src")
set(HTSLIB_INSTALL "${CMAKE_CURRENT_BINARY_DIR}/_deps/htslib-install")
set(HTSLIB_ARCHIVE "${HTSLIB_INSTALL}/lib/libhts.a")

if(NOT EXISTS "${HTSLIB_SOURCE}/htslib/hts.h"
   OR NOT EXISTS "${HTSLIB_SOURCE}/htscodecs/htscodecs/htscodecs.h")
    message(
        FATAL_ERROR
        "The pinned HTSlib submodule is missing. Run: "
        "git submodule update --init --recursive"
    )
endif()

file(
    STRINGS "${HTSLIB_SOURCE}/version.sh"
    HTSLIB_VERSION
    REGEX "^VERSION=1\\.24$"
)
file(
    STRINGS "${HTSLIB_SOURCE}/htscodecs/configure.ac"
    HTSCODECS_VERSION
    REGEX "^AC_INIT\\(htscodecs, 1\\.6\\.7\\)$"
)
if(NOT HTSLIB_VERSION OR NOT HTSCODECS_VERSION)
    message(FATAL_ERROR "The pinned HTSlib source has an unexpected version")
endif()

file(MAKE_DIRECTORY "${HTSLIB_INSTALL}/include" "${HTSLIB_INSTALL}/lib")
file(MAKE_DIRECTORY "${HTSLIB_BUILD}")
file(
    COPY "${HTSLIB_SOURCE}/"
    DESTINATION "${HTSLIB_BUILD}"
    PATTERN ".git" EXCLUDE
)
file(
    WRITE "${HTSLIB_BUILD}/htscodecs/htscodecs/version.h"
    "#define HTSCODECS_VERSION_TEXT \"1.6.7\"\n"
)

ExternalProject_Add(
    bsreadsim_htslib_build
    DOWNLOAD_COMMAND ""
    UPDATE_COMMAND ""
    SOURCE_DIR "${HTSLIB_BUILD}"
    BUILD_IN_SOURCE TRUE
    CONFIGURE_COMMAND
        "${BSREADSIM_AUTORECONF}" -i
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
    BUILD_COMMAND "${BSREADSIM_MAKE}" lib-static
    INSTALL_COMMAND
        "${CMAKE_COMMAND}" -E copy
        <SOURCE_DIR>/libhts.a
        "${HTSLIB_ARCHIVE}"
    COMMAND
        "${CMAKE_COMMAND}" -E copy_directory
        <SOURCE_DIR>/htslib
        "${HTSLIB_INSTALL}/include/htslib"
    BUILD_BYPRODUCTS "${HTSLIB_ARCHIVE}"
)

add_library(bsreadsim_htslib STATIC IMPORTED GLOBAL)
set_target_properties(
    bsreadsim_htslib
    PROPERTIES
        IMPORTED_LOCATION "${HTSLIB_ARCHIVE}"
        INTERFACE_INCLUDE_DIRECTORIES "${HTSLIB_INSTALL}/include"
        INTERFACE_LINK_LIBRARIES
            "ZLIB::ZLIB;Threads::Threads;${CMAKE_DL_LIBS};m"
)
add_dependencies(bsreadsim_htslib bsreadsim_htslib_build)

add_executable(
    htsim
    htsim/src/struct.cpp
    htsim/src/utility.cpp
    htsim/src/haplo.cpp
    htsim/src/methdb.cpp
    htsim/src/profile.cpp
    htsim/src/mode.cpp
    htsim/src/rrcut.cpp
    htsim/src/stream.cpp
    htsim/src/reads.cpp
    htsim/src/option.cpp
    htsim/src/htsim.cpp
)
target_compile_definitions(
    htsim PRIVATE BSREADSIM_VERSION="${BSREADSIM_VERSION}"
)
target_link_libraries(
    htsim PRIVATE bsreadsim_htslib ZLIB::ZLIB Threads::Threads
)
set_target_properties(
    htsim
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
)
install(TARGETS htsim RUNTIME DESTINATION bin)

set(TEST_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/tests")

if(BUILD_TESTING AND EXISTS "${TEST_ROOT}")
    find_package(Python3 REQUIRED COMPONENTS Interpreter)
    set(PYTHON_TEST_ROOT "${CMAKE_CURRENT_BINARY_DIR}/python-test")
    set(PYTHON_TEST_PACKAGE "${PYTHON_TEST_ROOT}/src")
    add_custom_target(
        bsreadsim_python_test_package ALL
        COMMAND "${CMAKE_COMMAND}" -E rm -rf "${PYTHON_TEST_ROOT}"
        COMMAND "${CMAKE_COMMAND}" -E make_directory "${PYTHON_TEST_PACKAGE}"
        COMMAND
            "${CMAKE_COMMAND}" -E copy_directory
            "${CMAKE_CURRENT_SOURCE_DIR}/src/bsreadsim"
            "${PYTHON_TEST_PACKAGE}/bsreadsim"
        COMMAND
            "${CMAKE_COMMAND}" -E copy_directory
            "${CMAKE_CURRENT_SOURCE_DIR}/data"
            "${PYTHON_TEST_PACKAGE}/bsreadsim/data"
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        COMMENT "Building the Python test package"
        VERBATIM
    )

    function(bsreadsim_add_python_test test_name script_name)
        add_test(
            NAME "${test_name}"
            COMMAND
                "${Python3_EXECUTABLE}"
                "${TEST_ROOT}/integration/${script_name}"
                "$<TARGET_FILE:htsim>"
        )
        set_tests_properties(
            "${test_name}"
            PROPERTIES ENVIRONMENT "PYTHONPATH=${PYTHON_TEST_PACKAGE}"
        )
    endfunction()

    add_test(
        NAME python_unit
        COMMAND
            "${Python3_EXECUTABLE}" -m unittest discover
            -s "${TEST_ROOT}/unit"
            -t "${CMAKE_CURRENT_SOURCE_DIR}"
    )
    set_tests_properties(
        python_unit
        PROPERTIES ENVIRONMENT "PYTHONPATH=${PYTHON_TEST_PACKAGE}"
    )

    add_test(
        NAME semantics
        COMMAND
            "${Python3_EXECUTABLE}" -m unittest discover
            -s "${TEST_ROOT}/semantic"
            -t "${CMAKE_CURRENT_SOURCE_DIR}"
    )
    set_tests_properties(
        semantics
        PROPERTIES
            ENVIRONMENT
                "PYTHONPATH=${PYTHON_TEST_PACKAGE};BSREADSIM_HTSIM=$<TARGET_FILE:htsim>"
            TIMEOUT 1800
    )

    bsreadsim_add_python_test(core_cross_language check_core_stream.py)
    bsreadsim_add_python_test(pipeline_cross_language check_pipeline_e2e.py)
    bsreadsim_add_python_test(bam_cross_language check_bam.py)
    bsreadsim_add_python_test(methdb_cross_language check_methdb.py)
    bsreadsim_add_python_test(
        variant_catalog_cross_language check_variant_catalog.py
    )
    bsreadsim_add_python_test(
        technology_outputs_cross_language check_technology_outputs.py
    )
    bsreadsim_add_python_test(
        example_data_cross_language check_example_data.py
    )
    bsreadsim_add_python_test(
        validate_cross_language check_validate.py
    )
endif()

if(
    BUILD_TESTING
    AND NOT EXISTS "${HTSIM_TEST_ROOT}"
    AND NOT EXISTS "${TEST_ROOT}"
)
    message(STATUS "Tests are absent; configuring production targets only")
endif()
