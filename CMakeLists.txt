cmake_minimum_required(VERSION 3.10)

project(DG_ExProj
    VERSION 0.1.0
    LANGUAGES CXX
)

set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)

# Add optimization flags explicitly
if (CMAKE_SYSTEM_PROCESSOR MATCHES x86_64)
    add_compile_options(-O3 -march=native -ffast-math -DNDEBUG)
elseif (CMAKE_SYSTEM_PROCESSOR MATCHES aarch64)
    add_compile_options(-O3 -march=armv8-a -ffast-math -DNDEBUG)
else()
    add_compile_options(-O3 -march=native -ffast-math -DNDEBUG)
endif()

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(DG_Sparse REQUIRED NO_MODULE)

add_library(DG_ExProj SHARED
    src/exproj_utils.cpp
    src/exproj_setup.cpp
    src/exproj_solve.cpp
)

target_compile_features(DG_ExProj
    PUBLIC
        cxx_std_11
)

target_include_directories(DG_ExProj
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(DG_ExProj
    PUBLIC
        Eigen3::Eigen
        DG_Sparse::DG_Sparse
)

set_target_properties(DG_ExProj PROPERTIES
    OUTPUT_NAME "dg_exproj"
    POSITION_INDEPENDENT_CODE ON
)

install(
    DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(
    TARGETS DG_ExProj
    EXPORT DG_ExProjTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(
    EXPORT DG_ExProjTargets
    FILE DG_ExProjTargets.cmake
    NAMESPACE DG_ExProj::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/DG_ExProj
)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/DG_ExProjConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/DG_ExProjConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/DG_ExProj
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/DG_ExProjConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(
    FILES
        ${CMAKE_CURRENT_BINARY_DIR}/DG_ExProjConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/DG_ExProjConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/DG_ExProj
)


# cmake_minimum_required(VERSION 3.10)

# project(run LANGUAGES CXX)

# # Make sure we're building with Release flags
# set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)

# # Set default C++ standard
# set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)

# # Add optimization flags explicitly
# if (CMAKE_SYSTEM_PROCESSOR MATCHES x86_64)
#     add_compile_options(-O3 -march=native -ffast-math -DNDEBUG)
# elseif (CMAKE_SYSTEM_PROCESSOR MATCHES aarch64)
#     add_compile_options(-O3 -march=armv8-a -ffast-math -DNDEBUG)
# else()
#     add_compile_options(-O3 -march=native -ffast-math -DNDEBUG)
# endif()

# # Find Eigen
# find_package(Eigen3 3.4 REQUIRED NO_MODULE)
# find_package(DG_Sparse REQUIRED NO_MODULE)

# add_executable( run main.cpp )

# include_directories(
#     ${PROJECT_SOURCE_DIR}/inc
#     ${PROJECT_SOURCE_DIR}/../utils/inc
#     # ${PROJECT_SOURCE_DIR}/../sparse_matrix/inc
# )

# add_library( solver SHARED
#     ${PROJECT_SOURCE_DIR}/src/exproj.cpp
#     ${PROJECT_SOURCE_DIR}/src/solver_setup.cpp
#     ${PROJECT_SOURCE_DIR}/src/solver_solve.cpp
# )

# add_library( utils SHARED
#     ${PROJECT_SOURCE_DIR}/../utils/src/utils.cpp
#     ${PROJECT_SOURCE_DIR}/../utils/src/linalg_utils.cpp
# )

# target_link_libraries( run Eigen3::Eigen utils solver ) 

# target_link_libraries( solver Eigen3::Eigen DG_Sparse::DG_Sparse utils )
# # target_link_libraries( solver Eigen3::Eigen utils )
# target_link_libraries( utils  Eigen3::Eigen )

# (Optional) Print flags at configure-time for debugging
message(STATUS "CMAKE_BUILD_TYPE  : ${CMAKE_BUILD_TYPE}")
message(STATUS "CMAKE_CXX_COMPILER: ${CMAKE_CXX_COMPILER}")
message(STATUS "Compile flags     : ${CMAKE_CXX_FLAGS}")
message(STATUS "HOST PROCESSOR    : ${CMAKE_HOST_SYSTEM_PROCESSOR}")
message(STATUS "TARGET PROCESSOR  : ${CMAKE_SYSTEM_PROCESSOR}")