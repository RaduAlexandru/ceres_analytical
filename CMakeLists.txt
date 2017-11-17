cmake_minimum_required(VERSION 2.6)

project (opt)


LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_modules)



find_package(Eigen REQUIRED)
find_package(LAPACK QUIET)
find_package(SuiteSparse)
#find_package(CXSparse)
find_package(Gflags)
find_package(Glog)
find_package(OpenMP QUIET)





include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-std=c++11" COMPILER_HAS_CXX11_FLAG)
if (CXX11 AND COMPILER_HAS_CXX11_FLAG)
  # Update CMAKE_REQUIRED_FLAGS used by CheckCXXSourceCompiles to include
  # -std=c++11 s/t we will detect the C++11 versions of unordered_map &
  # shared_ptr if they exist.
  set(CMAKE_REQUIRED_FLAGS -std=c++11)
endif (CXX11 AND COMPILER_HAS_CXX11_FLAG)


# Set the Ceres compile definitions for the unordered_map configuration.
include(FindUnorderedMap)
find_unordered_map()
if (UNORDERED_MAP_FOUND)
  if (HAVE_UNORDERED_MAP_IN_STD_NAMESPACE)
    list(APPEND CERES_COMPILE_OPTIONS CERES_STD_UNORDERED_MAP)
  endif(HAVE_UNORDERED_MAP_IN_STD_NAMESPACE)
  if (HAVE_UNORDERED_MAP_IN_TR1_NAMESPACE)
    list(APPEND CERES_COMPILE_OPTIONS CERES_STD_UNORDERED_MAP_IN_TR1_NAMESPACE)
  endif(HAVE_UNORDERED_MAP_IN_TR1_NAMESPACE)
  if (HAVE_TR1_UNORDERED_MAP_IN_TR1_NAMESPACE)
    list(APPEND CERES_COMPILE_OPTIONS CERES_TR1_UNORDERED_MAP)
  endif(HAVE_TR1_UNORDERED_MAP_IN_TR1_NAMESPACE)
else (UNORDERED_MAP_FOUND)
  message("-- Replacing unordered_map/set with map/set (warning: slower!), "
    "try enabling CXX11 option if you expect C++11 to be available.")
  list(APPEND CERES_COMPILE_OPTIONS CERES_NO_UNORDERED_MAP)
endif()

# Set the Ceres compile definitions for the shared_ptr configuration.
include(FindSharedPtr)
find_shared_ptr()
if (SHARED_PTR_FOUND)
  if (SHARED_PTR_TR1_MEMORY_HEADER)
    list(APPEND CERES_COMPILE_OPTIONS CERES_TR1_MEMORY_HEADER)
  endif (SHARED_PTR_TR1_MEMORY_HEADER)
  if (SHARED_PTR_TR1_NAMESPACE)
    list(APPEND CERES_COMPILE_OPTIONS CERES_TR1_SHARED_PTR)
  endif (SHARED_PTR_TR1_NAMESPACE)
else (SHARED_PTR_FOUND)
  message(FATAL_ERROR "Unable to find shared_ptr, try enabling CXX11 option "
    "if you expect C++11 to be available.")
endif (SHARED_PTR_FOUND)

# To ensure that CXX11 accurately reflects whether we are using C++11,
# check if it is required given where the potentially C++11 features Ceres
# uses were found, and disable it if C++11 is not being used.
if (CXX11)
  if (NOT HAVE_SHARED_PTR_IN_STD_NAMESPACE AND
      NOT HAVE_UNORDERED_MAP_IN_STD_NAMESPACE)
    message("-- Failed to find C++11 components in C++11 locations & "
      "namespaces, disabling CXX11.")
    update_cache_variable(CXX11 OFF)
  else()
    message("   ==============================================================")
    message("   Compiling Ceres using C++11.  This will result in a version ")
    message("   of Ceres that will require the use of C++11 in client code.")
    message("   ==============================================================")
    list(APPEND CERES_COMPILE_OPTIONS CERES_USE_CXX11)
    if (COMPILER_HAS_CXX11_FLAG)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    endif()
  endif()
endif(CXX11)






include_directories(
  include
  internal
  internal/ceres
  ${GLOG_INCLUDE_DIRS})

include_directories(SYSTEM ${EIGEN_INCLUDE_DIRS})
include_directories(${SUITESPARSE_INCLUDE_DIRS})
include_directories(${CXSPARSE_INCLUDE_DIRS})
include_directories(${GFLAGS_INCLUDE_DIRS})



SET(CMAKE_CXX_FLAGS "-std=c++0x")				# New C11
set(PROJECT_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)			# Source folder
set(PROJECT_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/src)		# Header folder

set(CMAKE_BUILD_TYPE release)

include_directories(${PROJECT_INCLUDE_DIR})
include_directories(${PROJECT_SOURCE_DIR})

include_directories(${CMAKE_SOURCE_DIR})



find_package(Ceres REQUIRED)
include_directories(${CERES_INCLUDE_DIR})
include_directories(${GFLAGS_INCLUDE_DIRS})
message( ${CERES_INCLUDE_DIR})
message ( ${GFLAGS_INCLUDE_DIRS})


add_definitions(-DCERES_GFLAGS_NAMESPACE=${GFLAGS_NAMESPACE})



add_executable( mpopt
               ${PROJECT_SOURCE_DIR}/mp_opt.cpp
               ${PROJECT_SOURCE_DIR}/ba_problem.cpp
               ${PROJECT_SOURCE_DIR}/custom_pose_parametrization.cpp)
target_link_libraries(mpopt ceres ${GFLAGS_LIBRARIES})
add_executable(opt
               ${PROJECT_SOURCE_DIR}/opt.cpp
               ${PROJECT_SOURCE_DIR}/bal_problem.cpp
               ${PROJECT_SOURCE_DIR}/custom_pose_parametrization.cpp)
target_link_libraries(opt ceres ${GFLAGS_LIBRARIES})
