#project(utopiae-orbit-dependencies
#        LANGUAGES CXX
#        VERSION 1.0.0)

# Include external libraries
include(ExternalProject)

# When gsl is required
if (ENABLE_GSL)
    find_package(GSL REQUIRED)
    list(APPEND MANDATORY_LIBRARIES "${GSL_LIBRARIES}")
    include_directories(${GSL_INCLUDE_DIRS})
endif()

# Smart math is required
include(cmake/external/smart-math.cmake)
list (APPEND MANDATORY_LIBRARIES "${SMART_MATH_STATIC_LIBRARY}")
list (APPEND MANDATORY_LIBRARIES "${SMART_MATH_LIBRARY}")
include_directories (${SMART_MATH_INCLUDE_DIR})

# Smart uq is optional
if (ENABLE_UQ)

    # Add compiler flag
    add_definitions(-D__USE_SMARTUQ)

    include(cmake/external/smart-uq.cmake)
    list (APPEND MANDATORY_LIBRARIES "${SMART_UQ_STATIC_LIBRARY}")
    list (APPEND MANDATORY_LIBRARIES "${SMART_UQ_LIBRARY}")
    include_directories (${SMART_UQ_INCLUDE_DIR})
endif()

# Pagmo is optional
if (ENABLE_PAGMO)

    # Add compiler flag
    add_definitions(-D__USE_PAGMO)

    include(cmake/external/pagmo.cmake)
    list (APPEND MANDATORY_LIBRARIES "${PAGMO_STATIC_LIBRARY}")
    list (APPEND MANDATORY_LIBRARIES "${PAGMO_LIBRARY}")

    # Eigen is requested to be installed in the system
    set(EIGEN_INCLUDE_DIR /usr/include/eigen3/)

    include_directories (${PAGMO_INCLUDE_DIR} ${EIGEN_INCLUDE_DIR})

endif()


