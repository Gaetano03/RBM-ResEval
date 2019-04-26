# Include external libraries
include(ExternalProject)

include(cmake/external/smart-math.cmake)
if (BUILD_STATIC)
    list(APPEND MANDATORY_LIBRARIES "${SMART_MATH_STATIC_LIBRARY}")
else ()
    list(APPEND MANDATORY_LIBRARIES "${SMART_MATH_LIBRARY}")
endif ()
include_directories("${SMART_MATH_INCLUDE_DIR}")

include(cmake/external/smart-uq.cmake)
if (BUILD_STATIC)
    list(APPEND MANDATORY_LIBRARIES "${SMART_UQ_STATIC_LIBRARY}")
else ()
    list(APPEND MANDATORY_LIBRARIES "${SMART_UQ_LIBRARY}")
endif ()
include_directories("${SMART_UQ_INCLUDE_DIR}")
