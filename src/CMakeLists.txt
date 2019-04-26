add_executable (SingleMODES SingleMODES.cpp)
target_link_libraries (SingleMODES MODES ${MANDATORY_LIBRARIES})

add_executable (Adaptive_Energy_Mode Adaptive_Energy_Mode.cpp)
target_link_libraries (Adaptive_Energy_Mode MODES ${MANDATORY_LIBRARIES})

add_executable (AdaptiveMODES AdaptiveMODES.cpp)
target_link_libraries (AdaptiveMODES MODES ${MANDATORY_LIBRARIES})

add_executable (Adaptive_ResEval Adaptive_ResEval.cpp)
target_link_libraries (Adaptive_ResEval MODES ${MANDATORY_LIBRARIES})

