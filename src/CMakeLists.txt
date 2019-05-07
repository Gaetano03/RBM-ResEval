add_executable (SingleMODES SingleMODES.cpp)
target_link_libraries (SingleMODES MODES ${MANDATORY_LIBRARIES})

add_executable (Adaptive_Energy_Mode Adaptive_Energy_Mode.cpp)
target_link_libraries (Adaptive_Energy_Mode MODES ${MANDATORY_LIBRARIES})

add_executable (AdaptiveMODES AdaptiveMODES.cpp)
target_link_libraries (AdaptiveMODES MODES ${MANDATORY_LIBRARIES})

add_executable (Adaptive_ResEval Adaptive_ResEval.cpp)
target_link_libraries (Adaptive_ResEval MODES ${MANDATORY_LIBRARIES})

add_executable (Adaptive_Rec_DTR Adaptive_Rec_DTR.cpp)
target_link_libraries (Adaptive_Rec_DTR MODES ${MANDATORY_LIBRARIES})

add_executable (POD_Rec POD_Rec.cpp)
target_link_libraries (POD_Rec MODES ${MANDATORY_LIBRARIES})

