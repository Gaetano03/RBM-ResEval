#add_executable (SingleMODES SingleMODES.cpp)
#target_link_libraries (SingleMODES MODES ${MANDATORY_LIBRARIES})

#add_executable (FLUID_MODES_DNS FLUID_MODES_DNS.cpp)
#target_link_libraries (FLUID_MODES_DNS MODES ${MANDATORY_LIBRARIES})

#add_executable (Adaptive_Energy_Mode Adaptive_Energy_Mode.cpp)
#target_link_libraries (Adaptive_Energy_Mode MODES ${MANDATORY_LIBRARIES})

#add_executable (AdaptiveMODES AdaptiveMODES.cpp)
#target_link_libraries (AdaptiveMODES MODES ${MANDATORY_LIBRARIES})

#add_executable (Adaptive_MODESConservative Adaptive_MODESConservative.cpp)
#target_link_libraries (Adaptive_MODESConservative MODES ${MANDATORY_LIBRARIES})

#add_executable (Adaptive_ResEval Adaptive_ResEval.cpp)
#target_link_libraries (Adaptive_ResEval MODES ${MANDATORY_LIBRARIES})

#add_executable (Adaptive_Rec_DTR Adaptive_Rec_DTR.cpp)
#target_link_libraries (Adaptive_Rec_DTR MODES ${MANDATORY_LIBRARIES})

#add_executable (POD_Rec POD_Rec.cpp)
#target_link_libraries (POD_Rec MODES ${MANDATORY_LIBRARIES})

#add_executable (pagmo_test pagmo_test.cpp)
#target_link_libraries (pagmo_test MODES ${MANDATORY_LIBRARIES})

function( ADD_RBM_TESTCASE arg1)

    add_executable(${arg1} ${arg1}.cpp)
    target_link_libraries(${arg1} MODES ${MANDATORY_LIBRARIES})

endfunction()


#ADD_RBM_TESTCASE(Adaptive_Energy_Mode)
#ADD_RBM_TESTCASE(AdaptiveMODES)
#ADD_RBM_TESTCASE(Adaptive_MODESConservative)
#ADD_RBM_TESTCASE(Adaptive_Rec_DTR)
ADD_RBM_TESTCASE(Adaptive_ResEval)
#ADD_RBM_TESTCASE(FLUID_MODES_DNS)
#ADD_RBM_TESTCASE(pagmo_test)
ADD_RBM_TESTCASE(ResEval_DtRes)
#ADD_RBM_TESTCASE(test)
#ADD_RBM_TESTCASE(SingleMODES)

