#CMake file for Adding executables

function( ADD_RBM_TESTCASE arg1)

    add_executable(${arg1} ${arg1}.cpp)
    target_link_libraries(${arg1} MODES ${MANDATORY_LIBRARIES})

endfunction()


#ADD_RBM_TESTCASE(Adaptive_Energy_Mode)
#ADD_RBM_TESTCASE(AdaptiveMODES)
#ADD_RBM_TESTCASE(Adaptive_MODESConservative)
ADD_RBM_TESTCASE(Adaptive_Rec_DTR)
ADD_RBM_TESTCASE(Adaptive_ResEval)
#ADD_RBM_TESTCASE(FLUID_MODES_DNS)
#ADD_RBM_TESTCASE(pagmo_test)
ADD_RBM_TESTCASE(POD_Rec)
ADD_RBM_TESTCASE(ResEval_DtRes)
ADD_RBM_TESTCASE(Stesttest)
#ADD_RBM_TESTCASE(test)
#ADD_RBM_TESTCASE(SingleMODES)

if(ENABLE_PAGMO)
    ADD_RBM_TESTCASE(Adaptive_Sampling)
    ADD_RBM_TESTCASE(AdaptSample_ResEval)
endif()

