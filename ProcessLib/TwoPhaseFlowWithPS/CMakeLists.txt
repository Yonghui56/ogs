APPEND_SOURCE_FILES(SOURCES)

add_library(TwoPhaseFlowWithPS ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS TwoPhaseFlowWithPS LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_link_libraries(TwoPhaseFlowWithPS
    PUBLIC ProcessLib
    PRIVATE ParameterLib
)

