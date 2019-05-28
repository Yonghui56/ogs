APPEND_SOURCE_FILES(SOURCES)

add_library(FinesTransport ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS FinesTransport LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_link_libraries(FinesTransport
    PUBLIC ProcessLib
    PRIVATE ParameterLib
)

