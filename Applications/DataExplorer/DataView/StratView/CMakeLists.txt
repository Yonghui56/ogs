set(SOURCES
    StratBar.cpp
    StratScene.cpp
    StratView.cpp
    StratWindow.cpp
)

set(HEADERS
    StratBar.h
    StratScene.h
    StratView.h
    StratWindow.h
)

set(SOURCE_DIR_REL ${CMAKE_CURRENT_SOURCE_DIR}/../../../..)
set(GUI_SOURCE_DIR_REL ${CMAKE_CURRENT_SOURCE_DIR}/../..)
include_directories(
    ${SOURCE_DIR_REL}/BaseLib
    ${SOURCE_DIR_REL}/GeoLib
    ${SOURCE_DIR_REL}/MathLib
    ${GUI_SOURCE_DIR_REL}/Base
)

file(GLOB_RECURSE UI_FILES *.ui)
source_group("UI Files" FILES ${UI_FILES})

add_library(QtStratView
    ${SOURCES}
    ${HEADERS}
    ${UIS}
)

target_link_libraries(QtStratView Qt4::QtGui QtBase)

set_property(TARGET QtStratView PROPERTY FOLDER "DataExplorer")
