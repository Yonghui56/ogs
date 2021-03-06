macro(InstallXmlSchemaFiles GLOB_EXPRESSION)
    file(GLOB XSD_FILES . ${GLOB_EXPRESSION})
    if(APPLE AND OGS_BUILD_GUI)
        install(FILES ${XSD_FILES} DESTINATION ${CMAKE_BINARY_DIR}/_CPack_Packages/Darwin/DragNDrop/${CPACK_PACKAGE_FILE_NAME}/ALL_IN_ONE/DataExplorer.app/Contents/MacOS COMPONENT ogs_gui)
    else()
        install(FILES ${XSD_FILES} DESTINATION bin COMPONENT ogs_cli)
    endif()
endmacro()
