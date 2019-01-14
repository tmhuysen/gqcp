# Provide a function that takes the name of target driver, and includes/links the appropriate libraries to it


function(configure_driver DRIVER_NAME)

    configure_executable(${DRIVER_NAME})


    # Install the driver in a separate location
    install(TARGETS ${DRIVER_NAME} RUNTIME DESTINATION ${BIN_INSTALL_DIR})

endfunction(configure_driver)