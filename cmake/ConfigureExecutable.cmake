# In this file, we will provide a function that takes the name of an executable, and includes/links the appropriate libraries to it.


function(configure_executable EXECUTABLE_NAME)

    # Include this project
    target_include_directories(${EXECUTABLE_NAME} PRIVATE ${PROJECT_INCLUDE_FOLDER})
    target_link_libraries(${EXECUTABLE_NAME} PRIVATE ${LIBRARY_NAME})


    # Include boost
    target_include_directories(${EXECUTABLE_NAME} PUBLIC ${Boost_INCLUDE_DIRS})
    target_link_libraries(${EXECUTABLE_NAME} PUBLIC ${Boost_LIBRARIES})


    # Include Eigen
    target_link_libraries(${EXECUTABLE_NAME} PUBLIC Eigen3::Eigen)


    # Include cpputil
    target_include_directories(${EXECUTABLE_NAME} PRIVATE ${cpputil_INCLUDE_DIRS})
    target_link_libraries(${EXECUTABLE_NAME} PRIVATE cpputil)


    # Include numopt
    target_include_directories(${EXECUTABLE_NAME} PUBLIC ${numopt_INCLUDE_DIRS})
    target_link_libraries(${EXECUTABLE_NAME} PUBLIC numopt)



endfunction(configure_executable)
