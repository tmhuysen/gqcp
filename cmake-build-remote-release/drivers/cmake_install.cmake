# Install script for directory: /home/tmhuysen/GQCG/drivers

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/usr/local/gqcp/bin/fci_lowdin" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/local/gqcp/bin/fci_lowdin")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/usr/local/gqcp/bin/fci_lowdin"
         RPATH "/home/tmhuysen/miniconda/lib")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/gqcp/bin/fci_lowdin")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/gqcp/bin" TYPE EXECUTABLE FILES "/home/tmhuysen/GQCG/cmake-build-remote-release/drivers/fci_lowdin")
  if(EXISTS "$ENV{DESTDIR}/usr/local/gqcp/bin/fci_lowdin" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/local/gqcp/bin/fci_lowdin")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/usr/local/gqcp/bin/fci_lowdin"
         OLD_RPATH "/home/tmhuysen/GQCG/cmake-build-remote-release/src:/home/tmhuysen/miniconda/lib:"
         NEW_RPATH "/home/tmhuysen/miniconda/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/usr/local/gqcp/bin/fci_lowdin")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/usr/local/gqcp/bin/hubbard" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/local/gqcp/bin/hubbard")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/usr/local/gqcp/bin/hubbard"
         RPATH "/home/tmhuysen/miniconda/lib")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/gqcp/bin/hubbard")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/gqcp/bin" TYPE EXECUTABLE FILES "/home/tmhuysen/GQCG/cmake-build-remote-release/drivers/hubbard")
  if(EXISTS "$ENV{DESTDIR}/usr/local/gqcp/bin/hubbard" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/local/gqcp/bin/hubbard")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/usr/local/gqcp/bin/hubbard"
         OLD_RPATH "/home/tmhuysen/GQCG/cmake-build-remote-release/src:/home/tmhuysen/miniconda/lib:"
         NEW_RPATH "/home/tmhuysen/miniconda/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/usr/local/gqcp/bin/hubbard")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/usr/local/gqcp/bin/oo_doci" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/local/gqcp/bin/oo_doci")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/usr/local/gqcp/bin/oo_doci"
         RPATH "/home/tmhuysen/miniconda/lib")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/gqcp/bin/oo_doci")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/gqcp/bin" TYPE EXECUTABLE FILES "/home/tmhuysen/GQCG/cmake-build-remote-release/drivers/oo_doci")
  if(EXISTS "$ENV{DESTDIR}/usr/local/gqcp/bin/oo_doci" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/local/gqcp/bin/oo_doci")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/usr/local/gqcp/bin/oo_doci"
         OLD_RPATH "/home/tmhuysen/GQCG/cmake-build-remote-release/src:/home/tmhuysen/miniconda/lib:"
         NEW_RPATH "/home/tmhuysen/miniconda/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/usr/local/gqcp/bin/oo_doci")
    endif()
  endif()
endif()

