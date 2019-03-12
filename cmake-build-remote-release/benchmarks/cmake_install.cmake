# Install script for directory: /home/tmhuysen/GQCG/benchmarks

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
  if(EXISTS "$ENV{DESTDIR}/usr/local/gqcp/bin/doci_case" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/local/gqcp/bin/doci_case")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/usr/local/gqcp/bin/doci_case"
         RPATH "/home/tmhuysen/miniconda/lib")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/gqcp/bin/doci_case")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/gqcp/bin" TYPE EXECUTABLE FILES "/home/tmhuysen/GQCG/cmake-build-remote-release/benchmarks/doci_case")
  if(EXISTS "$ENV{DESTDIR}/usr/local/gqcp/bin/doci_case" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/local/gqcp/bin/doci_case")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/usr/local/gqcp/bin/doci_case"
         OLD_RPATH "/home/tmhuysen/GQCG/cmake-build-remote-release/src:/home/tmhuysen/miniconda/lib:"
         NEW_RPATH "/home/tmhuysen/miniconda/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/usr/local/gqcp/bin/doci_case")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/usr/local/gqcp/bin/doci_matrix" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/local/gqcp/bin/doci_matrix")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/usr/local/gqcp/bin/doci_matrix"
         RPATH "/home/tmhuysen/miniconda/lib")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/gqcp/bin/doci_matrix")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/gqcp/bin" TYPE EXECUTABLE FILES "/home/tmhuysen/GQCG/cmake-build-remote-release/benchmarks/doci_matrix")
  if(EXISTS "$ENV{DESTDIR}/usr/local/gqcp/bin/doci_matrix" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/local/gqcp/bin/doci_matrix")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/usr/local/gqcp/bin/doci_matrix"
         OLD_RPATH "/home/tmhuysen/GQCG/cmake-build-remote-release/src:/home/tmhuysen/miniconda/lib:"
         NEW_RPATH "/home/tmhuysen/miniconda/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/usr/local/gqcp/bin/doci_matrix")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/usr/local/gqcp/bin/doci_matvec" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/local/gqcp/bin/doci_matvec")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/usr/local/gqcp/bin/doci_matvec"
         RPATH "/home/tmhuysen/miniconda/lib")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/gqcp/bin/doci_matvec")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/gqcp/bin" TYPE EXECUTABLE FILES "/home/tmhuysen/GQCG/cmake-build-remote-release/benchmarks/doci_matvec")
  if(EXISTS "$ENV{DESTDIR}/usr/local/gqcp/bin/doci_matvec" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/local/gqcp/bin/doci_matvec")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/usr/local/gqcp/bin/doci_matvec"
         OLD_RPATH "/home/tmhuysen/GQCG/cmake-build-remote-release/src:/home/tmhuysen/miniconda/lib:"
         NEW_RPATH "/home/tmhuysen/miniconda/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/usr/local/gqcp/bin/doci_matvec")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/usr/local/gqcp/bin/fci_hchain" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/local/gqcp/bin/fci_hchain")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/usr/local/gqcp/bin/fci_hchain"
         RPATH "/home/tmhuysen/miniconda/lib")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/gqcp/bin/fci_hchain")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/gqcp/bin" TYPE EXECUTABLE FILES "/home/tmhuysen/GQCG/cmake-build-remote-release/benchmarks/fci_hchain")
  if(EXISTS "$ENV{DESTDIR}/usr/local/gqcp/bin/fci_hchain" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/local/gqcp/bin/fci_hchain")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/usr/local/gqcp/bin/fci_hchain"
         OLD_RPATH "/home/tmhuysen/GQCG/cmake-build-remote-release/src:/home/tmhuysen/miniconda/lib:"
         NEW_RPATH "/home/tmhuysen/miniconda/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/usr/local/gqcp/bin/fci_hchain")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/usr/local/gqcp/bin/fci_matrix" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/local/gqcp/bin/fci_matrix")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/usr/local/gqcp/bin/fci_matrix"
         RPATH "/home/tmhuysen/miniconda/lib")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/gqcp/bin/fci_matrix")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/gqcp/bin" TYPE EXECUTABLE FILES "/home/tmhuysen/GQCG/cmake-build-remote-release/benchmarks/fci_matrix")
  if(EXISTS "$ENV{DESTDIR}/usr/local/gqcp/bin/fci_matrix" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/local/gqcp/bin/fci_matrix")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/usr/local/gqcp/bin/fci_matrix"
         OLD_RPATH "/home/tmhuysen/GQCG/cmake-build-remote-release/src:/home/tmhuysen/miniconda/lib:"
         NEW_RPATH "/home/tmhuysen/miniconda/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/usr/local/gqcp/bin/fci_matrix")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/usr/local/gqcp/bin/fci_matvec" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/local/gqcp/bin/fci_matvec")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/usr/local/gqcp/bin/fci_matvec"
         RPATH "/home/tmhuysen/miniconda/lib")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/gqcp/bin/fci_matvec")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/gqcp/bin" TYPE EXECUTABLE FILES "/home/tmhuysen/GQCG/cmake-build-remote-release/benchmarks/fci_matvec")
  if(EXISTS "$ENV{DESTDIR}/usr/local/gqcp/bin/fci_matvec" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/local/gqcp/bin/fci_matvec")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/usr/local/gqcp/bin/fci_matvec"
         OLD_RPATH "/home/tmhuysen/GQCG/cmake-build-remote-release/src:/home/tmhuysen/miniconda/lib:"
         NEW_RPATH "/home/tmhuysen/miniconda/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/usr/local/gqcp/bin/fci_matvec")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/usr/local/gqcp/bin/hubbard_diagonalization" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/local/gqcp/bin/hubbard_diagonalization")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/usr/local/gqcp/bin/hubbard_diagonalization"
         RPATH "/home/tmhuysen/miniconda/lib")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/gqcp/bin/hubbard_diagonalization")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/gqcp/bin" TYPE EXECUTABLE FILES "/home/tmhuysen/GQCG/cmake-build-remote-release/benchmarks/hubbard_diagonalization")
  if(EXISTS "$ENV{DESTDIR}/usr/local/gqcp/bin/hubbard_diagonalization" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/local/gqcp/bin/hubbard_diagonalization")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/usr/local/gqcp/bin/hubbard_diagonalization"
         OLD_RPATH "/home/tmhuysen/GQCG/cmake-build-remote-release/src:/home/tmhuysen/miniconda/lib:"
         NEW_RPATH "/home/tmhuysen/miniconda/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/usr/local/gqcp/bin/hubbard_diagonalization")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/usr/local/gqcp/bin/hubbard_matrix" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/local/gqcp/bin/hubbard_matrix")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/usr/local/gqcp/bin/hubbard_matrix"
         RPATH "/home/tmhuysen/miniconda/lib")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/gqcp/bin/hubbard_matrix")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/gqcp/bin" TYPE EXECUTABLE FILES "/home/tmhuysen/GQCG/cmake-build-remote-release/benchmarks/hubbard_matrix")
  if(EXISTS "$ENV{DESTDIR}/usr/local/gqcp/bin/hubbard_matrix" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/local/gqcp/bin/hubbard_matrix")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/usr/local/gqcp/bin/hubbard_matrix"
         OLD_RPATH "/home/tmhuysen/GQCG/cmake-build-remote-release/src:/home/tmhuysen/miniconda/lib:"
         NEW_RPATH "/home/tmhuysen/miniconda/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/usr/local/gqcp/bin/hubbard_matrix")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/usr/local/gqcp/bin/hubbard_matvec" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/local/gqcp/bin/hubbard_matvec")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/usr/local/gqcp/bin/hubbard_matvec"
         RPATH "/home/tmhuysen/miniconda/lib")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/gqcp/bin/hubbard_matvec")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/gqcp/bin" TYPE EXECUTABLE FILES "/home/tmhuysen/GQCG/cmake-build-remote-release/benchmarks/hubbard_matvec")
  if(EXISTS "$ENV{DESTDIR}/usr/local/gqcp/bin/hubbard_matvec" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/local/gqcp/bin/hubbard_matvec")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/usr/local/gqcp/bin/hubbard_matvec"
         OLD_RPATH "/home/tmhuysen/GQCG/cmake-build-remote-release/src:/home/tmhuysen/miniconda/lib:"
         NEW_RPATH "/home/tmhuysen/miniconda/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/usr/local/gqcp/bin/hubbard_matvec")
    endif()
  endif()
endif()

