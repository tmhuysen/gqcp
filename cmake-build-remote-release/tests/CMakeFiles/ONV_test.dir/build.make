# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/tmhuysen/GQCG

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/tmhuysen/GQCG/cmake-build-remote-release

# Include any dependencies generated for this target.
include tests/CMakeFiles/ONV_test.dir/depend.make

# Include the progress variables for this target.
include tests/CMakeFiles/ONV_test.dir/progress.make

# Include the compile flags for this target's objects.
include tests/CMakeFiles/ONV_test.dir/flags.make

tests/CMakeFiles/ONV_test.dir/FockSpace/ONV_test.cpp.o: tests/CMakeFiles/ONV_test.dir/flags.make
tests/CMakeFiles/ONV_test.dir/FockSpace/ONV_test.cpp.o: ../tests/FockSpace/ONV_test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/tmhuysen/GQCG/cmake-build-remote-release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tests/CMakeFiles/ONV_test.dir/FockSpace/ONV_test.cpp.o"
	cd /home/tmhuysen/GQCG/cmake-build-remote-release/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ONV_test.dir/FockSpace/ONV_test.cpp.o -c /home/tmhuysen/GQCG/tests/FockSpace/ONV_test.cpp

tests/CMakeFiles/ONV_test.dir/FockSpace/ONV_test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ONV_test.dir/FockSpace/ONV_test.cpp.i"
	cd /home/tmhuysen/GQCG/cmake-build-remote-release/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/tmhuysen/GQCG/tests/FockSpace/ONV_test.cpp > CMakeFiles/ONV_test.dir/FockSpace/ONV_test.cpp.i

tests/CMakeFiles/ONV_test.dir/FockSpace/ONV_test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ONV_test.dir/FockSpace/ONV_test.cpp.s"
	cd /home/tmhuysen/GQCG/cmake-build-remote-release/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/tmhuysen/GQCG/tests/FockSpace/ONV_test.cpp -o CMakeFiles/ONV_test.dir/FockSpace/ONV_test.cpp.s

tests/CMakeFiles/ONV_test.dir/FockSpace/ONV_test.cpp.o.requires:

.PHONY : tests/CMakeFiles/ONV_test.dir/FockSpace/ONV_test.cpp.o.requires

tests/CMakeFiles/ONV_test.dir/FockSpace/ONV_test.cpp.o.provides: tests/CMakeFiles/ONV_test.dir/FockSpace/ONV_test.cpp.o.requires
	$(MAKE) -f tests/CMakeFiles/ONV_test.dir/build.make tests/CMakeFiles/ONV_test.dir/FockSpace/ONV_test.cpp.o.provides.build
.PHONY : tests/CMakeFiles/ONV_test.dir/FockSpace/ONV_test.cpp.o.provides

tests/CMakeFiles/ONV_test.dir/FockSpace/ONV_test.cpp.o.provides.build: tests/CMakeFiles/ONV_test.dir/FockSpace/ONV_test.cpp.o


# Object files for target ONV_test
ONV_test_OBJECTS = \
"CMakeFiles/ONV_test.dir/FockSpace/ONV_test.cpp.o"

# External object files for target ONV_test
ONV_test_EXTERNAL_OBJECTS =

tests/ONV_test: tests/CMakeFiles/ONV_test.dir/FockSpace/ONV_test.cpp.o
tests/ONV_test: tests/CMakeFiles/ONV_test.dir/build.make
tests/ONV_test: src/libgqcp.so
tests/ONV_test: /home/tmhuysen/miniconda/lib/libboost_program_options.so
tests/ONV_test: /home/tmhuysen/miniconda/lib/libint2.so
tests/ONV_test: tests/CMakeFiles/ONV_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/tmhuysen/GQCG/cmake-build-remote-release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ONV_test"
	cd /home/tmhuysen/GQCG/cmake-build-remote-release/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ONV_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/CMakeFiles/ONV_test.dir/build: tests/ONV_test

.PHONY : tests/CMakeFiles/ONV_test.dir/build

tests/CMakeFiles/ONV_test.dir/requires: tests/CMakeFiles/ONV_test.dir/FockSpace/ONV_test.cpp.o.requires

.PHONY : tests/CMakeFiles/ONV_test.dir/requires

tests/CMakeFiles/ONV_test.dir/clean:
	cd /home/tmhuysen/GQCG/cmake-build-remote-release/tests && $(CMAKE_COMMAND) -P CMakeFiles/ONV_test.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/ONV_test.dir/clean

tests/CMakeFiles/ONV_test.dir/depend:
	cd /home/tmhuysen/GQCG/cmake-build-remote-release && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/tmhuysen/GQCG /home/tmhuysen/GQCG/tests /home/tmhuysen/GQCG/cmake-build-remote-release /home/tmhuysen/GQCG/cmake-build-remote-release/tests /home/tmhuysen/GQCG/cmake-build-remote-release/tests/CMakeFiles/ONV_test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/CMakeFiles/ONV_test.dir/depend

