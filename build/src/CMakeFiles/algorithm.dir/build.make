# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_SOURCE_DIR = /root/Zsj/SelfFile/WorkToSubmit

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /root/Zsj/SelfFile/WorkToSubmit/build

# Include any dependencies generated for this target.
include src/CMakeFiles/algorithm.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/algorithm.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/algorithm.dir/flags.make

src/CMakeFiles/algorithm.dir/Algorithm.cpp.o: src/CMakeFiles/algorithm.dir/flags.make
src/CMakeFiles/algorithm.dir/Algorithm.cpp.o: ../src/Algorithm.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/root/Zsj/SelfFile/WorkToSubmit/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/algorithm.dir/Algorithm.cpp.o"
	cd /root/Zsj/SelfFile/WorkToSubmit/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/algorithm.dir/Algorithm.cpp.o -c /root/Zsj/SelfFile/WorkToSubmit/src/Algorithm.cpp

src/CMakeFiles/algorithm.dir/Algorithm.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/algorithm.dir/Algorithm.cpp.i"
	cd /root/Zsj/SelfFile/WorkToSubmit/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /root/Zsj/SelfFile/WorkToSubmit/src/Algorithm.cpp > CMakeFiles/algorithm.dir/Algorithm.cpp.i

src/CMakeFiles/algorithm.dir/Algorithm.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/algorithm.dir/Algorithm.cpp.s"
	cd /root/Zsj/SelfFile/WorkToSubmit/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /root/Zsj/SelfFile/WorkToSubmit/src/Algorithm.cpp -o CMakeFiles/algorithm.dir/Algorithm.cpp.s

# Object files for target algorithm
algorithm_OBJECTS = \
"CMakeFiles/algorithm.dir/Algorithm.cpp.o"

# External object files for target algorithm
algorithm_EXTERNAL_OBJECTS =

../lib/libalgorithm.so: src/CMakeFiles/algorithm.dir/Algorithm.cpp.o
../lib/libalgorithm.so: src/CMakeFiles/algorithm.dir/build.make
../lib/libalgorithm.so: ../lib/libvecmat.so
../lib/libalgorithm.so: src/CMakeFiles/algorithm.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/root/Zsj/SelfFile/WorkToSubmit/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared library ../../lib/libalgorithm.so"
	cd /root/Zsj/SelfFile/WorkToSubmit/build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/algorithm.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/algorithm.dir/build: ../lib/libalgorithm.so

.PHONY : src/CMakeFiles/algorithm.dir/build

src/CMakeFiles/algorithm.dir/clean:
	cd /root/Zsj/SelfFile/WorkToSubmit/build/src && $(CMAKE_COMMAND) -P CMakeFiles/algorithm.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/algorithm.dir/clean

src/CMakeFiles/algorithm.dir/depend:
	cd /root/Zsj/SelfFile/WorkToSubmit/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /root/Zsj/SelfFile/WorkToSubmit /root/Zsj/SelfFile/WorkToSubmit/src /root/Zsj/SelfFile/WorkToSubmit/build /root/Zsj/SelfFile/WorkToSubmit/build/src /root/Zsj/SelfFile/WorkToSubmit/build/src/CMakeFiles/algorithm.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/algorithm.dir/depend

