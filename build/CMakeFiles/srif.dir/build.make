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
include CMakeFiles/srif.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/srif.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/srif.dir/flags.make

CMakeFiles/srif.dir/Srif.cpp.o: CMakeFiles/srif.dir/flags.make
CMakeFiles/srif.dir/Srif.cpp.o: ../Srif.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/root/Zsj/SelfFile/WorkToSubmit/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/srif.dir/Srif.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/srif.dir/Srif.cpp.o -c /root/Zsj/SelfFile/WorkToSubmit/Srif.cpp

CMakeFiles/srif.dir/Srif.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/srif.dir/Srif.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /root/Zsj/SelfFile/WorkToSubmit/Srif.cpp > CMakeFiles/srif.dir/Srif.cpp.i

CMakeFiles/srif.dir/Srif.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/srif.dir/Srif.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /root/Zsj/SelfFile/WorkToSubmit/Srif.cpp -o CMakeFiles/srif.dir/Srif.cpp.s

# Object files for target srif
srif_OBJECTS = \
"CMakeFiles/srif.dir/Srif.cpp.o"

# External object files for target srif
srif_EXTERNAL_OBJECTS =

../bin/srif: CMakeFiles/srif.dir/Srif.cpp.o
../bin/srif: CMakeFiles/srif.dir/build.make
../bin/srif: ../lib/libgnss.so
../bin/srif: ../lib/libalgorithm.so
../bin/srif: ../lib/libvecmat.so
../bin/srif: CMakeFiles/srif.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/root/Zsj/SelfFile/WorkToSubmit/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../bin/srif"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/srif.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/srif.dir/build: ../bin/srif

.PHONY : CMakeFiles/srif.dir/build

CMakeFiles/srif.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/srif.dir/cmake_clean.cmake
.PHONY : CMakeFiles/srif.dir/clean

CMakeFiles/srif.dir/depend:
	cd /root/Zsj/SelfFile/WorkToSubmit/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /root/Zsj/SelfFile/WorkToSubmit /root/Zsj/SelfFile/WorkToSubmit /root/Zsj/SelfFile/WorkToSubmit/build /root/Zsj/SelfFile/WorkToSubmit/build /root/Zsj/SelfFile/WorkToSubmit/build/CMakeFiles/srif.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/srif.dir/depend

