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
CMAKE_SOURCE_DIR = /home/subraman/Documents/LeMonADE

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/subraman/Documents/LeMonADE/build

# Include any dependencies generated for this target.
include projects/Examples/example4/CMakeFiles/Example4.dir/depend.make

# Include the progress variables for this target.
include projects/Examples/example4/CMakeFiles/Example4.dir/progress.make

# Include the compile flags for this target's objects.
include projects/Examples/example4/CMakeFiles/Example4.dir/flags.make

projects/Examples/example4/CMakeFiles/Example4.dir/ex4main.cpp.o: projects/Examples/example4/CMakeFiles/Example4.dir/flags.make
projects/Examples/example4/CMakeFiles/Example4.dir/ex4main.cpp.o: ../projects/Examples/example4/ex4main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/subraman/Documents/LeMonADE/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object projects/Examples/example4/CMakeFiles/Example4.dir/ex4main.cpp.o"
	cd /home/subraman/Documents/LeMonADE/build/projects/Examples/example4 && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Example4.dir/ex4main.cpp.o -c /home/subraman/Documents/LeMonADE/projects/Examples/example4/ex4main.cpp

projects/Examples/example4/CMakeFiles/Example4.dir/ex4main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Example4.dir/ex4main.cpp.i"
	cd /home/subraman/Documents/LeMonADE/build/projects/Examples/example4 && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/subraman/Documents/LeMonADE/projects/Examples/example4/ex4main.cpp > CMakeFiles/Example4.dir/ex4main.cpp.i

projects/Examples/example4/CMakeFiles/Example4.dir/ex4main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Example4.dir/ex4main.cpp.s"
	cd /home/subraman/Documents/LeMonADE/build/projects/Examples/example4 && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/subraman/Documents/LeMonADE/projects/Examples/example4/ex4main.cpp -o CMakeFiles/Example4.dir/ex4main.cpp.s

projects/Examples/example4/CMakeFiles/Example4.dir/ex4main.cpp.o.requires:

.PHONY : projects/Examples/example4/CMakeFiles/Example4.dir/ex4main.cpp.o.requires

projects/Examples/example4/CMakeFiles/Example4.dir/ex4main.cpp.o.provides: projects/Examples/example4/CMakeFiles/Example4.dir/ex4main.cpp.o.requires
	$(MAKE) -f projects/Examples/example4/CMakeFiles/Example4.dir/build.make projects/Examples/example4/CMakeFiles/Example4.dir/ex4main.cpp.o.provides.build
.PHONY : projects/Examples/example4/CMakeFiles/Example4.dir/ex4main.cpp.o.provides

projects/Examples/example4/CMakeFiles/Example4.dir/ex4main.cpp.o.provides.build: projects/Examples/example4/CMakeFiles/Example4.dir/ex4main.cpp.o


# Object files for target Example4
Example4_OBJECTS = \
"CMakeFiles/Example4.dir/ex4main.cpp.o"

# External object files for target Example4
Example4_EXTERNAL_OBJECTS =

bin/Example4: projects/Examples/example4/CMakeFiles/Example4.dir/ex4main.cpp.o
bin/Example4: projects/Examples/example4/CMakeFiles/Example4.dir/build.make
bin/Example4: lib/libLeMonADE.a
bin/Example4: projects/Examples/example4/CMakeFiles/Example4.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/subraman/Documents/LeMonADE/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../../bin/Example4"
	cd /home/subraman/Documents/LeMonADE/build/projects/Examples/example4 && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Example4.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
projects/Examples/example4/CMakeFiles/Example4.dir/build: bin/Example4

.PHONY : projects/Examples/example4/CMakeFiles/Example4.dir/build

projects/Examples/example4/CMakeFiles/Example4.dir/requires: projects/Examples/example4/CMakeFiles/Example4.dir/ex4main.cpp.o.requires

.PHONY : projects/Examples/example4/CMakeFiles/Example4.dir/requires

projects/Examples/example4/CMakeFiles/Example4.dir/clean:
	cd /home/subraman/Documents/LeMonADE/build/projects/Examples/example4 && $(CMAKE_COMMAND) -P CMakeFiles/Example4.dir/cmake_clean.cmake
.PHONY : projects/Examples/example4/CMakeFiles/Example4.dir/clean

projects/Examples/example4/CMakeFiles/Example4.dir/depend:
	cd /home/subraman/Documents/LeMonADE/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/subraman/Documents/LeMonADE /home/subraman/Documents/LeMonADE/projects/Examples/example4 /home/subraman/Documents/LeMonADE/build /home/subraman/Documents/LeMonADE/build/projects/Examples/example4 /home/subraman/Documents/LeMonADE/build/projects/Examples/example4/CMakeFiles/Example4.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : projects/Examples/example4/CMakeFiles/Example4.dir/depend

