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
CMAKE_SOURCE_DIR = "/mnt/c/Users/MIDUL JACOB/Documents/Projects/set-cover/tools/WiFi-Tools/tt-mjacob/transmission-time/transmission-time"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/mnt/c/Users/MIDUL JACOB/Documents/Projects/set-cover/tools/WiFi-Tools/tt-mjacob/transmission-time/transmission-time/cmake-build-debug"

# Include any dependencies generated for this target.
include CMakeFiles/efficiency.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/efficiency.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/efficiency.dir/flags.make

CMakeFiles/efficiency.dir/efficiency.cpp.o: CMakeFiles/efficiency.dir/flags.make
CMakeFiles/efficiency.dir/efficiency.cpp.o: ../efficiency.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/mnt/c/Users/MIDUL JACOB/Documents/Projects/set-cover/tools/WiFi-Tools/tt-mjacob/transmission-time/transmission-time/cmake-build-debug/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/efficiency.dir/efficiency.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/efficiency.dir/efficiency.cpp.o -c "/mnt/c/Users/MIDUL JACOB/Documents/Projects/set-cover/tools/WiFi-Tools/tt-mjacob/transmission-time/transmission-time/efficiency.cpp"

CMakeFiles/efficiency.dir/efficiency.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/efficiency.dir/efficiency.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/mnt/c/Users/MIDUL JACOB/Documents/Projects/set-cover/tools/WiFi-Tools/tt-mjacob/transmission-time/transmission-time/efficiency.cpp" > CMakeFiles/efficiency.dir/efficiency.cpp.i

CMakeFiles/efficiency.dir/efficiency.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/efficiency.dir/efficiency.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/mnt/c/Users/MIDUL JACOB/Documents/Projects/set-cover/tools/WiFi-Tools/tt-mjacob/transmission-time/transmission-time/efficiency.cpp" -o CMakeFiles/efficiency.dir/efficiency.cpp.s

CMakeFiles/efficiency.dir/efficiency.cpp.o.requires:

.PHONY : CMakeFiles/efficiency.dir/efficiency.cpp.o.requires

CMakeFiles/efficiency.dir/efficiency.cpp.o.provides: CMakeFiles/efficiency.dir/efficiency.cpp.o.requires
	$(MAKE) -f CMakeFiles/efficiency.dir/build.make CMakeFiles/efficiency.dir/efficiency.cpp.o.provides.build
.PHONY : CMakeFiles/efficiency.dir/efficiency.cpp.o.provides

CMakeFiles/efficiency.dir/efficiency.cpp.o.provides.build: CMakeFiles/efficiency.dir/efficiency.cpp.o


# Object files for target efficiency
efficiency_OBJECTS = \
"CMakeFiles/efficiency.dir/efficiency.cpp.o"

# External object files for target efficiency
efficiency_EXTERNAL_OBJECTS =

efficiency: CMakeFiles/efficiency.dir/efficiency.cpp.o
efficiency: CMakeFiles/efficiency.dir/build.make
efficiency: CMakeFiles/efficiency.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/mnt/c/Users/MIDUL JACOB/Documents/Projects/set-cover/tools/WiFi-Tools/tt-mjacob/transmission-time/transmission-time/cmake-build-debug/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable efficiency"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/efficiency.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/efficiency.dir/build: efficiency

.PHONY : CMakeFiles/efficiency.dir/build

CMakeFiles/efficiency.dir/requires: CMakeFiles/efficiency.dir/efficiency.cpp.o.requires

.PHONY : CMakeFiles/efficiency.dir/requires

CMakeFiles/efficiency.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/efficiency.dir/cmake_clean.cmake
.PHONY : CMakeFiles/efficiency.dir/clean

CMakeFiles/efficiency.dir/depend:
	cd "/mnt/c/Users/MIDUL JACOB/Documents/Projects/set-cover/tools/WiFi-Tools/tt-mjacob/transmission-time/transmission-time/cmake-build-debug" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/mnt/c/Users/MIDUL JACOB/Documents/Projects/set-cover/tools/WiFi-Tools/tt-mjacob/transmission-time/transmission-time" "/mnt/c/Users/MIDUL JACOB/Documents/Projects/set-cover/tools/WiFi-Tools/tt-mjacob/transmission-time/transmission-time" "/mnt/c/Users/MIDUL JACOB/Documents/Projects/set-cover/tools/WiFi-Tools/tt-mjacob/transmission-time/transmission-time/cmake-build-debug" "/mnt/c/Users/MIDUL JACOB/Documents/Projects/set-cover/tools/WiFi-Tools/tt-mjacob/transmission-time/transmission-time/cmake-build-debug" "/mnt/c/Users/MIDUL JACOB/Documents/Projects/set-cover/tools/WiFi-Tools/tt-mjacob/transmission-time/transmission-time/cmake-build-debug/CMakeFiles/efficiency.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/efficiency.dir/depend
