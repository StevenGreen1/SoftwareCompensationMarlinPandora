# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canoncical targets will work.
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
CMAKE_COMMAND = /afs/desy.de/project/ilcsoft/sw/x86_64_gcc48_sl6/CMake/2.8.5/bin/cmake

# The command to remove a file.
RM = /afs/desy.de/project/ilcsoft/sw/x86_64_gcc48_sl6/CMake/2.8.5/bin/cmake -E remove -f

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /afs/desy.de/project/ilcsoft/sw/x86_64_gcc48_sl6/CMake/2.8.5/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /usera/sg568/ilcsoft_v01_17_09/geant4Analysis/Analysis/MultiPlot

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /usera/sg568/ilcsoft_v01_17_09/geant4Analysis/Analysis/MultiPlot/build

# Utility rule file for NightlyTest.

CMakeFiles/NightlyTest:
	/afs/desy.de/project/ilcsoft/sw/x86_64_gcc48_sl6/CMake/2.8.5/bin/ctest -D NightlyTest

NightlyTest: CMakeFiles/NightlyTest
NightlyTest: CMakeFiles/NightlyTest.dir/build.make
.PHONY : NightlyTest

# Rule to build all files generated by this target.
CMakeFiles/NightlyTest.dir/build: NightlyTest
.PHONY : CMakeFiles/NightlyTest.dir/build

CMakeFiles/NightlyTest.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/NightlyTest.dir/cmake_clean.cmake
.PHONY : CMakeFiles/NightlyTest.dir/clean

CMakeFiles/NightlyTest.dir/depend:
	cd /usera/sg568/ilcsoft_v01_17_09/geant4Analysis/Analysis/MultiPlot/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /usera/sg568/ilcsoft_v01_17_09/geant4Analysis/Analysis/MultiPlot /usera/sg568/ilcsoft_v01_17_09/geant4Analysis/Analysis/MultiPlot /usera/sg568/ilcsoft_v01_17_09/geant4Analysis/Analysis/MultiPlot/build /usera/sg568/ilcsoft_v01_17_09/geant4Analysis/Analysis/MultiPlot/build /usera/sg568/ilcsoft_v01_17_09/geant4Analysis/Analysis/MultiPlot/build/CMakeFiles/NightlyTest.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/NightlyTest.dir/depend

