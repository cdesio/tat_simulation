# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/cdesio/TAT/tat_simu_2steps/photon_simulation

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/cdesio/TAT/tat_simu_2steps/photon_simulation/build

# Utility rule file for AlwaysCheckGit.

# Include any custom commands dependencies for this target.
include CMakeFiles/AlwaysCheckGit.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/AlwaysCheckGit.dir/progress.make

CMakeFiles/AlwaysCheckGit:
	/usr/local/bin/cmake -DRUN_CHECK_GIT_VERSION=1 -Dpre_configure_dir=/home/cdesio/TAT/tat_simu_2steps/photon_simulation/cmake -Dpost_configure_file=/home/cdesio/TAT/tat_simu_2steps/photon_simulation/build/generated -DGIT_HASH_CACHE= -P /home/cdesio/TAT/tat_simu_2steps/photon_simulation/cmake/CheckGit.cmake

AlwaysCheckGit: CMakeFiles/AlwaysCheckGit
AlwaysCheckGit: CMakeFiles/AlwaysCheckGit.dir/build.make
.PHONY : AlwaysCheckGit

# Rule to build all files generated by this target.
CMakeFiles/AlwaysCheckGit.dir/build: AlwaysCheckGit
.PHONY : CMakeFiles/AlwaysCheckGit.dir/build

CMakeFiles/AlwaysCheckGit.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/AlwaysCheckGit.dir/cmake_clean.cmake
.PHONY : CMakeFiles/AlwaysCheckGit.dir/clean

CMakeFiles/AlwaysCheckGit.dir/depend:
	cd /home/cdesio/TAT/tat_simu_2steps/photon_simulation/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/cdesio/TAT/tat_simu_2steps/photon_simulation /home/cdesio/TAT/tat_simu_2steps/photon_simulation /home/cdesio/TAT/tat_simu_2steps/photon_simulation/build /home/cdesio/TAT/tat_simu_2steps/photon_simulation/build /home/cdesio/TAT/tat_simu_2steps/photon_simulation/build/CMakeFiles/AlwaysCheckGit.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/AlwaysCheckGit.dir/depend
