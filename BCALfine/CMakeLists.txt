#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(BCALfine) #takes the name of the folder and has to be different from the main executable

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

##### added by Aram to find ROOT ###
###===================================###
## This allows CMake to find your module
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

### Find ROOT
find_package(ROOT)

### Add ROOT header paths
include_directories(${ROOT_INCLUDE_DIR})

### Setup generation of dictionaries
set(MCEvent_SRCS ${PROJECT_SOURCE_DIR}/src/MCEvent.cc)

set(MCEvent_HEADERS ${PROJECT_SOURCE_DIR}/include/MCEvent.hh)

set(MCEvent_LINKDEF ${PROJECT_SOURCE_DIR}/include/MCEventLinkDef.h )
set(MCEvent_DICTIONARY ${PROJECT_SOURCE_DIR}/src/MCEventDict.cc) 

execute_process(COMMAND rootcint -f ${MCEvent_DICTIONARY} -c ${MCEvent_HEADERS} ${MCEvent_LINKDEF})

set(MCEvent_SRCS ${MCEvent_SRCS} ${MCEvent_DICTIONARY}) 
##set(sources ${souces} ${MCEvent_SRCS})   ###sometimes need to coment this out
set(headers ${headers} ${MCEvent_HEADERS})

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(finestructure finestructure.cc ${sources} ${headers})
###target_link_libraries(finestructure ${Geant4_LIBRARIES} )
target_link_libraries(finestructure ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build EXAMPLE. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(MY_SCRIPTS
    run.mac 
    macro/vis.mac 
    energy_spectra/pinacle_6MV.txt
    energy_spectra/pinacle_10MV.txt
    energy_spectra/pinacle_18MV.txt
  )

foreach(_script ${MY_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
add_custom_target(BCALfine DEPENDS finestructure)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS finestructure DESTINATION bin)












# add the dictionary to the list of source files

###SET(MCEvent_SRCS ${MCEvent_SRCS} ${MCEvent_DICTIONARY}) 

####add_library(MCEvent ${MCEvent_SRCS})
###target_link_libraries(MCEvent ${ROOT_LIBRARIES})

#add_custom_target(MCEvent DEPENDS ${MCEvent_SRCS})
#SET(sources ${souces} DEPENDS MCEvent)
