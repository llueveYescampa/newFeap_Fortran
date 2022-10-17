# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
get_filename_component(CurrDir ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  
#Create a variable called linalg_SOURCES containing all .cpp files:

set(pcdept_Src 
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pcdate.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pcomp.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pctime.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pintio.f90
    PARENT_SCOPE
    )

#    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pclear.f90

