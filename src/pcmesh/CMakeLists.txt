# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
get_filename_component(CurrDir ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  
#Create a variable called linalg_SOURCES containing all .cpp files:

set(pcmsh1_Src 
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/blkgen.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/genvec.f90  
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pbcin.f90  
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pedges.f90  
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pelin.f90  
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pmatin.f90  
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pmesh.f90  
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/polar.f90  
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/sblk.f90  
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/sethis.f90
    PARENT_SCOPE
    )


