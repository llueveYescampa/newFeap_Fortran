# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
get_filename_component(CurrDir ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  
#Create a variable called linalg_SOURCES containing all .cpp files:

set(pcelm2_Src 
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/bmat02.f90  
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/elmt02.f90  
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/elpl02.f90  
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/gvc02.f90  
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/matl02.f90  
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/modl02.f90  
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/stcn02.f90  
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/strn02.f90
    PARENT_SCOPE
    )


