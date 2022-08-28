# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
get_filename_component(CurrDir ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  
#Create a variable called linalg_SOURCES containing all .cpp files:

set(pcelm5_Src 
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/bmat05.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/elmt05.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/matl05.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/modl05.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/tran05.f90
    PARENT_SCOPE
    )


