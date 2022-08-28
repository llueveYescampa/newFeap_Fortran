# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
get_filename_component(CurrDir ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  
#Create a variable called linalg_SOURCES containing all .cpp files:

set(pcelm6_Src 
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/coord.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/elmt06.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/flo06.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/stcn06.f90
    PARENT_SCOPE
    )


