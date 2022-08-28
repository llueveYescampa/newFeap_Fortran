# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
  get_filename_component(CurrDir ${CMAKE_CURRENT_SOURCE_DIR} NAME)

#Create a variable called linalg_SOURCES containing all .cpp files:

set(pcelm4_Src 
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/elmt04.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/modl04.f90
    PARENT_SCOPE
    )


