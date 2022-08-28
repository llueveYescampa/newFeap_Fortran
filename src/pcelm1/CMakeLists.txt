# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
  get_filename_component(CurrDir ${CMAKE_CURRENT_SOURCE_DIR} NAME)

#Create a variable called linalg_SOURCES containing all .cpp files:

set(pcelm1_Src 
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/elmt01.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/stcn01.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/ster01.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/stre01.f90
    PARENT_SCOPE
    )


