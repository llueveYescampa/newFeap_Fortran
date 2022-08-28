# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
get_filename_component(CurrDir ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  
#Create a variable called linalg_SOURCES containing all .cpp files:

set(pcplt1_Src 
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/frame.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pdefm.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pdevcl.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pdevop.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pline.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/plotl.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pplotf.f90
    PARENT_SCOPE
    )


