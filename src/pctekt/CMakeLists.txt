# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
get_filename_component(CurrDir ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  
#Create a variable called linalg_SOURCES containing all .cpp files:

set(pctekt_Src 
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/brk4.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/dplot.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/vclrwk.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/vclswk.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/vencur.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/vopnwk.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/vrqstr.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/vslcol.f90
    PARENT_SCOPE
    )


