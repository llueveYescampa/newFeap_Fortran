# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
get_filename_component(CurrDir ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  
#Create a variable called linalg_SOURCES containing all .cpp files:

set(pceign_Src 
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/chlbac.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/chlfac.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/chlfwd.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/colbac.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/eisql.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/geig.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/numass.f90  
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/scalev.f90  
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/sproja.f90  
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/sprojb.f90  
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/subsp.f90 
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/wprojm.f90
    PARENT_SCOPE
    )


