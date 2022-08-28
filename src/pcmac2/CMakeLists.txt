# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
get_filename_component(CurrDir ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  
#Create a variable called linalg_SOURCES containing all .cpp files:

set(pcmac2_Src 
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/addvec.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/just.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/modify.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pangl.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/param.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pform.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pload.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/propld.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/prtdis.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/prterr.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/prtrea.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/prtstr.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/ptrans.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/setci.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/update.f90
    PARENT_SCOPE
    )


