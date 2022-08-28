# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
get_filename_component(CurrDir ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  
#Create a variable called linalg_SOURCES containing all .cpp files:

set(pcmac1_Src 
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/formfe.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/gama.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/inaccl.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pcinitc.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/phstio.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pmacio.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pmacr1.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pmacr2.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pmacr.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/serchl.f90
    PARENT_SCOPE
    )


