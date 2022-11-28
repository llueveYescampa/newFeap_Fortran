# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
get_filename_component(CurrDir ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  
#Create a variable called linalg_SOURCES containing all .cpp files:

set(pasolv_Src 
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/dasbly.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/dasol.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/datest.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/datri.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/delfrt.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/dredu.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/profil.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/psolve.f90
    PARENT_SCOPE
    )




