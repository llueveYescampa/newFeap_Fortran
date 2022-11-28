# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
get_filename_component(CurrDir ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  
#Create a variable called linalg_SOURCES containing all .cpp files:

set(pfsolv_Src 
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/delfrt.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pbuff.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pfrtas.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pfrtb.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pfrtbk.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pfrtd.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pfrtf.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pfrtfw.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/prefrt.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/psolve.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/profil.f90
    PARENT_SCOPE
    )


