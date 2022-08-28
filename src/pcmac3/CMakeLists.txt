# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
get_filename_component(CurrDir ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  
#Create a variable called linalg_SOURCES containing all .cpp files:

set(pcmac3_Src 
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/acheck.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/ckisop.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/dot.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/elmlib.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pconsd.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pconsi.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pconsr.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pdisk.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pend.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pgauss.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/phelp.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pperror.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/psetm.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pstres.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/saxpb.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/shap2.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/shape.f90
    PARENT_SCOPE
    )


