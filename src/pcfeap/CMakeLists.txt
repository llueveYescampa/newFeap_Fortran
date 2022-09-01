## Include the directory itself as a path to include directories
#set(CMAKE_INCLUDE_CURRENT_DIR ON)
get_filename_component(CurrDir ${CMAKE_CURRENT_SOURCE_DIR} NAME)
#MESSAGE( STATUS "CurrDir:  "            ${CurrDir} )

  
#Create a variable called linalg_SOURCES containing all .cpp files:

set(pcfeap_Src 
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/main.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/filnam.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pcontr.f90  
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pdefil.f90  
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/pltstr.f90
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/prthed.f90
    PARENT_SCOPE
    )


