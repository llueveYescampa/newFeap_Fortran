# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

#Create a variable called linalg_SOURCES containing all .cpp files:

add_subdirectory(pcfeap)
add_subdirectory(pasolv)
add_subdirectory(pcdept)
add_subdirectory(pceign)
add_subdirectory(pcelm1)
add_subdirectory(pcelm2)
add_subdirectory(pcelm3)
add_subdirectory(pcelm4)
add_subdirectory(pcelm5)
add_subdirectory(pcelm6)
add_subdirectory(pcelm7)
add_subdirectory(pcfor5)
add_subdirectory(pcmac1)
add_subdirectory(pcmac2)
add_subdirectory(pcmac3)
add_subdirectory(pcmsh1)
add_subdirectory(pcplt1)
add_subdirectory(pctekt)
add_subdirectory(pfsolv)


set(newFeap_Src 
    ${pcfeap_Src}
    ${pasolv_Src}
    ${pcdept_Src}
    ${pceign_Src}
    ${pcelm1_Src}
    ${pcelm2_Src}
    ${pcelm3_Src}
    ${pcelm4_Src}
    ${pcelm5_Src}
    ${pcelm6_Src}
    ${pcelm7_Src}
    ${pcmac1_Src}
    ${pcmac2_Src}
    ${pcmac3_Src}
    ${pcmsh1_Src}
    ${pcplt1_Src}
    ${pctekt_Src}
    #${pfsolv_Src}
    #${pcfor5_Src}
    PARENT_SCOPE
    )


