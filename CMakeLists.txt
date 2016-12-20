

#find_package(MKL REQUIRED)
message("Pardiso include dir:   " ${PARDISO_INCLUDE_DIR})




ADD_SUBDIRECTORY(1dExamples)
ADD_SUBDIRECTORY(2dExamples)
ADD_SUBDIRECTORY(3dExamples)

#ADD_SUBDIRECTORY(GradientEnhancedDamageProject)



  
