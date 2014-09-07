PROJECT(DeformableRegistration2)

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
    INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

ADD_EXECUTABLE(DeformableRegistration2 DeformableRegistration2.cxx)

TARGET_LINK_LIBRARIES(DeformableRegistration2 ${ITK_LIBRARIES})