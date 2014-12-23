PROJECT(registerSliceToAtlas_filterVersion)

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
    INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

ADD_EXECUTABLE(registerSliceToAtlas_filterVersion registerSliceToAtlas.cxx)

TARGET_LINK_LIBRARIES(registerSliceToAtlas_filterVersion ${ITK_LIBRARIES})