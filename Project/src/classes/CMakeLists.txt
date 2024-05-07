list(APPEND DFN_headers_classes ${CMAKE_CURRENT_SOURCE_DIR}/Fracture.hpp)
list(APPEND DFN_sources_classes ${CMAKE_CURRENT_SOURCE_DIR}/Fracture.cpp)
list(APPEND DFN_headers_classes ${CMAKE_CURRENT_SOURCE_DIR}/TracesMesh.hpp)
list(APPEND DFN_headers_classes ${CMAKE_CURRENT_SOURCE_DIR}/PolygonalMesh.hpp)
list(APPEND DFN_includes_classes ${CMAKE_CURRENT_SOURCE_DIR})

set(DFN_sources_classes ${DFN_sources_classes} PARENT_SCOPE)
set(DFN_headers_classes ${DFN_headers_classes} PARENT_SCOPE)
set(DFN_includes_classes ${DFN_includes_classes} PARENT_SCOPE)