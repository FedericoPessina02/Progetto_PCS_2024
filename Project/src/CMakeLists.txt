list(APPEND DFN_headers_src ${CMAKE_CURRENT_SOURCE_DIR}/Utils.hpp)
list(APPEND DFN_headers_src ${CMAKE_CURRENT_SOURCE_DIR}/Algorithms.hpp)
list(APPEND DFN_headers_src ${CMAKE_CURRENT_SOURCE_DIR}/Tests.hpp)
list(APPEND DFN_sources_src ${CMAKE_CURRENT_SOURCE_DIR}/Utils.cpp)
list(APPEND DFN_sources_src ${CMAKE_CURRENT_SOURCE_DIR}/Algorithms.cpp)
list(APPEND DFN_includes_src ${CMAKE_CURRENT_SOURCE_DIR})

set(DFN_sources_src ${DFN_sources_src} PARENT_SCOPE)
set(DFN_headers_src ${DFN_headers_src} PARENT_SCOPE)
set(DFN_includes_src ${DFN_includes_src} PARENT_SCOPE)
