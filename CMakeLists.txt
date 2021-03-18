cmake_minimum_required(VERSION 3.10.2)
project(skbuild_pycrtm)
enable_language(Fortran)
list( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )
find_package(F2PY REQUIRED)
find_package(PythonLibs REQUIRED)
find_package(Python3 REQUIRED COMPONENTS NumPy)
find_package( OpenMP COMPONENTS Fortran )
find_package( NetCDF REQUIRED COMPONENTS Fortran )
include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${_Python3_NumPy_INCLUDE_DIR})
include_directories(${NetCDF_INCLUDE_DIRS})
#set(crtm_INCLUDE "/discover/nobackup/bkarpowi/github/JCSDA_crtm/crtm-bundle/crtm_v2.4.0/module/crtm/GNU/9.2.0/")
#set(crtm_LIB "/discover/nobackup/bkarpowi/github/JCSDA_crtm/crtm-bundle/crtm_v2.4.0/lib/")
set(crtm_INCLUDE "/discover/nobackup/bkarpowi/github/JCSDA_crtm/REL-2.4.0/crtm_v2.4.0-alpha/include")
set(crtm_LIB "/discover/nobackup/bkarpowi/github/JCSDA_crtm/REL-2.4.0/crtm_v2.4.0-alpha/lib")
set(fortran_src_file "${CMAKE_CURRENT_SOURCE_DIR}/pycrtm.f90")
set(f2py_module_name "pycrtm")
set(generated_module_file ${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}${PYTHON_EXTENSION_MODULE_SUFFIX})
add_compile_options("-lcrtm")
link_libraries("-lcrtm")
add_compile_options("-lnetcdf")
link_libraries("-lnetcdf")
add_compile_options("-lnetcdff")
link_libraries("-lnetcdff")

add_custom_target(${f2py_module_name} ALL
  DEPENDS ${generated_module_file}
  )
set(CDF_PATH "")
get_filename_component(CDF_PATH ${NetCDF_Fortran_LIBRARY} DIRECTORY)

#set(CMAKE_INSTALL_RPATH ${crtm_LIB})
#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
#set(CMAKE_INSTALL_RPATH ${CDF_PATH})
#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
#list(GET NetCDF_LIBRARIES 0 c_cdf)
#list(GET NetCDF_LIBRARIES 2 f_cdf)

add_custom_command(
  OUTPUT ${generated_module_file}
  COMMAND ${F2PY_EXECUTABLE}
    -m ${f2py_module_name}
    -h ${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}.pyf ${fortran_src_file}
    --overwrite-signature
  COMMAND ${F2PY_EXECUTABLE}
    -m ${f2py_module_name}
    -c
    ${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}.pyf 
    ${fortran_src_file}
    -I${crtm_INCLUDE}
    -I${NetCDF_INCLUDE_DIRS}
    -L${CDF_PATH}
    -L${NetCDF_INCLUDE_DIRS}
    -L${crtm_LIB}
    -lcrtm
    -lnetcdff
    -lgomp
    only: 
    wrap_forward 
    wrap_k_matrix
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
list(GET _Python3_INTERPRETER_PROPERTIES 6 version)
set(produced_so ${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}.${version}.so)

install(FILES ${produced_so} DESTINATION pycrtm)

#install(TARGETS pycrtm LIBRARY DESTINATION pycrtm)
