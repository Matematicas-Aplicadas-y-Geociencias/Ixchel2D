# Version minima de CMake requerido
cmake_minimum_required(VERSION 3.22)

# Nombre del proyecto y lenguaje de programacion a utilizar
project(IXCHEL2D LANGUAGES Fortran)

# Configuracion para compilador
option(GF  "Enable GFORTRAN use"  OFF)
option(NVF "Enable NVFORTRAN use" OFF)

# Configuracion para usar gpu
option(GPU "Enable GPU use" OFF)

# Selecciona el compilador

if( NVF )
  # Especificar el compilador de fortran a utilizar
  set(CMAKE_Fortran_COMPILER nvfortran)

  if(GPU)
    message(STATUS "Using OpenAcc GPU")
    
    # Arquitectura de GPU a utilizar
    set(CUDA_ARCHITECTURE ${TARGET})
    message(STATUS "CUDA Architecture: ${CUDA_ARCHITECTURE}")
    
    # Definir el nombre del ejecutable
    set(EXEC_NAME IXCHEL2D-gpu)
    
    # Banderas para utilizar OpenAcc con GPU
    set(ACC_GPU "-mcmodel=medium -fast -acc=gpu -gpu=${CUDA_ARCHITECTURE} -Minfo=accel")
  else()
    message(STATUS "Using OpenAcc MULTICORE")
    
    # Definir el nombre del ejecutable
    set(EXEC_NAME IXCHEL2D-multicore)
    
    # Banderas para utilizar OpenAcc con Multicore
    set(ACC_MULTICORE "-mcmodel=medium -fast -acc=multicore -Minfo=accel")
    
  endif() # cierra opciones de nvf

elseif( GF )
  
  # Especificar el compilador de fortran a utilizar
  set(CMAKE_Fortran_COMPILER gfortran)
    
  if(GPU)
    message(STATUS "Using OpenMP GPU")
      
    # Arquitectura de GPU a utilizar
    set(CUDA_ARCHITECTURE ${TARGET})
    message(STATUS "CUDA Architecture: ${CUDA_ARCHITECTURE}")
    
    # Definir el nombre del ejecutable
    set(EXEC_NAME IXCHEL2D-gpu)
    
    # Banderas para utilizar OpenAcc con GPU
    set(OMP_GPU "-fopenmp -foffload=nvptx-none -O3")
  else()
    message(STATUS "Using OpenMP MULTICORE")
    
    # Definir el nombre del ejecutable
    set(EXEC_NAME IXCHEL2D-multicore)
    
    # Banderas para utilizar OpenAcc con Multicore
    set(OMP_MULTICORE "-fopenmp -O3")
    
  endif() # cierra opciones de GF
  
  
endif()

# Definir las banderas que va a utilizar el compilador
if( NVF )
  
  if(GPU)
    set(CMAKE_Fortran_FLAGS ${ACC_GPU})
  else()
    set(CMAKE_Fortran_FLAGS ${ACC_MULTICORE})
  endif()

elseif( GF )

  if( GPU )
    set(CMAKE_Fortran_FLAGS ${OMP_GPU})
  else()
    set(CMAKE_Fortran_FLAGS ${OMP_MULTICORE})
  endif()
endif()

# Variable para listar los archivos fuentes
set(SOURCE_FILES
    mod_cond_frontera.f90
    mod_ec_energia.f90
    mod_ec_momento.f90
    mod_front_inmersa.f90
    mod_malla.f90
    mod_postproceso.f90
    mod_residuos.f90
    mod_solucionador.f90
    ixchel2D.f90)

# Creacion del ejecutable
add_executable(${EXEC_NAME} ${SOURCE_FILES})
