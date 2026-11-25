cmake_minimum_required(VERSION 3.16)
project(meosar-rx VERSION 0.1.4 LANGUAGES C)

set(CMAKE_C_STANDARD 11)
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_C_FLAGS_RELEASE "-O2")
option(NATIVE "Optimise for this CPU (-march=native)" ON)
if(NATIVE)
	add_compile_options(-march=native)
endif()
add_compile_options(-Wall -Wextra -Wno-unused-parameter)

find_package(PkgConfig REQUIRED)
find_package(Threads REQUIRED)
pkg_check_modules(FFTW3F REQUIRED fftw3f)
pkg_check_modules(LIBAIRSPY libairspy)
pkg_check_modules(LIBHACKRF libhackrf)

add_library(meosar STATIC src/dsp.c src/fgb.c src/fgb_demod.c src/fgb_repeat.c src/detect.c src/mids.c src/sgb.c src/sgb_rx.c)
target_include_directories(meosar PUBLIC src ${FFTW3F_INCLUDE_DIRS})
target_link_directories(meosar PUBLIC ${FFTW3F_LIBRARY_DIRS})
target_link_libraries(meosar PUBLIC ${FFTW3F_LIBRARIES} m Threads::Threads)

add_executable(meosar-rx src/main.c src/out.c src/input_file.c src/input_live.c)
target_link_libraries(meosar-rx meosar Threads::Threads)
# version.h from git describe, refreshed on every build (cmake/version.cmake)
add_custom_target(version_h
	COMMAND ${CMAKE_COMMAND} -DSRC=${CMAKE_SOURCE_DIR} -DOUT=${CMAKE_BINARY_DIR}/version.h
	        -DFALLBACK=${PROJECT_VERSION} -P ${CMAKE_SOURCE_DIR}/cmake/version.cmake
	BYPRODUCTS ${CMAKE_BINARY_DIR}/version.h)
add_dependencies(meosar-rx version_h)
target_include_directories(meosar-rx PRIVATE ${CMAKE_BINARY_DIR})

add_executable(test_fgb tests/test_fgb.c)
target_link_libraries(test_fgb meosar)
add_executable(test_repeat tests/test_repeat.c)
target_link_libraries(test_repeat meosar)
add_executable(test_sgb tests/test_sgb.c)
target_link_libraries(test_sgb meosar)

foreach(lib AIRSPY HACKRF)
	if(LIB${lib}_FOUND)
		target_compile_definitions(meosar-rx PRIVATE HAVE_${lib})
		target_include_directories(meosar-rx PRIVATE ${LIB${lib}_INCLUDE_DIRS})
		target_link_directories(meosar-rx PRIVATE ${LIB${lib}_LIBRARY_DIRS})
		target_link_libraries(meosar-rx ${LIB${lib}_LIBRARIES})
		message(STATUS "lib${lib} found: live input on")
	endif()
endforeach()

install(TARGETS meosar-rx RUNTIME DESTINATION bin)
