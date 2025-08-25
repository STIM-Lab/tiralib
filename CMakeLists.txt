 #Specify the version being used aswell as the language
cmake_minimum_required(VERSION 3.19)

project(tiratest)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(OpenGL REQUIRED)
find_package(glfw3 CONFIG REQUIRED)
find_package(glm CONFIG REQUIRED)
find_package(imgui CONFIG REQUIRED core glfw-binding opengl3-binding)
find_package(GLEW REQUIRED)
find_package(TIRA REQUIRED)
find_package(MKL REQUIRED)

configure_file(test/tiratest.py
				tiratest.py COPYONLY)

file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/toolchain.txt" ${CMAKE_TOOLCHAIN_FILE})