# минимальная версия CMake, которая нужна для этого файла
cmake_minimum_required(VERSION 3.16)

# название проекта и язык (CXX означает C++)
project(fusion_sim CXX)

# используем стандарт C++17 (например, запись 1'000'000 появилась в C++14)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# программа этапа 1: из файла src/main.cpp сделать программу fusion_sim
add_executable(fusion_sim src/main.cpp)

# программа этапа 2: из файла src/tunneling.cpp сделать программу tunneling
add_executable(tunneling src/tunneling.cpp)

# программа этапа 3: из файла src/gamow.cpp сделать программу gamow
add_executable(gamow src/gamow.cpp)

# программа этапа 4: зависимость скорости реакций от температуры
add_executable(rate src/rate.cpp)
