
include_directories(${Boost_INCLUDE_DIR} ${GTSAM_INCLUDE_DIR})
include (GenerateExportHeader) #WIN32 export symbols

file(GLOB subdir_srcs "*.cpp")
file(GLOB subdir_headers "*.h")

add_executable(dynamics ${subdir_srcs} ${subdir_headers})
target_link_libraries(dynamics ${Boost_LIBRARIES} gtsam)

