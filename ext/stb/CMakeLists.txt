### stb ###                
add_library(stb INTERFACE stb_image.h stb_image_write.h)
target_include_directories(stb INTERFACE .)