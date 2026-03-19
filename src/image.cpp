#include <image.h>

#include <core/md_log.h>
#include <core/md_os.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

static void write_func(void* context, void* data, int size) {
    ASSERT(context);
    md_file_t* file = (md_file_t*)context;
    md_file_write(*file, data, size);
}

static inline md_file_t open_file(str_t filename) {
    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
        MD_LOG_ERROR("Failed to open file");
    }
    return file;
}

bool image_write_jpg(str_t filename, const void* rgba, int width, int height, int quality) {
    md_file_t file = open_file(filename);
    bool result = false;
    if (md_file_valid(file)) {
        result = stbi_write_jpg_to_func(write_func, &file, width, height, 4, rgba, quality) != 0;
        md_file_close(&file);
    }
    return result;
}

bool image_write_png(str_t filename, const void* rgba, int width, int height) {
    md_file_t file = open_file(filename);
    bool result = false;
    if (md_file_valid(file)) {
        result = stbi_write_png_to_func(write_func, &file, width, height, 4, rgba, width * sizeof(uint32_t)) != 0;
        md_file_close(&file);
    }
    return result;
}

bool image_write_bmp(str_t filename, const void* rgba, int width, int height) {
    md_file_t file = open_file(filename);
    bool result = false;
    if (md_file_valid(file)) {
        result = stbi_write_bmp_to_func(write_func, &file, width, height, 4, rgba) != 0;
        md_file_close(&file);
    }
    return result;
}
