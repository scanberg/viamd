### MD4C: Markdown parser (https://github.com/mity/md4c, MIT) ###
# Vendored: md4c.c, md4c.h, LICENSE.md are unmodified copies of upstream.
project(md4c LANGUAGES C)

add_library(md4c STATIC md4c.c md4c.h)

target_include_directories(md4c
	PUBLIC
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# Text is UTF-8 (md4c's default on all platforms except when MD4C_USE_UTF16 is defined, which we do not).
target_compile_features(md4c PRIVATE c_std_99)
