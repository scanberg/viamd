//#define NOMINMAX

#include "ramachandran.h"
#include <core/gl.h>
#include <core/log.h>
#include <gfx/gl_utils.h>

namespace ramachandran {

// Segmentation texture data
constexpr int seg_width = 36;
constexpr int seg_height = 36;
constexpr GLenum seg_data_format = GL_RGB;
constexpr unsigned int seg_data[] = {
    0xCC997ACC, 0x5500CC55, 0x00CC5500, 0xCC5500CC, 0x5500CC55, 0x00CC5500, 0xCC5500CC, 0x5500CC55, 0x00CC997A, 0xCC997ACC, 0x997AD8BF, 0xADD8BFAD,
    0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFC8, 0xCCB7C8CC, 0xB7C8CCB7, 0xC8CCB7C8, 0xCCB7FFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFD8BFAD,
    0xD8BFADD8, 0xBFADD8BF, 0xADCC997A, 0xCC997ACC, 0x5500CC55, 0x00CC5500, 0xCC5500CC, 0x5500CC55, 0x00CC5500, 0xCC5500CC, 0x5500CC55, 0x00CC5500,
    0xCC997ACC, 0x997ACC99, 0x7AD8BFAD, 0xD8BFADFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF,
    0xFFFFFFFF, 0xFFFFFFFF, 0xFFD8BFAD, 0xD8BFADD8, 0xBFADCC99, 0x7ACC997A, 0xCC997ACC, 0x5500CC55, 0x00CC5500, 0xCC5500CC, 0x5500CC55, 0x00CC5500,
    0xCC5500CC, 0x5500CC55, 0x00CC5500, 0xCC997ACC, 0x997ACC99, 0x7AD8BFAD, 0xD8BFADD8, 0xBFADD8BF, 0xADFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF,
    0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFD8BFAD, 0xD8BFADCC, 0x997ACC99, 0x7ACC997A, 0xCC997ACC, 0x997ACC55, 0x00CC5500,
    0xCC5500CC, 0x5500CC55, 0x00CC5500, 0xCC5500CC, 0x5500CC55, 0x00CC5500, 0xCC5500CC, 0x997ACC99, 0x7ACC997A, 0xCC997AD8, 0xBFADD8BF, 0xADFFFFFF,
    0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFD8BFAD, 0xD8BFADD8, 0xBFADCC99, 0x7ACC997A,
    0xCC997ACC, 0x997ACC55, 0x00CC5500, 0xCC5500CC, 0x5500CC55, 0x00CC5500, 0xCC5500CC, 0x5500CC55, 0x00CC5500, 0xCC5500CC, 0x5500CC99, 0x7ACC997A,
    0xD8BFADD8, 0xBFADD8BF, 0xADFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFD8BFAD,
    0xD8BFADD8, 0xBFADD8BF, 0xADCC997A, 0xCC997ACC, 0x997ACC55, 0x00CC5500, 0xCC5500CC, 0x5500CC55, 0x00CC5500, 0xCC5500CC, 0x5500CC55, 0x00CC5500,
    0xCC5500CC, 0x997ACC99, 0x7ACC997A, 0xD8BFADD8, 0xBFADD8BF, 0xADFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF,
    0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xD8BFADD8, 0xBFADD8BF, 0xADCC997A, 0xCC997ACC, 0x997ACC99, 0x7ACC5500, 0xCC5500CC, 0x5500CC55, 0x00CC5500,
    0xCC5500CC, 0x5500CC55, 0x00CC5500, 0xCC997ACC, 0x997ACC99, 0x7ACC997A, 0xD8BFADD8, 0xBFADD8BF, 0xADFFFFFF, 0xFFFFFFFF, 0xFFFFB7C2, 0xCCB7C2CC,
    0xB7C2CCB7, 0xC2CCB7C2, 0xCCB7C2CC, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFD8, 0xBFADD8BF, 0xADCC997A, 0xCC997ACC, 0x997ACC99, 0x7ACC997A,
    0xCC5500CC, 0x5500CC55, 0x00CC5500, 0xCC5500CC, 0x5500CC55, 0x00CC997A, 0xCC997ACC, 0x997ACC99, 0x7ACC997A, 0xCC997AD8, 0xBFADD8BF, 0xADB7C2CC,
    0xB7C2CCB7, 0xC2CCB7C2, 0xCCB7C2CC, 0xB7C2CCB7, 0xC2CCB7C2, 0xCCB7C2CC, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFD8, 0xBFADD8BF, 0xADCC997A,
    0xCC997ACC, 0x997ACC99, 0x7ACC997A, 0xCC997ACC, 0x5500CC99, 0x7ACC997A, 0xCC5500CC, 0x997ACC99, 0x7ACC997A, 0xCC997ACC, 0x997ACC99, 0x7ACC997A,
    0xCC997AD8, 0xBFADD8BF, 0xADB7C2CC, 0xB7C2CCB7, 0xC2CCB7C2, 0xCCB7C2CC, 0x669ECC66, 0x9ECCB7C2, 0xCCB7C2CC, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF,
    0xFFFFFFD8, 0xBFADD8BF, 0xADCC997A, 0xCC997ACC, 0x997ACC99, 0x7ACC997A, 0xCC997ACC, 0x997ACC99, 0x7ACC997A, 0xCC997ACC, 0x997ACC99, 0x7ACC997A,
    0xCC997ACC, 0x997AD8BF, 0xADD8BFAD, 0xD8BFADD8, 0xBFADD8BF, 0xADB7C2CC, 0xB7C2CC66, 0x9ECC669E, 0xCC669ECC, 0x669ECCB7, 0xC2CCB7C2, 0xCCB7C2CC,
    0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFD8, 0xBFADD8BF, 0xADD8BFAD, 0xCC997ACC, 0x997ACC99, 0x7ACC997A, 0xCC997ACC, 0x997ACC99, 0x7ACC997A,
    0xCC997ACC, 0x997ACC99, 0x7ACC997A, 0xD8BFADD8, 0xBFADD8BF, 0xADD8BFAD, 0xD8BFADD8, 0xBFADD8BF, 0xADB7C2CC, 0xB7C2CC66, 0x9ECC669E, 0xCC669ECC,
    0x669ECCB7, 0xC2CCB7C2, 0xCCB7C2CC, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFD8, 0xBFADD8BF, 0xADD8BFAD, 0xD8BFADCC, 0x997ACC99, 0x7ACC997A,
    0xCC997ACC, 0x997ACC99, 0x7ACC997A, 0xCC997ACC, 0x997ACC99, 0x7ACC997A, 0xCC997AD8, 0xBFADD8BF, 0xADFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFB7C2CC,
    0xB7C2CCB7, 0xC2CC669E, 0xCC669ECC, 0x669ECCB7, 0xC2CCB7C2, 0xCCB7C2CC, 0xB7C2CCFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFD8, 0xBFADD8BF, 0xADD8BFAD,
    0xD8BFADD8, 0xBFADCC99, 0x7ACC997A, 0xCC997ACC, 0x997ACC99, 0x7ACC997A, 0xCC997ACC, 0x997ACC99, 0x7ACC997A, 0xCC997AD8, 0xBFADD8BF, 0xADFFFFFF,
    0xFFFFFFFF, 0xFFFFFFFF, 0xFFB7C2CC, 0xB7C2CC66, 0x9ECC669E, 0xCC669ECC, 0x669ECC66, 0x9ECC669E, 0xCCB7C2CC, 0xB7C2CCFF, 0xFFFFFFFF, 0xFFFFFFFF,
    0xFFFFFFD8, 0xBFADD8BF, 0xADD8BFAD, 0xB7CCBEB7, 0xCCBE7ACC, 0x957ACC95, 0x7ACC957A, 0xCC957ACC, 0x957ACC95, 0x7ACC957A, 0xCC957ACC, 0x957ACC95,
    0xB7CCBEB7, 0xCCBEB7CC, 0xBEFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFB7C2CC, 0xB7C2CC66, 0x9ECC669E, 0xCC0070CC, 0x669ECC66, 0x9ECCB7C2, 0xCCB7C2CC,
    0xB7C2CCFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFB7, 0xCCBEB7CC, 0xBEB7CCBE, 0xB7CCBE7A, 0xCC957ACC, 0x957ACC95, 0x7ACC957A, 0xCC957ACC, 0x957ACC95,
    0x7ACC957A, 0xCC957ACC, 0x957ACC95, 0xB7CCBEB7, 0xCCBEB7CC, 0xBEFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFB7C2CC, 0xB7C2CCB7, 0xC2CC669E, 0xCC0070CC,
    0x0070CC66, 0x9ECCB7C2, 0xCCB7C2CC, 0xB7C2CCFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFB7, 0xCCBEB7CC, 0xBEB7CCBE, 0xB7CCBE7A, 0xCC957ACC, 0x957ACC95,
    0x7ACC957A, 0xCC957ACC, 0x9500CC44, 0x7ACC957A, 0xCC957ACC, 0x957ACC95, 0x7ACC95B7, 0xCCBEB7CC, 0xBEB7CCBE, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFB7C2CC,
    0xB7C2CCB7, 0xC2CC669E, 0xCC669ECC, 0x669ECC66, 0x9ECC669E, 0xCCB7C2CC, 0xB7C2CCFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFB7, 0xCCBEB7CC, 0xBEB7CCBE,
    0xB7CCBE7A, 0xCC957ACC, 0x957ACC95, 0x7ACC957A, 0xCC9500CC, 0x4400CC44, 0x00CC4400, 0xCC447ACC, 0x957ACC95, 0x7ACC95B7, 0xCCBEB7CC, 0xBEB7CCBE,
    0xB7CCBEFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xB7C2CCB7, 0xC2CCB7C2, 0xCC669ECC, 0x669ECC66, 0x9ECCB7C2, 0xCCB7C2CC, 0xB7C2CCFF, 0xFFFFFFFF, 0xFFFFFFFF,
    0xFFFFFFB7, 0xCCBEB7CC, 0xBEB7CCBE, 0xB7CCBE7A, 0xCC957ACC, 0x957ACC95, 0x7ACC957A, 0xCC9500CC, 0x4400CC44, 0x00CC4400, 0xCC4400CC, 0x447ACC95,
    0x7ACC957A, 0xCC95B7CC, 0xBEB7CCBE, 0xB7CCBEFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xB7C2CCB7, 0xC2CCB7C2, 0xCCB7C2CC, 0x669ECC66, 0x9ECC669E, 0xCCB7C2CC,
    0xB7C2CCFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFB7, 0xCCBEB7CC, 0xBEB7CCBE, 0xB7CCBE7A, 0xCC957ACC, 0x957ACC95, 0x7ACC9500, 0xCC4400CC, 0x4400CC44,
    0x00CC4400, 0xCC4400CC, 0x4400CC44, 0x7ACC957A, 0xCC957ACC, 0x95B7CCBE, 0xB7CCBEB7, 0xCCBEFFFF, 0xFFFFFFFF, 0xB7C2CCB7, 0xC2CCB7C2, 0xCC669ECC,
    0xB7C2CC66, 0x9ECC669E, 0xCCB7C2CC, 0xB7C2CCFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFB7, 0xCCBEB7CC, 0xBEB7CCBE, 0xB7CCBE7A, 0xCC957ACC, 0x957ACC95,
    0x7ACC957A, 0xCC9500CC, 0x4400CC44, 0x00CC4400, 0xCC4400CC, 0x4400CC44, 0x00CC447A, 0xCC957ACC, 0x95B7CCBE, 0xB7CCBEB7, 0xCCBEB7CC, 0xBEFFFFFF,
    0xB7C2CCB7, 0xC2CCB7C2, 0xCCB7C2CC, 0xB7C2CCB7, 0xC2CCB7C2, 0xCCB7C2CC, 0xB7C2CCFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFB7, 0xCCBEB7CC, 0xBEB7CCBE,
    0xB7CCBEB7, 0xCCBE7ACC, 0x957ACC95, 0x7ACC957A, 0xCC9500CC, 0x4400CC44, 0x00CC4400, 0xCC4400CC, 0x4400CC44, 0x00CC4400, 0xCC447ACC, 0x957ACC95,
    0xB7CCBEB7, 0xCCBEB7CC, 0xBEFFFFFF, 0xB7C2CCB7, 0xC2CCB7C2, 0xCCB7C2CC, 0xB7C2CCB7, 0xC2CCB7C2, 0xCCB7C2CC, 0xB7C2CCFF, 0xFFFFFFFF, 0xFFFFFFFF,
    0xFFFFFFB7, 0xCCBEB7CC, 0xBEB7CCBE, 0xB7CCBEB7, 0xCCBE7ACC, 0x957ACC95, 0x7ACC957A, 0xCC957ACC, 0x9500CC44, 0x00CC4400, 0xCC4400CC, 0x4400CC44,
    0x00CC4400, 0xCC447ACC, 0x957ACC95, 0x7ACC95B7, 0xCCBEB7CC, 0xBEB7CCBE, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF,
    0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFB7, 0xCCBEB7CC, 0xBEB7CCBE, 0xB7CCBEB7, 0xCCBEB7CC, 0xBE7ACC95, 0x7ACC957A, 0xCC957ACC, 0x957ACC95,
    0x00CC4400, 0xCC4400CC, 0x4400CC44, 0x00CC4400, 0xCC4400CC, 0x447ACC95, 0x7ACC95B7, 0xCCBEB7CC, 0xBEB7CCBE, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF,
    0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFB7, 0xCCBEB7CC, 0xBEB7CCBE, 0xB7CCBE7A, 0xCC957ACC, 0x957ACC95,
    0x7ACC957A, 0xCC957ACC, 0x957ACC95, 0x7ACC9500, 0xCC4400CC, 0x4400CC44, 0x00CC4400, 0xCC4400CC, 0x447ACC95, 0x7ACC957A, 0xCC95B7CC, 0xBEB7CCBE,
    0xB7CCBEFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFB7, 0xCCBEB7CC, 0xBEB7CCBE,
    0x7ACC957A, 0xCC95B7CC, 0xBE7ACC95, 0x7ACC957A, 0xCC957ACC, 0x957ACC95, 0x7ACC957A, 0xCC957ACC, 0x9500CC44, 0x00CC4400, 0xCC4400CC, 0x447ACC95,
    0x7ACC957A, 0xCC95B7CC, 0xBEB7CCBE, 0xB7CCBEFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF,
    0xFFFFFFB7, 0xCCBEB7CC, 0xBEB7CCBE, 0xB7CCBEB7, 0xCCBEB7CC, 0xBE7ACC95, 0x7ACC957A, 0xCC957ACC, 0x957ACC95, 0x7ACC957A, 0xCC957ACC, 0x957ACC95,
    0x7ACC957A, 0xCC957ACC, 0x957ACC95, 0x7ACC957A, 0xCC957ACC, 0x95B7CCBE, 0xB7CCBEFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF,
    0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFB7, 0xCCBEB7CC, 0xBEB7CCBE, 0xB7CCBEB7, 0xCCBEB7CC, 0xBEB7CCBE, 0xB7CCBE7A, 0xCC957ACC, 0x957ACC95,
    0x7ACC957A, 0xCC957ACC, 0x957ACC95, 0x7ACC957A, 0xCC957ACC, 0x957ACC95, 0x7ACC957A, 0xCC95B7CC, 0xBEB7CCBE, 0xB7CCBEFF, 0xFFFFFFFF, 0xFFFFFFFF,
    0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFB7, 0xCCBEB7CC, 0xBEB7CCBE, 0xD8BFADB7, 0xCCBEB7CC, 0xBEB7CCBE,
    0xB7CCBE7A, 0xCC957ACC, 0x957ACC95, 0x7ACC957A, 0xCC957ACC, 0x957ACC95, 0x7ACC957A, 0xCC957ACC, 0x95B7CCBE, 0xB7CCBEB7, 0xCCBEB7CC, 0xBEB7CCBE,
    0xB7CCBEFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF,
    0xD8BFADD8, 0xBFADD8BF, 0xADB7CCBE, 0xB7CCBEB7, 0xCCBED8BF, 0xADD8BFAD, 0xD8BFADD8, 0xBFADB7CC, 0xBE7ACC95, 0xB7CCBEB7, 0xCCBEB7CC, 0xBEB7CCBE,
    0xB7CCBEB7, 0xCCBEB7CC, 0xBEB7CCBE, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF,
    0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xD8BFADD8, 0xBFADD8BF, 0xADD8BFAD, 0xD8BFADD8, 0xBFADD8BF, 0xADCC997A, 0xCC997AD8, 0xBFADB7CC, 0xBEB7CCBE,
    0xB7CCBEB7, 0xCCBEB7CC, 0xBEB7CCBE, 0xB7CCBEFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF,
    0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xD8BFADD8, 0xBFADD8BF, 0xADD8BFAD, 0xD8BFADCC, 0x997ACC99, 0x7AD8BFAD,
    0xD8BFADD8, 0xBFADD8BF, 0xADD8BFAD, 0xD8BFADD8, 0xBFADFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFC8, 0xCCB7C8CC, 0xB7C8CCB7,
    0xC8CCB7C8, 0xCCB7FFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xD8BFADD8, 0xBFADD8BF, 0xADD8BFAD,
    0xD8BFADD8, 0xBFADCC99, 0x7ACC997A, 0xCC997ACC, 0x997AD8BF, 0xADD8BFAD, 0xD8BFADD8, 0xBFADFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF,
    0xFFFFFFC8, 0xCCB7C8CC, 0xB7C8CCB7, 0xC8CCB7C8, 0xCCB7FFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF,
    0xD8BFADCC, 0x997ACC99, 0x7ACC997A, 0xD8BFADCC, 0x997ACC99, 0x7ACC997A, 0xCC997ACC, 0x997AD8BF, 0xADD8BFAD, 0xD8BFADD8, 0xBFADD8BF, 0xADFFFFFF,
    0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFC8, 0xCCB7C8CC, 0xB7B7CC51, 0xC8CCB7C8, 0xCCB7C8CC, 0xB7FFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF,
    0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xD8BFADCC, 0x997ACC99, 0x7ACC997A, 0xCC997ACC, 0x997ACC99, 0x7ACC997A, 0xCC997ACC, 0x997ACC99, 0x7AD8BFAD,
    0xD8BFADD8, 0xBFADD8BF, 0xADFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFC8, 0xCCB7C8CC, 0xB7B7CC51, 0xC8CCB7C8, 0xCCB7C8CC, 0xB7FFFFFF,
    0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFD8, 0xBFADD8BF, 0xADD8BFAD, 0xD8BFADCC, 0x997ACC99, 0x7ACC997A, 0xCC997ACC, 0x997ACC99, 0x7ACC997A,
    0xCC997ACC, 0x997ACC99, 0x7ACC997A, 0xCC997AD8, 0xBFADD8BF, 0xADFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFC8, 0xCCB7C8CC, 0xB7B7CC51,
    0xB7CC51C8, 0xCCB7C8CC, 0xB7FFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFD8, 0xBFADD8BF, 0xADD8BFAD, 0xCC997ACC, 0x997ACC99, 0x7ACC997A,
    0xCC997ACC, 0x997ACC99, 0x7ACC997A, 0xCC997ACC, 0x997ACC99, 0x7ACC997A, 0xCC997AD8, 0xBFADD8BF, 0xADFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF,
    0xFFFFFFC8, 0xCCB7C8CC, 0xB7C8CCB7, 0xC8CCB7C8, 0xCCB7C8CC, 0xB7FFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFD8, 0xBFADD8BF, 0xADCC997A};

// Accumulation texture data
constexpr int acc_width = 1024;
constexpr int acc_height = 1024;

static GLuint segmentation_tex = 0;
static GLuint accumulation_tex = 0;
static GLuint coord_tex = 0;
static GLuint coord_buf = 0;
static GLuint fbo = 0;
static GLuint program = 0;

static GLuint vao = 0;
static GLuint ibo = 0;
static GLuint vbo = 0;

static GLint uniform_loc_coord_tex = -1;
static GLint uniform_loc_instance_offset = -1;
static GLint uniform_loc_radius = -1;
static GLint uniform_loc_color = -1;
static GLint uniform_loc_outline = -1;

GLuint get_accumulation_texture() { return accumulation_tex; }
GLuint get_segmentation_texture() { return segmentation_tex; }

// @NOTE: This should generate a quad with a certain size in texture coordinates
constexpr const char* v_shader_src = R"(
#version 150 core

uniform int u_instance_offset = 0;
uniform samplerBuffer u_tex_coord;
uniform float u_radius;
out vec2 uv;

void main() {
	int VID = gl_VertexID;
	int IID = gl_InstanceID + u_instance_offset;

	vec2 coord = texelFetch(u_tex_coord, IID).xy;
	uv = vec2(VID / 2, VID % 2) * 2.0 - 1.0; 

	gl_Position = vec4(coord * 2.0 - 1.0 + uv * u_radius, 0, 1);
}
)";

// @NOTE: Do some radial falloff based on uv coordinate
constexpr const char* f_shader_src = R"(
#version 150 core

uniform vec4 u_color;
uniform float u_outline;
in vec2 uv;
out vec4 out_frag;

float step_dist(float edge, float dist) {
	float factor = 1.0; // <-- value can be played around with a bit
	float mask = step(edge, dist);
	float step_w = factor * length(vec2(dFdx(mask), dFdy(mask)));
	return smoothstep(-step_w/2.0, step_w/2.0, mask);
}

void main() {
	float dist = sqrt(dot(uv, uv));
	if (u_outline > 0) {
		vec4 inner_color = u_color.rgba;
		vec4 rim_color = vec4(0,0,0,1);
		vec4 outer_color = vec4(0,0,0,0);
		out_frag = mix(inner_color, mix(rim_color, outer_color, step_dist(1.0, dist)), step_dist(1.0 - u_outline, dist));
	} else {
		float falloff = max(0, 1.0 - dist);
		out_frag = vec4(u_color.rgb, u_color.a * falloff);	
	}
}
)";

void initialize() {
    constexpr int BUFFER_SIZE = 1024;
    char buffer[BUFFER_SIZE];

    GLuint v_shader = glCreateShader(GL_VERTEX_SHADER);
    GLuint f_shader = glCreateShader(GL_FRAGMENT_SHADER);
    defer {
        glDeleteShader(v_shader);
        glDeleteShader(f_shader);
    };
    glShaderSource(v_shader, 1, &v_shader_src, 0);
    glShaderSource(f_shader, 1, &f_shader_src, 0);

    glCompileShader(v_shader);
    if (gl::get_shader_compile_error(buffer, BUFFER_SIZE, v_shader)) {
        LOG_ERROR("Compiling ramachandran vertex shader:\n%s\n", buffer);
    }
    glCompileShader(f_shader);
    if (gl::get_shader_compile_error(buffer, BUFFER_SIZE, f_shader)) {
        LOG_ERROR("Compiling ramachandran fragment shader:\n%s\n", buffer);
    }

    program = glCreateProgram();
    glAttachShader(program, v_shader);
    glAttachShader(program, f_shader);
    glLinkProgram(program);
    if (gl::get_program_link_error(buffer, BUFFER_SIZE, program)) {
        LOG_ERROR("Linking ramachandran program:\n%s\n", buffer);
    }

    glDetachShader(program, v_shader);
    glDetachShader(program, f_shader);

    uniform_loc_coord_tex = glGetUniformLocation(program, "u_coord_tex");
    uniform_loc_instance_offset = glGetUniformLocation(program, "u_instance_offset");
    uniform_loc_radius = glGetUniformLocation(program, "u_radius");
    uniform_loc_color = glGetUniformLocation(program, "u_color");
    uniform_loc_outline = glGetUniformLocation(program, "u_outline");

    if (!segmentation_tex) {
        glGenTextures(1, &segmentation_tex);
        glBindTexture(GL_TEXTURE_2D, segmentation_tex);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, seg_width, seg_height, 0, seg_data_format, GL_UNSIGNED_BYTE, seg_data);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glBindTexture(GL_TEXTURE_2D, 0);
    }

    if (!accumulation_tex) {
        glGenTextures(1, &accumulation_tex);
        glBindTexture(GL_TEXTURE_2D, accumulation_tex);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, acc_width, acc_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glBindTexture(GL_TEXTURE_2D, 0);
    }

    if (!coord_buf) {
        glGenBuffers(1, &coord_buf);
    }

    if (!coord_tex) {
        glGenTextures(1, &coord_tex);
    }

    if (!fbo) {
        glGenFramebuffers(1, &fbo);
        glBindFramebuffer(GL_FRAMEBUFFER, fbo);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, accumulation_tex, 0);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
    }

    if (!vbo) {
        glGenBuffers(1, &vbo);
    }

    if (!ibo) {
        constexpr unsigned char data[4] = {0, 1, 2, 3};
        glGenBuffers(1, &ibo);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, 4, data, GL_STATIC_DRAW);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    }
}

void shutdown() {
    if (segmentation_tex) glDeleteTextures(1, &segmentation_tex);
    if (accumulation_tex) glDeleteTextures(1, &accumulation_tex);
    if (coord_buf) glDeleteBuffers(1, &coord_buf);
    if (coord_tex) glDeleteTextures(1, &coord_tex);
    if (fbo) glDeleteFramebuffers(1, &fbo);
    if (vbo) glDeleteBuffers(1, &vbo);
}

void clear_accumulation_texture() {
    GLint last_viewport[4];
    glGetIntegerv(GL_VIEWPORT, last_viewport);

    glViewport(0, 0, acc_width, acc_height);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fbo);
    glDrawBuffer(GL_COLOR_ATTACHMENT0);

    glClearColor(0, 0, 0, 0);
    glClear(GL_COLOR_BUFFER_BIT);

    glViewport(last_viewport[0], last_viewport[1], (GLsizei)last_viewport[2], (GLsizei)last_viewport[3]);
}

void compute_accumulation_texture(Array<const BackboneAngles> angles, vec4 color, float radius, float outline) {
    constexpr float ONE_OVER_TWO_PI = 1.f / (2.f * math::PI);

    struct Coord {
        unsigned short x, y;
    };

    // Use fast scratch memory here
    Coord* coords = (Coord*)TMP_MALLOC((angles.count) * sizeof(Coord));
    defer { TMP_FREE(coords); };

    int32 count = 0;
    for (const auto& angle : angles) {
        if (angle.phi == 0 || angle.psi == 0) continue;
        vec2 coord = vec2(angle.phi, angle.psi) * ONE_OVER_TWO_PI + 0.5f;  // [-PI, PI] -> [0, 1]
        coord.y = 1.f - coord.y;
        coords[count].x = (unsigned short)(coord.x * 0xffff);
        coords[count].y = (unsigned short)(coord.y * 0xffff);
        count++;
    }

    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, count * sizeof(Coord), coords, GL_STREAM_DRAW);

    // Backup GL state
    GLint last_polygon_mode[2];
    glGetIntegerv(GL_POLYGON_MODE, last_polygon_mode);
    GLint last_viewport[4];
    glGetIntegerv(GL_VIEWPORT, last_viewport);
    GLenum last_blend_src_rgb;
    glGetIntegerv(GL_BLEND_SRC_RGB, (GLint*)&last_blend_src_rgb);
    GLenum last_blend_dst_rgb;
    glGetIntegerv(GL_BLEND_DST_RGB, (GLint*)&last_blend_dst_rgb);
    GLenum last_blend_src_alpha;
    glGetIntegerv(GL_BLEND_SRC_ALPHA, (GLint*)&last_blend_src_alpha);
    GLenum last_blend_dst_alpha;
    glGetIntegerv(GL_BLEND_DST_ALPHA, (GLint*)&last_blend_dst_alpha);
    GLenum last_blend_equation_rgb;
    glGetIntegerv(GL_BLEND_EQUATION_RGB, (GLint*)&last_blend_equation_rgb);
    GLenum last_blend_equation_alpha;
    glGetIntegerv(GL_BLEND_EQUATION_ALPHA, (GLint*)&last_blend_equation_alpha);
    GLboolean last_enable_blend = glIsEnabled(GL_BLEND);
    GLboolean last_enable_cull_face = glIsEnabled(GL_CULL_FACE);
    GLboolean last_enable_depth_test = glIsEnabled(GL_DEPTH_TEST);
    GLboolean last_enable_scissor_test = glIsEnabled(GL_SCISSOR_TEST);

    // RENDER TO ACCUMULATION TEXTURE

    glViewport(0, 0, acc_width, acc_height);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fbo);
    glDrawBuffer(GL_COLOR_ATTACHMENT0);

    glEnable(GL_BLEND);
    glBlendEquation(GL_FUNC_ADD);
    glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
    glDisable(GL_CULL_FACE);
    glDisable(GL_DEPTH_TEST);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    // Texture 0
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_BUFFER, coord_tex);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RG16, vbo);

    glUseProgram(program);
    glUniform1i(uniform_loc_coord_tex, 0);

    // Draw
    glUniform1f(uniform_loc_radius, radius * 0.01f);
    glUniform1i(uniform_loc_instance_offset, 0);
    glUniform4fv(uniform_loc_color, 1, &color[0]);
    glUniform1f(uniform_loc_outline, outline);

    glBindVertexArray(vao);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
    glDrawElementsInstanced(GL_TRIANGLE_STRIP, 4, GL_UNSIGNED_BYTE, 0, count);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    glUseProgram(0);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);

    // Restore modified GL state
    glBlendEquationSeparate(last_blend_equation_rgb, last_blend_equation_alpha);
    glBlendFuncSeparate(last_blend_src_rgb, last_blend_dst_rgb, last_blend_src_alpha, last_blend_dst_alpha);
    if (last_enable_blend)
        glEnable(GL_BLEND);
    else
        glDisable(GL_BLEND);
    if (last_enable_cull_face)
        glEnable(GL_CULL_FACE);
    else
        glDisable(GL_CULL_FACE);
    if (last_enable_depth_test)
        glEnable(GL_DEPTH_TEST);
    else
        glDisable(GL_DEPTH_TEST);
    if (last_enable_scissor_test)
        glEnable(GL_SCISSOR_TEST);
    else
        glDisable(GL_SCISSOR_TEST);
    glPolygonMode(GL_FRONT_AND_BACK, last_polygon_mode[0]);
    glViewport(last_viewport[0], last_viewport[1], (GLsizei)last_viewport[2], (GLsizei)last_viewport[3]);
}

}  // namespace ramachandran
