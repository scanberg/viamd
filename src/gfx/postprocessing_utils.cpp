/*-----------------------------------------------------------------------
  Copyright (c) 2014, NVIDIA. All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:
   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.
   * Neither the name of its contributors may be used to endorse
     or promote products derived from this software without specific
     prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
  PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
  OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
-----------------------------------------------------------------------*/

// Shaders for HBAO are based on nVidias examples and are copyright protected as stated above

#include <gfx/postprocessing_utils.h>

#include <core/md_str.h>
#include <core/md_log.h>
#include <core/md_hash.h>

#include <gfx/gl_utils.h>

#include <float.h>

#include <shaders.inl>

#define PUSH_GPU_SECTION(lbl)                                                                       \
    {                                                                                               \
        if (glPushDebugGroup) glPushDebugGroup(GL_DEBUG_SOURCE_APPLICATION, GL_KHR_debug, -1, lbl); \
    }
#define POP_GPU_SECTION()                       \
    {                                           \
        if (glPopDebugGroup) glPopDebugGroup(); \
    }

namespace postprocessing {

// @TODO: Use half-res render targets for SSAO
// @TODO: Use shared textures for all postprocessing operations
// @TODO: Use some kind of unified pipeline for all post processing operations

static const unsigned short bluenoise_64x64_data[] = {
    0x0ee0,0x6e90,0xcb80,0x35a0,0x6310,0x29b0,0x6cf0,0xccd0,0x8b00,0xf890,0x7230,0x0b00,0xc560,0x4280,0xa2a0,0xedd0,0x96d0,0xb5d0,0x3130,0xd870,0xa120,0x1000,0xaee0,0xc810,0xdc50,0x2c20,0x63c0,0x0e70,0x7410,0x3070,0x68c0,0xf5f0,0xc7f0,0x41a0,0x5480,0x18d0,0xf200,0x31c0,0x0500,0x7160,0x1bb0,0xe7e0,0xa000,0x1920,0x44c0,0x9330,0x7490,0x4d90,0xc5c0,0x20d0,0x54f0,0xa4d0,0x0a60,0x4ae0,0xb5c0,0x9750,0xdd10,0xaec0,0x33b0,0x6880,0xc980,0x2060,0x6170,0xace0,
    0x9da0,0x50a0,0xf9d0,0x7d50,0xa3f0,0xf050,0xb470,0x1980,0x3160,0x47f0,0xb840,0x54b0,0xdaf0,0x68d0,0x1380,0x3710,0xdd80,0x0580,0x65c0,0xc120,0x75b0,0x5680,0xf6b0,0x1f40,0x7780,0x4150,0xa2d0,0xf080,0xba90,0x1b20,0xe3a0,0xa180,0x2480,0xec40,0x6650,0x8f80,0x77b0,0xa580,0x97f0,0xbe60,0x5510,0xd6b0,0x6150,0x33c0,0xa8f0,0xf2d0,0xe150,0xb7e0,0x3e10,0x7e60,0xfa00,0xd710,0x9280,0xe530,0x7200,0x2b60,0xcd60,0x11d0,0xe8c0,0x9370,0xa540,0xdac0,0x88f0,0xedf0,
    0x3b40,0x2460,0xb060,0x05d0,0x4b20,0xd8e0,0x40a0,0x7800,0xd3b0,0xe750,0x9930,0x2680,0x85d0,0xfac0,0x79b0,0x5c90,0xc440,0x7fd0,0xf430,0x2630,0xe120,0x85c0,0x3390,0x9ae0,0xb6a0,0x0940,0x8240,0x5060,0xdab0,0x8dc0,0x4650,0xb460,0x3840,0x0be0,0xaa50,0xc220,0x2250,0x5d10,0xdcf0,0x8780,0x2900,0xc9c0,0x78e0,0xb4f0,0x0740,0x68b0,0x2510,0x13c0,0xef10,0x6550,0xb270,0x2db0,0x5dc0,0xc450,0x0e80,0x64a0,0x89d0,0x7690,0x5090,0x2780,0x7c00,0x4400,0x1460,0xc370,
    0xd670,0x6760,0xe790,0x90a0,0xc4c0,0x1150,0x8510,0x5d30,0xab80,0x0640,0x6530,0x3ad0,0xa760,0x0d10,0xce30,0xad40,0x4630,0x1af0,0xa660,0x3cb0,0x95a0,0x0070,0x4980,0xd500,0x5c40,0xfba0,0xcc20,0x3670,0xaca0,0x0260,0x5f90,0xcd70,0x74f0,0x8660,0xd610,0x4760,0xe620,0x1270,0x3a90,0xf0d0,0x4590,0x0ef0,0xf750,0x4fd0,0xc470,0x87b0,0xd3a0,0xa410,0x8e80,0x0930,0x9b50,0x1be0,0x83a0,0x3b70,0xf740,0xaa60,0x4520,0xbd40,0x1a80,0xd0c0,0xb560,0xf880,0x5720,0x74e0,
    0xb990,0x8320,0x1b90,0x5bb0,0x2fd0,0xf3f0,0x9d00,0x23a0,0xff10,0x8f50,0xc750,0xe1e0,0xbb90,0x4cf0,0x9360,0x2410,0xe280,0x8ee0,0x50b0,0xb9c0,0x6060,0xe8f0,0xc270,0x6b90,0x2a90,0x14b0,0x9210,0x6e30,0x2440,0x7e50,0xf910,0x9670,0x13d0,0x59e0,0xfe90,0x3050,0x7fe0,0xb340,0xc640,0x6800,0xaaf0,0x9bd0,0x2210,0x9110,0xdde0,0x4120,0x5c70,0x2ff0,0x5110,0xcf10,0x43a0,0xe0e0,0xbb50,0x52e0,0xd640,0x1ee0,0xec00,0x9c90,0xf4d0,0x5fd0,0x36d0,0x0080,0x9a90,0x2d40,
    0x1080,0x46f0,0x98f0,0xe0c0,0x7470,0xb780,0x4f70,0xc010,0x3730,0x1260,0x7d00,0x58d0,0x1c90,0x70e0,0xef30,0x3430,0x6c40,0xfe10,0x0b50,0xcd50,0x72f0,0x15a0,0xad10,0x8b80,0xed00,0xa460,0xe290,0x42c0,0xc460,0xeb10,0x51e0,0x2b90,0xe030,0xb6d0,0x9d70,0x0710,0x6e80,0x9460,0x5100,0x0110,0x81b0,0xe1f0,0x6f40,0x3820,0x1600,0x7ef0,0xfc00,0x72a0,0xc230,0xebc0,0x7910,0x6910,0x01b0,0xa2e0,0x9000,0x3450,0x59a0,0x07f0,0x8450,0x6f00,0xabd0,0xe250,0x8af0,0xeb80,
    0xa4c0,0xfcf0,0xd230,0x3cf0,0xa810,0x0240,0x69a0,0xdc00,0x71f0,0x4940,0xeaa0,0x2dd0,0x9fa0,0xd600,0x0280,0xb540,0x5740,0x7d40,0xd7b0,0x2de0,0x9de0,0xf7e0,0x3910,0x1fa0,0x7b90,0x5310,0xb3c0,0x62e0,0x0e20,0x9f60,0xbbb0,0x3bd0,0x6d30,0x1fe0,0x4cb0,0xcb00,0xed80,0x25b0,0xd530,0xfbf0,0x3230,0x5bd0,0xbcc0,0xed60,0x9f20,0xad70,0x0350,0xb640,0x1e00,0x3680,0xa840,0xfec0,0x28b0,0xccf0,0x6e20,0x7c80,0xdfb0,0xc590,0x3fd0,0xd920,0x1760,0x4ce0,0xc700,0x5e00,
    0x70b0,0x09b0,0x54d0,0x25c0,0xca50,0x8a50,0x1d30,0xee90,0x96e0,0xcba0,0xb160,0x8940,0xf9c0,0x4340,0x8210,0xc920,0x9b80,0x1710,0xb090,0x4420,0x82f0,0x5860,0xc950,0x47a0,0xd8b0,0x0420,0x3110,0xcea0,0x7530,0x1ce0,0x8790,0xd480,0xa9b0,0xf390,0x8dd0,0x3590,0x60f0,0xa1c0,0x75e0,0x1b50,0xb0c0,0xcfb0,0x0bc0,0x53a0,0xd6e0,0x63d0,0x2a50,0xe4c0,0x9730,0x59f0,0x0ce0,0x8ba0,0x4ac0,0xe670,0x1480,0xb630,0x2710,0xa7e0,0x9350,0x2e30,0xbe40,0x7d80,0x2070,0x39a0,
    0xb190,0x8ca0,0xbe80,0x7d20,0xf770,0x5fa0,0x4270,0xa440,0x2980,0x07e0,0x62f0,0x1820,0xc070,0x67a0,0x2620,0xe830,0x3900,0x6440,0xefd0,0x2270,0xe040,0x0950,0x9200,0x6980,0xb950,0xfef0,0x85b0,0x98a0,0xf150,0xdce0,0x4820,0x0480,0x6280,0x7ae0,0xc290,0x0e60,0xdca0,0xb820,0x3d80,0x4ad0,0x88c0,0x95c0,0x27a0,0x4450,0x7850,0xc680,0x8a10,0x49e0,0xda20,0x8060,0xb3d0,0xc720,0x5f60,0x9ac0,0x4190,0xf120,0x5050,0x65e0,0x0f10,0xfe00,0xa100,0x69e0,0xf320,0xcda0,
    0xdda0,0x17b0,0x67b0,0x3090,0xe520,0x0cd0,0xb5f0,0xd120,0x57c0,0x78b0,0xd9b0,0x3bb0,0x52a0,0x94d0,0xab20,0x1050,0x4f30,0xc320,0x8e30,0x7000,0xbea0,0xa630,0xeb00,0x2bf0,0xa070,0x1a20,0x5d70,0x3b00,0x2760,0x55c0,0xb1f0,0x9490,0x1850,0xe240,0x4260,0x56d0,0x8730,0x1680,0xe980,0xc140,0xf5b0,0x6030,0xa7b0,0xf180,0x1b00,0x3570,0xf8d0,0x12d0,0x6dc0,0x3db0,0xf3a0,0x1910,0x3200,0xada0,0x8350,0x04f0,0xd420,0x76c0,0xe800,0x45c0,0x5790,0x0680,0x9500,0x2a10,
    0x4eb0,0xefe0,0xa020,0xadf0,0x4a90,0x9440,0x8140,0x3520,0xfc80,0x8d60,0xa710,0xf250,0x0a50,0xe2f0,0xd1a0,0x7580,0xf680,0xa260,0x05b0,0x4b60,0x33e0,0x6100,0x1230,0xd320,0x7370,0x4d20,0xe420,0xaa80,0xc790,0x6b20,0xf8f0,0x32f0,0xbce0,0x9fb0,0x29f0,0xfcc0,0xa960,0x6950,0x7ce0,0x2d70,0x08b0,0x6e60,0xe0a0,0x82e0,0xb6e0,0x9bc0,0x57e0,0xa570,0xc0c0,0x2650,0xd150,0x76a0,0xdc90,0x6a90,0xfa50,0xbe10,0x90f0,0x1d10,0xc820,0x8640,0xb810,0xe2e0,0x3ef0,0x8300,
    0xc400,0x04a0,0x3ba0,0xd520,0x1f00,0x6bb0,0xdf70,0x16d0,0xc300,0x4880,0x2330,0xb720,0x6e50,0x3360,0x87d0,0x5cd0,0x2c30,0xdc10,0x7e10,0xce70,0xfa60,0xb1e0,0x8770,0x4020,0xf370,0x8ea0,0xbdf0,0x1590,0x7bc0,0x08f0,0x83d0,0xe7d0,0x5ba0,0x7320,0xcc40,0x91a0,0x03d0,0xd660,0x5220,0x9d30,0xc8d0,0x1520,0x3ed0,0xceb0,0x0450,0x6ad0,0xd490,0xec10,0x0720,0x91d0,0x5200,0x9f30,0x0dc0,0x22d0,0x5700,0x3c50,0x2cd0,0x5ef0,0xae70,0x3560,0x1250,0xd4d0,0xaa20,0x6140,
    0x7290,0x97b0,0x5ad0,0x7a00,0xc620,0xf060,0x5520,0x9e50,0x02c0,0x6930,0xec20,0x7f50,0x9c00,0x1e20,0xbd00,0x4700,0xaf80,0x1890,0x3e90,0x96f0,0x1f50,0x51d0,0x7940,0xc430,0x23b0,0x0ab0,0x65d0,0xda40,0x44d0,0xa490,0xce10,0x1e80,0x3f20,0x0d70,0xf010,0x4a40,0x21d0,0xb310,0x3640,0xde50,0x8d50,0xae50,0x5890,0x2670,0x9320,0x4e50,0x2eb0,0x7a40,0xb150,0x6160,0x44b0,0xbc30,0xe990,0x8de0,0xcdd0,0xa4e0,0xe470,0x99d0,0xf670,0x6c20,0x24f0,0x7860,0xfa30,0x1b70,
    0xb010,0xe110,0xf700,0x2ba0,0x0fa0,0xba00,0x3e00,0x75a0,0xb020,0x30c0,0xd4c0,0x5930,0x4000,0xcd10,0xfeb0,0x0150,0x9260,0xecc0,0x6840,0xba40,0xe6d0,0x0c80,0xdec0,0x5b20,0xac00,0x3810,0x9810,0xf6a0,0x2550,0xba10,0x5010,0x8f10,0xd9a0,0xafd0,0x81c0,0x6b30,0xc7d0,0x5ca0,0xf870,0x1a60,0x4740,0xeba0,0x7960,0xffe0,0xba80,0xe2b0,0x4220,0x1fb0,0x8560,0xde60,0xfc50,0x3620,0x7f70,0xb370,0x48f0,0x1330,0x7b40,0x0980,0x4dd0,0xcd40,0x8fd0,0x9ec0,0x4900,0x3340,
    0x8970,0x0c10,0x43e0,0x8230,0xa640,0x8db0,0x24d0,0xf690,0xca60,0x9530,0x18b0,0xe3f0,0x0e50,0xa4b0,0x6370,0xd990,0x77a0,0xc520,0x5440,0x2af0,0xa4a0,0x8b60,0x3010,0x9c10,0xd0b0,0xeae0,0x8220,0x5750,0x7110,0x3480,0x9cc0,0xfb20,0x6400,0x2c10,0xa290,0x3950,0xe730,0x9910,0x7380,0x8390,0xbdb0,0x6570,0x31a0,0xa370,0x0d90,0x8890,0xc3d0,0xf450,0x0c40,0xa470,0x2ad0,0x7170,0x00a0,0x5df0,0xf040,0x6ed0,0xdb80,0xc060,0x4090,0xed10,0x5a60,0x01a0,0xe880,0xc9f0,
    0x5630,0xb7b0,0xd080,0x66f0,0xde20,0x5db0,0x4df0,0xd770,0x84a0,0x62b0,0x4b50,0xc020,0x8c50,0x7180,0x4e90,0x24c0,0x3660,0x8680,0x0990,0xd3d0,0x7190,0x4480,0xfc70,0x6900,0x0030,0x4b90,0xb210,0x1030,0xe490,0xd560,0x03f0,0x7720,0x1780,0xc2e0,0xe140,0x0750,0xbb60,0x1210,0x2b70,0xa730,0x0020,0xcb30,0x1e70,0xd570,0x70c0,0x5f70,0x9c60,0x52d0,0x6870,0xcab0,0x16a0,0xd5e0,0x9cb0,0xc500,0x2120,0x3240,0x8b90,0xb100,0x1bd0,0x8170,0xa8c0,0xbcb0,0x67f0,0x2820,
    0xfd60,0x1ec0,0x9d20,0x34f0,0xec80,0x05f0,0xaaa0,0x1cb0,0x38c0,0x0970,0xf940,0xac90,0x2da0,0xe9f0,0xb7a0,0x99e0,0xf070,0xa860,0x5d60,0xf570,0x1470,0xb4c0,0xc710,0x1c60,0x7810,0xbf70,0x2cb0,0x8f70,0xc490,0x5f30,0x42b0,0xb5a0,0x8a40,0x56c0,0x4710,0x7c40,0x8f90,0x5170,0xedb0,0x40d0,0xdcd0,0x5420,0x9720,0x3cd0,0xf160,0x15e0,0x33f0,0xda10,0xb410,0x3d40,0x8920,0x55a0,0xe200,0x3f60,0xa980,0xff30,0x5320,0x6730,0x2ac0,0xd7e0,0x37d0,0x14a0,0xdb00,0x7db0,
    0xa550,0x6f60,0x4ba0,0x1370,0xc210,0x9400,0x7150,0xbb30,0xe8b0,0x9e90,0x7760,0x2040,0x8280,0x3d30,0x0b80,0xc910,0x1810,0x4780,0xddb0,0x93b0,0x38e0,0x80c0,0x50d0,0xdbf0,0x9420,0x4060,0xf8b0,0x1a00,0xa230,0x80a0,0xf1e0,0x2f10,0xea20,0xcec0,0x2400,0xfe50,0x6330,0xd1d0,0xb070,0x6c10,0x8c00,0xfa10,0x7e20,0xb550,0x4b00,0xac10,0xe6c0,0x2600,0x74c0,0x93c0,0xbf00,0xf4e0,0x6620,0x82c0,0x15b0,0x95e0,0xd1c0,0x06b0,0x9dc0,0xf920,0x7260,0x5070,0x95b0,0x3f30,
    0xe6b0,0xc8f0,0x8b10,0xf290,0x7950,0x28d0,0xdae0,0x4490,0x6990,0xc900,0x5610,0xdf10,0xd340,0x5e40,0xfb00,0x6a00,0x7e30,0xb230,0x21c0,0xbfe0,0x64e0,0xa210,0x2720,0xef60,0xa9e0,0x61c0,0xd360,0x6c70,0x5140,0x2280,0xa880,0x0ed0,0x9870,0x7060,0xab70,0x1450,0x9f50,0x3260,0x1d60,0xc1f0,0x0fc0,0x5ee0,0x27c0,0x09d0,0xc580,0x8e90,0x7fb0,0x02f0,0xfae0,0x49c0,0x1cf0,0x2e20,0x08e0,0xb770,0x4b30,0x7300,0xe400,0xb8b0,0x5cf0,0xc800,0x8cd0,0xf110,0xb260,0x0870,
    0x1960,0x3000,0x5bf0,0xacc0,0x3c40,0x56f0,0xff80,0x0eb0,0x8fb0,0x3100,0xb610,0x46d0,0x14f0,0xa690,0x9170,0x5120,0xe370,0x3310,0x74a0,0x02e0,0xcf90,0xe5d0,0x10b0,0x5a20,0x3400,0x0690,0x7dc0,0xe8a0,0x3920,0xc8e0,0xda80,0x4c20,0x65f0,0x3a40,0xddd0,0xc550,0x4360,0xf280,0x82b0,0xe320,0x3770,0xa5f0,0xd6c0,0xea50,0x6ab0,0x56b0,0xce20,0xa350,0x6020,0xd2b0,0xab90,0x79a0,0x9b10,0xed30,0xcb50,0x2830,0x3a70,0x7cb0,0x4670,0x0dd0,0x3060,0x2160,0xc2d0,0x6220,
    0x7b60,0x9b70,0xdd20,0x00d0,0xcd90,0xb450,0xa080,0x7dd0,0x22c0,0xeed0,0x0370,0x9960,0x6c50,0xc2f0,0x2ab0,0x06a0,0xd210,0x9e70,0x5850,0xfdf0,0x8e70,0x3fb0,0x7050,0x88b0,0xc970,0xb250,0x97a0,0x1410,0xb960,0x8d10,0x7830,0xf6c0,0xbf30,0x01e0,0x86a0,0x5870,0x77e0,0x07d0,0x9890,0x4ea0,0x7240,0xbc70,0x4540,0x9cd0,0x2f50,0x1f20,0x4030,0xbb00,0x12a0,0x3750,0xe4b0,0x51f0,0xd960,0x1160,0x5c10,0x8a80,0xf530,0x1ed0,0xac30,0xea70,0x6a50,0x84e0,0x49d0,0xd680,
    0xf490,0xb800,0x2590,0x6be0,0x8500,0x1750,0x4d80,0xd4b0,0x63e0,0xa9c0,0x8440,0xe570,0x36a0,0xf5c0,0x8720,0xbb70,0x41e0,0xeda0,0x8370,0x2b50,0x4c40,0xa850,0xbc10,0xf550,0x2010,0xd900,0x4610,0xff40,0x5990,0x08d0,0x32e0,0x1da0,0xafb0,0x92d0,0x29c0,0xe780,0xb930,0x66c0,0xd970,0x2380,0xf7b0,0x9050,0x17f0,0x7b80,0xfdc0,0xdd70,0x73d0,0xf210,0x9860,0x6cd0,0x84f0,0xc3b0,0x3e80,0x6b00,0xb1a0,0xa310,0x02a0,0xc130,0x90d0,0xd260,0x9bf0,0xe1b0,0xa9a0,0x38b0,
    0x9030,0x52b0,0x4530,0xf9f0,0x9630,0xe330,0x3630,0xbd20,0xf4c0,0x3f40,0x5a40,0xcc60,0x7890,0x4d60,0x1f10,0xac60,0x6340,0x18f0,0xc6f0,0xb110,0x1c20,0xe090,0x0a90,0x7a60,0x4f50,0x6790,0x28e0,0xa5d0,0x6fe0,0xe450,0x9d10,0x60d0,0xcbd0,0x51a0,0xf980,0x18e0,0xa430,0x3040,0xca00,0xacd0,0x5940,0x0270,0xcfe0,0x6290,0xaf10,0x0810,0x89c0,0x5080,0x2c70,0xeb50,0x0510,0x2430,0x8fa0,0xf310,0x3150,0x4840,0xdcc0,0x6300,0x50f0,0x3ce0,0x1190,0x5970,0x0460,0x7210,
    0x16c0,0xe810,0x0bf0,0xbfd0,0x5d90,0x29a0,0x7670,0x0800,0x9220,0x1a50,0xb900,0x27b0,0x0b70,0xd6a0,0xe860,0x7330,0x9830,0x0df0,0xdd00,0x6aa0,0x9580,0x5f00,0x3ae0,0x9e20,0xe7c0,0x8fe0,0xc2b0,0x3c10,0x8400,0xd000,0x47c0,0xeec0,0x7de0,0x3df0,0x71b0,0xd180,0x4a70,0x8df0,0x0ec0,0x3fe0,0xebd0,0x8530,0x35d0,0xb910,0x48a0,0x9590,0xc4e0,0x1ba0,0xae00,0xccb0,0x9eb0,0x5920,0xb9d0,0x1ad0,0xcdb0,0x82d0,0x6fa0,0xfcb0,0x2660,0x76e0,0xefa0,0xbd50,0x2bb0,0xcf80,
    0xadc0,0x7ea0,0x9ef0,0xd380,0x3d00,0xab50,0xea40,0xca80,0x5210,0x6f70,0xe0b0,0x9e80,0xb1b0,0x9120,0x5ac0,0x3940,0xfa80,0x51c0,0x7bd0,0x3440,0xf660,0xc410,0xd2c0,0x2d90,0x1500,0xb080,0x0160,0xdbd0,0x19e0,0xb4d0,0x25e0,0x0af0,0xa650,0xe070,0x10a0,0xb480,0x6050,0xf0e0,0x7e00,0x6d80,0xc090,0xa170,0xdba0,0x2180,0xee50,0x5b50,0x3aa0,0xd8a0,0x6180,0x7b10,0x4690,0xfde0,0xa7a0,0x77f0,0xe480,0x0b60,0x9b20,0x1510,0xb130,0xcbe0,0xa340,0x81a0,0xfad0,0x6350,
    0xc4b0,0x3180,0x6e00,0x1c10,0x8bc0,0x67c0,0x1220,0xa220,0x8200,0x2e50,0xfd10,0x48c0,0x68e0,0x12f0,0xc390,0x2c40,0xa3a0,0xcf40,0xbb20,0x4570,0x03c0,0x8620,0x6ee0,0x5540,0xfb50,0x7f80,0x5f50,0xf3e0,0x52c0,0x9230,0x6680,0xc110,0x2e90,0x8900,0x9900,0x2200,0x3980,0x9e10,0xdff0,0x2a20,0x1a90,0x4d70,0x6780,0x78d0,0x0f90,0xa4f0,0x8030,0xf7d0,0x0a40,0x3510,0xdd90,0x1120,0x6600,0x3b80,0x5290,0xbd80,0x2f20,0xd800,0x89f0,0x36c0,0x1c80,0x4bc0,0x94a0,0x3ff0,
    0x21b0,0x55d0,0xdb30,0xf0f0,0x4be0,0xb740,0xf8c0,0x4320,0xd8c0,0xbee0,0x0180,0x8960,0x3be0,0xf2e0,0x7fc0,0xd980,0x0630,0x8b40,0x25a0,0xabc0,0xeb70,0x1f90,0xb660,0xa590,0x4050,0xc740,0x7430,0x3190,0xa2b0,0xe600,0x7750,0xf520,0x57f0,0xd5b0,0x6a70,0xfe70,0xc7b0,0x0430,0x5650,0xd0a0,0xb570,0x97d0,0xf9e0,0x3080,0xe340,0xcdc0,0x2530,0x6b80,0xb440,0x88e0,0xbf50,0x96b0,0x2920,0xeea0,0x8d70,0xabb0,0x5f40,0x4390,0xed20,0x5880,0x6bf0,0xb8d0,0x0ae0,0xe060,
    0xf540,0x8590,0xb120,0x0490,0xc770,0x26d0,0x7ca0,0x5af0,0x2130,0x96a0,0x6000,0xd240,0xa790,0x1dd0,0xb670,0x4b80,0x6450,0xef40,0x74d0,0x5c80,0x9470,0x4c10,0xe460,0x0c20,0xd750,0x97e0,0x2360,0x4a60,0xce80,0x07c0,0x39f0,0xae10,0x4550,0x0da0,0xb940,0x4e20,0x77c0,0xaab0,0x8c10,0xf4b0,0x42d0,0x0890,0xc650,0x8e00,0x53f0,0xbc40,0x4290,0x9fc0,0xecd0,0x1550,0x5500,0x7340,0xc840,0x0380,0xd460,0x1eb0,0xf840,0x7ba0,0xc260,0x0120,0xe720,0xd3c0,0x75f0,0xa090,
    0x60c0,0x0f70,0x9680,0x7630,0x34c0,0x9c20,0xe270,0x0aa0,0xad00,0x3650,0xeb40,0x7a50,0x5430,0xe3b0,0x6ff0,0x9710,0x1630,0x3f10,0xdfa0,0x1090,0xc1a0,0x7c70,0x2b00,0x6660,0x8950,0x1640,0xf100,0xb710,0x84d0,0x1970,0x8d20,0xc860,0x1f60,0xe4d0,0x8470,0x2880,0xea30,0x13f0,0x34d0,0x70d0,0x8160,0x5e60,0x1660,0xab30,0x6ea0,0x02d0,0x9450,0x4ee0,0x2e80,0xcfd0,0xe560,0x40e0,0xae20,0x83f0,0x4a10,0x9f00,0x6ec0,0x1130,0x95f0,0xa6f0,0x2a70,0x8ce0,0x3540,0x4890,
    0xd110,0xc160,0x3eb0,0xe590,0x5410,0x6480,0xd090,0x8f20,0xf240,0x6d10,0xc600,0x0e90,0x2970,0xcb40,0x33a0,0xffa0,0xaeb0,0xc830,0x9ee0,0x34a0,0xfab0,0xd170,0xa0f0,0x39b0,0xea90,0x5950,0xa6e0,0x6d50,0xdaa0,0x61d0,0xfbe0,0x9f10,0x72d0,0x5e90,0x9620,0x3b50,0xc280,0x64c0,0xde00,0x2080,0xb220,0xeef0,0xd540,0x3c00,0x8670,0xffd0,0xdb40,0x1b60,0x77d0,0x6380,0xa1e0,0x20b0,0xfbb0,0x6a60,0xdf40,0x2c80,0xcde0,0x3a10,0xe380,0x4dc0,0x6590,0xffc0,0x1620,0xaad0,
    0x2800,0x6bd0,0x1db0,0xfce0,0xa680,0x11b0,0xb9a0,0x3ee0,0x4ff0,0x1950,0xb4b0,0x45d0,0x9d60,0x8ab0,0x0ad0,0x58a0,0x82a0,0x2770,0x6a80,0x8d00,0x5660,0x06e0,0x7270,0xbaf0,0x48d0,0xc570,0x0400,0x2bd0,0x3fc0,0xbfc0,0x50c0,0x3300,0xdc20,0x0040,0xf590,0xb0b0,0xa110,0x4990,0x9130,0xcb60,0x9db0,0x2cc0,0x4c30,0xe7f0,0x2290,0x6080,0xc1b0,0xac40,0xf380,0x8ec0,0x06d0,0xb890,0x33d0,0x5b80,0x1400,0xba50,0x5530,0xb000,0x8570,0x1de0,0xb750,0xc780,0x8080,0xe9d0,
    0x4d50,0xb620,0x8bd0,0x7bb0,0xc9e0,0x2d80,0x7360,0x86d0,0xd700,0xa330,0x8020,0xf730,0xdbb0,0x62a0,0xbed0,0xe7b0,0x4850,0xd860,0xb870,0x1d00,0x41f0,0xad30,0xdea0,0x2110,0x92a0,0xfea0,0x7700,0xe170,0x9950,0x0ff0,0xac50,0x2690,0x7eb0,0xcd00,0x53c0,0x1b10,0xd810,0x0790,0xfb30,0x5570,0x10d0,0x6ac0,0x7cd0,0xb920,0xa2f0,0x34b0,0x0d40,0x8100,0x3a30,0xc5d0,0x4870,0xd820,0x9660,0x7a80,0xf340,0x8cf0,0x0620,0xee80,0x7450,0xdbe0,0x41c0,0x5ab0,0x0590,0x9480,
    0xe180,0x58f0,0xda00,0x0010,0x4910,0xe960,0x1ff0,0xfa70,0x06c0,0x3170,0x5a80,0x20c0,0x7220,0x3bc0,0xa800,0x1b40,0x7510,0x0000,0xf6f0,0xca90,0x8070,0xf260,0x5e30,0x1290,0x8340,0x32b0,0xb200,0x5490,0x8000,0xf090,0x69c0,0xe950,0xbc90,0x4180,0x8ac0,0x68f0,0x7900,0x2ed0,0x8360,0x3c80,0xbcd0,0xe230,0x9560,0x1800,0xcce0,0x7480,0x50e0,0xe940,0x26b0,0x5a50,0x7100,0x1020,0xe500,0xa610,0x3f90,0xc6e0,0x6470,0x9990,0x32d0,0x1060,0x9fe0,0x6de0,0xf2f0,0x37e0,
    0xa3c0,0x1580,0x3210,0x9ea0,0xbe30,0x5e50,0xaf30,0x99b0,0x6940,0xc4d0,0xe9b0,0x9160,0x04d0,0xcdf0,0xf1d0,0x94c0,0x31b0,0xa280,0x6070,0x4d00,0x9880,0x2a60,0xc2c0,0xa510,0xd400,0x64d0,0x0d50,0xc7a0,0x1fc0,0xd0f0,0x4830,0x91b0,0x14e0,0xa700,0x2310,0xb760,0xebb0,0xc4f0,0x5fc0,0xa6d0,0xf300,0x01f0,0x4460,0x5a00,0xf510,0x8ff0,0xd8d0,0xb3a0,0x9be0,0xd0e0,0xf6e0,0x8810,0x51b0,0x27e0,0x19d0,0xd5f0,0x4960,0x2470,0xfa90,0xbd70,0xcf60,0x21e0,0xb140,0x7560,
    0xc250,0x8520,0xef90,0x69d0,0x8f00,0x3b30,0xcf30,0x53d0,0xdcb0,0x4210,0xbae0,0xac70,0x4d30,0x2840,0xb7c0,0x5560,0x8820,0xe430,0xb2c0,0x0d60,0xeb30,0x3870,0x6f10,0x4f80,0xe900,0x43d0,0xf5e0,0x8c80,0x3b10,0xa320,0x0820,0x7590,0x5f20,0xd450,0xf830,0x39c0,0x1200,0x92c0,0xdad0,0x2220,0x7310,0xd060,0xafa0,0x29e0,0x6890,0x0670,0x3f00,0x14c0,0x66b0,0x1d40,0xaa30,0x31f0,0xc4a0,0xb400,0x6ce0,0x7fa0,0xe9e0,0xa9f0,0x5d00,0x7aa0,0x8e50,0x46e0,0xd690,0x0b10,
    0x6090,0x4350,0xcbc0,0x25d0,0xf760,0x0910,0x7f10,0x16b0,0x28f0,0x78f0,0x1110,0x8580,0xfcd0,0x6750,0x7af0,0x0fb0,0xd130,0x3e30,0x6df0,0x7d10,0x1990,0xcf00,0x8ad0,0x0290,0x9ab0,0x2750,0xa930,0x7010,0x5cb0,0xbad0,0xfd20,0x3410,0xded0,0x4f10,0x9980,0x72b0,0x5770,0xae80,0x0a20,0x4c80,0x8b70,0x3530,0x8180,0xdeb0,0xc330,0xa7c0,0xfbc0,0x8420,0xbfb0,0x4410,0x76b0,0x0220,0x6200,0xfdd0,0x9bb0,0x0a10,0xbb40,0x8980,0x3d70,0x0200,0xe680,0x53e0,0x2d10,0xfb40,
    0x9390,0x1b80,0xb2f0,0x7390,0x4fe0,0xaa70,0xe390,0x94e0,0xf5d0,0xa5a0,0x5e10,0x3500,0xd760,0x9d40,0x4440,0xee10,0x1e50,0xc420,0x27d0,0xd940,0xa8e0,0x5960,0xb8a0,0xdf50,0x7840,0xbf40,0x16f0,0xd7f0,0xe540,0x28a0,0x97c0,0x8490,0xc1c0,0x0340,0x2d30,0xcb20,0xe640,0x41d0,0xfe80,0xa190,0x63b0,0xf020,0x1ac0,0x9c70,0x49a0,0x78c0,0x30f0,0x5550,0xe5a0,0x9780,0xefb0,0xd740,0x9040,0x3c60,0xde80,0x2ea0,0x5280,0x1740,0xdb70,0xc480,0x6670,0x9aa0,0xabf0,0x7e40,
    0xe2d0,0x5710,0xec30,0x0e30,0xd850,0x3330,0xbbc0,0x6d90,0x4af0,0xcbb0,0xe700,0x1f70,0xc080,0x0920,0xdef0,0xa6a0,0x5b70,0x9920,0xfe40,0x4a20,0x9270,0xf140,0x41b0,0x1e40,0x3470,0xf2c0,0x5270,0x8110,0x0cc0,0x4bb0,0xaff0,0x1dc0,0xec90,0x67d0,0xa530,0x8870,0x1ea0,0x7f00,0xbeb0,0x2c50,0xca40,0xb490,0x54e0,0x0cf0,0xeac0,0x20a0,0x93d0,0xd4e0,0x0b90,0x2700,0xb520,0x4ca0,0x7df0,0x1e60,0x5c00,0xca70,0x7520,0xf600,0xa060,0x3420,0x1d70,0xee00,0x0f60,0x3c20,
    0xd280,0x2b20,0xa360,0x7b70,0x9a70,0x6260,0x20e0,0x8be0,0x3a80,0x03a0,0x98d0,0x7500,0x52f0,0x8bb0,0x2f70,0x6d20,0xbd10,0x7f90,0x05e0,0x6580,0x2fe0,0x0bd0,0x6b50,0x8460,0xb2e0,0x61a0,0x93f0,0xca10,0x39d0,0x6b70,0xcf70,0x5900,0x7ac0,0x4370,0xb590,0x0de0,0xd270,0x5d50,0x6d70,0x1670,0x91c0,0x3dd0,0xd580,0x6e70,0xba30,0x6190,0xc9a0,0xae90,0x7090,0x5f80,0x36b0,0xcc90,0x1170,0xa160,0xe8d0,0xaf00,0x9180,0x47b0,0x6b60,0xb330,0x86e0,0xcaf0,0x71a0,0xbcf0,
    0x0300,0x8930,0x3890,0xc6a0,0x44f0,0xff70,0x0b30,0xc3a0,0xef50,0xafe0,0x6710,0xd430,0x3f50,0xb530,0xf950,0x1610,0x4da0,0x37b0,0xae30,0xe580,0xbef0,0xcae0,0xa1a0,0xfaf0,0xd2d0,0x4470,0x0560,0xed40,0xa620,0x8b30,0xf7c0,0x1040,0x3760,0xd9f0,0xf440,0x4e30,0x9cf0,0x35c0,0xdf20,0xf610,0x0470,0x7920,0xe310,0x2b30,0x87e0,0xf0c0,0x3a00,0x17a0,0x8cc0,0xe190,0xa5c0,0xf820,0x6e40,0xb8e0,0x42e0,0x0ca0,0x26f0,0xd5a0,0x0730,0x5810,0xfe60,0x2570,0x4c90,0x62d0,
    0xb650,0xf8a0,0x6cb0,0x1720,0xb1d0,0x8550,0xd3e0,0x5b00,0x7c30,0x2d60,0x1530,0xf3c0,0x80f0,0x2500,0x95d0,0xcee0,0xead0,0x8840,0xd5c0,0x14d0,0x7730,0x8d40,0x4cd0,0x2810,0x11c0,0x9af0,0x7680,0x2d50,0xb8c0,0x1a70,0xe130,0x9b60,0xbe00,0x8d80,0x2420,0x7130,0xe920,0xaf20,0x47e0,0x85f0,0xa740,0x5b40,0x9970,0x4600,0xa940,0x00b0,0x5190,0xfed0,0x4510,0x7f60,0x0550,0x5600,0x2260,0x8850,0xf170,0x6560,0x8190,0xe300,0xc360,0x7930,0x4040,0xa560,0xde10,0x9840,
    0x1ca0,0x4750,0xd8f0,0xe6a0,0x53b0,0x2960,0xebe0,0x3e40,0xa3e0,0x4e10,0xde40,0xc630,0x5a70,0xa890,0x63f0,0x0330,0x71d0,0x23c0,0x9f70,0x5f10,0x3c70,0x1c70,0xe3d0,0x5b60,0xad80,0xe770,0xc530,0x6520,0x5640,0x46c0,0x72e0,0x2b80,0xaa90,0x6240,0x05c0,0xc6b0,0x9250,0x13a0,0x2730,0xbbe0,0xcf20,0x1f30,0xf360,0x1360,0xc670,0x7d70,0xa0e0,0xd7d0,0xbf60,0x2d20,0x98b0,0xc6c0,0xdb20,0x30e0,0xc0b0,0x9b30,0x3930,0xa900,0x1aa0,0x94b0,0xe740,0x1180,0x3350,0x76f0,
    0xab40,0x5dd0,0x92e0,0x0f00,0xa670,0x67e0,0x9770,0x1c30,0xb9f0,0x8d30,0x0c70,0x9c30,0x1b30,0x3860,0xe210,0xc350,0x4500,0xb3b0,0xf780,0x5160,0xed90,0xb8f0,0xd830,0x7040,0x81e0,0x3800,0x2030,0xff20,0xd4f0,0x8600,0x0c60,0xccc0,0x5260,0x8090,0xfd70,0x3290,0x5a10,0x7a10,0xdb90,0x6820,0x5040,0x3740,0xb180,0x6610,0xe760,0x3120,0x7280,0x5e20,0x1ef0,0xb2b0,0xee70,0x62c0,0x4810,0x79d0,0x1340,0x5130,0xfa40,0x5e70,0x2e40,0xb970,0x4f40,0x8b20,0xcfa0,0xf400,
    0x26a0,0x7ff0,0xcc70,0x32c0,0x7a70,0xc190,0x00e0,0x7460,0xcc50,0xfc10,0x6ca0,0x43c0,0xb7f0,0x7880,0xf640,0x54a0,0x8fc0,0x2e00,0x0a00,0x7da0,0x9510,0x31e0,0xa420,0x0060,0x4800,0xbbf0,0x90e0,0x0860,0xa140,0xb280,0xee20,0x38f0,0xe440,0x1ae0,0xd620,0xb360,0x4230,0xef80,0xa250,0x0770,0xf860,0x7610,0x8a70,0xd330,0x22a0,0xb730,0x9100,0x0e00,0xe360,0x6d40,0x3bf0,0x0a80,0x9290,0xadb0,0xe510,0xcc00,0x0230,0x7350,0xd1f0,0xf1f0,0x64f0,0x0830,0xc0f0,0x5580,
    0xdc80,0x0540,0xba60,0x4240,0xef70,0xd6f0,0x4a80,0xe220,0x5da0,0x2f90,0x85a0,0xea60,0xd7a0,0x2b10,0x0f20,0x8310,0xa300,0xde70,0xbd30,0xce50,0x6410,0x1100,0x8760,0xca30,0xf5a0,0x5400,0xdd30,0x7820,0x4140,0x23d0,0x6720,0xc040,0x9060,0x48b0,0x6eb0,0x9a10,0x0ea0,0xc000,0x8330,0x2fa0,0x9690,0xc310,0x0f50,0x4a50,0x5830,0xf900,0x3e60,0xcbf0,0x4ed0,0xaa00,0x8910,0xd2a0,0xfca0,0x1d20,0x6af0,0xa380,0x89e0,0x4680,0xaf90,0x2300,0x7f40,0x9f80,0x6db0,0x3de0,
    0x8800,0x9a50,0xfd50,0x6390,0x19a0,0x8c20,0xa990,0x38a0,0x21a0,0xb0f0,0x0700,0x5240,0x93a0,0xaa10,0x6960,0xcc10,0x1650,0x3d90,0x6da0,0x2190,0xfd90,0x42a0,0xb390,0x29d0,0x6130,0x1700,0xab10,0x2f00,0xf1b0,0x5b30,0x99a0,0x7a30,0x0130,0xa720,0xecb0,0x24e0,0x6420,0xd350,0x5390,0x1c40,0xe840,0x4110,0xa5b0,0xe160,0x9a30,0x0660,0x8410,0x9fd0,0xf2b0,0x1440,0x7640,0x2870,0xbd90,0x5a90,0x3da0,0x2b40,0xdc70,0x1570,0x92b0,0x3a50,0xd950,0x18a0,0xb1c0,0xec50,
    0x12b0,0x4ec0,0x70f0,0xa130,0x2ca0,0x5910,0xf800,0x8130,0x9e40,0xee40,0xd040,0x6270,0x2020,0xc0e0,0x4920,0xe7a0,0x5c20,0xf3d0,0xabe0,0x4c60,0x9c50,0xd9c0,0x6f30,0xec70,0x98e0,0xd1b0,0x6c80,0x8a20,0xc5f0,0x0fe0,0xd9e0,0xfb10,0x2f80,0x58c0,0xbc50,0x89b0,0x3830,0xfaa0,0xa870,0x6f90,0xce00,0x61b0,0x2990,0x7d30,0x6d60,0xd880,0xc170,0x2580,0x5c30,0x3580,0xeb20,0x9c80,0x4c50,0x8380,0xb680,0xeeb0,0x7870,0xc540,0xe870,0x5be0,0xf990,0x4a30,0xc990,0x2f30,
    0xa8b0,0xe100,0x2100,0xcac0,0xb350,0x08a0,0xbe70,0x1240,0x6ae0,0x4170,0xbba0,0x7a90,0x39e0,0xfee0,0x0190,0x99c0,0x3320,0x7be0,0x8c90,0x0840,0xc3c0,0x5690,0x1bc0,0x3720,0x7d60,0x0b40,0x45b0,0xe710,0xb580,0x5000,0x1cc0,0xacb0,0x4200,0xcb10,0x15c0,0xe000,0x7740,0x0530,0x4620,0x8da0,0xb7d0,0x0a30,0xfe20,0xb040,0x19b0,0x31d0,0x4860,0x7980,0xdc40,0xb4a0,0x6640,0xc7c0,0x0900,0xd650,0x12e0,0x9a00,0x6250,0x4e60,0x0610,0xa400,0x2a30,0x8e20,0x78a0,0x6110,
    0xc200,0x3ab0,0x7ed0,0xea10,0x49b0,0xd4a0,0x73f0,0xde30,0x5030,0x8f40,0x17d0,0xa1b0,0xd930,0x8880,0x71c0,0xaf50,0xdc60,0x19f0,0xd300,0x2d00,0x7570,0xe6e0,0x9140,0xa8d0,0xbca0,0xfc20,0x24a0,0x9f40,0x37a0,0x72c0,0x86f0,0xd100,0x6850,0x80e0,0x9ed0,0x4e40,0xaf70,0xc690,0x22e0,0xdb10,0x56e0,0x8430,0x3b60,0xc940,0xea00,0x6460,0xa7d0,0xf720,0x0360,0x8cb0,0x1e90,0x4330,0xe5f0,0x7030,0xf930,0x35f0,0x2140,0xb380,0x83e0,0xbd60,0x6bc0,0xd3f0,0x08c0,0xf270,
    0x55f0,0x9430,0x0c90,0x6740,0x89a0,0x3e20,0x9850,0x2f40,0xc6d0,0xf580,0x28c0,0xe4e0,0x0e40,0x4fa0,0xc890,0x2340,0x4310,0xbbd0,0x65b0,0xefc0,0xb0a0,0x3e50,0x04e0,0xdfd0,0x4e70,0x5bc0,0xcca0,0x6500,0x0410,0xe020,0x9570,0x26c0,0xe650,0x0960,0xf4a0,0x2bc0,0x5eb0,0xee30,0x9b00,0x3460,0xf1c0,0x1770,0x9e00,0x4cc0,0x8e10,0xb980,0x1280,0x9740,0x5470,0xd050,0x7cf0,0xadd0,0x2ec0,0x9240,0x5250,0xab00,0xce40,0xf4f0,0x40f0,0xe1d0,0x15d0,0x9ba0,0x42f0,0xb5e0,
    0x1a40,0xd7c0,0x3250,0xfa20,0xad60,0x1540,0xf0a0,0x5ff0,0x0440,0xa820,0x59c0,0x6e10,0xb700,0x30b0,0x6120,0xe9c0,0x90b0,0x54c0,0xa010,0x10e0,0x8690,0x5ed0,0xc9b0,0x2930,0x8990,0x19c0,0xb030,0x7f20,0xf7a0,0xbf20,0x4560,0x5d40,0xb600,0x3a20,0x71e0,0xd5d0,0x9010,0x10c0,0x7bf0,0x6770,0xab60,0xd0d0,0x73e0,0x5d80,0x2230,0xd590,0x4130,0x6cc0,0xc050,0x3960,0xff50,0xa150,0x5fe0,0xc380,0x0090,0x8a60,0x66e0,0x0d80,0x75c0,0x32a0,0x57a0,0xfd00,0x2560,0x8480,
    0xa500,0x74b0,0xc660,0x59d0,0x9dd0,0x25f0,0x79e0,0xb510,0x8750,0xd550,0x4660,0x8250,0x9610,0xf130,0xa520,0x16e0,0x7ee0,0xfc40,0x3610,0xc0a0,0x48e0,0xf810,0x79f0,0xa3b0,0x6a30,0xef00,0x3d60,0x91e0,0x2dc0,0x1300,0xa8a0,0xf000,0x1790,0x9a20,0xc240,0x4790,0x1e30,0xb850,0x3fa0,0xe610,0x02b0,0x2cf0,0xbf90,0xf6d0,0x09a0,0x8050,0xedc0,0x2aa0,0xdf90,0x0c00,0x4de0,0x15f0,0xd910,0x2320,0x7a20,0xe2a0,0x47d0,0xa040,0xd730,0x92f0,0xae40,0xc340,0x68a0,0xdfe0,
    0xed50,0x4770,0xb2d0,0x0250,0xe1c0,0x4e00,0xc2a0,0xe630,0x35b0,0x1f80,0xfbd0,0xc1e0,0x3b20,0x09c0,0xcd30,0x7420,0xb3e0,0x03e0,0xddc0,0x6f20,0x21f0,0xd410,0x1490,0x3370,0xbe20,0xd9d0,0x1010,0xd030,0x5620,0x6c30,0xc760,0x79c0,0x5230,0x8630,0xff60,0xa770,0x6b40,0xde90,0xc510,0x5360,0x88a0,0x45e0,0x9310,0xe1a0,0x37c0,0xaed0,0x57b0,0x8a30,0x9d80,0xb430,0x7400,0x86c0,0xede0,0xb240,0x37f0,0xf350,0xba20,0x2a00,0x1a10,0xece0,0x7e90,0x03b0,0x4f00,0x2fc0,
    0x1140,0x6210,0x2890,0x7f30,0xd020,0x6a20,0x1930,0x43b0,0x7080,0x9d50,0x1310,0x6540,0xe010,0x2610,0x5820,0x3f80,0xd780,0x2ae0,0x4f60,0x9a40,0xaae0,0x8e40,0xe9a0,0x44a0,0x9800,0x7440,0x4c00,0xa480,0xe3c0,0x8860,0x3600,0x22f0,0xdd50,0x01d0,0x2e60,0x5760,0x0cb0,0x8270,0x2520,0xa200,0xfc60,0xb3f0,0x1350,0x6810,0xa0a0,0xc880,0x75d0,0x1c00,0xf470,0x6360,0x2740,0xc960,0x4250,0x5a30,0x9b90,0x0f80,0x5590,0x7070,0xc8c0,0x6040,0x3a60,0xd250,0xa0b0,0x8ef0,
    0xcd20,0xf7f0,0x9700,0x3c90,0x8a90,0xff00,0xa5e0,0x9380,0xd140,0xb2a0,0x5380,0x7b50,0xaac0,0x8a00,0xf630,0x9790,0x6630,0x8700,0xc5b0,0xee60,0x0a70,0x55b0,0xb690,0x6320,0x0310,0xfdb0,0x2640,0xb860,0x0760,0xf650,0x9e30,0xcd80,0xaea0,0x65a0,0xe5c0,0xd1e0,0x9540,0xf2a0,0x35e0,0x7250,0x5d20,0xd2f0,0x7cc0,0x2450,0x5020,0xdbc0,0x0210,0x46b0,0xcfc0,0x3700,0xe5e0,0x91f0,0x0570,0x6c00,0xd200,0x8260,0xdd60,0xa9d0,0x8e60,0x4aa0,0x2150,0xf560,0x6dd0,0xb790,
    0x7ad0,0x1870,0xda60,0xbc80,0x0fd0,0x30a0,0x57d0,0x0ba0,0xf330,0x2850,0xe910,0x00c0,0xd6d0,0x4b40,0xbab0,0x0bb0,0xa450,0xe550,0x1840,0x60b0,0x3ca0,0x80d0,0x27f0,0xdf30,0xca20,0x87a0,0x5c50,0x7b20,0x4010,0x61e0,0x1ab0,0x4ab0,0x8f60,0x3d10,0xba70,0x7710,0x43f0,0xaf60,0xcb70,0x07b0,0x1cd0,0xe3e0,0x3dc0,0xbdd0,0xf190,0x2ee0,0x9820,0xbb10,0xa830,0x1420,0x5180,0xbf10,0xa3d0,0xf970,0x2e70,0x1730,0x4070,0xfe30,0x0600,0xbde0,0xe0d0,0x87c0,0x0d20,0x40c0,
    0x5ae0,0xa920,0x4d40,0x6690,0xe8e0,0xaef0,0xcad0,0x7c50,0x5fb0,0x3990,0x9090,0xc5a0,0x2fb0,0x6f50,0x1d90,0xc930,0x3550,0x4730,0x7650,0xbc60,0xfb80,0xd190,0x7140,0x9ff0,0x1bf0,0x3780,0xac20,0xeb60,0x9600,0xbdc0,0xda50,0x7120,0xf220,0x12c0,0xa050,0x2790,0x1860,0x6230,0xe930,0x9b40,0x4d10,0xa780,0x8f30,0x0e10,0x70a0,0x8610,0x5cc0,0xfd30,0x6a40,0x8120,0xdb60,0x7620,0x2050,0xb0d0,0x8c30,0xc180,0x69b0,0x9a60,0x3280,0x76d0,0xb0e0,0x5450,0x2c60,0xe5b0,
    0xc0d0,0x34e0,0x2240,0x9f90,0x7600,0x4160,0x20f0,0x8830,0xdf60,0xb9e0,0x69f0,0x4580,0xa1d0,0xeaf0,0x5c60,0x7c10,0xf8e0,0x22b0,0xade0,0x9550,0x30d0,0x1430,0xb170,0x4e80,0xf460,0xc150,0x11a0,0xd470,0x2be0,0x0c50,0x8150,0x3030,0x59b0,0xc5e0,0x83b0,0xfb60,0x5330,0x8c70,0xbff0,0x7e70,0x2860,0xf9a0,0x63a0,0xcef0,0xad50,0xe4f0,0x2000,0x3ec0,0x90c0,0x0b20,0xeee0,0x3380,0x4720,0x60e0,0xe850,0x4fb0,0x2390,0xcc30,0xe970,0x1320,0x64b0,0xa1f0,0xcb90,0x9520,
    0x7020,0xeff0,0x8650,0xdee0,0x0520,0xf1a0,0xbe50,0x4970,0x0850,0xa910,0x1830,0xfda0,0x84b0,0x0db0,0xb420,0xd220,0x8d90,0xdc30,0x5730,0x01c0,0x85e0,0xe350,0x6700,0x4080,0x8ed0,0x7770,0x5460,0x6a10,0x44e0,0xfc30,0xa600,0xb500,0xe820,0x04b0,0x6b10,0xcff0,0xdf80,0x0d30,0x38d0,0x6ba0,0xdb50,0xb830,0x3270,0x5670,0x4430,0x0880,0xc3e0,0xd790,0x2ce0,0xb320,0x5800,0xa030,0xcf50,0x07a0,0x7c60,0xda30,0xa6b0,0x5b90,0x81d0,0x4300,0xd890,0x1d50,0xfb70,0x00f0,
    0xaf40,0x1390,0xc9d0,0x5300,0x93e0,0x60a0,0xd510,0x9ce0,0xf620,0x55e0,0x7790,0xd840,0x26e0,0x98c0,0x3d50,0x4fc0,0x11f0,0x6920,0xc100,0xf230,0x49f0,0xc8b0,0x9ad0,0x0c30,0xea80,0x24b0,0xa240,0xe2c0,0x9070,0xc870,0x5150,0x2170,0x99f0,0x3ac0,0x4930,0xa750,0x2df0,0xb290,0xf480,0x45f0,0x0050,0x8740,0x17e0,0x9760,0xf500,0x7b30,0xa0d0,0x4f90,0x73b0,0xc7e0,0x1c50,0x83c0,0xf710,0xbb80,0x13e0,0x3970,0x9150,0x09e0,0xf790,0xbac0,0x9080,0x36f0,0x7c20,0x4950,
    0xd630,0x6860,0x2a80,0xb880,0x3880,0x1880,0x6f80,0x3020,0x23e0,0x9190,0xcc80,0x3790,0xbfa0,0x6430,0xf420,0xa950,0x2c90,0x9ca0,0x3af0,0x73a0,0xa6c0,0x1e10,0x5980,0xd720,0x3220,0xb6f0,0xce90,0x0140,0x7c90,0x1690,0x6510,0x7660,0xd2e0,0x8b50,0xbda0,0x1df0,0x7990,0x5780,0x9340,0xd370,0xa2c0,0x7550,0xc850,0xdf00,0x23f0,0x6490,0xb6b0,0x1560,0xf960,0x9640,0xe0f0,0x66a0,0x3e70,0x9940,0x6fb0,0xf030,0xb300,0x4b10,0x6c90,0x2950,0xa970,0xeb90,0x5e80,0x8c60,
    0x3f70,0xf410,0xa270,0x7ab0,0xfd40,0xac80,0x81f0,0xebf0,0xb5b0,0x66d0,0x10f0,0x4bf0,0xe6f0,0x0390,0x7540,0xe050,0x8290,0xcaa0,0x0ac0,0xe690,0x2910,0x7e80,0xbc00,0xff90,0x6fd0,0x84c0,0x5ea0,0x3cc0,0xb050,0xecf0,0x3690,0xdfc0,0x11e0,0xf9b0,0x5de0,0xeab0,0x9df0,0xc730,0x13b0,0x2540,0x5ce0,0xef20,0x4db0,0xaa40,0x3b90,0x8ae0,0xeca0,0x3490,0x5ec0,0x45a0,0x0100,0x2940,0xaba0,0x5340,0x2090,0xc8a0,0x2e10,0xe260,0x1a30,0xd440,0x5370,0x1070,0xc610,0x2370,
    0x9650,0x0780,0x5b10,0x4640,0xd290,0x0f30,0x4f20,0xc3f0,0x3d20,0xddf0,0x7ec0,0xae60,0x94f0,0x58b0,0xb9b0,0x1900,0x46a0,0xfc90,0xafc0,0x6010,0xd070,0x9300,0x3850,0x05a0,0xa7f0,0x18c0,0x4a00,0xf850,0x2350,0x9e60,0xc030,0x4b70,0xad20,0x2c00,0x6fc0,0x09f0,0x3c30,0x8540,0x6830,0xe080,0xbaa0,0x2ef0,0x0d00,0x6c60,0xbe90,0x04c0,0xda90,0x8010,0xad90,0xd010,0x8c40,0xbec0,0xec60,0xd310,0x8710,0x61f0,0xa0c0,0x7b00,0xbf80,0x86b0,0x9c40,0x73c0,0xb6c0,0xe410,
    0x80b0,0xdd40,0xbc20,0x1d80,0x8eb0,0xe4a0,0x9a80,0x0170,0x5aa0,0xa390,0x1fd0,0xf0b0,0x2f60,0x88d0,0xd160,0x2490,0x6ef0,0x5350,0x8aa0,0x17c0,0x4380,0xed70,0x6970,0x4ef0,0x8bf0,0xe660,0xc1d0,0x9410,0xd390,0x56a0,0x87f0,0x0650,0x7d90,0x96c0,0xda70,0xb4e0,0xce60,0x4c70,0xfb90,0xacf0,0x40b0,0x9020,0x8040,0xfff0,0xced0,0x58e0,0x2a40,0x9d90,0x0f40,0xe890,0x6d00,0x36e0,0x7970,0x1940,0xfd80,0x4100,0x0320,0x5840,0xf3b0,0x3ea0,0x06f0,0xffb0,0x4bd0,0x3140,
};

static inline bool is_orthographic_proj_matrix(const float P[4][4]) { return P[2][3] == 0.0f; }

static struct {
    GLuint vao = 0;
    GLuint v_shader_fs_quad = 0;
    GLuint tex_width = 0;
    GLuint tex_height = 0;

    struct {
        GLuint fbo = 0;
        GLuint tex_rgba8 = 0;
    } tmp;

    struct {
        GLuint fbo = 0;
        GLuint tex_color[2] = {0, 0};
        GLuint tex_temporal_buffer[2] = {0, 0};  // These are dedicated and cannot be use as intermediate buffers by other shaders
    } targets;

    struct {
        GLuint fbo = 0;
        GLuint tex_tilemax = 0;
        GLuint tex_neighbormax = 0;
        int32_t tex_width = 0;
        int32_t tex_height = 0;
    } velocity;

    struct {
        GLuint fbo_0 = 0;
        GLuint fbo_1 = 0;
        GLuint fbo_2 = 0;
        GLuint fbo_3 = 0;
        GLuint texture = 0;
        struct {
            GLuint program_persp = 0;
            GLuint program_ortho = 0;
        } linearize;
        struct {
            GLuint program;
        } downsample;
    } linear_depth;

    struct {
        GLuint tex_random = 0;
        GLuint ubo_hbao_data = 0;
        GLuint fbo = 0;
        GLuint tex[2] = {};

        struct {
            GLuint program_persp = 0;
            GLuint program_ortho = 0;
        } hbao;

        struct {
            GLuint program = 0;
        } blur;
    } ssao;

    struct {
        GLuint program = 0;
        struct {
            GLint tex_half_res = -1;
            GLint tex_color = -1;
            GLint tex_depth = -1;
            GLint pixel_size = -1;
            GLint focus_point = -1;
            GLint focus_scale = -1;
            GLint time = -1;
        } uniform_loc;

        struct {
            GLuint fbo = 0;
            GLuint program = 0;
            struct {
                GLuint color_coc = 0;
            } tex;
            struct {
                GLint tex_depth = -1;
                GLint tex_color = -1;
                GLint focus_point = -1;
                GLint focus_scale = -1;
            } uniform_loc;
        } half_res;
    } bokeh_dof;

    struct {
        GLuint program = 0;
        struct {
            GLint mode = -1;
            GLint tex_color = -1;
        } uniform_loc;
    } tonemapping;

    struct {
        GLuint program = 0;
        struct {
            GLint tex_rgba = -1;
        } uniform_loc;
    } luma;

    struct {
        GLuint program = 0;
        struct {
            GLint tex_rgbl = -1;
            GLint rcp_res  = -1;
            GLint tc_scl   = -1;
        } uniform_loc;
    } fxaa;

    struct {
        struct {
            GLuint program = 0;
            struct {
                GLint tex_linear_depth = -1;
                GLint tex_main = -1;
                GLint tex_prev = -1;
                GLint tex_vel = -1;
                GLint tex_vel_neighbormax = -1;
                GLint texel_size = -1;
                GLint time = -1;
                GLint feedback_min = -1;
                GLint feedback_max = -1;
                GLint motion_scale = -1;
                GLint jitter_uv = -1;
            } uniform_loc;
        } with_motion_blur;
        struct {
            GLuint program = 0;
            struct {
                GLint tex_linear_depth = -1;
                GLint tex_main = -1;
                GLint tex_prev = -1;
                GLint tex_vel = -1;
                GLint texel_size = -1;
                GLint time = -1;
                GLint feedback_min = -1;
                GLint feedback_max = -1;
                GLint motion_scale = -1;
                GLint jitter_uv = -1;
            } uniform_loc;
        } no_motion_blur;
    } temporal;
} gl;

static constexpr str_t v_shader_src_fs_quad = STR_LIT(
R"(
#version 150 core

out vec2 tc;

uniform vec2 u_tc_scl = vec2(1,1);

void main() {
	uint idx = uint(gl_VertexID) % 3U;
	gl_Position = vec4(
		(float( idx     &1U)) * 4.0 - 1.0,
		(float((idx>>1U)&1U)) * 4.0 - 1.0,
		0, 1.0);
	tc = (gl_Position.xy * 0.5 + 0.5) * u_tc_scl;
}
)");

static constexpr str_t f_shader_src_linearize_depth = STR_LIT(
R"(
#ifndef PERSPECTIVE
#define PERSPECTIVE 1
#endif

// z_n * z_f,  z_n - z_f,  z_f, *not used*
uniform vec4 u_clip_info;
uniform sampler2D u_tex_depth;

float ReconstructCSZ(float d, vec4 clip_info) {
#if PERSPECTIVE
    return (clip_info[0] / (d*clip_info[1] + clip_info[2]));
#else
    return (clip_info[1] + clip_info[2] - d*clip_info[1]);
#endif
}

out vec4 out_frag;

void main() {
  float d = texelFetch(u_tex_depth, ivec2(gl_FragCoord.xy), 0).x;
  out_frag = vec4(ReconstructCSZ(d, u_clip_info));
}
)");

static constexpr str_t f_shader_src_downsample_depth = STR_LIT(R"(
#version 150 core

uniform sampler2D u_tex_linear_depth;
uniform int u_src_lod;

out vec4 out_frag;

#define OP max

void main() {
    ivec2 base_coord = ivec2(gl_FragCoord.xy) * 2;

	float d00 = texelFetch(u_tex_linear_depth, base_coord + ivec2(0,0), u_src_lod).x;
	float d01 = texelFetch(u_tex_linear_depth, base_coord + ivec2(0,1), u_src_lod).x;
	float d10 = texelFetch(u_tex_linear_depth, base_coord + ivec2(1,0), u_src_lod).x;
	float d11 = texelFetch(u_tex_linear_depth, base_coord + ivec2(1,1), u_src_lod).x;

	out_frag = vec4(OP(OP(d00, d01), OP(d10, d11)));
}
)");

static GLuint setup_program_from_source(str_t name, str_t f_shader_src, str_t defines = {}) {
    GLuint f_shader = gl::compile_shader_from_source(f_shader_src, GL_FRAGMENT_SHADER, defines);
    GLuint program = 0;

    if (f_shader) {
        char buffer[1024];
        program = glCreateProgram();

        glAttachShader(program, gl.v_shader_fs_quad);
        glAttachShader(program, f_shader);
        glLinkProgram(program);
        if (gl::get_program_link_error(buffer, sizeof(buffer), program)) {
            MD_LOG_ERROR("Error while linking %.*s program:\n%s", (int)name.len, name.ptr, buffer);
            glDeleteProgram(program);
            return 0;
        }

        glDetachShader(program, gl.v_shader_fs_quad);
        glDetachShader(program, f_shader);
        glDeleteShader(f_shader);
    }

    return program;
}

namespace ssao {
#ifndef AO_RANDOM_TEX_SIZE
#define AO_RANDOM_TEX_SIZE 4
#endif

struct HBAOData {
    float radius_to_screen;
    float neg_inv_r2;
    float n_dot_v_bias;
    float z_max;

    vec2_t inv_full_res;
    float ao_multiplier;
    float pow_exponent;

    vec4_t proj_info;

    vec4_t sample_pattern[32];
};

void setup_ubo_hbao_data(HBAOData* data, int width, int height, const float proj_mat[4][4], float intensity, float radius, float bias) {
    ASSERT(data);

    // From intel ASSAO
    static constexpr float SAMPLE_PATTERN[] = {
        0.78488064,  0.56661671,  1.500000, -0.126083,
        0.26022232,  -0.29575172, 1.500000, -1.064030,
        0.10459357,  0.08372527,  1.110000, -2.730563,
        -0.68286800, 0.04963045,  1.090000, -0.498827,
        -0.13570161, -0.64190155, 1.250000, -0.532765,
        -0.26193795, -0.08205118, 0.670000, -1.783245,
        -0.61177456, 0.66664219,  0.710000, -0.044234,
        0.43675563,  0.25119025,  0.610000, -1.167283,
        0.07884444,  0.86618668,  0.640000, -0.459002,
        -0.12790935, -0.29869005, 0.600000, -1.729424,
        -0.04031125, 0.02413622,  0.600000, -4.792042,
        0.16201244,  -0.52851415, 0.790000, -1.067055,
        -0.70991218, 0.47301072,  0.640000, -0.335236,
        0.03277707,  -0.22349690, 0.600000, -1.982384,
        0.68921727,  0.36800742,  0.630000, -0.266718,
        0.29251814,  0.37775412,  0.610000, -1.422520,
        -0.12224089, 0.96582592,  0.600000, -0.426142,
        0.11071457,  -0.16131058, 0.600000, -2.165947,
        0.46562141,  -0.59747696, 0.600000, -0.189760,
        -0.51548797, 0.11804193,  0.600000, -1.246800,
        0.89141309,  -0.42090443, 0.600000,  0.028192, 
        -0.32402530, -0.01591529, 0.600000, -1.543018,
        0.60771245,  0.41635221,  0.600000, -0.605411,
        0.02379565,  -0.08239821, 0.600000, -3.809046,
        0.48951152,  -0.23657045, 0.600000, -1.189011,
        -0.17611565, -0.81696892, 0.600000, -0.513724,
        -0.33930185, -0.20732205, 0.600000, -1.698047,
        -0.91974425, 0.05403209,  0.600000,  0.062246,
        -0.15064627, -0.14949332, 0.600000, -1.896062,
        0.53180975,  -0.35210401, 0.600000, -0.758838,
        0.41487166,  0.81442589,  0.600000, -0.505648,
        -0.24106961, -0.32721516, 0.600000, -1.665244};
    constexpr float METERS_TO_VIEWSPACE = 1.f;

    vec4_t proj_info;
    float proj_scl;
    float z_max;

    if (!is_orthographic_proj_matrix(proj_mat)) {
        // Persp
        proj_info = {
            2.0f / (proj_mat[0][0]),                    // (x) * (R - L)/N
            2.0f / (proj_mat[1][1]),                    // (y) * (T - B)/N
            -(1.0f - proj_mat[2][0]) / proj_mat[0][0],  // L/N
            -(1.0f + proj_mat[2][1]) / proj_mat[1][1]   // B/N
        };

        // proj_scl = float(height) / (math::tan(fovy * 0.5f) * 2.0f);
        proj_scl = float(height) * proj_mat[1][1] * 0.5f;
        z_max = (float)(proj_mat[3][2] / (proj_mat[2][2] + 1.0));
    } else {
        // Ortho
        proj_info = {
            2.0f / (proj_mat[0][0]),                    // ((x) * R - L)
            2.0f / (proj_mat[1][1]),                    // ((y) * T - B)
            -(1.0f + proj_mat[3][0]) / proj_mat[0][0],  // L
            -(1.0f - proj_mat[3][1]) / proj_mat[1][1]   // B
        };
        proj_scl = float(height) / proj_info.y;
        z_max = (float)((-2.0 + proj_mat[3][2]) / proj_mat[2][2]);
    }

    float r = radius * METERS_TO_VIEWSPACE;

    data->radius_to_screen = r * 0.5f * proj_scl;
    data->neg_inv_r2 = -1.f / (r * r);
    data->n_dot_v_bias = CLAMP(bias, 0.f, 1.f - FLT_EPSILON);
    data->z_max = z_max * 0.99f;
    data->inv_full_res = vec2_set(1.f / (float)width, 1.f / (float)height);
    data->ao_multiplier = 1.f / (1.f - data->n_dot_v_bias);
    data->pow_exponent = MAX(intensity, 0.f);
    data->proj_info = proj_info;
    MEMCPY(&data->sample_pattern, SAMPLE_PATTERN, sizeof(SAMPLE_PATTERN));

}

void initialize_rnd_tex(GLuint rnd_tex) {
    const int buffer_size = AO_RANDOM_TEX_SIZE * AO_RANDOM_TEX_SIZE;
    int16_t buffer[buffer_size * 4];

    for (int i = 0; i < buffer_size; i++) {
        double rand1 = bluenoise_64x64_data[i * 3 + 0] / (double)UINT16_MAX;
        double rand2 = bluenoise_64x64_data[i * 3 + 1] / (double)UINT16_MAX;
        double rand3 = bluenoise_64x64_data[i * 3 + 2] / (double)UINT16_MAX;
        double angle = TWO_PI * rand1;

        buffer[i * 4 + 0] = (int16_t)(INT16_MAX * cos(angle));
        buffer[i * 4 + 1] = (int16_t)(INT16_MAX * sin(angle));
        buffer[i * 4 + 2] = (int16_t)(INT16_MAX * rand2);
        buffer[i * 4 + 3] = (int16_t)(INT16_MAX * rand3);
    }

    glBindTexture(GL_TEXTURE_2D, rnd_tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16_SNORM, AO_RANDOM_TEX_SIZE, AO_RANDOM_TEX_SIZE, 0, GL_RGBA, GL_SHORT, buffer);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glBindTexture(GL_TEXTURE_2D, 0);
}

float compute_sharpness(float radius) { return 4.0f / sqrtf(radius); }

void initialize(int width, int height) {
    gl.ssao.hbao.program_persp = setup_program_from_source(STR_LIT("ssao persp"), {(const char*)ssao_frag, ssao_frag_size}, STR_LIT("#define AO_PERSPECTIVE 1"));
    gl.ssao.hbao.program_ortho = setup_program_from_source(STR_LIT("ssao ortho"), {(const char*)ssao_frag, ssao_frag_size}, STR_LIT("#define AO_PERSPECTIVE 0"));
    gl.ssao.blur.program       = setup_program_from_source(STR_LIT("ssao blur"),  {(const char*)blur_frag, blur_frag_size});
    
    if (!gl.ssao.fbo) glGenFramebuffers(1, &gl.ssao.fbo);

    if (!gl.ssao.tex_random) glGenTextures(1, &gl.ssao.tex_random);
    if (!gl.ssao.tex[0])     glGenTextures(2, gl.ssao.tex);

    if (!gl.ssao.ubo_hbao_data) glGenBuffers(1, &gl.ssao.ubo_hbao_data);

    initialize_rnd_tex(gl.ssao.tex_random);

    glBindTexture(GL_TEXTURE_2D, gl.ssao.tex[0]);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R8, width, height, 0, GL_RED, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindTexture(GL_TEXTURE_2D, gl.ssao.tex[1]);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R8, width, height, 0, GL_RED, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindTexture(GL_TEXTURE_2D, 0);

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.ssao.fbo);
    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.ssao.tex[0], 0);
    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, gl.ssao.tex[1], 0);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);

    glBindBuffer(GL_UNIFORM_BUFFER, gl.ssao.ubo_hbao_data);
    glBufferData(GL_UNIFORM_BUFFER, sizeof(HBAOData), nullptr, GL_DYNAMIC_DRAW);
}

void shutdown() {
    if (gl.ssao.fbo) glDeleteFramebuffers(1, &gl.ssao.fbo);
    if (gl.ssao.tex_random) glDeleteTextures(1, &gl.ssao.tex_random);
    if (gl.ssao.tex[0]) glDeleteTextures(2, gl.ssao.tex);
    if (gl.ssao.ubo_hbao_data) glDeleteBuffers(1, &gl.ssao.ubo_hbao_data);
    if (gl.ssao.hbao.program_persp) glDeleteProgram(gl.ssao.hbao.program_persp);
    if (gl.ssao.hbao.program_ortho) glDeleteProgram(gl.ssao.hbao.program_ortho);
    if (gl.ssao.blur.program) glDeleteProgram(gl.ssao.blur.program);
}

}  // namespace ssao

namespace fxaa {
void initialize() {
    gl.luma.program = setup_program_from_source(STR_LIT("luma"), {(const char*)luma_frag, luma_frag_size});
    gl.luma.uniform_loc.tex_rgba = glGetUniformLocation(gl.luma.program, "u_tex_rgba");

    str_t defines = STR_LIT("#define FXAA_PC 1\n#define FXAA_GLSL_130 1\n#define FXAA_QUALITY__PRESET 12");
    gl.fxaa.program = setup_program_from_source(STR_LIT("fxaa"), {(const char*)fxaa_frag, fxaa_frag_size}, defines);
    gl.fxaa.uniform_loc.tex_rgbl = glGetUniformLocation(gl.fxaa.program, "u_tex_rgbl");
    gl.fxaa.uniform_loc.rcp_res  = glGetUniformLocation(gl.fxaa.program, "u_rcp_res");
    gl.fxaa.uniform_loc.tc_scl   = glGetUniformLocation(gl.fxaa.program, "u_tc_scl");

}

void shutdown() {
    if (gl.luma.program) glDeleteProgram(gl.luma.program);
    if (gl.fxaa.program) glDeleteProgram(gl.fxaa.program);
}
}

namespace highlight {

static struct {
    GLuint program = 0;
    GLuint selection_texture = 0;
    struct {
        GLint texture_atom_idx = -1;
        GLint buffer_selection = -1;
        GLint highlight = -1;
        GLint selection = -1;
        GLint outline = -1;
    } uniform_loc;
} highlight;

void initialize() {
    highlight.program = setup_program_from_source(STR_LIT("highlight"), {(const char*)highlight_frag, highlight_frag_size});
    if (!highlight.selection_texture) glGenTextures(1, &highlight.selection_texture);
    highlight.uniform_loc.texture_atom_idx = glGetUniformLocation(highlight.program, "u_texture_atom_idx");
    highlight.uniform_loc.buffer_selection = glGetUniformLocation(highlight.program, "u_buffer_selection");
    highlight.uniform_loc.highlight = glGetUniformLocation(highlight.program, "u_highlight");
    highlight.uniform_loc.selection = glGetUniformLocation(highlight.program, "u_selection");
    highlight.uniform_loc.outline = glGetUniformLocation(highlight.program, "u_outline");
}

void shutdown() {
    if (highlight.program) glDeleteProgram(highlight.program);
}
}  // namespace highlight

namespace hsv {

static struct {
    GLuint program = 0;
    struct {
        GLint texture_color = -1;
        GLint hsv_scale = -1;
    } uniform_loc;
} gl;

void initialize() {
    gl.program = setup_program_from_source(STR_LIT("scale hsv"), {(const char*)scale_hsv_frag, scale_hsv_frag_size});
    gl.uniform_loc.texture_color = glGetUniformLocation(gl.program, "u_texture_atom_color");
    gl.uniform_loc.hsv_scale = glGetUniformLocation(gl.program, "u_hsv_scale");
}

void shutdown() {
    if (gl.program) glDeleteProgram(gl.program);
}
}  // namespace hsv

namespace compose {

struct ubo_data_t {
    mat4_t inv_proj_mat;
    vec3_t bg_color;
    float  time;
    vec3_t env_radiance;
    float  roughness;
    vec3_t dir_radiance;
    float  F0;
    vec3_t light_dir;
};

static struct {
    GLuint program = 0;
    GLuint ubo = 0;
    struct {
        GLint uniform_data = -1;
        GLint texture_depth = -1;
        GLint texture_color = -1;
        GLint texture_normal = -1;
    } uniform_loc;
} compose;

void initialize() {
    compose.program = setup_program_from_source(STR_LIT("compose deferred"), {(const char*)compose_deferred_frag, compose_deferred_frag_size});
    
    if (compose.ubo == 0) {
        glGenBuffers(1, &compose.ubo);
        glBindBuffer(GL_UNIFORM_BUFFER, compose.ubo);
        glBufferData(GL_UNIFORM_BUFFER, sizeof(ubo_data_t), 0, GL_DYNAMIC_DRAW);
        glBindBuffer(GL_UNIFORM_BUFFER, 0);
    }

    compose.uniform_loc.texture_depth       = glGetUniformLocation(compose.program, "u_texture_depth");
    compose.uniform_loc.texture_color       = glGetUniformLocation(compose.program, "u_texture_color");
    compose.uniform_loc.texture_normal      = glGetUniformLocation(compose.program, "u_texture_normal");
    compose.uniform_loc.uniform_data        = glGetUniformBlockIndex(compose.program, "UniformData");
}

void shutdown() {
    if (compose.program)    glDeleteProgram(compose.program);
    if (compose.ubo)        glDeleteBuffers(1, &compose.ubo);
}
}  // namespace compose

namespace tonemapping {

static struct {
    GLuint program = 0;
    struct {
        GLint texture = -1;
    } uniform_loc;
} passthrough;

static struct {
    GLuint program = 0;
    struct {
        GLint texture = -1;
        GLint exposure = -1;
        GLint gamma = -1;
    } uniform_loc;
} exposure_gamma;

static struct {
    GLuint program = 0;
    struct {
        GLint texture = -1;
        GLint exposure = -1;
        GLint gamma = -1;
    } uniform_loc;
} filmic;

static struct {
    GLuint program = 0;
    struct {
        GLint texture = -1;
        GLint exposure = -1;
        GLint gamma = -1;
    } uniform_loc;
} aces;

static struct {
    GLuint program_forward = 0;
    GLuint program_inverse = 0;
    struct {
        GLint texture = -1;
    } uniform_loc;
} fast_reversible;

void initialize() {
    {
        // PASSTHROUGH
        passthrough.program = setup_program_from_source(STR_LIT("Passthrough"), {(const char*)passthrough_frag, passthrough_frag_size});
        passthrough.uniform_loc.texture = glGetUniformLocation(passthrough.program, "u_texture");
    }
    {
        // EXPOSURE GAMMA
        exposure_gamma.program = setup_program_from_source(STR_LIT("Exposure Gamma"), {(const char*)exposure_gamma_frag, exposure_gamma_frag_size});
        exposure_gamma.uniform_loc.texture = glGetUniformLocation(exposure_gamma.program, "u_texture");
        exposure_gamma.uniform_loc.exposure = glGetUniformLocation(exposure_gamma.program, "u_exposure");
        exposure_gamma.uniform_loc.gamma = glGetUniformLocation(exposure_gamma.program, "u_gamma");
    }
    {
        // UNCHARTED
        filmic.program = setup_program_from_source(STR_LIT("Filmic"), {(const char*)uncharted_frag, uncharted_frag_size});
        filmic.uniform_loc.texture = glGetUniformLocation(filmic.program, "u_texture");
        filmic.uniform_loc.exposure = glGetUniformLocation(filmic.program, "u_exposure");
        filmic.uniform_loc.gamma = glGetUniformLocation(filmic.program, "u_gamma");
    }
    {
        // ACES
        aces.program = setup_program_from_source(STR_LIT("ACES"), {(const char*)aces_frag, aces_frag_size});
        aces.uniform_loc.texture = glGetUniformLocation(filmic.program, "u_texture");
        aces.uniform_loc.exposure = glGetUniformLocation(filmic.program, "u_exposure");
        aces.uniform_loc.gamma = glGetUniformLocation(filmic.program, "u_gamma");
    }
    {
        // Fast Reversible (For AA) (Credits to Brian Karis: http://graphicrants.blogspot.com/2013/12/tone-mapping.html)
        fast_reversible.program_forward = setup_program_from_source(STR_LIT("Fast Reversible"), {(const char*)fast_reversible_frag, fast_reversible_frag_size}, STR_LIT("#define USE_INVERSE 0"));
        fast_reversible.program_inverse = setup_program_from_source(STR_LIT("Fast Reversible"), {(const char*)fast_reversible_frag, fast_reversible_frag_size}, STR_LIT("#define USE_INVERSE 1"));
        fast_reversible.uniform_loc.texture = glGetUniformLocation(fast_reversible.program_forward, "u_texture");
    }
}

void shutdown() {
    if (passthrough.program) glDeleteProgram(passthrough.program);
    if (exposure_gamma.program) glDeleteProgram(exposure_gamma.program);
    if (filmic.program) glDeleteProgram(filmic.program);
    if (aces.program) glDeleteProgram(aces.program);
    if (fast_reversible.program_forward) glDeleteProgram(fast_reversible.program_forward);
    if (fast_reversible.program_inverse) glDeleteProgram(fast_reversible.program_inverse);
}

}  // namespace tonemapping

namespace dof {
void initialize(int32_t width, int32_t height) {
    {
        gl.bokeh_dof.half_res.program = setup_program_from_source(STR_LIT("DOF prepass"), {(const char*)dof_half_res_prepass_frag, dof_half_res_prepass_frag_size});
        if (gl.bokeh_dof.half_res.program) {
            gl.bokeh_dof.half_res.uniform_loc.tex_depth   = glGetUniformLocation(gl.bokeh_dof.half_res.program, "u_tex_depth");
            gl.bokeh_dof.half_res.uniform_loc.tex_color   = glGetUniformLocation(gl.bokeh_dof.half_res.program, "u_tex_color");
            gl.bokeh_dof.half_res.uniform_loc.focus_point = glGetUniformLocation(gl.bokeh_dof.half_res.program, "u_focus_point");
            gl.bokeh_dof.half_res.uniform_loc.focus_scale = glGetUniformLocation(gl.bokeh_dof.half_res.program, "u_focus_scale");
        }
    }

    if (!gl.bokeh_dof.half_res.tex.color_coc) {
        glGenTextures(1, &gl.bokeh_dof.half_res.tex.color_coc);
    }
    glBindTexture(GL_TEXTURE_2D, gl.bokeh_dof.half_res.tex.color_coc);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width / 2, height / 2, 0, GL_RGBA, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D, 0);

    if (!gl.bokeh_dof.half_res.fbo) {
        glGenFramebuffers(1, &gl.bokeh_dof.half_res.fbo);
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.bokeh_dof.half_res.fbo);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.bokeh_dof.half_res.tex.color_coc, 0);
        GLenum status = glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER);
        if (status != GL_FRAMEBUFFER_COMPLETE) {
            MD_LOG_ERROR("Something went wrong when generating framebuffer for DOF");
        }
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
    }

    // DOF
    {
        gl.bokeh_dof.program = setup_program_from_source(STR_LIT("Bokeh DOF"), {(const char*)dof_frag, dof_frag_size});
        if (gl.bokeh_dof.program) {
            gl.bokeh_dof.uniform_loc.tex_color = glGetUniformLocation(gl.bokeh_dof.program, "u_half_res");
            gl.bokeh_dof.uniform_loc.tex_color = glGetUniformLocation(gl.bokeh_dof.program, "u_tex_color");
            gl.bokeh_dof.uniform_loc.tex_depth = glGetUniformLocation(gl.bokeh_dof.program, "u_tex_depth");
            gl.bokeh_dof.uniform_loc.pixel_size = glGetUniformLocation(gl.bokeh_dof.program, "u_texel_size");
            gl.bokeh_dof.uniform_loc.focus_point = glGetUniformLocation(gl.bokeh_dof.program, "u_focus_depth");
            gl.bokeh_dof.uniform_loc.focus_scale = glGetUniformLocation(gl.bokeh_dof.program, "u_focus_scale");
            gl.bokeh_dof.uniform_loc.time = glGetUniformLocation(gl.bokeh_dof.program, "u_time");
        }
    }
}

void shutdown() {}
}  // namespace dof

namespace blit {
static GLuint program_tex = 0;
static GLuint program_col = 0;
static GLint uniform_loc_texture = -1;
static GLint uniform_loc_color = -1;

constexpr str_t f_shader_src_tex = STR_LIT(R"(
#version 150 core

uniform sampler2D u_texture;

out vec4 out_frag;

void main() {
    out_frag = texelFetch(u_texture, ivec2(gl_FragCoord.xy), 0);
}
)");

constexpr str_t f_shader_src_col = STR_LIT(R"(
#version 150 core

uniform vec4 u_color;
out vec4 out_frag;

void main() {
	out_frag = u_color;
}
)");

void initialize() {
    program_tex = setup_program_from_source(STR_LIT("blit texture"), f_shader_src_tex);
    uniform_loc_texture = glGetUniformLocation(program_tex, "u_texture");

    program_col = setup_program_from_source(STR_LIT("blit color"), f_shader_src_col);
    uniform_loc_color = glGetUniformLocation(program_col, "u_color");
}

void shutdown() {
    if (program_tex) glDeleteProgram(program_tex);
    if (program_col) glDeleteProgram(program_col);
}
}  // namespace blit

namespace blur {
static GLuint program_gaussian = 0;
static GLuint program_box = 0;
static GLint uniform_loc_texture = -1;
static GLint uniform_loc_inv_res_dir = -1;

constexpr str_t f_shader_src_gaussian = STR_LIT(R"(
#version 150 core

#define KERNEL_RADIUS 5

uniform sampler2D u_texture;
uniform vec2      u_inv_res_dir;

in vec2 tc;
out vec4 out_frag;

float blur_weight(float r) {
    const float sigma = KERNEL_RADIUS * 0.5;
    const float falloff = 1.0 / (2.0*sigma*sigma);
    float w = exp2(-r*r*falloff);
    return w;
}

void main() {
    vec2 uv = tc;
    vec4  c_tot = texture(u_texture, uv);
    float w_tot = 1.0;

    for (float r = 1; r <= KERNEL_RADIUS; ++r) {
        float w = blur_weight(r);
        vec4  c = texture(u_texture, uv + u_inv_res_dir * r);
        c_tot += c * w;
        w_tot += w;
    }
    for (float r = 1; r <= KERNEL_RADIUS; ++r) {
        float w = blur_weight(r);
        vec4  c = texture(u_texture, uv - u_inv_res_dir * r);
        c_tot += c * w;
        w_tot += w;
    }

    out_frag = c_tot / w_tot;
}
)");

constexpr str_t f_shader_src_box = STR_LIT(R"(
#version 150 core

uniform sampler2D u_texture;
out vec4 out_frag;

void main() {
    vec4 c = vec4(0);
    c += texelFetch(u_texture, ivec2(gl_FragCoord.xy) + ivec2(-1, -1), 0);
    c += texelFetch(u_texture, ivec2(gl_FragCoord.xy) + ivec2( 0, -1), 0);
    c += texelFetch(u_texture, ivec2(gl_FragCoord.xy) + ivec2(+1, -1), 0);
    c += texelFetch(u_texture, ivec2(gl_FragCoord.xy) + ivec2(-1,  0), 0);
    c += texelFetch(u_texture, ivec2(gl_FragCoord.xy) + ivec2( 0,  0), 0);
    c += texelFetch(u_texture, ivec2(gl_FragCoord.xy) + ivec2(+1,  0), 0);
    c += texelFetch(u_texture, ivec2(gl_FragCoord.xy) + ivec2(-1, +1), 0);
    c += texelFetch(u_texture, ivec2(gl_FragCoord.xy) + ivec2( 0, +1), 0);
    c += texelFetch(u_texture, ivec2(gl_FragCoord.xy) + ivec2(+1, +1), 0);

    out_frag = c / 9.0;
}
)");

void initialize() {
    program_gaussian = setup_program_from_source(STR_LIT("gaussian blur"), f_shader_src_gaussian);
    uniform_loc_texture = glGetUniformLocation(program_gaussian, "u_texture");
    uniform_loc_inv_res_dir = glGetUniformLocation(program_gaussian, "u_inv_res_dir");

    program_box = setup_program_from_source(STR_LIT("box blur"), f_shader_src_box);
}

void shutdown() {
    if (program_gaussian) glDeleteProgram(program_gaussian);
    if (program_box) glDeleteProgram(program_box);
}
}  // namespace blit

namespace velocity {
#define VEL_TILE_SIZE 8

struct {
    GLuint program = 0;
    struct {
		GLint tex_depth = -1;
        GLint curr_clip_to_prev_clip_mat = -1;
        GLint jitter_uv = -1;
    } uniform_loc;
} blit_velocity;

struct {
    GLuint program = 0;
    struct {
        GLint tex_vel = -1;
        GLint tex_linear_depth = -1;
    } uniform_loc;
} blit_tilemax;

struct {
    GLuint program = 0;
    struct {
        GLint tex_vel = -1;
        GLint tex_linear_depth = -1;
        GLint tex_vel_texel_size = -1;
    } uniform_loc;
} blit_neighbormax;

struct {
    GLuint program = 0;
    struct {
        GLint tex_vel = -1;
        GLint texel_size = -1;
    } uniform_loc;
} blit_dilate;

void initialize(int32_t width, int32_t height) {
    {
        blit_velocity.program = setup_program_from_source(STR_LIT("screen-space velocity"), {(const char*)blit_velocity_frag, blit_velocity_frag_size});
		blit_velocity.uniform_loc.tex_depth = glGetUniformLocation(blit_velocity.program, "u_tex_depth");
        blit_velocity.uniform_loc.curr_clip_to_prev_clip_mat = glGetUniformLocation(blit_velocity.program, "u_curr_clip_to_prev_clip_mat");
        blit_velocity.uniform_loc.jitter_uv = glGetUniformLocation(blit_velocity.program, "u_jitter_uv");

    }
    {
        str_t defines = STR_LIT("#define TILE_SIZE " STRINGIFY_VAL(VEL_TILE_SIZE));
        blit_tilemax.program = setup_program_from_source(STR_LIT("tilemax"), {(const char*)blit_tilemax_frag, blit_tilemax_frag_size}, defines);
        blit_tilemax.uniform_loc.tex_vel = glGetUniformLocation(blit_tilemax.program, "u_tex_vel");
        blit_tilemax.uniform_loc.tex_linear_depth = glGetUniformLocation(blit_tilemax.program, "u_tex_linear_depth");
    }
    {
        blit_neighbormax.program = setup_program_from_source(STR_LIT("neighbormax"), {(const char*)blit_neighbormax_frag, blit_neighbormax_frag_size});
        blit_neighbormax.uniform_loc.tex_vel = glGetUniformLocation(blit_neighbormax.program, "u_tex_vel");
        blit_neighbormax.uniform_loc.tex_linear_depth = glGetUniformLocation(blit_neighbormax.program, "u_tex_linear_depth");
        blit_neighbormax.uniform_loc.tex_vel_texel_size = glGetUniformLocation(blit_neighbormax.program, "u_tex_vel_texel_size");
    }
    {
        blit_dilate.program = setup_program_from_source(STR_LIT("velocity dilate"), { (const char*)blit_velocity_dilate_frag, blit_velocity_dilate_frag_size });
        blit_dilate.uniform_loc.tex_vel = glGetUniformLocation(blit_dilate.program, "u_tex_vel");
        blit_dilate.uniform_loc.texel_size = glGetUniformLocation(blit_dilate.program, "u_texel_size");
    }

    if (!gl.velocity.tex_tilemax) {
        glGenTextures(1, &gl.velocity.tex_tilemax);
    }

    if (!gl.velocity.tex_neighbormax) {
        glGenTextures(1, &gl.velocity.tex_neighbormax);
    }

    gl.velocity.tex_width  = width  / VEL_TILE_SIZE;
    gl.velocity.tex_height = height / VEL_TILE_SIZE;

    glBindTexture(GL_TEXTURE_2D, gl.velocity.tex_tilemax);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RG16F, gl.velocity.tex_width, gl.velocity.tex_height, 0, GL_RG, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D, 0);

    glBindTexture(GL_TEXTURE_2D, gl.velocity.tex_neighbormax);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RG16F, gl.velocity.tex_width, gl.velocity.tex_height, 0, GL_RG, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D, 0);

    if (!gl.velocity.fbo) {
        glGenFramebuffers(1, &gl.velocity.fbo);
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.velocity.fbo);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.velocity.tex_tilemax, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, gl.velocity.tex_neighbormax, 0);
        GLenum status = glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER);
        if (status != GL_FRAMEBUFFER_COMPLETE) {
            MD_LOG_ERROR("Something went wrong in creating framebuffer for velocity");
        }
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
    }
}

void shutdown() {
    if (blit_velocity.program) glDeleteProgram(blit_velocity.program);
    if (blit_tilemax.program) glDeleteProgram(blit_tilemax.program);
    if (blit_neighbormax.program) glDeleteProgram(blit_neighbormax.program);
    if (blit_dilate.program) glDeleteProgram(blit_dilate.program);
    if (gl.velocity.tex_tilemax) glDeleteTextures(1, &gl.velocity.tex_tilemax);
    if (gl.velocity.tex_neighbormax) glDeleteTextures(1, &gl.velocity.tex_neighbormax);
    if (gl.velocity.fbo) glDeleteFramebuffers(1, &gl.velocity.fbo);
}
}  // namespace velocity

namespace temporal {
void initialize() {
    {
        gl.temporal.with_motion_blur.program = setup_program_from_source(STR_LIT("temporal aa + motion-blur"), {(const char*)temporal_frag, temporal_frag_size});
        gl.temporal.no_motion_blur.program   = setup_program_from_source(STR_LIT("temporal aa"), {(const char*)temporal_frag, temporal_frag_size}, STR_LIT("#define USE_MOTION_BLUR 0\n"));

        gl.temporal.with_motion_blur.uniform_loc.tex_linear_depth = glGetUniformLocation(gl.temporal.with_motion_blur.program, "u_tex_linear_depth");
        gl.temporal.with_motion_blur.uniform_loc.tex_main = glGetUniformLocation(gl.temporal.with_motion_blur.program, "u_tex_main");
        gl.temporal.with_motion_blur.uniform_loc.tex_prev = glGetUniformLocation(gl.temporal.with_motion_blur.program, "u_tex_prev");
        gl.temporal.with_motion_blur.uniform_loc.tex_vel = glGetUniformLocation(gl.temporal.with_motion_blur.program, "u_tex_vel");
        gl.temporal.with_motion_blur.uniform_loc.tex_vel_neighbormax = glGetUniformLocation(gl.temporal.with_motion_blur.program, "u_tex_vel_neighbormax");
        gl.temporal.with_motion_blur.uniform_loc.texel_size = glGetUniformLocation(gl.temporal.with_motion_blur.program, "u_texel_size");
        gl.temporal.with_motion_blur.uniform_loc.jitter_uv = glGetUniformLocation(gl.temporal.with_motion_blur.program, "u_jitter_uv");
        gl.temporal.with_motion_blur.uniform_loc.time = glGetUniformLocation(gl.temporal.with_motion_blur.program, "u_time");
        gl.temporal.with_motion_blur.uniform_loc.feedback_min = glGetUniformLocation(gl.temporal.with_motion_blur.program, "u_feedback_min");
        gl.temporal.with_motion_blur.uniform_loc.feedback_max = glGetUniformLocation(gl.temporal.with_motion_blur.program, "u_feedback_max");
        gl.temporal.with_motion_blur.uniform_loc.motion_scale = glGetUniformLocation(gl.temporal.with_motion_blur.program, "u_motion_scale");

        gl.temporal.no_motion_blur.uniform_loc.tex_linear_depth = glGetUniformLocation(gl.temporal.no_motion_blur.program, "u_tex_linear_depth");
        gl.temporal.no_motion_blur.uniform_loc.tex_main = glGetUniformLocation(gl.temporal.no_motion_blur.program, "u_tex_main");
        gl.temporal.no_motion_blur.uniform_loc.tex_prev = glGetUniformLocation(gl.temporal.no_motion_blur.program, "u_tex_prev");
        gl.temporal.no_motion_blur.uniform_loc.tex_vel = glGetUniformLocation(gl.temporal.no_motion_blur.program, "u_tex_vel");
        gl.temporal.no_motion_blur.uniform_loc.texel_size = glGetUniformLocation(gl.temporal.no_motion_blur.program, "u_texel_size");
        gl.temporal.no_motion_blur.uniform_loc.jitter_uv = glGetUniformLocation(gl.temporal.no_motion_blur.program, "u_jitter_uv");
        gl.temporal.no_motion_blur.uniform_loc.time = glGetUniformLocation(gl.temporal.no_motion_blur.program, "u_time");
        gl.temporal.no_motion_blur.uniform_loc.feedback_min = glGetUniformLocation(gl.temporal.no_motion_blur.program, "u_feedback_min");
        gl.temporal.no_motion_blur.uniform_loc.feedback_max = glGetUniformLocation(gl.temporal.no_motion_blur.program, "u_feedback_max");
        gl.temporal.no_motion_blur.uniform_loc.motion_scale = glGetUniformLocation(gl.temporal.no_motion_blur.program, "u_motion_scale");
    }
}

void shutdown() {}
}  // namespace temporal

namespace sharpen {
static GLuint program = 0;
void initialize() {
    constexpr str_t f_shader_src_sharpen = STR_LIT(
 R"(#version 150 core

    uniform sampler2D u_tex;
    uniform float u_weight;
    out vec4 out_frag;

    void main() {
        vec3 cc = texelFetch(u_tex, ivec2(gl_FragCoord.xy), 0).rgb;
        vec3 cl = texelFetch(u_tex, ivec2(gl_FragCoord.xy) + ivec2(-1, 0), 0).rgb;
        vec3 ct = texelFetch(u_tex, ivec2(gl_FragCoord.xy) + ivec2( 0, 1), 0).rgb;
        vec3 cr = texelFetch(u_tex, ivec2(gl_FragCoord.xy) + ivec2( 1, 0), 0).rgb;
        vec3 cb = texelFetch(u_tex, ivec2(gl_FragCoord.xy) + ivec2( 0,-1), 0).rgb;

        float pos_weight = 1.0 + u_weight;
        float neg_weight = -u_weight * 0.25;
        out_frag = vec4(vec3(pos_weight * cc + neg_weight * (cl + ct + cr + cb)), 1.0);
    })");
    program = setup_program_from_source(STR_LIT("sharpen"), f_shader_src_sharpen);
}

void sharpen(GLuint in_texture, float weight = 1.0f) {
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, in_texture);

    glUseProgram(program);
    glUniform1i(glGetUniformLocation(sharpen::program, "u_tex"), 0);
    glUniform1f(glGetUniformLocation(sharpen::program, "u_weight"), weight);

    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);

    glUseProgram(0);
}

void shutdown() {
    if (program) glDeleteProgram(program);
}
}

void initialize(int width, int height) {
    if (!gl.vao) glGenVertexArrays(1, &gl.vao);

    gl.v_shader_fs_quad = gl::compile_shader_from_source(v_shader_src_fs_quad, GL_VERTEX_SHADER);

    // LINEARIZE DEPTH

    gl.linear_depth.linearize.program_persp = setup_program_from_source(STR_LIT("linearize depth persp"), f_shader_src_linearize_depth, STR_LIT("#version 150 core\n#define PERSPECTIVE 1"));
    gl.linear_depth.linearize.program_ortho = setup_program_from_source(STR_LIT("linearize depth ortho"), f_shader_src_linearize_depth, STR_LIT("#version 150 core\n#define PERSPECTIVE 0"));
    gl.linear_depth.downsample.program = setup_program_from_source(STR_LIT("linear depth downsample"), f_shader_src_downsample_depth);

    if (gl.linear_depth.texture)
        glDeleteTextures(1, &gl.linear_depth.texture);
    
    {
        glGenTextures(1, &gl.linear_depth.texture);
        glBindTexture(GL_TEXTURE_2D, gl.linear_depth.texture);
        glTexStorage2D(GL_TEXTURE_2D, 4, GL_R16F, width, height);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0);
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 3);
        glBindTexture(GL_TEXTURE_2D, 0);
    }

    {
        GLenum status = 0;
        if (!gl.linear_depth.fbo_0)
            glGenFramebuffers(1, &gl.linear_depth.fbo_0);

        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.linear_depth.fbo_0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.linear_depth.texture, 0);
        status = glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER);
        if (status != GL_FRAMEBUFFER_COMPLETE) {
            MD_LOG_ERROR("Something went wrong in creating framebuffer for depth linearization");
        }

        if (!gl.linear_depth.fbo_1)
            glGenFramebuffers(1, &gl.linear_depth.fbo_1);

        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.linear_depth.fbo_1);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.linear_depth.texture, 1);
        status = glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER);
        if (status != GL_FRAMEBUFFER_COMPLETE) {
            MD_LOG_ERROR("Something went wrong in creating framebuffer for half linear depth");
        }

        if (!gl.linear_depth.fbo_2)
            glGenFramebuffers(1, &gl.linear_depth.fbo_2);

        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.linear_depth.fbo_2);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.linear_depth.texture, 2);
        status = glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER);
        if (status != GL_FRAMEBUFFER_COMPLETE) {
            MD_LOG_ERROR("Something went wrong in creating framebuffer for quarter linear depth");
        }

        if (!gl.linear_depth.fbo_3)
            glGenFramebuffers(1, &gl.linear_depth.fbo_3);

        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.linear_depth.fbo_3);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.linear_depth.texture, 3);
        status = glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER);
        if (status != GL_FRAMEBUFFER_COMPLETE) {
            MD_LOG_ERROR("Something went wrong in creating framebuffer for eighth linear depth");
        }

        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
    }

    // COLOR
    if (!gl.targets.tex_color[0]) glGenTextures(2, gl.targets.tex_color);
    glBindTexture(GL_TEXTURE_2D, gl.targets.tex_color[0]);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R11F_G11F_B10F, width, height, 0, GL_RGB, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D, gl.targets.tex_color[1]);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R11F_G11F_B10F, width, height, 0, GL_RGB, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D, 0);

    if (!gl.targets.tex_temporal_buffer[0]) glGenTextures(2, gl.targets.tex_temporal_buffer);
    glBindTexture(GL_TEXTURE_2D, gl.targets.tex_temporal_buffer[0]);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB16, width, height, 0, GL_RGB, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D, gl.targets.tex_temporal_buffer[1]);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB16, width, height, 0, GL_RGB, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D, 0);

    if (!gl.targets.fbo) {
        glGenFramebuffers(1, &gl.targets.fbo);
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.targets.fbo);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.targets.tex_color[0], 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, gl.targets.tex_color[1], 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT2, GL_TEXTURE_2D, gl.targets.tex_temporal_buffer[0], 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT3, GL_TEXTURE_2D, gl.targets.tex_temporal_buffer[1], 0);

        GLenum status = glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER);
        if (status != GL_FRAMEBUFFER_COMPLETE) {
            MD_LOG_ERROR("Something went wrong in creating framebuffer for targets");
        }

        GLenum buffers[] = {GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2, GL_COLOR_ATTACHMENT3};
        glDrawBuffers(4, buffers);
        glClear(GL_COLOR_BUFFER_BIT);
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
    }

    if (!gl.tmp.tex_rgba8) glGenTextures(1, &gl.tmp.tex_rgba8);
    glBindTexture(GL_TEXTURE_2D, gl.tmp.tex_rgba8);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D, 0);

    if (!gl.tmp.fbo) {
        glGenFramebuffers(1, &gl.tmp.fbo);
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.tmp.fbo);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.tmp.tex_rgba8, 0);
    }

    gl.tex_width = width;
    gl.tex_height = height;

    ssao::initialize(width, height);
    dof::initialize(width, height);
    velocity::initialize(width, height);
    highlight::initialize();
    hsv::initialize();
    tonemapping::initialize();
    temporal::initialize();
    blit::initialize();
    blur::initialize();
    sharpen::initialize();
    compose::initialize();
    fxaa::initialize();
}

void shutdown() {
    ssao::shutdown();
    dof::shutdown();
    velocity::shutdown();
    highlight::shutdown();
    hsv::shutdown();
    tonemapping::shutdown();
    temporal::shutdown();
    blit::shutdown();
    blur::shutdown();
    sharpen::shutdown();
    compose::shutdown();
    fxaa::shutdown();

    if (gl.vao) glDeleteVertexArrays(1, &gl.vao);
    //if (gl.vbo) glDeleteBuffers(1, &gl.vbo);
    if (gl.v_shader_fs_quad) glDeleteShader(gl.v_shader_fs_quad);
    if (gl.tmp.fbo) glDeleteFramebuffers(1, &gl.tmp.fbo);
    if (gl.tmp.tex_rgba8) glDeleteTextures(1, &gl.tmp.tex_rgba8);
}

void compute_linear_depth(GLuint depth_tex, float near_plane, float far_plane, bool orthographic = false) {
    const vec4_t clip_info {near_plane * far_plane, near_plane - far_plane, far_plane, 0};

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, depth_tex);

    GLuint program = orthographic ? gl.linear_depth.linearize.program_ortho : gl.linear_depth.linearize.program_persp;
    glUseProgram(program);
    glUniform1i(glGetUniformLocation (program, "u_tex_depth"), 0);
    glUniform4fv(glGetUniformLocation(program, "u_clip_info"), 1, &clip_info.x);

    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
}

void downsample_depth(GLuint linear_depth_tex, int src_lod) {
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, linear_depth_tex);

    glUseProgram(gl.linear_depth.downsample.program);
    glUniform1i(glGetUniformLocation(gl.linear_depth.downsample.program, "u_tex_linear_depth"), 0);
    glUniform1i(glGetUniformLocation(gl.linear_depth.downsample.program, "u_src_lod"), src_lod);

    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
}

void compute_ssao(GLuint linear_depth_tex, GLuint normal_tex, const float proj_mat[4][4], float intensity, float radius, float bias) {
    ASSERT(glIsTexture(linear_depth_tex));
    ASSERT(glIsTexture(normal_tex));

    GLint last_fbo;
    GLint last_viewport[4];
    GLint last_draw_buffer;
    GLint last_scissor_box[4];
    glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &last_fbo);
    glGetIntegerv(GL_VIEWPORT, last_viewport);
    glGetIntegerv(GL_SCISSOR_BOX, last_scissor_box);
    glGetIntegerv(GL_DRAW_BUFFER, &last_draw_buffer);

    int width  = last_viewport[2];
    int height = last_viewport[3];

    const bool ortho = is_orthographic_proj_matrix(proj_mat);
    const float sharpness = ssao::compute_sharpness(radius);
    const vec2_t inv_res = vec2_t{ 1.f / (float)width, 1.f / (float)height };

    glBindVertexArray(gl.vao);

    ssao::HBAOData ubo_data = {};
    ssao::setup_ubo_hbao_data(&ubo_data, width, height, proj_mat, intensity, radius, bias);
    glBindBuffer(GL_UNIFORM_BUFFER, gl.ssao.ubo_hbao_data);
    glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(ssao::HBAOData), &ubo_data);
    glBindBuffer(GL_UNIFORM_BUFFER, 0);

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.ssao.fbo);
    glDrawBuffer(GL_COLOR_ATTACHMENT0);
    glViewport(0, 0, width, height);
    glScissor(0, 0, width, height);
    glClearColor(1,1,1,1);
    glClear(GL_COLOR_BUFFER_BIT);

    GLuint program = ortho ? gl.ssao.hbao.program_ortho : gl.ssao.hbao.program_persp;

    PUSH_GPU_SECTION("HBAO")
    glUseProgram(program);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, linear_depth_tex);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, normal_tex);
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, gl.ssao.tex_random);

    glBindBufferBase(GL_UNIFORM_BUFFER, 0, gl.ssao.ubo_hbao_data);
    glUniformBlockBinding(program, glGetUniformBlockIndex(program, "u_control_buffer"), 0);
    glUniform1i(glGetUniformLocation(program, "u_tex_linear_depth"), 0);
    glUniform1i(glGetUniformLocation(program, "u_tex_normal"), 1);
    glUniform1i(glGetUniformLocation(program, "u_tex_random"), 2);

    glDrawArrays(GL_TRIANGLES, 0, 3);
    POP_GPU_SECTION()

    PUSH_GPU_SECTION("BLUR");
    glUseProgram(gl.ssao.blur.program);

    glUniform1i(glGetUniformLocation(gl.ssao.blur.program, "u_tex_linear_depth"), 0);
    glUniform1i(glGetUniformLocation(gl.ssao.blur.program, "u_tex_ao"), 1);
    glUniform1f(glGetUniformLocation(gl.ssao.blur.program, "u_sharpness"), sharpness);
    glUniform1f(glGetUniformLocation(gl.ssao.blur.program, "u_zmax"), ubo_data.z_max);
    glUniform2f(glGetUniformLocation(gl.ssao.blur.program, "u_inv_res_dir"), inv_res.x, 0);

    glActiveTexture(GL_TEXTURE1);

    // BLUR FIRST
    PUSH_GPU_SECTION("1st")
    glDrawBuffer(GL_COLOR_ATTACHMENT1);
    glBindTexture(GL_TEXTURE_2D, gl.ssao.tex[0]);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    POP_GPU_SECTION()

    glUniform2f(glGetUniformLocation(gl.ssao.blur.program, "u_inv_res_dir"), 0, inv_res.y);

    glEnable(GL_BLEND);
    glBlendFunc(GL_ZERO, GL_SRC_COLOR);

    // BLUR SECOND AND BLEND RESULT
    PUSH_GPU_SECTION("2nd")
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, last_fbo);
    glDrawBuffer(last_draw_buffer);
    glViewport(last_viewport[0], last_viewport[1], last_viewport[2], last_viewport[3]);
    glScissor(last_scissor_box[0], last_scissor_box[1], last_scissor_box[2], last_scissor_box[3]);
    glBindTexture(GL_TEXTURE_2D, gl.ssao.tex[1]);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    POP_GPU_SECTION()

    glDisable(GL_BLEND);

    glBindVertexArray(0);
    POP_GPU_SECTION()
}

static void compose_deferred(GLuint depth_tex, GLuint color_tex, GLuint normal_tex, const mat4_t& inv_proj_matrix, const vec3_t bg_color, float time) {
    ASSERT(glIsTexture(depth_tex));
    ASSERT(glIsTexture(color_tex));
    ASSERT(glIsTexture(normal_tex));

    const vec3_t env_radiance = bg_color * 0.25f;
    const vec3_t dir_radiance = {10, 10, 10};
    const float roughness = 0.4f;
    const float F0 = 0.04f;
    const vec3_t L = {0.57735026918962576451f, 0.57735026918962576451f, 0.57735026918962576451f}; // 1.0 / sqrt(3)

    compose::ubo_data_t data = {
        .inv_proj_mat = inv_proj_matrix,
        .bg_color = bg_color,
        .time = time,
        .env_radiance = env_radiance,
        .roughness = roughness,
        .dir_radiance = dir_radiance,
        .F0 = F0,
        .light_dir = L,
    };

    glBindBuffer(GL_UNIFORM_BUFFER, compose::compose.ubo);
    glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(data), &data);
    glBindBuffer(GL_UNIFORM_BUFFER, 0);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, depth_tex);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, color_tex);
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, normal_tex);

    GLuint program = compose::compose.program;

    glUseProgram(program);

    glBindBufferBase(GL_UNIFORM_BUFFER, 0, compose::compose.ubo);
    glUniformBlockBinding(program, compose::compose.uniform_loc.uniform_data, 0);
    glUniform1i (compose::compose.uniform_loc.texture_depth, 0);
    glUniform1i (compose::compose.uniform_loc.texture_color, 1);
    glUniform1i (compose::compose.uniform_loc.texture_normal, 2);

    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);

    glUseProgram(0);
}

void highlight_selection(GLuint atom_idx_tex, GLuint selection_buffer, const vec3_t& highlight, const vec3_t& selection, const vec3_t& outline) {
    ASSERT(glIsTexture(atom_idx_tex));
    ASSERT(glIsBuffer(selection_buffer));

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, atom_idx_tex);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_BUFFER, highlight::highlight.selection_texture);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_R8UI, selection_buffer);

    glUseProgram(highlight::highlight.program);
    glUniform1i(highlight::highlight.uniform_loc.texture_atom_idx, 0);
    glUniform1i(highlight::highlight.uniform_loc.buffer_selection, 1);
    glUniform3fv(highlight::highlight.uniform_loc.highlight, 1, &highlight.x);
    glUniform3fv(highlight::highlight.uniform_loc.selection, 1, &selection.x);
    glUniform3fv(highlight::highlight.uniform_loc.outline, 1, &outline.x);
    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
}

void half_res_color_coc(GLuint linear_depth_tex, GLuint color_tex, float focus_point, float focus_scale) {
    PUSH_GPU_SECTION("DOF Prepass");
    GLint last_viewport[4];
    glGetIntegerv(GL_VIEWPORT, last_viewport);
    glViewport(0, 0, gl.tex_width / 2, gl.tex_height / 2);

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.bokeh_dof.half_res.fbo);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, linear_depth_tex);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, color_tex);

    glUseProgram(gl.bokeh_dof.half_res.program);

    glUniform1i(gl.bokeh_dof.half_res.uniform_loc.tex_depth, 0);
    glUniform1i(gl.bokeh_dof.half_res.uniform_loc.tex_color, 1);
    glUniform1f(gl.bokeh_dof.half_res.uniform_loc.focus_point, focus_point);
    glUniform1f(gl.bokeh_dof.half_res.uniform_loc.focus_scale, focus_scale);

    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);

    glViewport(last_viewport[0], last_viewport[1], last_viewport[2], last_viewport[3]);
    POP_GPU_SECTION();
}

void apply_dof(GLuint linear_depth_tex, GLuint color_tex, float focus_point, float focus_scale, float time) {
    ASSERT(glIsTexture(linear_depth_tex));
    ASSERT(glIsTexture(color_tex));

    const vec2_t pixel_size = vec2_t{1.f / gl.tex_width, 1.f / gl.tex_height};

    GLint last_fbo;
    glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &last_fbo);

    half_res_color_coc(linear_depth_tex, color_tex, focus_point, focus_scale);

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, last_fbo);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, gl.bokeh_dof.half_res.tex.color_coc);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, linear_depth_tex);
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, color_tex);

    glUseProgram(gl.bokeh_dof.program);
    glUniform1i(gl.bokeh_dof.uniform_loc.tex_half_res, 0);
    glUniform1i(gl.bokeh_dof.uniform_loc.tex_depth, 1);
    glUniform1i(gl.bokeh_dof.uniform_loc.tex_color, 2);
    glUniform2f(gl.bokeh_dof.uniform_loc.pixel_size, pixel_size.x, pixel_size.y);
    glUniform1f(gl.bokeh_dof.uniform_loc.focus_point, focus_point);
    glUniform1f(gl.bokeh_dof.uniform_loc.focus_scale, focus_scale);
    glUniform1f(gl.bokeh_dof.uniform_loc.time, time);

    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
    glActiveTexture(GL_TEXTURE0);
}

void apply_tonemapping(GLuint color_tex, Tonemapping tonemapping, float exposure, float gamma) {
    ASSERT(glIsTexture(color_tex));

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, color_tex);

    switch (tonemapping) {
        case Tonemapping_ExposureGamma:
            glUseProgram(tonemapping::exposure_gamma.program);
            glUniform1i(tonemapping::exposure_gamma.uniform_loc.texture, 0);
            glUniform1f(tonemapping::exposure_gamma.uniform_loc.exposure, exposure);
            glUniform1f(tonemapping::exposure_gamma.uniform_loc.gamma, gamma);
            break;
        case Tonemapping_Filmic:
            glUseProgram(tonemapping::filmic.program);
            glUniform1i(tonemapping::filmic.uniform_loc.texture, 0);
            glUniform1f(tonemapping::filmic.uniform_loc.exposure, exposure);
            glUniform1f(tonemapping::filmic.uniform_loc.gamma, gamma);
            break;
        case Tonemapping_ACES:
            glUseProgram(tonemapping::aces.program);
            glUniform1i(tonemapping::aces.uniform_loc.texture, 0);
            glUniform1f(tonemapping::aces.uniform_loc.exposure, exposure);
            glUniform1f(tonemapping::aces.uniform_loc.gamma, gamma);
            break;
        case Tonemapping_Passthrough:
        default:
            glUseProgram(tonemapping::passthrough.program);
            glUniform1i(tonemapping::passthrough.uniform_loc.texture, 0);
            break;
    }

    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
}

void apply_aa_tonemapping(GLuint color_tex) {
    ASSERT(glIsTexture(color_tex));

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, color_tex);

    glUseProgram(tonemapping::fast_reversible.program_forward);
    glUniform1i(tonemapping::fast_reversible.uniform_loc.texture, 0);

    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
}

void apply_inverse_aa_tonemapping(GLuint color_tex) {
    ASSERT(glIsTexture(color_tex));

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, color_tex);

    glUseProgram(tonemapping::fast_reversible.program_inverse);
    glUniform1i(tonemapping::fast_reversible.uniform_loc.texture, 0);

    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
}

void blit_static_velocity(GLuint depth_tex, const ViewParam& view_param) {

    //mat4_t curr_clip_to_prev_clip_mat = view_param.matrix.previous.view_proj * view_param.matrix.inverse.view_proj;
    mat4_t curr_clip_to_prev_clip_mat = view_param.matrix.prev.proj * view_param.matrix.prev.view * view_param.matrix.inv.view * view_param.matrix.inv.proj;

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, depth_tex);

    const vec2_t res = view_param.resolution;
    const vec2_t jitter_uv_cur = view_param.jitter.curr / res;
    const vec2_t jitter_uv_prev = view_param.jitter.prev / res;
    vec4_t jitter_uv = {jitter_uv_cur.x, jitter_uv_cur.y, jitter_uv_prev.x, jitter_uv_prev.y};

    glUseProgram(velocity::blit_velocity.program);
	glUniform1i(velocity::blit_velocity.uniform_loc.tex_depth, 0);
    glUniformMatrix4fv(velocity::blit_velocity.uniform_loc.curr_clip_to_prev_clip_mat, 1, GL_FALSE, &curr_clip_to_prev_clip_mat.elem[0][0]);
    glUniform4fv(velocity::blit_velocity.uniform_loc.jitter_uv, 1, jitter_uv.elem);
    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
}

void blit_tilemax(GLuint velocity_tex, GLuint linear_depth_tex) {
    ASSERT(glIsTexture(velocity_tex));
    ASSERT(glIsTexture(linear_depth_tex));

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, velocity_tex);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, linear_depth_tex);

    glUseProgram(velocity::blit_tilemax.program);
    glUniform1i(velocity::blit_tilemax.uniform_loc.tex_vel, 0);
    glUniform1i(velocity::blit_tilemax.uniform_loc.tex_linear_depth, 1);
    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
}

void blit_neighbormax(GLuint velocity_tex, GLuint linear_depth_tex, int tex_width, int tex_height) {
    ASSERT(glIsTexture(velocity_tex));
    ASSERT(glIsTexture(linear_depth_tex));
    const vec2_t texel_size = {1.f / tex_width, 1.f / tex_height};

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, velocity_tex);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, linear_depth_tex);

    glUseProgram(velocity::blit_neighbormax.program);
    glUniform1i(velocity::blit_neighbormax.uniform_loc.tex_vel, 0);
    glUniform1i(velocity::blit_neighbormax.uniform_loc.tex_linear_depth, 1);
    glUniform2fv(velocity::blit_neighbormax.uniform_loc.tex_vel_texel_size, 1, &texel_size.x);
    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
}

void blit_velocity_dilate(GLuint velocity_tex, int tex_width, int tex_height) {
    ASSERT(glIsTexture(velocity_tex));
    const vec2_t texel_size = {1.f / tex_width, 1.f / tex_height};

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, velocity_tex);

    glUseProgram(velocity::blit_dilate.program);
    glUniform1i(velocity::blit_dilate.uniform_loc.tex_vel, 0);
    glUniform2fv(velocity::blit_dilate.uniform_loc.texel_size, 1, &texel_size.x);
    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
}

void apply_temporal_aa(GLuint linear_depth_tex, GLuint color_tex, GLuint velocity_tex, GLuint velocity_neighbormax_tex, const vec2_t& curr_jitter, const vec2_t& prev_jitter, float feedback_min,
                       float feedback_max, float motion_scale, float time) {
    ASSERT(glIsTexture(linear_depth_tex));
    ASSERT(glIsTexture(color_tex));
    ASSERT(glIsTexture(velocity_tex));
    ASSERT(glIsTexture(velocity_neighbormax_tex));

    static int target = 0;
    target = (target + 1) % 2;

    const int dst_buf = target;
    const int src_buf = (target + 1) % 2;

    const vec2_t res = {(float)gl.tex_width, (float)gl.tex_height};
    const vec2_t inv_res = 1.0f / res;
    const vec4_t texel_size = vec4_t{inv_res.x, inv_res.y, res.x, res.y};
    const vec2_t jitter_uv_curr = curr_jitter / res;
    const vec2_t jitter_uv_prev = prev_jitter / res;
    const vec4_t jitter_uv = vec4_t{jitter_uv_curr.x, jitter_uv_curr.y, jitter_uv_prev.x, jitter_uv_prev.y};

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, linear_depth_tex);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, color_tex);

    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, gl.targets.tex_temporal_buffer[src_buf]);

    glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_2D, velocity_tex);

    glActiveTexture(GL_TEXTURE4);
    glBindTexture(GL_TEXTURE_2D, velocity_neighbormax_tex);

    GLint bound_buffer;
    glGetIntegerv(GL_DRAW_BUFFER, &bound_buffer);

    GLenum draw_buffers[2];
    draw_buffers[0] = GL_COLOR_ATTACHMENT2 + dst_buf;  // tex_temporal_buffer[0 or 1]
    draw_buffers[1] = bound_buffer;                    // assume that this is part of the same gbuffer

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.targets.fbo);
    glViewport(0, 0, gl.tex_width, gl.tex_height);
    glDrawBuffers(2, draw_buffers);

    if (motion_scale != 0.f) {
        glUseProgram(gl.temporal.with_motion_blur.program);

        glUniform1i(gl.temporal.with_motion_blur.uniform_loc.tex_linear_depth, 0);
        glUniform1i(gl.temporal.with_motion_blur.uniform_loc.tex_main, 1);
        glUniform1i(gl.temporal.with_motion_blur.uniform_loc.tex_prev, 2);
        glUniform1i(gl.temporal.with_motion_blur.uniform_loc.tex_vel, 3);
        glUniform1i(gl.temporal.with_motion_blur.uniform_loc.tex_vel_neighbormax, 4);

        glUniform4fv(gl.temporal.with_motion_blur.uniform_loc.texel_size, 1, &texel_size.x);
        glUniform4fv(gl.temporal.with_motion_blur.uniform_loc.jitter_uv, 1, &jitter_uv.x);
        glUniform1f(gl.temporal.with_motion_blur.uniform_loc.time, time);
        glUniform1f(gl.temporal.with_motion_blur.uniform_loc.feedback_min, feedback_min);
        glUniform1f(gl.temporal.with_motion_blur.uniform_loc.feedback_max, feedback_max);
        glUniform1f(gl.temporal.with_motion_blur.uniform_loc.motion_scale, motion_scale);
    } else {
        glUseProgram(gl.temporal.no_motion_blur.program);

        glUniform1i(gl.temporal.no_motion_blur.uniform_loc.tex_linear_depth, 0);
        glUniform1i(gl.temporal.no_motion_blur.uniform_loc.tex_main, 1);
        glUniform1i(gl.temporal.no_motion_blur.uniform_loc.tex_prev, 2);
        glUniform1i(gl.temporal.no_motion_blur.uniform_loc.tex_vel, 3);

        glUniform4fv(gl.temporal.no_motion_blur.uniform_loc.texel_size, 1, &texel_size.x);
        glUniform4fv(gl.temporal.no_motion_blur.uniform_loc.jitter_uv, 1, &jitter_uv.x);
        glUniform1f(gl.temporal.no_motion_blur.uniform_loc.time, time);
        glUniform1f(gl.temporal.no_motion_blur.uniform_loc.feedback_min, feedback_min);
        glUniform1f(gl.temporal.no_motion_blur.uniform_loc.feedback_max, feedback_max);
        glUniform1f(gl.temporal.no_motion_blur.uniform_loc.motion_scale, motion_scale);
    }

    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
}

void scale_hsv(GLuint color_tex, vec3_t hsv_scale) {
    GLint last_fbo;
    GLint last_viewport[4];
    GLint last_scissor_box[4];
    GLint last_draw_buffer;
    glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &last_fbo);
    glGetIntegerv(GL_VIEWPORT, last_viewport);
    glGetIntegerv(GL_SCISSOR_BOX, last_scissor_box);
    glGetIntegerv(GL_DRAW_BUFFER, &last_draw_buffer);

    GLint w, h;

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, color_tex);

    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH,  &w);
    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &h);

    glBindVertexArray(gl.vao);

    glViewport(0, 0, w, h);
    glScissor(0, 0, w, h);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.tmp.fbo);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.tmp.tex_rgba8, 0);
    glDrawBuffer(GL_COLOR_ATTACHMENT0);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, color_tex);

    glUseProgram(hsv::gl.program);
    glUniform1i(hsv::gl.uniform_loc.texture_color, 0);
    glUniform3fv(hsv::gl.uniform_loc.hsv_scale, 1, &hsv_scale.x);
    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);

    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, color_tex, 0);
    glDrawBuffer(GL_COLOR_ATTACHMENT1);
    blit_texture(gl.tmp.tex_rgba8);

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, last_fbo);
    glViewport(last_viewport[0], last_viewport[1], last_viewport[2], last_viewport[3]);
    glScissor(last_scissor_box[0], last_scissor_box[1], last_scissor_box[2], last_scissor_box[3]);
    glDrawBuffer(last_draw_buffer);
}

void blit_texture(GLuint tex) {
    ASSERT(glIsTexture(tex));
    glUseProgram(blit::program_tex);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, tex);
    glUniform1i(blit::uniform_loc_texture, 0);
    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
}

void blit_color(vec4_t color) {
    glUseProgram(blit::program_col);
    glUniform4fv(blit::uniform_loc_color, 1, &color.x);
    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
}

void blur_texture_gaussian(GLuint tex, int num_passes) {
    ASSERT(glIsTexture(tex));
    ASSERT(num_passes > 0);

    GLint last_fbo;
    GLint last_viewport[4];
    GLint last_draw_buffer[8];
    glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &last_fbo);
    glGetIntegerv(GL_VIEWPORT, last_viewport);
    for (int i = 0; i < 8; ++i) glGetIntegerv(GL_DRAW_BUFFER0 + i, &last_draw_buffer[i]);

    GLint w, h;

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, tex);

    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, &w);
    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &h);

    glBindVertexArray(gl.vao);

    glUseProgram(blur::program_gaussian);
    glUniform1i(blur::uniform_loc_texture, 0);

    glViewport(0, 0, w, h);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.tmp.fbo);

    for (int i = 0; i < num_passes; ++i) {
        glBindTexture(GL_TEXTURE_2D, tex);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.tmp.tex_rgba8, 0);
        glDrawBuffer(GL_COLOR_ATTACHMENT0);
        glUniform2f(blur::uniform_loc_inv_res_dir, 1.0f / w, 0.0f);
        glDrawArrays(GL_TRIANGLES, 0, 3);

        glBindTexture(GL_TEXTURE_2D, gl.tmp.tex_rgba8);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex, 0);
        glDrawBuffer(GL_COLOR_ATTACHMENT0);
        glUniform2f(blur::uniform_loc_inv_res_dir, 0.0f, 1.0f / h);
        glDrawArrays(GL_TRIANGLES, 0, 3);
    }

    glUseProgram(0);
    glBindVertexArray(0);

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, last_fbo);
    glViewport(last_viewport[0], last_viewport[1], last_viewport[2], last_viewport[3]);
    for (int i = 0; i < 8; ++i) glDrawBuffers(8, (GLenum*)last_draw_buffer);
}

void blur_texture_box(GLuint tex, int num_passes) {
    ASSERT(glIsTexture(tex));
    ASSERT(num_passes > 0);

    GLint last_fbo;
    GLint last_viewport[4];
    GLint last_draw_buffer[8];
    glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &last_fbo);
    glGetIntegerv(GL_VIEWPORT, last_viewport);
    for (int i = 0; i < 8; ++i) glGetIntegerv(GL_DRAW_BUFFER0 + i, &last_draw_buffer[i]);

    GLint w, h;

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, tex);

    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, &w);
    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &h);

    glBindVertexArray(gl.vao);

    glUseProgram(blur::program_box);

    glViewport(0, 0, w, h);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.tmp.fbo);

    for (int i = 0; i < num_passes; ++i) {
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.tmp.tex_rgba8, 0);
        glDrawBuffer(GL_COLOR_ATTACHMENT0);
        glDrawArrays(GL_TRIANGLES, 0, 3);

        glBindTexture(GL_TEXTURE_2D, gl.tmp.tex_rgba8);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex, 0);
        glDrawBuffer(GL_COLOR_ATTACHMENT0);
        glDrawArrays(GL_TRIANGLES, 0, 3);
    }

    glUseProgram(0);
    glBindVertexArray(0);

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, last_fbo);
    glViewport(last_viewport[0], last_viewport[1], last_viewport[2], last_viewport[3]);
    for (int i = 0; i < 8; ++i) glDrawBuffers(8, (GLenum*)last_draw_buffer);
}

static void compute_luma(GLuint tex) {
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, tex);

    glBindVertexArray(gl.vao);
    glUseProgram(gl.luma.program);
    glUniform1i(gl.luma.uniform_loc.tex_rgba, 0);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glUseProgram(0);
    glBindVertexArray(0);
}

static void compute_fxaa(GLuint tex_rgbl, int width, int height) {
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, tex_rgbl);

    int w, h;
    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, &w);
    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &h);

    float scl_x = (float)width  / (float)w;
    float scl_y = (float)height / (float)h;

    glBindVertexArray(gl.vao);
    glUseProgram(gl.fxaa.program);
    glUniform2f(gl.fxaa.uniform_loc.tc_scl, scl_x, scl_y);
    glUniform1i(gl.fxaa.uniform_loc.tex_rgbl, 0);
    glUniform2f(gl.fxaa.uniform_loc.rcp_res, 1.0f / (float)w, 1.0f / (float)h);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glUseProgram(0);
    glBindVertexArray(0);
}

void shade_and_postprocess(const Descriptor& desc, const ViewParam& view_param) {
    ASSERT(glIsTexture(desc.input_textures.depth));
    ASSERT(glIsTexture(desc.input_textures.color));
    ASSERT(glIsTexture(desc.input_textures.normal));
    if (desc.temporal_aa.enabled) {
        ASSERT(glIsTexture(desc.input_textures.velocity));
    }

    // For seeding noise
    static float time = 0.f;
    time = time + 0.01f;
    if (time > 100.f) time -= 100.f;
    //time = 0.0f;
    //static unsigned int frame = 0;
    //frame = frame + 1;

    const auto near_dist = view_param.clip_planes.near;
    const auto far_dist = view_param.clip_planes.far;
    const auto ortho = is_orthographic_proj_matrix(view_param.matrix.curr.proj.elem);

    GLint last_fbo;
    GLint last_viewport[4];
    GLint last_draw_buffer;
    GLint last_scissor_box[4];
    glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &last_fbo);
    glGetIntegerv(GL_VIEWPORT, last_viewport);
    glGetIntegerv(GL_SCISSOR_BOX, last_scissor_box);
    glGetIntegerv(GL_DRAW_BUFFER, &last_draw_buffer);

    int width = last_viewport[2];
    int height = last_viewport[3];

    if (width > (int)gl.tex_width || height > (int)gl.tex_height) {
        initialize(width, height);
    }

    //glViewport(0, 0, gl.tex_width, gl.tex_height);
    glEnable(GL_SCISSOR_TEST);
    glScissor(0, 0, width, height);
    glViewport(0, 0, width, height);
    glBindVertexArray(gl.vao);

    PUSH_GPU_SECTION("Linearize Depth") {
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.linear_depth.fbo_0);
        glClearColor(far_dist,0,0,0);
        glClear(GL_COLOR_BUFFER_BIT);
        compute_linear_depth(desc.input_textures.depth, near_dist, far_dist, ortho);
    }
    POP_GPU_SECTION()

    if (desc.ambient_occlusion.enabled) {
        PUSH_GPU_SECTION("Generate Linear Depth Mipmaps") {
            glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.linear_depth.fbo_1);
            glViewport(0, 0, gl.tex_width / 2, gl.tex_height / 2);
            downsample_depth(gl.linear_depth.texture, 0);
            glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.linear_depth.fbo_2);
            glViewport(0, 0, gl.tex_width / 4, gl.tex_height / 4);
            downsample_depth(gl.linear_depth.texture, 1);
            glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.linear_depth.fbo_3);
            glViewport(0, 0, gl.tex_width / 8, gl.tex_height / 8);
            downsample_depth(gl.linear_depth.texture, 2);
        }
        POP_GPU_SECTION()
    }

    if (desc.temporal_aa.enabled) {
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.velocity.fbo);
        glViewport(0, 0, gl.velocity.tex_width, gl.velocity.tex_height);
        glScissor(0, 0, gl.velocity.tex_width, gl.velocity.tex_height);
        glClearColor(0,0,0,0);
        glClear(GL_COLOR_BUFFER_BIT);

        PUSH_GPU_SECTION("Velocity: Tilemax") {
            glDrawBuffer(GL_COLOR_ATTACHMENT0);
            blit_tilemax(desc.input_textures.velocity, gl.linear_depth.texture);
        }
        POP_GPU_SECTION()

        PUSH_GPU_SECTION("Velocity: Neighbormax") {
            glDrawBuffer(GL_COLOR_ATTACHMENT1);
            blit_neighbormax(gl.velocity.tex_tilemax, gl.linear_depth.texture, gl.velocity.tex_width, gl.velocity.tex_height);
        }
        POP_GPU_SECTION()
    }

    const GLenum draw_buffers[2] = {
        GL_COLOR_ATTACHMENT0,
        GL_COLOR_ATTACHMENT1,
    };
     
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.targets.fbo);
    glViewport(0, 0, width, height);
    glScissor(0, 0, width, height);
    glDrawBuffers(2, draw_buffers);
    glClearColor(0,0,0,0);
    glClear(GL_COLOR_BUFFER_BIT);

    GLenum dst_buffer = GL_COLOR_ATTACHMENT1;
    GLuint src_texture = gl.targets.tex_color[0];

    auto swap_target = [&dst_buffer, &src_texture]() {
        dst_buffer = dst_buffer == GL_COLOR_ATTACHMENT0 ? GL_COLOR_ATTACHMENT1 : GL_COLOR_ATTACHMENT0;
        src_texture = src_texture == gl.targets.tex_color[0] ? gl.targets.tex_color[1] : gl.targets.tex_color[0];
    };
    glDrawBuffer(dst_buffer);

    PUSH_GPU_SECTION("Compose")
    compose_deferred(desc.input_textures.depth, desc.input_textures.color, desc.input_textures.normal, view_param.matrix.inv.proj, desc.background.color, time);
    POP_GPU_SECTION()

    if (desc.ambient_occlusion.enabled) {
        PUSH_GPU_SECTION("SSAO")
        compute_ssao(gl.linear_depth.texture, desc.input_textures.normal, view_param.matrix.curr.proj.elem, desc.ambient_occlusion.intensity, desc.ambient_occlusion.radius, desc.ambient_occlusion.bias);
        POP_GPU_SECTION()
    }

    PUSH_GPU_SECTION("Tonemapping")
    swap_target();
    glDrawBuffer(dst_buffer);
    Tonemapping tonemapper = desc.tonemapping.enabled ? desc.tonemapping.mode : Tonemapping_Passthrough;
    apply_tonemapping(src_texture, tonemapper, desc.tonemapping.exposure, desc.tonemapping.gamma);
    POP_GPU_SECTION()

    if (desc.depth_of_field.enabled) {
        swap_target();
        glDrawBuffer(dst_buffer);
        PUSH_GPU_SECTION("DOF")
        apply_dof(gl.linear_depth.texture, src_texture, desc.depth_of_field.focus_depth, desc.depth_of_field.focus_scale, time);
        POP_GPU_SECTION()
    }

    if (desc.input_textures.transparency) {
        PUSH_GPU_SECTION("Add Transparency")
            glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        blit_texture(desc.input_textures.transparency);
        glDisable(GL_BLEND);
        POP_GPU_SECTION()
    }

    if (desc.fxaa.enabled) {
        swap_target();
        PUSH_GPU_SECTION("Luma")
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.tmp.fbo);
        glDrawBuffer(GL_COLOR_ATTACHMENT0);
        compute_luma(src_texture);
        POP_GPU_SECTION()
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.targets.fbo);
        glDrawBuffer(dst_buffer);
        PUSH_GPU_SECTION("FXAA")
        compute_fxaa(gl.tmp.tex_rgba8, width, height);
        POP_GPU_SECTION()
    }

    if (desc.temporal_aa.enabled) {
        swap_target();
        glDrawBuffer(dst_buffer);
        const float feedback_min = desc.temporal_aa.feedback_min;
        const float feedback_max = desc.temporal_aa.feedback_max;
        const float motion_scale = desc.temporal_aa.motion_blur.enabled ? desc.temporal_aa.motion_blur.motion_scale : 0.f;
        if (motion_scale != 0.f)
            PUSH_GPU_SECTION("Temporal AA + Motion Blur")
        else
            PUSH_GPU_SECTION("Temporal AA")

        apply_temporal_aa(gl.linear_depth.texture, src_texture, desc.input_textures.velocity, gl.velocity.tex_neighbormax, view_param.jitter.curr, view_param.jitter.prev, feedback_min, feedback_max, motion_scale, time);
        POP_GPU_SECTION()
    }
     
    if (desc.sharpen.enabled) {
        PUSH_GPU_SECTION("Sharpen")
        swap_target();
        glDrawBuffer(dst_buffer);
        sharpen::sharpen(src_texture, desc.sharpen.weight);
        POP_GPU_SECTION()
    }

    // Activate backbuffer or whatever was bound before
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, last_fbo);
    glDrawBuffer(last_draw_buffer);
    glViewport(last_viewport[0], last_viewport[1], last_viewport[2], last_viewport[3]);
    glScissor(last_scissor_box[0], last_scissor_box[1], last_scissor_box[2], last_scissor_box[3]);
    glDisable(GL_SCISSOR_TEST);

    swap_target();
    glDepthMask(0);
    blit_texture(src_texture);

    glDepthMask(1);
    glColorMask(1, 1, 1, 1);
}

}  // namespace postprocessing

// #gbuffer
void init_gbuffer(GBuffer* gbuf, int width, int height) {
    ASSERT(gbuf);

    bool attach_textures_deferred = false;
    if (!gbuf->fbo) {
        glGenFramebuffers(1, &gbuf->fbo);
        attach_textures_deferred = true;
    }

    if (!gbuf->tex.depth) glGenTextures(1, &gbuf->tex.depth);
    if (!gbuf->tex.color) glGenTextures(1, &gbuf->tex.color);
    if (!gbuf->tex.normal) glGenTextures(1, &gbuf->tex.normal);
    if (!gbuf->tex.velocity) glGenTextures(1, &gbuf->tex.velocity);
    if (!gbuf->tex.transparency) glGenTextures(1, &gbuf->tex.transparency);
    if (!gbuf->tex.picking) glGenTextures(1, &gbuf->tex.picking);
    if (!gbuf->pbo_picking.color[0]) glGenBuffers((int)ARRAY_SIZE(gbuf->pbo_picking.color), gbuf->pbo_picking.color);
    if (!gbuf->pbo_picking.depth[0]) glGenBuffers((int)ARRAY_SIZE(gbuf->pbo_picking.depth), gbuf->pbo_picking.depth);

    glBindTexture(GL_TEXTURE_2D, gbuf->tex.depth);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH24_STENCIL8, width, height, 0, GL_DEPTH_COMPONENT, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindTexture(GL_TEXTURE_2D, gbuf->tex.color);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindTexture(GL_TEXTURE_2D, gbuf->tex.normal);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RG16, width, height, 0, GL_RG, GL_UNSIGNED_SHORT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindTexture(GL_TEXTURE_2D, gbuf->tex.velocity);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RG16F, width, height, 0, GL_RG, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindTexture(GL_TEXTURE_2D, gbuf->tex.transparency);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindTexture(GL_TEXTURE_2D, gbuf->tex.picking);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    /*
    glBindTexture(GL_TEXTURE_2D, gbuf->tex.temporal_accumulation[0]);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_BGRA, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindTexture(GL_TEXTURE_2D, gbuf->tex.temporal_accumulation[1]);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_BGRA, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    */

    for (uint32_t i = 0; i < ARRAY_SIZE(gbuf->pbo_picking.color); ++i) {
        glBindBuffer(GL_PIXEL_PACK_BUFFER, gbuf->pbo_picking.color[i]);
        glBufferData(GL_PIXEL_PACK_BUFFER, 4, NULL, GL_DYNAMIC_READ);
        glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
    }

    for (uint32_t i = 0; i < ARRAY_SIZE(gbuf->pbo_picking.depth); ++i) {
        glBindBuffer(GL_PIXEL_PACK_BUFFER, gbuf->pbo_picking.depth[i]);
        glBufferData(GL_PIXEL_PACK_BUFFER, 4, NULL, GL_DYNAMIC_READ);
        glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
    }

    glBindTexture(GL_TEXTURE_2D, 0);

    gbuf->width = width;
    gbuf->height = height;

    const GLenum draw_buffers[] = {GL_COLOR_ATTACHMENT_COLOR, GL_COLOR_ATTACHMENT_NORMAL, GL_COLOR_ATTACHMENT_VELOCITY,
        GL_COLOR_ATTACHMENT_TRANSPARENCY, GL_COLOR_ATTACHMENT_PICKING};

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gbuf->fbo);
    if (attach_textures_deferred) {
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, gbuf->tex.depth, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_STENCIL_ATTACHMENT, GL_TEXTURE_2D, gbuf->tex.depth, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT_COLOR, GL_TEXTURE_2D, gbuf->tex.color, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT_NORMAL, GL_TEXTURE_2D, gbuf->tex.normal, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT_VELOCITY, GL_TEXTURE_2D, gbuf->tex.velocity, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT_TRANSPARENCY, GL_TEXTURE_2D, gbuf->tex.transparency, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT_PICKING, GL_TEXTURE_2D, gbuf->tex.picking, 0);
    }
    ASSERT(glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE);
    glDrawBuffers((int)ARRAY_SIZE(draw_buffers), draw_buffers);
    glClearColor(0, 0, 0, 0);
    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
}

void clear_gbuffer(GBuffer* gbuffer) {
    const GLenum draw_buffers[] = {GL_COLOR_ATTACHMENT_COLOR, GL_COLOR_ATTACHMENT_NORMAL, GL_COLOR_ATTACHMENT_VELOCITY, GL_COLOR_ATTACHMENT_PICKING, GL_COLOR_ATTACHMENT_TRANSPARENCY};

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gbuffer->fbo);
    glViewport(0, 0, gbuffer->width, gbuffer->height);

    glDepthMask(1);
    glColorMask(1, 1, 1, 1);
    glStencilMask(0xFF);
    const vec4_t zero    = {0,0,0,0};
    const vec4_t picking = {1,1,1,1};

    // Setup gbuffer and clear textures
    PUSH_GPU_SECTION("Clear G-buffer") {
        // Clear color+alpha, normal, velocity, emissive, post_tonemap and depth+stencil
        glDrawBuffers((int)ARRAY_SIZE(draw_buffers), draw_buffers);
        glClearBufferfv(GL_COLOR, 0, zero.elem);
        glClearBufferfv(GL_COLOR, 1, zero.elem);
        glClearBufferfv(GL_COLOR, 2, zero.elem);
        glClearBufferfv(GL_COLOR, 3, picking.elem);
        glClearBufferfv(GL_COLOR, 4, zero.elem);
        glClearBufferfi(GL_DEPTH_STENCIL, 0, 1.0f, 0x01);
    }
    POP_GPU_SECTION()
}

void destroy_gbuffer(GBuffer* gbuf) {
    ASSERT(gbuf);
    if (gbuf->fbo) glDeleteFramebuffers(1, &gbuf->fbo);
    if (gbuf->tex.depth) glDeleteTextures(1, &gbuf->tex.depth);
    if (gbuf->tex.color) glDeleteTextures(1, &gbuf->tex.color);
    if (gbuf->tex.normal) glDeleteTextures(1, &gbuf->tex.normal);
    if (gbuf->tex.transparency) glDeleteTextures(1, &gbuf->tex.transparency);
    if (gbuf->tex.picking) glDeleteTextures(1, &gbuf->tex.picking);
    if (gbuf->tex.temporal_accumulation) glDeleteTextures((int)ARRAY_SIZE(gbuf->tex.temporal_accumulation), gbuf->tex.temporal_accumulation);

    if (gbuf->pbo_picking.color[0]) glDeleteBuffers((int)ARRAY_SIZE(gbuf->pbo_picking.color), gbuf->pbo_picking.color);
    if (gbuf->pbo_picking.depth[0]) glDeleteBuffers((int)ARRAY_SIZE(gbuf->pbo_picking.depth), gbuf->pbo_picking.depth);
}

// #picking
void extract_gbuffer_picking_idx_and_depth(uint32_t* out_idx, float* out_depth, GBuffer* gbuf, int x, int y) {
    uint32_t idx = 0;
    float depth = 0;

#if EXPERIMENTAL_GFX_API
    if (use_gfx) {
        idx = md_gfx_get_picking_idx();
        depth = md_gfx_get_picking_depth();
        md_gfx_query_picking((uint32_t)x, (uint32_t)y);
    }
    else {
#endif
        ASSERT(gbuf);
        const uint32_t N = (uint32_t)ARRAY_SIZE(gbuf->pbo_picking.color);
        uint32_t frame = gbuf->pbo_picking.frame++;
        uint32_t queue = (frame) % N;
        uint32_t read  = (frame + N-1) % N;

        uint8_t  color[4];

        PUSH_GPU_SECTION("READ PICKING DATA")
        glBindFramebuffer(GL_READ_FRAMEBUFFER, gbuf->fbo);
        glReadBuffer(GL_COLOR_ATTACHMENT_PICKING);

        // Queue async reads from current frame to pixel pack buffer
        glBindBuffer(GL_PIXEL_PACK_BUFFER, gbuf->pbo_picking.color[queue]);
        glReadPixels(x, y, 1, 1, GL_BGRA, GL_UNSIGNED_BYTE, 0);

        glBindBuffer(GL_PIXEL_PACK_BUFFER, gbuf->pbo_picking.depth[queue]);
        glReadPixels(x, y, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, 0);

        // Read values from previous frames pixel pack buffer
        glBindBuffer(GL_PIXEL_PACK_BUFFER, gbuf->pbo_picking.color[read]);
        glGetBufferSubData(GL_PIXEL_PACK_BUFFER, 0, sizeof(color), color);

        glBindBuffer(GL_PIXEL_PACK_BUFFER, gbuf->pbo_picking.depth[read]);
        glGetBufferSubData(GL_PIXEL_PACK_BUFFER, 0, sizeof(depth), &depth);

        glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
        glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);
        POP_GPU_SECTION()

        // BGRA
        idx = (color[0] << 16) | (color[1] << 8) | (color[2] << 0) | (color[3] << 24);

#if EXPERIMENTAL_GFX_API
    }
#endif

    if (out_idx) {
        *out_idx = idx;
    }
    if (out_depth) {
        *out_depth = depth;
    }
}
