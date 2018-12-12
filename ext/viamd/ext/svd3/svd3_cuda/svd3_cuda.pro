TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp
HEADERS  += svd3_cuda.h

CUDA_SOURCES += \
    svd3.cu

OTHER_FILES += \
    svd3.cu

CUDA_DIR = /contrib/projects/cuda5-toolkit
INCLUDEPATH += $$CUDA_DIR/include
INCLUDEPATH += $$CUDA_DIR/samples/common/inc
QMAKE_LIBDIR += $$CUDA_DIR/lib64

LIBS += -lcudart -lcuda

# GPU ARCH
# this gets passed as the gpu-architecture flag to nvcc compiler
# specifying particular architectures enable certain features, limited to the compute capability
# of the GPU. compute capabilities listed here http://en.wikipedia.org/wiki/CUDA
# MSLAB GeForce 460 seems to have compute capability 2.1
CUDA_ARCH = sm_21
#CUDA_ARCH = sm_35

# custom NVCC flags
NVCCFLAGS = --compiler-options -fno-strict-aliasing -use_fast_math --ptxas-options=-v

# Prepare the extra compiler configuration (taken from the nvidia forum - i'm not an expert in this part)
CUDA_INC = $$join(INCLUDEPATH,' -I','-I',' ') -I$$_PRO_FILE_PWD_

# compile CUDA kernels using nvcc
cuda.commands = $$CUDA_DIR/bin/nvcc -m64 -g -G -arch=$$CUDA_ARCH -c $$NVCCFLAGS $$CUDA_INC $$LIBS  ${QMAKE_FILE_NAME} -o ${QMAKE_FILE_OUT} \
    2>&1 | sed -r \"s/\\(([0-9]+)\\)/:\\1/g\" 1>&2
# Prepare the extra compiler configuration (taken from the nvidia forum - i'm not an expert in this part)
cuda.input = CUDA_SOURCES
cuda.output = ${OBJECTS_DIR}${QMAKE_FILE_BASE}_cuda.o # suffix needed for this to work?
# Tell Qt that we want add more stuff to the Makefile
QMAKE_EXTRA_COMPILERS += cuda
