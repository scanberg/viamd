#include <iostream>
#include <stdio.h>

extern "C" void runCudaPart();

int main()
{
    runCudaPart();
    return 0;
}

