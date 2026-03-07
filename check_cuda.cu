#include <cuda_runtime.h>
#include <stdio.h>

int main() {
    int count = 0;
    cudaError_t err = cudaGetDeviceCount(&count);
    if (err != cudaSuccess) {
        printf("CUDA Error Code: %d (%s)\n", (int)err, cudaGetErrorString(err));
        return 1;
    }
    printf("Device count: %d\n", count);
    return 0;
}
