// Minimal CUDA DCT round-trip test
// Verifies that DCT(IDCT(x)) = x and finds the correct normalization.
#include <cufft.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include <math.h>

#define N 8  // 1D test size

#define CHECK_CUDA(call) \
    { \
        cudaError_t err = call; \
        if (err != cudaSuccess) { \
            printf("CUDA Error: %s at line %d\n", cudaGetErrorString(err), __LINE__); \
            return 1; \
        } \
    }

#define CHECK_CUFFT(call) \
    { \
        cufftResult res = call; \
        if (res != CUFFT_SUCCESS) { \
            printf("cuFFT Error: %d at line %d\n", (int)res, __LINE__); \
            return 1; \
        } \
    }

int main()
{
    int w = N;
    int w2 = 2 * w;
    const double pi = 3.14159265358979323846;

    double h_in[N] = {1, 2, 3, 4, 5, 6, 7, 8};
    double h_ext[2*N];
    for (int i = 0; i < w; i++) h_ext[i] = h_in[i];
    for (int i = 0; i < w; i++) h_ext[2*w-1-i] = h_in[i];

    double *d_ext;
    cufftDoubleComplex *d_complex;
    CHECK_CUDA(cudaMalloc(&d_ext, w2 * sizeof(double)));
    CHECK_CUDA(cudaMalloc(&d_complex, (w + 1) * sizeof(cufftDoubleComplex)));
    CHECK_CUDA(cudaMemcpy(d_ext, h_ext, w2 * sizeof(double), cudaMemcpyHostToDevice));

    cufftHandle plan;
    CHECK_CUFFT(cufftPlan1d(&plan, w2, CUFFT_D2Z, 1));
    CHECK_CUFFT(cufftExecD2Z(plan, d_ext, d_complex));
    CHECK_CUDA(cudaDeviceSynchronize());

    cufftDoubleComplex h_complex[N+1];
    CHECK_CUDA(cudaMemcpy(h_complex, d_complex, (w + 1) * sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost));

    printf("Raw FFT Output:\n");
    for (int k = 0; k <= w; k++) {
        printf("  fft[%d] = %f + %fj\n", k, h_complex[k].x, h_complex[k].y);
    }

    double h_modes[N];
    double norm = 2.0 * sqrt(2.0 * w);
    printf("\nDCT-II modes (norm=%.4f):\n", norm);
    for (int k = 0; k < w; k++) {
        double theta = -pi * k / (2.0 * w);
        double raw = h_complex[k].x * cos(theta) - h_complex[k].y * sin(theta);
        h_modes[k] = raw / norm;
        printf("  mode[%d] = %f\n", k, h_modes[k]);
    }

    cufftDoubleComplex h_inv_in[N+1] = {};
    for (int k = 0; k < w; k++) {
        double theta = +pi * k / (2.0 * w);
        h_inv_in[k].x = h_modes[k] * cos(theta);
        h_inv_in[k].y = h_modes[k] * sin(theta);
    }

    cufftDoubleComplex *d_inv_in;
    double *d_inv_out;
    CHECK_CUDA(cudaMalloc(&d_inv_in, (w + 1) * sizeof(cufftDoubleComplex)));
    CHECK_CUDA(cudaMalloc(&d_inv_out, w2 * sizeof(double)));
    CHECK_CUDA(cudaMemcpy(d_inv_in, h_inv_in, (w + 1) * sizeof(cufftDoubleComplex), cudaMemcpyHostToDevice));

    cufftHandle plan2;
    CHECK_CUFFT(cufftPlan1d(&plan2, w2, CUFFT_Z2D, 1));
    CHECK_CUFFT(cufftExecZ2D(plan2, d_inv_in, d_inv_out));
    CHECK_CUDA(cudaDeviceSynchronize());

    double h_inv_out[2*N];
    CHECK_CUDA(cudaMemcpy(h_inv_out, d_inv_out, w2 * sizeof(double), cudaMemcpyDeviceToHost));

    printf("\nIDCT output (G[k] = mode[k]*exp(+j*pi*k/2N), C2R output):\n");
    for (int i = 0; i < w; i++) {
        printf("  reconstructed[%i] = %f   (original = %.1f)\n", i, h_inv_out[i], h_in[i]);
    }
    printf("\nScale factor should be: %f\n", h_inv_out[0] / h_in[0]);

    cufftDestroy(plan);
    cufftDestroy(plan2);
    cudaFree(d_ext);
    cudaFree(d_complex);
    cudaFree(d_inv_in);
    cudaFree(d_inv_out);
    return 0;
}
