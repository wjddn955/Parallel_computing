/**
 * Copyright 1993-2015 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 */

#include <stdio.h>
#include <assert.h>
#include <omp.h>
#include <sys/time.h>

__global__ void
matrixMulCUDAkernel(float *C, float *A, float *B, int n)
{
    float Cvalue = 0;
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    for ( int e = 0; e < n; ++e)
        Cvalue += A[row * n + e] * B[e * n + col];
    C[row * n + col] = Cvalue;
}

void matrixMulCPU_serial(float *C, float *A, float *B, int n)
{
    for (int i = 0; i < (int)(n*n); i++)
    {
	
        float temp = 0;
        int x = i % n;
        int y = i / n;
        for( int k = 0; k < (int)(n); k++)
        {
            temp += A[n * y + k] * B[n * k + x];
        }

        C[n*y + x] = temp;
    }
}

void matrixMulCPU_parallel(float *C, float *A, float *B, int n)
{
    omp_set_num_threads(16);
    #pragma omp parallel for
    for (int i = 0; i < (int)(n*n); i++)
    {
	
        float temp = 0;
        int x = i % n;
        int y = i / n;
        for( int k = 0; k < (int)(n); k++)
        {
            temp += A[n * y + k] * B[n * k + x];
        }

        C[n*y + x] = temp;
    }
}

void matrixMulCUDA( float *h_C, float *h_A, float *h_B, int block_size, int n)
{
    unsigned int mem_size = sizeof(float) * n * n;
    // Allocate device memory
    float *d_A, *d_B, *d_C;

    cudaMalloc((void **) &d_A, mem_size);
    cudaMalloc((void **) &d_B, mem_size);
    cudaMalloc((void **) &d_C, mem_size);

    // copy host memory to device
    cudaMemcpy(d_A, h_A, mem_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, h_B, mem_size, cudaMemcpyHostToDevice);

    // Setup execution parameters
    dim3 threads(block_size, block_size); // (block_size, block_size, 1) default setting 1
    dim3 grid(n / threads.x, n / threads.y);

    // Execute the kernel
    matrixMulCUDAkernel<<< grid, threads >>>(d_C, d_A, d_B, n);

    cudaDeviceSynchronize();

    // Copy result from device to host
    cudaMemcpy(h_C, d_C, mem_size, cudaMemcpyDeviceToHost);

    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);
}

__global__ void
matrixMulTiledCUDAkernel(float *C, float *A, float *B, int tile_size, int n) 
{
    // Shared memory for tiles
    __shared__ float ds_A[tile_size][tile_size];
    __shared__ float ds_B[tile_size][tile_size];

    int row = blockIdx.y * tile_size + threadIdx.y;
    int col = blockIdx.x * tile_size + threadIdx.x;

    float Cvalue = 0;

    // Need to proceed loading data for tile n/tilie_size times
    for (int e = 0; e < n/tile_size; ++e) {
        // Load from A & B to tile
        ds_A[threadIdx.y][threadIdx.x] = A[row*n + (e * tile_size) + threadIdx.x];
        ds_B[threadIdx.y][threadIdx.x] = B[(e * tile_size + threadIdx.y) * n + col];
        __syncthreads();

        // Partial computation
        for (int i = 0; i < tile_size; ++i) {
            Cvalue += ds_A[threadIdx.y][i] * ds_B[i][threadIdx.x];
        }
        __syncthreads();
    }
    C[row * n + col] = Cvalue;

}

void matrixMulCUDA_tiled( float *h_C, float *h_A, float *h_B, int block_size, int n)
{
    // TODO: implement the code for tiled matrix multiplication

    unsigned int mem_size = sizeof(float) * n * n;
    // Allocate device memory
    float *d_A, *d_B, *d_C;

    cudaMalloc((void **) &d_A, mem_size);
    cudaMalloc((void **) &d_B, mem_size);
    cudaMalloc((void **) &d_C, mem_size);

    // copy host memory to device
    cudaMemcpy(d_A, h_A, mem_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, h_B, mem_size, cudaMemcpyHostToDevice);

    // Setup execution parameters
    dim3 threads(block_size, block_size); // (block_size, block_size, 1) default setting 1
    dim3 grid(n / threads.x, n / threads.y);

    matrixMulTiledCUDAkernel<<<grid, threads>>>(d_C, d_A, d_B, block_size, n)

    cudaDeviceSynchronize();

    // Copy result from device to host
    cudaMemcpy(h_C, d_C, mem_size, cudaMemcpyDeviceToHost);

    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);
}

bool error_checking(float *reference, float *result, int n, double eps)
{
    printf("Checking for error with eps=%.1e\n",eps);

    bool correct = true;
    for (int i = 0; i < (int)(n * n); i++)
    {
        double abs_err = fabs(reference[i] - result[i]);
        double abs_val = fabs(result[i]);
        double rel_err = abs_err/abs_val;

        if(rel_err > eps){
            // Remove comment if you want detailed error report
            //int x = i % n;
            //int y = i / n;
            //printf("Error! Matrix[%d][%d]=%.8f, ref=%.8f\n", y,x, result[i], reference[i]);
	    correct = false;
        }
    }

    printf("%s\n", correct ? "Result = PASS" : "Result = FAIL");

    return correct;
}

void initMatrix(float *data, int size)
{
    for (int i = 0; i < size; ++i)
    {
        data[i] = ((float)random())/10000;
    }
}

void initMatrix0(float *data, int size)
{
    for (int i = 0; i < size; ++i)
    {
        data[i] = 0;
    }
}

float get_elapsed_time(struct timeval start_time, struct timeval stop_time)
{
    float elapsed_time = (stop_time.tv_sec - start_time.tv_sec) +
                (float)(stop_time.tv_usec - start_time.tv_usec)/1000000;

    return elapsed_time;
}

int matrixMultiply(int block_size, int n)
{
    // Allocate memory for matrices A and B
    unsigned int size = n * n;
    unsigned int mem_size = sizeof(float) * size;
    float *matrix_A = (float *)malloc(mem_size);
    float *matrix_B = (float *)malloc(mem_size);

    bool correctness = true;

    // Initialize matrix A, B
    srandom(0);
    initMatrix(matrix_A, size);
    initMatrix(matrix_B, size);

    // Allocate matrix C
    float *matrix_C = (float *) malloc(mem_size);

    float *reference = (float *) malloc(mem_size);
    initMatrix0(reference, size);

    struct timeval start_time;
    struct timeval stop_time;
    
    
    printf("\nStarting CPU serial\n");
    gettimeofday(&start_time, NULL);

    matrixMulCPU_serial(reference, matrix_A, matrix_B, n);

    gettimeofday(&stop_time, NULL);
    printf("Ended CPU serial\n");

    float elapsed_time = get_elapsed_time(start_time, stop_time);
    printf("CPU serial execution time:%f\n",elapsed_time);

    printf("Use CPU serial result as error checking reference\n");



    initMatrix0(matrix_C, size);

    printf("\nStarting CPU parallel\n");
    gettimeofday(&start_time, NULL);

    matrixMulCPU_parallel(matrix_C, matrix_A, matrix_B, n);

    gettimeofday(&stop_time, NULL);
    printf("Ended CPU parallel\n");

    elapsed_time = get_elapsed_time(start_time, stop_time);
    printf("CPU parallel execution time:%f\n", elapsed_time);

    correctness &= error_checking(reference, matrix_C, n, 0);



    initMatrix0(matrix_C, size);

    printf("\nStarting GPU baseline\n");
    gettimeofday(&start_time, NULL);

    matrixMulCUDA(matrix_C, matrix_A, matrix_B, block_size, n);    

    gettimeofday(&stop_time, NULL);
    printf("Ended GPU baseline\n");

    elapsed_time = get_elapsed_time(start_time, stop_time);
    printf("GPU baseline execution time:%f\n", elapsed_time);

    correctness &= error_checking(reference, matrix_C, n, 1.e-6);



    initMatrix0(matrix_C, size);

    printf("\nStarting GPU tiled\n");
    gettimeofday(&start_time, NULL);

    matrixMulCUDA_tiled(matrix_C, matrix_A, matrix_B, block_size, n);

    gettimeofday(&stop_time, NULL);
    printf("Ended GPU tiled\n");

    elapsed_time = get_elapsed_time(start_time, stop_time);
    printf("GPU tiled execution time:%f\n", elapsed_time);

    correctness &= error_checking(reference, matrix_C, n, 1.e-6);


  // Clean up memory
    free(matrix_A);
    free(matrix_B);
    free(matrix_C);
    free(reference);

    cudaDeviceReset();

    if(correctness==true)
        return EXIT_SUCCESS;
    else
        return EXIT_FAILURE;
}


/**
 * Program main
 */
int main(int argc, char **argv)
{
    printf("[Matrix Multiply] - Starting...\n");

    if (argc != 2) {
        printf("Usage: %s <num_row/column> (must be multiple of 32)",argv[0]);
        exit(EXIT_FAILURE);
    }
    int n = atoi(argv[1]);
    if (n < 32 || n%32 != 0){
        printf("<num_row/column> must be multiple of 32");
        exit(EXIT_FAILURE);
    }

    int block_size = 32;

    printf("Matrix(%d,%d)\n", n, n);

    int matrix_result = matrixMultiply(block_size, n);

    exit(matrix_result);
}
