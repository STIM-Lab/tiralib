#ifndef CUDA_CALLABLE

// define the CUDA_CALLABLE macro (will prefix all members)
#ifdef __CUDACC__
#define CUDA_CALLABLE __host__ __device__ inline
#else
#define CUDA_CALLABLE
#endif

#ifdef __CUDACC__
#define CUDA_UNCALLABLE __host__ inline
#else
#define CUDA_UNCALLABLE
#endif

#endif
