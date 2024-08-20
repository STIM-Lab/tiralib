#ifndef STIM_CUDA_ERROR_H
#define STIM_CUDA_ERROR_H

#include <stdio.h>
#include <iostream>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cufft.h"
#include "cublas_v2.h"

//handle error macro
static void cuHandleError( cudaError_t err, const char *file,  int line ) {
   	if (err != cudaSuccess) {
            printf("%s in %s at line %d\n", cudaGetErrorString( err ),  file, line );

   	}
}
#define HANDLE_ERROR( err ) (cuHandleError( err, __FILE__, __LINE__ ))
static void cufftHandleError( cufftResult err, const char*file, int line )
{
    if (err != CUFFT_SUCCESS)
    {
        if(err == CUFFT_INVALID_PLAN)
            std::cout<<"The plan parameter is not a valid handle."<<std::endl;
        else if(err == CUFFT_ALLOC_FAILED)
            std::cout<<"Allocation failed."<<std::endl;
        else if(err == CUFFT_INVALID_VALUE)
            std::cout<<"At least one of the parameters idata, odata, and direction is not valid."<<std::endl;
        else if(err == CUFFT_INTERNAL_ERROR)
            std::cout<<"An internal driver error was detected."<<std::endl;
        else if(err == CUFFT_EXEC_FAILED)
            std::cout<<"CUFFT failed to execute the transform on the GPU."<<std::endl;
        else if(err == CUFFT_SETUP_FAILED)
            std::cout<<"The CUFFT library failed to initialize."<<std::endl;
        else
            std::cout<<"Unknown error: "<<err<<std::endl;

    }
}
#define CUFFT_HANDLE_ERROR( err ) (cufftHandleError( err, __FILE__, __LINE__ ))

static void cublasHandleError( cublasStatus_t err, const char*file, int line ){
	if(err != CUBLAS_STATUS_SUCCESS){
		if(err == CUBLAS_STATUS_NOT_INITIALIZED)
			std::cout<<"CUBLAS_STATUS_NOT_INITIALIZED" <<" in file "<<file<<" line "<<std::endl;
		else if(err == CUBLAS_STATUS_ALLOC_FAILED)
			std::cout<<"CUBLAS_STATUS_ALLOC_FAILED" <<" in file "<<file<<" line "<<std::endl;
		else if(err == CUBLAS_STATUS_INVALID_VALUE)
			std::cout<<"CUBLAS_STATUS_INVALID_VALUE" <<" in file "<<file<<" line "<<std::endl;
		else if(err == CUBLAS_STATUS_ARCH_MISMATCH)
			std::cout<<"CUBLAS_STATUS_ARCH_MISMATCH" <<" in file "<<file<<" line "<<std::endl;
		else if(err == CUBLAS_STATUS_MAPPING_ERROR)
			std::cout<<"CUBLAS_STATUS_MAPPING_ERROR" <<" in file "<<file<<" line "<<std::endl;
		else if(err == CUBLAS_STATUS_EXECUTION_FAILED)
			std::cout<<"CUBLAS_STATUS_EXECUTION_FAILED" <<" in file "<<file<<" line "<<std::endl;
		else if(err == CUBLAS_STATUS_INTERNAL_ERROR)
			std::cout<<"CUBLAS_STATUS_INTERNAL_ERROR" <<" in file "<<file<<" line "<<std::endl;
		else
			std::cout<<"Unknown error"<<" in file "<<file<<" line "<<std::endl;
	}
}
#define CUBLAS_HANDLE_ERROR( err ) (cublasHandleError( err, __FILE__, __LINE__ ))


#endif
