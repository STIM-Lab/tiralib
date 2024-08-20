#ifndef STIM_COLORMAP_H
#define STIM_COLORMAP_H

#include <string>
#include <stdlib.h>
#include <cmath>

#ifdef _WIN32
	#include <float.h>
#endif

#ifdef __CUDACC__
#include "cublas_v2.h"
#include <tira/cuda/cudatools/error.h>
#endif

//saving an image to a file uses the CImg library
	//this currently throws a lot of "unreachable" warnings (as of GCC 4.8.2, nvcc 6.5.12)
#include <tira/image.h>


#define BREWER_CTRL_PTS 11

static float  BREWERCP[BREWER_CTRL_PTS*4] = {0.192157f, 0.211765f, 0.584314f, 1.0f,
                                      0.270588f, 0.458824f, 0.705882f, 1.0f,
                                      0.454902f, 0.678431f, 0.819608f, 1.0f,
                                      0.670588f, 0.85098f, 0.913725f, 1.0f,
                                      0.878431f, 0.952941f, 0.972549f, 1.0f,
                                      1.0f, 1.0f, 0.74902f, 1.0f,
                                      0.996078f, 0.878431f, 0.564706f, 1.0f,
                                      0.992157f, 0.682353f, 0.380392f, 1.0f,
                                      0.956863f, 0.427451f, 0.262745f, 1.0f,
                                      0.843137f, 0.188235f, 0.152941f, 1.0f,
                                      0.647059f, 0.0f, 0.14902f, 1.0f};

//static float  BREWERCP[BREWER_CTRL_PTS * 4] = { 0.192157f, 0.584314f, 0.211765f, 1.0f,
//										0.270588f, 0.705882f, 0.458824f, 1.0f,
//										0.454902f, 0.819608f, 0.678431f, 1.0f,
//										0.670588f, 0.913725f, 0.85098f, 1.0f,
//										0.878431f, 0.972549f, 0.952941f, 1.0f,
//										1.0f, 0.74902f, 1.0f, 1.0f,
//										0.996078f, 0.878431f, 0.564706f, 1.0f,
//										0.992157f, 0.682353f, 0.380392f, 1.0f,
//										0.956863f, 0.427451f, 0.262745f, 1.0f,
//										0.843137f, 0.188235f, 0.152941f, 1.0f,
//										0.647059f, 0.0f, 0.14902f, 1.0f };

//static float  BREWERCP[BREWER_CTRL_PTS * 4] = { 0.0f, 0.407843f, 0.215686f, 1.0f,
//										0.101960f, 0.596078f, 0.313725f, 1.0f,
//										0.4f, 0.741176f, 0.388235f, 1.0f,
//										0.650980f, 0.850980f, 0.415686f, 1.0f,
//										0.850980f, 0.937254f, 0.545098f, 1.0f,
//										1.0f, 1.0f, 0.749019f, 1.0f,
//										0.996078f, 0.878431f, 0.545098f, 1.0f,
//										0.992156f, 0.682352f, 0.380392f, 1.0f,
//										0.956862f, 0.427450f, 0.262745f, 1.0f,
//										0.843137f, 0.188235f, 0.152941f, 1.0f,
//										0.647058f, 0.0f, 0.149019f, 1.0f };


#ifdef __CUDACC__
texture<float4, cudaTextureType1D> cudaTexBrewer;
static cudaArray* gpuBrewer;
#endif

namespace tira {

	namespace colormap {
		enum colormapType { cmBrewer, cmGrayscale, cmRainbow };

		static void buffer2image(unsigned char* buffer, std::string filename, size_t width, size_t height)
		{
			/*unsigned char* non_interleaved = (unsigned char*)malloc(x_size * y_size * 3);
			unsigned int S = x_size * y_size;

			for(unsigned int i = 0; i < S; i++){
				non_interleaved[i + 0 * S] = buffer[i * 3 + 0];
				non_interleaved[i + 1 * S] = buffer[i * 3 + 1];
				non_interleaved[i + 2 * S] = buffer[i * 3 + 2];
			}*/

			//create an image object
			//cimg_library::CImg<unsigned char> image(non_interleaved, x_size, y_size, 1, 3);
			//image.save(filename.c_str());
			image<unsigned char> I(buffer, width, height, 3);
			//I.set_interleaved_rgb(buffer, width, height);
			I.save(filename);
		}

#ifdef __CUDACC__
		static void initBrewer()
		{
			//initialize the Brewer colormap

			//allocate CPU space
			float4 cpuColorMap[BREWER_CTRL_PTS];

			//define control rtsPoints
			cpuColorMap[0] = make_float4(0.192157f, 0.211765f, 0.584314f, 1.0f);
			cpuColorMap[1] = make_float4(0.270588f, 0.458824f, 0.705882f, 1.0f);
			cpuColorMap[2] = make_float4(0.454902f, 0.678431f, 0.819608f, 1.0f);
			cpuColorMap[3] = make_float4(0.670588f, 0.85098f, 0.913725f, 1.0f);
			cpuColorMap[4] = make_float4(0.878431f, 0.952941f, 0.972549f, 1.0f);
			cpuColorMap[5] = make_float4(1.0f, 1.0f, 0.74902f, 1.0f);
			cpuColorMap[6] = make_float4(0.996078f, 0.878431f, 0.564706f, 1.0f);
			cpuColorMap[7] = make_float4(0.992157f, 0.682353f, 0.380392f, 1.0f);
			cpuColorMap[8] = make_float4(0.956863f, 0.427451f, 0.262745f, 1.0f);
			cpuColorMap[9] = make_float4(0.843137f, 0.188235f, 0.152941f, 1.0f);
			cpuColorMap[10] = make_float4(0.647059f, 0.0f, 0.14902f, 1.0f);


			int width = BREWER_CTRL_PTS;
			int height = 0;


			// allocate array and copy colormap data
			cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32, 32, 32, 32, cudaChannelFormatKindFloat);

			HANDLE_ERROR(cudaMallocArray(&gpuBrewer, &channelDesc, width, height));

			HANDLE_ERROR(cudaMemcpyToArray(gpuBrewer, 0, 0, cpuColorMap, sizeof(float4) * width, cudaMemcpyHostToDevice));

			// set texture parameters
			cudaTexBrewer.addressMode[0] = cudaAddressModeClamp;
			//texBrewer.addressMode[1] = cudaAddressModeClamp;
			cudaTexBrewer.filterMode = cudaFilterModeLinear;
			cudaTexBrewer.normalized = true;  // access with normalized texture coordinates

			// Bind the array to the texture
			HANDLE_ERROR(cudaBindTextureToArray(cudaTexBrewer, gpuBrewer, channelDesc));

		}

		static void destroyBrewer()
		{
			HANDLE_ERROR(cudaFreeArray(gpuBrewer));
		}

		template<class T>
		__global__ static void applyBrewer(T* gpuSource, unsigned char* gpuDest, unsigned int N, T minVal = 0, T maxVal = 1)
		{

			int i = blockIdx.y * gridDim.x * blockDim.x + blockIdx.x * blockDim.x + threadIdx.x;
			if (i >= N) return;

			//compute the normalized value on [minVal maxVal]
			float a = (gpuSource[i] - minVal) / (maxVal - minVal);

			//compensate for the additional space at the edges
			a *= (T)(BREWER_CTRL_PTS - 1) / (T)(BREWER_CTRL_PTS);

			//lookup the color
			float shift = (T)1 / (2 * BREWER_CTRL_PTS);
			float4 color = tex1D(cudaTexBrewer, a + shift);
			//float4 color = tex1D(cudaTexBrewer, a);

			gpuDest[i * 3 + 0] = 255 * color.x;
			gpuDest[i * 3 + 1] = 255 * color.y;
			gpuDest[i * 3 + 2] = 255 * color.z;
		}

		template<class T>
		__global__ static void applyGrayscale(T* gpuSource, unsigned char* gpuDest, unsigned int N, T minVal = 0, T maxVal = 1)
		{
			int i = blockIdx.y * gridDim.x * blockDim.x + blockIdx.x * blockDim.x + threadIdx.x;
			if (i >= N) return;

			//compute the normalized value on [minVal maxVal]
			float a = (gpuSource[i] - minVal) / (maxVal - minVal);

			//threshold
			if (a > 1)
				a = 1;
			if (a < 0)
				a = 0;

			gpuDest[i * 3 + 0] = 255 * a;
			gpuDest[i * 3 + 1] = 255 * a;
			gpuDest[i * 3 + 2] = 255 * a;
		}

		template<class T>
		static void gpu2gpu(T* gpuSource, unsigned char* gpuDest, unsigned int nVals, T minVal = 0, T maxVal = 1, colormapType cm = cmGrayscale, int blockDim = 128)
		{
			//This function converts a scalar field on the GPU to a color image on the GPU
			int gridX = (nVals + blockDim - 1) / blockDim;
			int gridY = 1;
			if (gridX > 65535)
			{
				gridY = (gridX + 65535 - 1) / 65535;
				gridX = 65535;
			}
			dim3 dimGrid(gridX, gridY);
			if (cm == cmGrayscale)
				applyGrayscale << <dimGrid, blockDim >> > (gpuSource, gpuDest, nVals, minVal, maxVal);
			else if (cm == cmBrewer)
			{
				initBrewer();
				applyBrewer << <dimGrid, blockDim >> > (gpuSource, gpuDest, nVals, minVal, maxVal);
				destroyBrewer();
			}

		}

		template<class T>
		static void gpu2cpu(T* gpuSource, unsigned char* cpuDest, unsigned int nVals, T minVal, T maxVal, colormapType cm = cmGrayscale)
		{
			//this function converts a scalar field on the GPU to a color image on the CPU

			//first create the color image on the GPU

			//allocate GPU memory for the color image
			unsigned char* gpuDest;
			HANDLE_ERROR(cudaMalloc((void**)&gpuDest, sizeof(unsigned char) * nVals * 3));

			//create the image on the gpu
			gpu2gpu(gpuSource, gpuDest, nVals, minVal, maxVal, cm);

			//copy the image from the GPU to the CPU
			HANDLE_ERROR(cudaMemcpy(cpuDest, gpuDest, sizeof(unsigned char) * nVals * 3, cudaMemcpyDeviceToHost));

			HANDLE_ERROR(cudaFree(gpuDest));

		}

		template<typename T>
		static void gpu2image(T* gpuSource, std::string fileDest, unsigned int x_size, unsigned int y_size, T valMin, T valMax, colormapType cm = cmGrayscale)
		{
			//allocate a color buffer
			unsigned char* cpuBuffer = NULL;
			cpuBuffer = (unsigned char*)malloc(sizeof(unsigned char) * 3 * x_size * y_size);

			//do the mapping
			gpu2cpu<T>(gpuSource, cpuBuffer, x_size * y_size, valMin, valMax, cm);

			//copy the buffer to an image
			buffer2image(cpuBuffer, fileDest, x_size, y_size);

			free(cpuBuffer);
		}

		/// save a GPU image to a file using automatic scaling
		template<typename T>
		static void gpu2image(T* gpuSource, std::string fileDest, unsigned int x_size, unsigned int y_size, colormapType cm = cmGrayscale) {
			size_t N = x_size * y_size;								//calculate the total number of elements in the image

			cublasStatus_t stat;
			cublasHandle_t handle;

			stat = cublasCreate(&handle);							//create a cuBLAS handle
			if (stat != CUBLAS_STATUS_SUCCESS) {						//test for failure
				printf("CUBLAS initialization failed\n");
				exit(1);
			}

			int i_min, i_max;
			stat = cublasIsamin(handle, (int)N, gpuSource, 1, &i_min);
			if (stat != CUBLAS_STATUS_SUCCESS) {						//test for failure
				printf("CUBLAS Error: failed to calculate minimum r value.\n");
				exit(1);
			}
			stat = cublasIsamax(handle, (int)N, gpuSource, 1, &i_max);
			if (stat != CUBLAS_STATUS_SUCCESS) {						//test for failure
				printf("CUBLAS Error: failed to calculate maximum r value.\n");
				exit(1);
			}
			cublasDestroy(handle);

			i_min--;				//cuBLAS uses 1-based indexing for Fortran compatibility
			i_max--;
			T v_min, v_max;											//allocate space to store the minimum and maximum values
			HANDLE_ERROR(cudaMemcpy(&v_min, gpuSource + i_min, sizeof(T), cudaMemcpyDeviceToHost));		//copy the min and max values from the device to the CPU
			HANDLE_ERROR(cudaMemcpy(&v_max, gpuSource + i_max, sizeof(T), cudaMemcpyDeviceToHost));



			gpu2image<T>(gpuSource, fileDest, x_size, y_size, min(v_min, v_max), max(v_min, v_max), cm);
		}

#endif

		template<class T>
		static void cpuApplyBrewer(T* cpuSource, unsigned char* cpuDest, size_t N, T minVal = 0, T maxVal = 1)
		{
			for (size_t i = 0; i < N; i++)
			{
				//compute the normalized value on [minVal maxVal]
				float a;
				if (minVal != maxVal)
					a = (cpuSource[i] - minVal) / (maxVal - minVal);
				else
					a = 0.5;
#ifdef _WIN32
				if (!_finite(a)) a = 1;							//deal with infinite and NaN values (return maximum in all cases)
#else
				if (!std::isfinite(a)) a = 1;
#endif
				else if (a < 0) a = 0;
				else if (a > 1) a = 1;

				float c = a * (float)(BREWER_CTRL_PTS - 1);
				int ptLow = (int)c;
				float m = c - (float)ptLow;
				//std::cout<<m<<std::endl;

				float r, g, b;
				if (ptLow == BREWER_CTRL_PTS - 1)
				{
					r = BREWERCP[ptLow * 4 + 0];
					g = BREWERCP[ptLow * 4 + 1];
					b = BREWERCP[ptLow * 4 + 2];
				}
				else
				{
					r = BREWERCP[ptLow * 4 + 0] * (1 - m) + BREWERCP[(ptLow + 1) * 4 + 0] * m;
					g = BREWERCP[ptLow * 4 + 1] * (1 - m) + BREWERCP[(ptLow + 1) * 4 + 1] * m;
					b = BREWERCP[ptLow * 4 + 2] * (1 - m) + BREWERCP[(ptLow + 1) * 4 + 2] * m;
				}


				cpuDest[i * 3 + 0] = (unsigned char)(255 * r);
				cpuDest[i * 3 + 1] = (unsigned char)(255 * g);
				cpuDest[i * 3 + 2] = (unsigned char)(255 * b);

			}
		}

		template<class T>
		static void cpu2cpu(T* cpuSource, unsigned char* cpuDest, size_t nVals, T valMin, T valMax, colormapType cm = cmGrayscale)
		{

			if (cm == cmBrewer)
				cpuApplyBrewer(cpuSource, cpuDest, nVals, valMin, valMax);
			else if (cm == cmGrayscale)
			{
				int i;
				float a;
				float range = valMax - valMin;

				for (i = 0; i < nVals; i++)
				{
					//normalize to the range [valMin valMax]
					if (range != 0)
						a = (cpuSource[i] - valMin) / range;
					else
						a = 0.5;

					if (a < 0) a = 0;
					if (a > 1) a = 1;

					cpuDest[i * 3 + 0] = (unsigned char)(255 * a);
					cpuDest[i * 3 + 1] = (unsigned char)(255 * a);
					cpuDest[i * 3 + 2] = (unsigned char)(255 * a);
				}
			}
		}

		template<class T>
		static void cpu2cpu(T* cpuSource, unsigned char* cpuDest, unsigned long long nVals, colormapType cm = cmGrayscale)
		{
			//computes the max and min range automatically

			//find the largest magnitude value
			T maxVal = cpuSource[0];
			T minVal = cpuSource[0];
			for (int i = 1; i < nVals; i++)
			{
				if (cpuSource[i] > maxVal)
					maxVal = cpuSource[i];
				if (cpuSource[i] < minVal)
					minVal = cpuSource[i];
			}

			cpu2cpu(cpuSource, cpuDest, nVals, minVal, maxVal, cm);

		}



		template<typename T>
		static void cpu2image(T* cpuSource, std::string fileDest, size_t x_size, size_t y_size, T valMin, T valMax, colormapType cm = cmGrayscale)
		{
			//allocate a color buffer
			unsigned char* cpuBuffer = (unsigned char*)malloc(sizeof(unsigned char) * 3 * x_size * y_size);

			//do the mapping
			cpu2cpu<T>(cpuSource, cpuBuffer, x_size * y_size, valMin, valMax, cm);

			//copy the buffer to an image
			buffer2image(cpuBuffer, fileDest, x_size, y_size);

			free(cpuBuffer);

		}

		template<typename T>
		static void cpu2image(T* cpuSource, std::string fileDest, size_t x_size, size_t y_size, colormapType cm = cmGrayscale)
		{
			//allocate a color buffer
			unsigned char* cpuBuffer = (unsigned char*)malloc(sizeof(unsigned char) * 3 * x_size * y_size);

			//do the mapping
			cpu2cpu<T>(cpuSource, cpuBuffer, x_size * y_size, cm);

			//copy the buffer to an image
			buffer2image(cpuBuffer, fileDest, x_size, y_size);

			free(cpuBuffer);

		}
	}	// end namespace colormap
}	//end namespace tira
#endif

