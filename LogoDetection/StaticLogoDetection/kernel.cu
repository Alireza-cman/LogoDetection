#include "device_atomic_functions.hpp"
#include "device_functions.hpp"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "ImgProcess.h"
#include <stdio.h>
#define Pi  3.14159265359
__global__  void general2final_kernel(int iw, int ih, float *source, unsigned char *dest)
{
	int x = blockDim.x*blockIdx.x + threadIdx.x;
	int y = blockDim.y*blockIdx.y + threadIdx.y;
	dest[iw*y + x] = (unsigned char)source[iw*y + x];
}
__global__  void treshold_kernel(int iw, int ih, int binary_treshold, unsigned char *source, unsigned char *dest)
{
	int x = blockDim.x*blockIdx.x + threadIdx.x;
	int y = blockDim.y*blockIdx.y + threadIdx.y;

	if (y < 10)
		dest[iw*y + x] = 0;
	if (y > ih - 10)
		dest[iw*y + x] = 0;
	if (x < 10)
		dest[iw*y + x] = 0;
	if (x > iw - 10)
		dest[iw*y + x] = 0;

	if ((unsigned char)dest[iw*y + x] > 60)
		dest[iw*y + x] = 255;
	else
		dest[iw*y + x] = 0;


	__syncthreads();

}
__global__  void Profile_kernel(int iw, int ih, unsigned char *source, double *dest)
{
	int x = blockDim.x*blockIdx.x + threadIdx.x;
	int y = blockDim.y*blockIdx.y + threadIdx.y;
	if (x >= 0 && x < iw  && y >= 0 && y < ih)
	{
		dest[x] += source[iw*y + x];
		//atomicAdd(&dest[x], source[iw*y + x]); // so better but I dont know why it doesnt declared :| 
	}


}

__global__ void sinc_kernel(int iw, int ih, double a1, double a2, unsigned char *source, unsigned char *dest)
{
	int x = blockDim.x*blockIdx.x + threadIdx.x;
	int y = blockDim.y*blockIdx.y + threadIdx.y;
	int Landa = 800, Betta = 100;
	double brightness;
	int offset_cols = iw;
	int offset_rows = ih;
	double P_a1 = a1*a1;
	double P_a2 = a2*a2;
	brightness = -1 * Landa*sin(Pi*pow(P_a1*(x - offset_cols / 2)*(x - offset_cols / 2)*1.0 + P_a2*(y - offset_rows / 2)*(y - offset_rows / 2)*1.0, 0.5)) / (Pi*pow(P_a1*(x - offset_cols / 2)*(x - offset_cols / 2)*1.0 + P_a2*(y - offset_rows / 2)*(y - offset_rows / 2)*1.0, 0.5)) + Betta; // Y must be more than X in rectangular image when cols is more than rows

	if (brightness > 255)
		brightness = 255;
	if (brightness < 0)

		brightness = 0;
	if (brightness < 50)
		dest[iw*y + x] = (unsigned char)brightness;
}
__global__ void generalgradient_kernel(int iw, int ih, int frameCount, unsigned char *source, float *dest)
{
	int x = blockDim.x*blockIdx.x + threadIdx.x;
	int y = blockDim.y*blockIdx.y + threadIdx.y;
	float temp;
	if (x > 0 && x < iw - 1 && y > 0 && y < ih - 1)
	{
		temp = dest[iw*y + x];
		dest[iw*y + x] = (float)(1.0*((frameCount - 1)*temp + source[iw*y + x]) / frameCount);
	}

}


__global__ void boxfilter_kernel(int iw, int ih, unsigned char *source, unsigned char *dest, int bw, int bh)
{
	int x = blockDim.x*blockIdx.x + threadIdx.x;
	int y = blockDim.y*blockIdx.y + threadIdx.y;

	int count = 0;
	float sum = 0.0;

	for (int j = -(bh / 2); j <= (bh / 2); j++)
		for (int i = -(bw / 2); i <= (bw / 2); i++)
		{
			if ((x + i) < iw && (x + i) >= 0 && (y + j) < ih && (y + j) >= 0)
			{
				sum += (float)source[((y + j)*iw) + (x + i)];
				count++;
			}
		}
	sum /= (float)count * 2;
	dest[(y*iw) + x] = (unsigned char)sum;
}
__global__ void sobelfilter_kernel(int iw, int ih, unsigned char *source, unsigned char *dest)
{
	int x = blockDim.x*blockIdx.x + threadIdx.x;
	int y = blockDim.y*blockIdx.y + threadIdx.y;
	if (x > 0 && x < iw - 1 && y > 0 && y < ih - 1)
	{
		int gx = -1 * source[iw*(y - 1) + (x - 1)] + source[iw*(y - 1) + (x + 1)] +
			-2 * source[iw*y + (x - 1)] + 2 * source[iw*y + (x + 1)] +
			-1 * source[iw*(y + 1) + (x - 1)] + source[iw*(y + 1) + (x + 1)];
		int gy = -source[iw*(y - 1) + (x - 1)] - 2 * source[iw*(y - 1) + x]
			- source[iw*(y - 1) + (x + 1)] +
			source[iw*(y + 1) + (x - 1)] + 2 * source[iw*(y + 1) + x] +
			source[iw*(y + 1) + (x + 1)];
		dest[iw*y + x] = (unsigned char)sqrt((float)gx*(float)gx + (float)gy*float(gy));

	}
}

extern "C" void boxfilter(int iw, int ih, unsigned char *source, unsigned char *dest, int bw, int bh)
{
	unsigned char *dev_source, *dev_dest;
	cudaHostGetDevicePointer(&dev_source, source, 0);
	cudaHostGetDevicePointer(&dev_dest, dest, 0);

	dim3 blocks(iw / 16, ih / 16);
	dim3 threads(16, 16);
	boxfilter_kernel << <blocks, threads >> >(iw, ih, dev_source, dev_dest, bw, bh);
	cudaThreadSynchronize();
}
extern "C" void sobelfilter(int iw, int ih, unsigned char *source, unsigned char *dest)
{
	// allocate memory for bitmap 
	unsigned char *dev_source, *dev_dest;
	cudaHostGetDevicePointer(&dev_source, source, 0);
	cudaHostGetDevicePointer(&dev_dest, dest, 0);
	dim3	block(iw / 16, ih / 16);
	dim3	threads(16, 16);
	sobelfilter_kernel << <block, threads >> >(iw, ih, dev_source, dev_dest);
	cudaThreadSynchronize();

}
extern "C" unsigned char* createImageBuffer(unsigned int bytes)
{
	unsigned char *ptr = NULL;
	cudaSetDeviceFlags(cudaDeviceMapHost);
	cudaHostAlloc(&ptr, bytes, cudaHostAllocMapped);
	return ptr;

}
extern "C" float * createImageBufferFloat(unsigned int Bytes)
{
	float *ptr = NULL;
	cudaSetDeviceFlags(cudaDeviceMapHost);
	cudaHostAlloc(&ptr, Bytes, cudaHostAllocMapped);
	return ptr;

}
void desetroyImageBuffer(unsigned char* bytes)
{
	cudaFreeHost(bytes);
}
void sinc(int iw, int ih, double a1, double a2, unsigned char *source, unsigned char *dest)
{
	unsigned char *dev_source, *dev_dest;
	cudaHostGetDevicePointer(&dev_source, source, 0);
	cudaHostGetDevicePointer(&dev_dest, dest, 0);

	dim3 blocks(iw / 16, ih / 16);
	dim3 threads(16, 16);
	sinc_kernel << <blocks, threads >> >(iw, ih, a1, a2, dev_source, dev_dest);
	cudaThreadSynchronize();
}
void generalgradient(int iw, int ih, int frameCount, unsigned char *source, float  *dest)
{
	unsigned char *dev_source;
	float *dev_dest;
	cudaHostGetDevicePointer(&dev_source, source, 0);
	cudaHostGetDevicePointer(&dev_dest, dest, 0);

	dim3 blocks(iw / 16, ih / 16);
	dim3 threads(16, 16);
	generalgradient_kernel << <blocks, threads >> >(iw, ih, frameCount, dev_source, dev_dest);
	cudaThreadSynchronize();
}
extern "C" void treshold(int iw, int ih, int binary_treshold, unsigned char *source, unsigned char *dest)
{
	unsigned char *dev_source, *dev_dest;
	cudaHostGetDevicePointer(&dev_source, source, 0);
	cudaHostGetDevicePointer(&dev_dest, dest, 0);

	dim3 blocks(iw / 16, ih / 16);
	dim3 threads(16, 16);
	treshold_kernel << <blocks, threads >> >(iw, ih, binary_treshold, dev_source, dev_dest);
	cudaThreadSynchronize();
}
extern "C" void profile(int iw, int ih, unsigned char  *img, double *myarray)
{
	unsigned char *dev_source;
	double *dev_dest;
	cudaHostGetDevicePointer(&dev_source, img, 0);
	cudaHostGetDevicePointer(&dev_dest, myarray, 0);

	dim3 blocks(iw / 16, ih / 16);
	dim3 threads(16, 16);
	Profile_kernel << <blocks, threads >> >(iw, ih, dev_source, dev_dest);
	cudaThreadSynchronize();
}
extern "C" double * createdouble(double Bytes)
{
	double *ptr = NULL;
	cudaSetDeviceFlags(cudaDeviceMapHost);
	cudaHostAlloc(&ptr, Bytes, cudaHostAllocMapped);
	return ptr;

}
extern "C" void general2final(int iw, int ih, float  *source, unsigned char *dest)
{
	float *dev_source;
	unsigned char *dev_dest;
	cudaHostGetDevicePointer(&dev_source, source, 0);
	cudaHostGetDevicePointer(&dev_dest, dest, 0);

	dim3 blocks(iw / 16, ih / 16);
	dim3 threads(16, 16);
	general2final_kernel << <blocks, threads >> >(iw, ih, dev_source, dev_dest);
	cudaThreadSynchronize();
}


