#include <iostream>

extern "C" void boxfilter(int iw, int ih, unsigned char *source, unsigned char *dest, int bw, int bh); // function for bluring the images
extern "C" void sobelfilter(int iw, int ih, unsigned char *source, unsigned char *dest); // to detecting the edges of the picture 
extern "C" void treshold(int iw, int ih, int binary_treshold, unsigned char *source, unsigned char *dest); // function to set the pixels to saturated level 
extern "C" unsigned char * createImageBuffer(unsigned int Bytes); // allocate memory in ram as pinned memory
extern "C" double * createdouble(double Bytes); // pinned memory doubled
extern "C" float * createImageBufferFloat(unsigned int Bytes); //allocate memory for picture which is declared as float for each of its p
extern "C" void	 desetroyImageBuffer(unsigned char* bytes); // destroy the pinned memory 
extern "C" void sinc(int iw, int ih, double a1, double a2, unsigned char *source, unsigned char *dest); // eliminate the pixels which is placed in center of the picture. 
extern "C" void generalgradient(int iw, int ih, int frameCount, unsigned char  *source, float *dest); //  make gradient from the picture
extern "C" void profile(int iw, int ih, unsigned char  *img, double *myarray); // calculate the histogram of the picture
extern "C" void general2final(int iw, int ih, float  *source, unsigned char *dest); // transfer the data from generalized frame to final picture. 