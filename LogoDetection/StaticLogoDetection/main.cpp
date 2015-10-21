#include <fstream>
#include <opencv2\core\core.hpp>
#include <opencv2\imgproc\imgproc.hpp>
#include <opencv2\highgui\highgui.hpp>
#include <string>
#include <ctime>
#include "ImgProcess.h"

using namespace std;
using namespace cv;
///////////////////////////////////////////////////////////////////////////////////
#define Pi  3.14159265359
////////////////////////////////////////////////////////////////////////////////

bool isAPlaceBlank(Mat *I){
	for (int i = 0; i < (*I).rows; i++){
		for (int j = 0; j < (*I).cols; j++){
			if ((*I).at<uchar>(i, j) != 0)
				return true;
		}
	}
	return false;
}
void differential(double arrey[], int size)
{
	for (int i = 1; i<size; i++)
	{
		arrey[i - 1] = arrey[i] - arrey[i - 1];
	}
	arrey[size - 1] = 0;

}

void profile(Mat image, double arrey[])
{
	for (int i = 0; i < image.cols; i++)
	{
		for (int j = 0; j < image.rows; j++)
			arrey[i] += image.at<uchar>(j, i);
	}

}
// in this function, the variance of the an array will be calculated 
double variance(double arrey[], double size) 
{
	double sum = 0;
	double average = 0, var = 0;
	for (int i = 0; i<size; i++)
		sum += arrey[i];
	average = sum / size;
	for (int i = 0; i<size; i++)
		var += (arrey[i] - average)*(arrey[i] - average);
	var = var / size;
	return var;
}
// hisogram the logos on the side of picture to detrive the coordinate of the picture. 
int projection(Mat final, int final_coord[2][2][10]) 
{

	int *colsy;
	colsy = (int *)malloc(final.cols*sizeof(int));

	for (int i = 0; i < final.cols; i++)
	{
		colsy[i] = 0;

	}
	int key = 0;
	int counter = 0;
	int Coord_arrey[10][2] = { 0 };
	int sum = 0;
	int min_sum = 50000;
	int tmp_coord;
	for (int i = 0; i < (int)(final.cols); i++)
	{
		for (int j = 0; j < (int)final.rows; j++)
		{
			colsy[i] += int(final.at<uchar>(j, i));

		}
	}

	//////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////
	for (int i = 0; i < final.cols; i++)
	{
		if (colsy[i] != 0 && key == 0)
		{
			key = 1;
			tmp_coord = i;
			sum += colsy[i];

		}
		if (colsy[i] != 0 && key == 1)
			sum += colsy[i];

		if ((final.cols - 1) == i && colsy[i] != 0 && sum >min_sum && key == 1)
		{
			Coord_arrey[counter][0] = tmp_coord;
			Coord_arrey[counter][1] = i;
			counter++;
			sum = 0;
			key = 0;
		}

		if (colsy[i] == 0 && key == 1)
		{
			key = 0;
			for (int p = 0; p<10; p++)
			{
				if (i + p<final.cols){
					if (colsy[i + p] != 0)
					{
						key = 1;
						break;
					}
				}
			}

			if (key == 0 && sum>min_sum)
			{
				Coord_arrey[counter][0] = tmp_coord;
				Coord_arrey[counter][1] = i - 1;
				counter++;
				sum = 0;
			}
		}

	}


	////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////
	int *rowsy;
	rowsy = (int *)malloc(final.rows*sizeof(int));
	int counter2[10] = { 0 };
	int Coord_arrey_2[10][2][10] = { 0 };
	sum = 0;
	for (int c = 0; c<counter; c++)
	{
		key = 0;
		for (int i = 0; i < final.rows; i++)
		{
			rowsy[i] = 0;

		}

		for (int i = 0; i < final.rows; i++)
		{
			for (int j = Coord_arrey[c][0]; j < Coord_arrey[c][1]; j++)
			{
				rowsy[i] += (int)(final.at<uchar>(i, j));

			}
		}
		///////////////////////////////////////////////////////////////////////////////
		for (int i = 0; i < final.rows; i++)
		{
			if (rowsy[i] != 0 && key == 0)
			{

				key = 1;
				tmp_coord = i;
				sum += rowsy[i];
			}

			if (rowsy[i] != 0 && key == 1)
				sum += rowsy[i];

			if ((final.rows - 1) == i && rowsy[i] != 0 && sum >min_sum && key == 1)
			{
				Coord_arrey_2[counter2[c]][0][c] = tmp_coord;
				Coord_arrey_2[counter2[c]][1][c] = i;
				counter2[c]++;
				sum = 0;
				key = 0;
			}

			if (rowsy[i] == 0 && key == 1)
			{
				key = 0;
				for (int p = 0; p<10 && i + p<final.rows; p++)
				{
					if (i + p<final.rows){
						if (rowsy[i + p] != 0)
						{
							key = 1;
							break;
						}
					}
				}

				if (key == 0 && sum>min_sum)
				{
					Coord_arrey_2[counter2[c]][0][c] = tmp_coord;
					Coord_arrey_2[counter2[c]][1][c] = i - 1;
					counter2[c]++;
					sum = 0;
				}
			}
		}
	}
	free(rowsy);
	/////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////
	int f_counter = 0;
	key = 0;
	sum = 0;
	for (int c = 0; c<counter; c++)
	{
		if (counter2[c]>1)
		{


			for (int d = 0; d<counter2[c]; d++)
			{

				final_coord[1][0][f_counter] = Coord_arrey_2[d][0][c];
				final_coord[1][1][f_counter] = Coord_arrey_2[d][1][c];

				for (int i = 0; i < final.cols; i++)
				{
					colsy[i] = 0;

				}

				for (int i = Coord_arrey[c][0]; i <= Coord_arrey[c][1]; i++)
				{
					for (int j = Coord_arrey_2[d][0][c]; j <= Coord_arrey_2[d][1][c]; j++)
					{
						colsy[i] += (int)(final.at<uchar>(j, i));

					}

				}

				//////////////////////////////////////////////////////////////////////////////////
				for (int i = Coord_arrey[c][0]; i <= Coord_arrey[c][1]; i++)
				{
					if (colsy[i] != 0 && key == 0)
					{
						key = 1;
						tmp_coord = i;
						sum += colsy[i];
					}

					if (colsy[i] != 0 && key == 1)
						sum += colsy[i];

					if (Coord_arrey[c][1] == i && colsy[i] != 0 && sum >min_sum && key == 1)
					{
						final_coord[0][0][f_counter] = tmp_coord;
						final_coord[0][1][f_counter] = i;
						f_counter++;
						key = 0;
						sum = 0;

					}
					if (colsy[i] == 0 && key == 1)
					{

						key = 0;
						for (int p = 0; p<10 && i + p <= Coord_arrey[1][c]; p++)
						{
							if (colsy[i + p] != 0)
							{
								key = 1;
								break;
							}
						}

						if (key == 0 && sum>min_sum)
						{
							final_coord[0][0][f_counter] = tmp_coord;
							final_coord[0][1][f_counter] = i - 1;
							f_counter++;
							sum = 0;
						}
					}

				}


			}

		}

		else
		{
			final_coord[0][0][f_counter] = Coord_arrey[c][0];
			final_coord[1][0][f_counter] = Coord_arrey_2[counter2[c] - 1][0][c];
			final_coord[0][1][f_counter] = Coord_arrey[c][1];
			final_coord[1][1][f_counter] = Coord_arrey_2[counter2[c] - 1][1][c];
			f_counter++;
		}
	}
	return f_counter;
}

/*void sinc(double x, double y, double  a1, double a2, Mat *image)
{
int Landa = 800, Betta = 100;
double brightness;
int offset_cols = (*image).cols;
int offset_rows = (*image).rows;
double P_a1 = a1*a1;
double P_a2 = a2*a2;
brightness = -1 * Landa*sin(Pi*pow(P_a1*(x - offset_cols / 2)*(x - offset_cols / 2)*1.0 + P_a2*(y - offset_rows / 2)*(y - offset_rows / 2)*1.0, 0.5)) / (Pi*pow(P_a1*(x - offset_cols / 2)*(x - offset_cols / 2)*1.0 + P_a2*(y - offset_rows / 2)*(y - offset_rows / 2)*1.0, 0.5)) + Betta; // Y must be more than X in rectangular image when cols is more than rows

if (brightness > 255)
brightness = 255;
if (brightness < 0)

brightness = 0;
if (brightness < 50)
(*image).at<uchar>(y, x) = (uchar)brightness;
}*/
double finding_Alpha(double M)
{
	double a = 6.1*sqrt(2.0) / (M*Pi);
	return a;
}



//////////////////////////////////////////////////////////////////////////////////
int main()
{
	//ofstream myfile; 
	//myfile.open("alireza.txt");
	double  binary_treshold = 60, var_dif_treshhold = 10000, min_frame_count = 100, frame_count = 1;
	double *colsy;
	//createdouble()
	double var = 10;
	cout << endl << "detecting the logo(s)..." << endl;

	string address = "C://Users//Alireza//Downloads//Video//h_holand_czech.avi";
	VideoCapture capture(address);
	Mat frame;
	capture >> frame;
	Mat img_org = frame.clone();
	/*
	camera.set(CV_CAP_PROP_FRAME_WIDTH, 640);
	camera.set(CV_CAP_PROP_FRAME_HEIGHT, 480);*/
	//camera >> frame;
	Mat sGray(frame.size(), CV_8U, createImageBuffer(frame.size().width*frame.size().height));
	Mat dGray(frame.size(), CV_8U, createImageBuffer(frame.size().width*frame.size().height));
	Mat eGray(frame.size(), CV_8U, createImageBuffer(frame.size().width*frame.size().height));
	Mat final(frame.size(), CV_8U, createImageBuffer(frame.size().width*frame.size().height));
	Mat generalizedGradient(frame.size(), CV_32FC1, createImageBufferFloat(frame.size().width*frame.size().height));

	//Mat generalizedGradient(frame.rows, frame.cols, CV_32FC1);
	cvtColor(frame, dGray, CV_BGR2GRAY);
	cvtColor(frame, eGray, CV_BGR2GRAY);
	//colsy = (double *)malloc(dGray.cols*sizeof(double));
	colsy = createdouble(dGray.cols*sizeof(double));
	for (int i = 0; i < dGray.cols; i++)
		colsy[i] = 0;
	//Mat generalizedGradient(frame.rows, frame.cols, CV_32FC1);
	//Mat final(frame.rows, frame.cols, CV_8UC1);

	int flag_var = 0;

	clock_t begin = clock();
	for (int frame_count01; frame_count<500; frame_count++)
	{
		capture >> frame;
		
		cvtColor(frame, sGray, CV_BGR2GRAY);

		boxfilter(frame.size().width, frame.size().height, sGray.data, dGray.data, 5, 5);
		///// sobel filter, without scharre algorithm //////////////////
		
		sobelfilter(frame.size().width, frame.size().height, dGray.data, eGray.data);

		/* sinc : in this function I calculate the the focal point of semi-ellipse shape and then appied it on the image.
		mostly, the center of image doesn't consist of specific logos. this fact make an condition for us to decide to eliminate it
		to reduce the amount of processing in out alghorithm*/
		sinc(frame.size().width, frame.size().height, finding_Alpha(eGray.cols), finding_Alpha(eGray.rows), eGray.data, eGray.data);
		
		generalgradient(frame.size().width, frame.size().height, frame_count, eGray.data, (float *)generalizedGradient.data);

		general2final(frame.size().width, frame.size().height, (float*)generalizedGradient.data, final.data);

		/* make the pixels of the image binary in such a way that saturate the pixels to 255 or 0 */
		
		treshold(frame.size().width, frame.size().height, binary_treshold, final.data, final.data);

		profile(frame.size().width, frame.size().height, final.data, colsy);

		differential(colsy, final.cols);

		

		if ((abs(var - variance(colsy, final.cols - 1))> var_dif_treshhold))
		{
			flag_var = 0;
			var = variance(colsy, generalizedGradient.cols - 1);
		}
		else if (flag_var<10)
			flag_var++;
		else if (flag_var == 10 && frame_count>min_frame_count)
		{
			break;
		}

		for (int i = 0; i < dGray.cols; i++)
			colsy[i] = 0;

	}
	/////////////////////////////////////final steps////////////////////////////////////
	/* In this scope, I finalize the position of the Logos and eliminate it in the orginal picture. */
	cout << endl << endl << "the number of used frames : " << frame_count << endl;
	int final_coord[2][2][10] = { 0 };
	int x;
	x = projection(final, final_coord);
	cout << endl << "number of found logo(s): " << x << endl << endl << "the coordinations of logos: " << endl << endl;
	for (int i = 0; i < x; i++)
		cout << "row1 = " << final_coord[1][0][i] << "	col1 = " << final_coord[0][0][i] << "	row2 = " << final_coord[1][1][i] << "	col2 = " << final_coord[0][1][i] + 2 << endl << endl;
	
	for (int p = 0; p < x; p++)
		for (int i = final_coord[1][0][p]; i <= final_coord[1][1][p]; i++)
			for (int j = final_coord[0][0][p]; j <= final_coord[0][1][p]; j++)
				for (int k = 0; k < 3; k++)
				{
					img_org.at<Vec3b>(i, j)[k] = 0;
				}
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout << elapsed_secs << endl;

	imshow("Final picture", img_org);
	//imwrite("Gray_Image.jpg", img_org);
	waitKey(0);

	///////////////////////// destroy the memories /////////////////////////
	desetroyImageBuffer(sGray.data);
	desetroyImageBuffer(dGray.data);
	desetroyImageBuffer(eGray.data);
	desetroyImageBuffer(final.data);
	desetroyImageBuffer(generalizedGradient.data);
	///////////////////////////////////////////////////////////////////////
}
