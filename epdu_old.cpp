#define WANT_STREAM
#define WANT_MATH
#define WANT_TIME

#include "include.h"
#include "rng.h"
#include "main.h"
#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include "epdu.h"
#include <algorithm>
#include <math.h>

//#ifdef use_namespace
using namespace KW_RNG;
//#endif

epdu::epdu(std::string inputFile)
{
	double num;
	numbins = 0;
	std::ifstream inFile;
	std::string buffer;
	int i;

	outFile.open("output.txt", ios::out);
	if (outFile.fail()){
	    cout << "Could not open 'output.txt'" << endl;
		exit(1);
    }
	inFile.open(inputFile.c_str(), ios::in);
	if (inFile.fail()){
	    cout << "Could not open " << inputFile.c_str() << endl;
		exit(1);
    }

    inFile >> buffer;

	// load the bin counts into the vector a and count the number of bins.
	while (buffer!=";"){
		num = atof(buffer.c_str());
		a.push_back(num);
		inFile >> buffer;
		numbins++;
    }
	// because of the "curse of dimensionality", you need more 
	// iterations when you have more bins (= dimensions).
	B = 5000*ceil(numbins/10.0);
	inFile >> buffer;
	// load the smoothing parameter
	s = atof(buffer.c_str());
	inFile >> buffer; // this is a semi-colon
	inFile >> buffer; // this is either a number or the EOF
	// load the optional simultaneous credibility level
	if (! inFile.fail()){
		num = atof(buffer.c_str());
		alpha = num;
		inFile >> buffer; // this is a semi-colon
		inFile >> buffer; // this is either a number or the EOF
	} else {
		alpha = 0.05;
	}
	// load the optional prior distribution
	if (! inFile.fail()) {
		while (buffer!=";"){
			num = atof(buffer.c_str());
			prior.push_back(num);
			inFile >> buffer;
		}
	} else {
		for (int i0=0; i0<a.size(); i0++){
			prior.push_back(1);
		}
	}
    inFile.close();

	// initialize tempcol
	for (i=0; i<numbins; i++){
		temprow.push_back(0);
	}
	// initialize temprow
	for (i=0; i<B; i++){
		tempcol.push_back(0);
	}
	// initialize r
	for (i=0; i<numbins; i++){
		r.push_back(tempcol);		
	}

	cl = a;
	cu = a;

	// calculate the number of "grains" k
	k = 0;
	for (i=0; i<numbins; i++){
		k += a[i];
	}

	// initialize withinCI
	for (i=0; i<B; i++){
		withinCI.push_back(1);
	}
}

epdu::~epdu()
{
	outFile.close();
}

bool epdu::isSmooth(void)
{
	bool answer = (s!=0);
	return answer;
}

void epdu::sbcb()
{
	double CIalpha;
	std::vector<std::vector<double> > sortedR;

	dirrnd();
	CIalpha = findCBalpha();
	sortedR = sortCols(r);
	confPrctile(sortedR,CIalpha);
	printResults();
}

void epdu::ssbcb()
{
	double CIalpha;
	std::vector<std::vector<double> > sortedR;

	smoothdirrnd2();
	CIalpha = findCBalpha();
	sortedR = sortCols(r);
	confPrctile(sortedR,CIalpha);
	printResults();
}

void epdu::printResults()
{
	//  for formatting purposes:
	int width = 2+(int)(log10(k));
	int i,j,num,numwidth;

	// write to the console:
	cout << "  bin number |";
	for (i=0; i<numbins; i++){
		cout << setw(width) << i+1;
	}
	cout << endl;
	cout << " ------------|";
	for (i=0; i<numbins; i++){
		for (int j=0; j<width; j++){
			cout << "-";
		}
	}
	cout << endl;
	cout << " lower bound |";
	for (i=0; i<numbins; i++){
		cout << setw(width) << (int)(cl[i]*k);
	}
	cout << endl;
	cout << "    observed |";
	for (i=0; i<numbins; i++){
		cout << setw(width) << (int)a[i];
	}
	cout << endl;
	cout << " upper bound |";
	for (i=0; i<numbins; i++){
		cout << setw(width) << (int)(cu[i]*k);
	}
	cout << endl;

	// write the same stuff to the output file:
	outFile << "  bin number | ";
	for (i=0; i<numbins; i++){
		outFile << i+1;
		for (int j=0; j<(width-1-(int)(log10(i+1))); j++) {
			outFile << " ";
		}
	}
	outFile << "\n";
	outFile << " ------------|";
	for (i=0; i<numbins; i++){
		for (int j=0; j<width; j++){
			outFile << "-";
		}
	}
	outFile << "\n";
	outFile << " lower bound | ";
	for (i=0; i<numbins; i++){
		num = (int)(cl[i]*k);
		if (num==0) {
			numwidth = 0;
		} else {
			numwidth = (int)(log10((int)(num)));
		}
		outFile << num;
		for (j=0; j<(width-1-numwidth); j++) {
			outFile << " ";
		}
	}
	outFile << "\n";
	outFile << "    observed | ";
	for (i=0; i<numbins; i++){
		num = (int)(a[i]);
		if (num==0) {
			numwidth = 0;
		} else {
			numwidth = (int)(log10((int)(num)));
		}
		outFile << num;
		for (j=0; j<(width-1-numwidth); j++) {
			outFile << " ";
		}
	}
	outFile << "\n";
	outFile << " upper bound | ";
	for (i=0; i<numbins; i++){
		num = (int)(cu[i]*k);
		if (num==0) {
			numwidth = 0;
		} else {
			numwidth = (int)(log10((int)(num)));
		}
		outFile << num;
		for (j=0; j<(width-1-numwidth); j++) {
			outFile << " ";
		}
	}
	outFile << "\n";
}

double epdu::findCBalpha()
{
	int numiterations = 10;
	double CaL = 0, CaU = 1, CIalpha, measCBalpha;
	std::vector<std::vector<double> > sortedR;

	sortedR = sortCols(r);
	// in numiterations iterations, this algorithm does a binary search
	// for the best independent confidence level that corresponds to an alpha
	// level simultanous credibility level
	for (int i=0; i<numiterations; i++) {
		for (int j=0; j<B; j++){
			withinCI[j] = 1;
		}
		numFallWithin = B;
		CIalpha = (CaL+CaU)/2;
		confPrctile(sortedR,CIalpha);
		measCBalpha = recCB(r,cu,cl);
		if (measCBalpha > alpha){
			CaU = CIalpha;
		} else {
			CaL = CIalpha;
		}
	}
	return CIalpha;
}

double epdu::recCB(std::vector<std::vector<double> > &matrix, 
				   std::vector<double> &cUpper, std::vector<double> &cLower)
{
	double measCBalpha = 0.05;
	std::vector<std::vector<double> > subMatrix;
	std::vector<double> subCupper, subClower;
	bool fallsWithin;

	// count the number of values in a column that fall inside the
	// simultaneous credibility band given by cUpper and cLower
	for (int i=0; i<B; i++){
		if (withinCI[i]){
			fallsWithin = ((matrix[0][i]>=cLower[0]) && (matrix[0][i]<=cUpper[0]));
			if (!fallsWithin){
				withinCI[i] = 0;
				numFallWithin -= 1;
			}
		}
	}
	if (matrix.size() == 1){
		// calculate the observed simultaneous credibility level
		measCBalpha = 1.0 - (double)numFallWithin/(double)B;
	} else {
		subMatrix = matrix; subMatrix.erase(subMatrix.begin());
		subCupper = cUpper; subCupper.erase(subCupper.begin());
		subClower = cLower; subClower.erase(subClower.begin());
		measCBalpha = recCB(subMatrix,subCupper,subClower);
	}

	return measCBalpha;	
}

void epdu::confPrctile(std::vector<std::vector<double> > &sortedR, double CIalpha)
{
	int index;

	for (int i=0; i<numbins; i++){
		index = (int)(B-1)*CIalpha/2;
		cl[i] = sortedR[i][index];
		index = (int)(B-1)*(1-CIalpha/2);
		cu[i] = sortedR[i][index];
	}
}

std::vector<std::vector<double> > epdu::sortCols(std::vector<std::vector<double> > &matrix)
{
	std::vector<std::vector<double> > sorted;
	std::vector<double> temp;
	for (int i=0; i<numbins; i++){
		temp = matrix[i];
		std::sort(temp.begin(),temp.end());
		sorted.push_back(temp);
	}
	return sorted;
}

void epdu::dirrnd()
{
	double rowSum;
	int i,j;
	RNG x;

	// according to Devroye (1986), generate samples from the Dirichlet posterior 
	// by generating B Gamma distributed numbers for each of the numbins bins.
	for (i=0; i<numbins; i++){
		for (j=0; j<B; j++){
			r[i][j] = x.gamma(a[i]+prior[i],1);
		}
	}
	for (j=0; j<B; j++){
		rowSum = 0;
		for (i=0; i<numbins; i++){
			rowSum += r[i][j];
		}
		for (i=0; i<numbins; i++){
			r[i][j]/=rowSum;
		}
	 }
}

void epdu::smoothdirrnd1()
{
	double roughness, secderiv, randoval, rowSum;
	RNG x;
	int i,j,numbars=0;
	
	for (j=0; j<B; j++){
		while (true) {
			roughness = 0;
			rowSum = 0;
			for (i=0; i<numbins; i++){	
				r[i][j] = x.gamma(a[i]+prior[i],1);
				rowSum+=r[i][j];
			}
			for (i=0; i<numbins; i++){
				r[i][j]/=rowSum;
			}
			for (i=1; i<numbins-1; i++){
				secderiv = r[i-1][j]-2*r[i][j]+r[i+1][j];
				roughness += secderiv*secderiv;
			}
			roughness = exp(-s*roughness);
			randoval = x.uniform();
			if (randoval<roughness) {
				break;
			}
		}
		if (64.0*((double)j/B) > numbars) {
			cout << "|";
			cout << flush;
			numbars++;
		}
	}
	cout << endl;
	cout << endl;
}

void epdu::smoothdirrnd2()
{
	double sumrow;
	int i,j;

	// Smoothed credibility bands are implemented through filtering on 
	// a "sliding window" of 3 bins at a time.
	recSmoothDirRnd(1);
	for (i=0; i<B; i++) {
		sumrow = sumMatrixRow(r,i,1,numbins);
		for (j=0; j<numbins; j++) {
			r[j][i]/=sumrow;
		}
	}
	// the following lines finish up the progress bar.
	int missingbars = 64 - numbins*floor(64/numbins);
	for (i=0; i<missingbars; i++){
		cout << "|";
	}
	cout << flush;
	cout << endl;
	cout << endl;
}

void epdu::recSmoothDirRnd(int j)
{
	if (j>=numbins-1){
		return;
	}
	// add a new column that is smoothed based on the neighbouring 2 bins.
	fixCol(j);
	// update the progress bar.
	for (int i=0; i< floor(64/numbins); i++) {
		if ((j==1)||(j==numbins-2)){
			cout << "|";
		}
		cout << "|";
	}
	cout << flush;
	// move to the next bin.
	recSmoothDirRnd(j+1);
}

void epdu::fixCol(int j)
{
	double rr, rowSum, sumnew, roughness, temp1, temp2, temp3, normcst, randoval, sumfixed = 0;
	RNG x;

	for (int i1=0; i1<B; i1++) {
		if (j!=1) {
			sumfixed = sumMatrixRow(r,i1,0,j+1);
		}
		while (true) {
			sumnew = 0;
			// as in Devroye (1986), generate new columns from a Gamma distribution.
			for (int i2=j-1; i2<numbins; i2++){
				temprow[i2] = x.gamma(a[i2]+prior[i2],1);
				sumnew += temprow[i2];
			}
			if (j!=1) {
				sumnew -= temprow[j-1];
			}
			rowSum = sumfixed + sumnew;
			if (j==1) {
				temp1 = temprow[j-1]/rowSum;
			} else {
				temp1 = r[j-1][i1]/rowSum;
			}
			temp2 = temprow[j]/rowSum;
			temp3 = temprow[j+1]/rowSum;
			normcst = temp1 + temp2 + temp3;
			temp1 /= normcst;
			temp2 /= normcst;
			temp3 /= normcst;
			roughness = temp1 - 2*temp2 + temp3;
			roughness*=roughness;
			rr = exp(-s*roughness);
			// the next line generates a random number from an exponential distribution:
			randoval = x.uniform();
			// this random number "filters" the posterior distribution. The higher s,
			// the more selective the filter.
			if (randoval<rr) {
				if (j==1) {
					r[0][i1] = temprow[0];
					r[1][i1] = temprow[1];
					if (numbins==3) {
						r[2][i1] = temprow[2];
					}
				} else if (j==numbins-2) {
					r[numbins-2][i1] = temprow[numbins-2];
					r[numbins-1][i1] = temprow[numbins-1];
				} else {
					r[j][i1] = temprow[j];
				}
				break;
			}
		}
	}
}

double epdu::sumMatrixRow(std::vector<std::vector<double> > &matrix, int rownum, int begin, int end)
{
	double rowsum = 0;
	for (int i=begin; i<end; i++) {
		rowsum += matrix[i][rownum];
	}
	return rowsum;
}

double epdu::sumVector(std::vector<double> &vec, int begin, int end)
{
	double rowsum = 0;
	for (int i=begin; i<end; i++) {
		rowsum += vec[i];
	}
	return rowsum;
}