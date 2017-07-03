#define WANT_STREAM
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

#ifdef use_namespace
using namespace KW_RNG;
#endif

void preamble();

int main()
{
   time_lapse tl;      // measure program run time

   {
	  std::ofstream outFile;
	  std::string fileName = "input.txt";
	  epdu *myEpdu;

	  preamble();
	  myEpdu = new epdu(fileName);
	  if (myEpdu->isSmooth()){
		  myEpdu->ssbcb();
	  } else {
		  myEpdu->sbcb();
	  }

	  delete myEpdu;
   }

   return 0;
}

void preamble()
{
	cout << "an Estimator of Probability Density and its Uncertainties (EPDU)" << endl;
	cout << "----------------------------------------------------------------" << endl;
	cout << "Version 1.2, 8-16-2004, Pieter Vermeesch" << endl;
	cout << endl;
}

//************** elapsed time class ****************

time_lapse::time_lapse()
{
   start_time = ((double)clock())/(double)CLOCKS_PER_SEC;
}

time_lapse::~time_lapse()
{
   char dummy[2];
   double time = ((double)clock())/(double)CLOCKS_PER_SEC - start_time;
   cout << endl;
   cout << "Elapsed (processor) time = " << setw(10) << setprecision(4) << time << " seconds" << endl;
   cout << endl;
   cout << "Press [enter] to quit: ";
   cin.getline (dummy,2);
}




