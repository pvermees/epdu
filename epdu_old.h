#ifndef _epdu_h
#define _epdu_h

/* The epdu class was written by P.Vermeesch, last updated on 7-3-2004.
 * It can be downloaded from http://pangea.stanford.edu/research/noble/epdu
 * The program relies on the Newran02B random number generator library
 * written by R.B. Davies and freely available on the Web at
 * http://www.robertnz.net.
 */

class epdu {

public:

	/* constructor takes one argument: the name of the input file
	 */
	epdu(std::string inputFile);

	/* destructor. No dynamic memory is allocated, so this is an empty box.
	 */
	~epdu();

	/* isSmooth returns false if s=0 in input.txt,
	 *          returns true if s>0
	 */
	bool isSmooth(void);

	/* sbcb calculates Simultaneous Bayesian Credibility Bands
	 * for non-smoothed histograms.
	 */
	void sbcb(void);

	/* ssbcb calculates Smooth Simultaneous Bayesian Credibility Bands
	 */
	void ssbcb(void);

private:

	/* dirrnd generates a matrix r of size [B x numbins] with random
	 * samples from a Dirichlet distribution with parameters a + prior
	 * using the method described in Devroye (1986)
	 */
	void dirrnd(void);

	/* smoothdirrnd is similar to dirrnd, but with smoothing.
	 * This program is slower than dirrnd, especially for large 
	 * smoothing parameter s
	 */
	void smoothdirrnd(void);

	/* recSmoothDirRnd is a helper function for smoothdirrnd, recursively
	 * adds new columns to the matrix r of random samples of the smoothed posterior
	 */
	void recSmoothDirRnd(int j);
	
	/* findCBalpha is a helper function for dirrnd and smoothdirrnd.
	 * It find the independent confidence levels for each column of r
	 * that yield the desired simultaneous credibility levels of the
	 * entire matrix.
	 */
	double findCBalpha(void);

	/* fixCol is a helper function for the recSmoothDirRnd function.
	 * Performs the bulk of the work of the smoothdirrnd function.
	 */
	void fixCol(int j);

	/* sortCols sorts the r matrix column by column, so that the percentiles
	 * can be calculated. Returns a matrix of the same size as r.
	 */
	std::vector<std::vector<Real> > sortCols(std::vector<std::vector<Real> > &matrix);

	/* uses sortedR, which is the output from the sortCols function,
	 * and calculates the CIalpha percentiles of it.
	 */
	void confPrctile(std::vector<std::vector<Real> > &sortedR, double CIalpha);

	/* recCB is recursive helper function to findCBalpha that computes the
	 * simultaneous credibility level of matrix that corresponds to the independent
	 * confidence intervals given by the vectors cUpper and cLower.
	 */
	double recCB(std::vector<std::vector<Real> > &matrix, 
				   std::vector<Real> &cUpper, std::vector<Real> &cLower);

	/* sumMatrixRow computes the sum of the rownum-th row of matrix (which is a vector of
	 * column vectors) from begin until (but not including) end.
	 */
	Real sumMatrixRow(std::vector<std::vector<Real> > &matrix, int rownum, int begin, int end);

	/* sumVector sums all the elements of a vector of Real variables.
	 */
	Real sumVector(std::vector<Real> &vec, int begin, int end);

	/* pritResults prints and formats the results.
	 */
	void printResults();

	/* private data of the vector of Reals type:
	 * a = the observed bin counts (size [1 x numbins]).
	 * prior = prior distribution, default = uniform (size [1 x m])
	 * cu = upper credibility bound (size [1 x numbins]).
	 * cl = lower credibility bound (size [1 x numbins]).
	 * temprow = temporary vector of size [1 x numbins].
	 * tempcol = temporary vector of size [B x 1].
	 */
	std::vector<Real> a, prior, cu, cl, temprow, tempcol;

	/* withinCI is a vector of values that indicate whether a column value
	 * falls within (1) or outside (0) of a (independent) confidence interval.
	 */
	std::vector<int> withinCI;

	/* alpha = the simultaneous credibility level, default = 0.05,
	 * s = smoothing parameter: s>=0
	 */
	double alpha, s;

	/* numbins = the length of the a vector
	 * B = the number of samples that have to be taken of the posterior distribution
	 * numFallWithin tallies the number of values in a column that fall within the confidence
	 * interval, this is the sum of the vector returned by the withinCI function.
	 * k = the number of "grains", equals the sum of the values in the vector a.
	 */
	int numbins, B, numFallWithin, k;

	/* the matrix r contains the samples from the posterior distribution, size [B x numbins].
	 */
	std::vector<std::vector<Real> > r;

	/* output file, named "output.txt".
	 */
	std::ofstream outFile;

};

#endif