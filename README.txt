This  README  file describes  how  to install  and  run  the EPDU  (an
Estimator of Probability Density and its Uncertainties) program, which
calculates simultaneous Bayesian  credibility bands for histograms, as
described  in the  manuscript  submitted to  the  American Journal  of
Science.

To install  EPDU, copy it  into the same  directory as its  input file
"input.txt".  The  binary executable is named  epduWIN.exe for Windows
(version   98SE   or   higher)   and   epduMAC   for   the   Macintosh
(OS-X). "input.txt" should be  formatted as in the following examples,
which  calculate simultaneous  Bayesian  95% credibility  bands for  a
histogram with 10 bins:

5 10 5 10 5 0 6 ;
0 ;
0.05 ;
1 1 1 1 1 1 1 ;

or, equivalently:

5 10 5 10 5 0 6 ; 0 ; 0.05 ; 1 1 1 1 1 1 1 ;

or, equivalently:

5 10 5 10 5 0 6 ;
0 ;
0.05 ;

or, equivalently:

5 10 5 10 5 0 6 ; 0 ; 0.05 ;

or, equivalently:

5 10 5 10 5 0 6 ;
0 ;

The  first line  contains  the 10  histogram  bin counts  (n_j in  the
paper).   They must  be separated  by a  single space  and ended  by a
semi-colon, which itself is also  separated from the last bin count by
a single space.  The second  entry is the smoothing parameter s, which
must be a positive number. If  s=0, no smoothing is applied. If s>0, a
different function  is used  that does allow  smoothing. If  the input
histogram is  very rough and a  large smoothing parameter  is used, it
can  take a long  time before  the program  finishes.  A  progress bar
appears when s>0.   It is a good idea to start  with a small smoothing
factor (s=1) and gradually increase it over several program runs.  The
next entries are optional. They can be placed on separate lines, or on
the same line as the bin counts and/or smoothing parameter.  Following
the bin  counts is the confidence  level.  The default  value is 0.05.
Again, it must be followed by a single space and a semi-colon. Finally
comes  an   optional  row  of   values  corresponding  to   the  prior
distribution.  This row must contain  the same number of values as the
first row (in the above example:  7). The default is a flat prior (all
1's),  but any  other prior  can be  specified as  well.  For example,
Jeffrey's prior for 7 bins would be implemented as:

0.5 0.5 0.5 0.5 0.5 0.5 0.5 ;

If  the third  row  is omitted,  also  the prior  information must  be
omitted.   However,  it is  possible  to  specify  a confidence  level
different from  0.05, and still  omit the prior information,  in which
case the default uniform prior is used.










