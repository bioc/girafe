
using namespace std;

// the above stuff is pure C++, thus needs to be outside the 'extern "C"'


#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Error.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Utils.h> 
 
  /* MAIN PART */

extern "C" {


 /*  coverage for AlignedGenomeIntervalsObjects 
     Written by Joern Toedling, 2010.
  */

  SEXP
  girafe_coverage(SEXP starts, SEXP ends, SEXP val, SEXP cl) {
    /* starts and ends must be passed sorted in ascending order
       of starts, which is ensured in R; assumed to be on the same chr */
    int i, nprotect = 0;
    // i: current position
    int runlen = 0;
    int currentval = 0;
    // length of current run

    SEXP res;

    int * s = INTEGER(starts);
    int * e = INTEGER(ends);
    int chrend = INTEGER(cl)[0];
    
    int posstart =0;
    int posend   =0;

    int * v = INTEGER(val);

    int len = length(starts);
    // upper border for number of expected runs:
    int maxlen = 2 * len + 1;
 
    /* we keep two lists:
     xx growing list of summed up values
     yy growing list of run lengths 
    list<int> xx;
    list<int>::iterator xit, tmpxit;
    list<int> yy;
    list<int>::iterator yit, tmpyit; */
    // use arrays instead;
    int xx [maxlen];
    int yy [maxlen];
    int xxit =0;

    // initialize
    for (i = 0; i < maxlen; i++ ){
      xx[i] = 0;
      yy[i] = 0;
    }

    // go along chromosome
    for (i = 1; i <= chrend; i++ ) { // i is current position
      
      // if i == posstart: add current value, run length
      if ( (i == s[posstart]) || (i >  e[posend]) ){
	// write out current run
	xx[xxit] = currentval;
	yy[xxit] = runlen;	
	xxit++;
	runlen = 0;
	// drop the reads that are finished at this position
	while ( posend < len && i > e[posend] ) {
	  currentval = currentval - v[posend];
	  posend++;
	}
	//  are we done?
	if (posend >= len){
	  runlen = chrend - i + 1;
	  break; // leave loop 
	}
	// add reads that start at this position
	while ( posstart < len  && i == s[posstart] ){
	  currentval = currentval+v[posstart];
	  posstart++;
	}
      }
      runlen++;
      // allow break with STRG+C
      R_CheckUserInterrupt();
    }
    
    //add last run of zeros:
    xx[xxit] = currentval;
    yy[xxit] = runlen;
    xxit++;
    // now 'xxit' contains the real rumber of observed runs

    // prepare results:
    PROTECT( res = allocVector( INTSXP, xxit * 2) );
    nprotect++;
    // fill in result: 
    for (i = 0; i < xxit; i++ ) {
      INTEGER(res)[i]      = xx[i]; //values
      INTEGER(res)[xxit+i] = yy[i]; //run-lengths
    }
    /* set the dimensions of the result: it's an array with number of
     *  runs for nrows and two columns */
    SEXP dim;
    PROTECT( dim = allocVector(INTSXP, 2));
    nprotect++;
    INTEGER(dim)[0] = xxit;
    INTEGER(dim)[1] = 2;
    SET_DIM(res, dim);

    // finish:
    UNPROTECT(nprotect);
    return res;
  } // coverage

}// extern C
