\name{Additional reduce-methods}
\docType{methods}
\alias{reduce-methods}
\alias{reduce,CompressedIRangesList-method}
\alias{reduce,IRanges-method}
\alias{reduce,IntegerRanges-method}
\alias{reduce,IntegerRangesList-method}
\title{Auxiliary methods for Function reduce in Package 'girafe'}
\description{
  This methods were written to resurrect the functionality of the
  'reduce' method of package 'IRanges' for objects belonging to classes
  defined in 'IRanges'.
  This functions had been overwritten by the later import of package
  'intervals'. See the corresponding help pages of package IRanges for
  more details on these methods.
}
\section{Methods}{
\describe{
\item{\code{signature(x = "CompressedIRangesList")}}{see help page in
  package \sQuote{IRanges} }
\item{\code{signature(x = "IRanges")}}{see help page in
  package \sQuote{IRanges} }
\item{\code{signature(x = "IntegerRanges")}}{see help page in
  package \sQuote{IRanges} }
\item{\code{signature(x = "IntegerRangesList")}}{see help page in
  package \sQuote{IRanges} }
}}
\seealso{ \code{\link[IRanges]{IntegerRangesList-class}},
  \code{\link[IRanges]{IntegerRanges-class}},
  \code{\link[IRanges]{IRanges-class}}
}
\keyword{methods}
\keyword{internal}
