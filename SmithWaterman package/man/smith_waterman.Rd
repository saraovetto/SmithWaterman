% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\name{smith_waterman}
\alias{smith_waterman}
\title{Local Sequence Alignment}
\usage{
smith_waterman(sequenceA, sequenceB, match = 2, mismatch = -1, gap = -1)
}
\arguments{
\item{sequenceA}{First nucleotide sequence. Can be a string, DNAString or DNAStringSet object.}

\item{sequenceB}{Second nucleotide sequence. Can be a string, DNAString or DNAStringSet object.}

\item{match}{A positive number associated with a nucleotide match. Default: 2}

\item{mismatch}{A negative number associated with a nucleotide mismatch. Default: -1}

\item{gap}{A negative number associated with the gap penalty. Default: -1}
}
\value{
A list containing the optimal alignment (DNAString object), the resulting score and a matrix with the start and end positions.
If the two sequences are shorter than 20 nucleotides, the scoring matrix and the traceback matrix are added to the list.
}
\description{
This function performs pairwise local sequence alignment of two nucleotide sequences
given as input. It uses a scoring system and the Smith Waterman greedy algorithm to
compute the scoring and traceback matrices and returns one of the possible optimal alignments.
}
\examples{
library(Biostrings)
smith_waterman('GCATGCG', DNAString('GATTACA'))
}
\references{
\url{https://en.wikipedia.org/wiki/Smith–Waterman_algorithm}\cr
}
\seealso{
\code{\link{pairwise_score}}\cr
\code{\link{fill_matrices}}\cr
\code{\link{best_alignment}}\cr
}
\author{
Sara Ometto\cr Politecnico di Milano\cr Maintainer: Sara
Ometto\cr E-Mail: <sara.ometto@mail.polimi.it> or <sara.ometto@studenti.unimi.it>
}
