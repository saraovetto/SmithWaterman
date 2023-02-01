# Smith - Waterman Algorithm

#' Calculate Pairwise Score
#'
#' This function assigns to a pair of nucleotides a score, whether the two nucleotides match or not.
#' The value depends on the parameters match and mismatch, which are part of the scoring system decided a priori.
#'
#' @param nt1 Nucleotide of the first sequence, a character.
#' @param nt2 Nucletide of the second sequence, a character.
#' @param match A positive number associated with a nucleotide match.
#' @param mismatch A negative number associated with a nucleotide mismatch.
#' @return An integer representing the score associated to the comparison of the two nucleotides.

pairwise_score <- function(nt1, nt2, match, mismatch) {
  if (nt1 == nt2) return(match)
  else return(mismatch)
}


#' Fill the Scoring Matrix and Traceback Matrix
#'
#' For each pair of nucleotides, the function assigns a score using the Smith Waterman greedy algorithm.
#' It inserts the value in the scoring matrix and, depending on the score,
#' it stores the direction of the alignment in the traceback matrix.
#'
#' @param scoringMat Initialized empty Scoring Matrix of integers, with the first row and first column set to 0.
#' @param tracebackMat Initialized empty Traceback Matrix of characters, with the first row and column set to 0.
#' @param seqA A vector of the first nucleotide sequence with "-" as first element.
#' @param seqB A vector of the second nucleotide sequence with "-" as first element.
#' @param match A positive number associated with a nucleotide match.
#' @param mismatch A negative number associated with a nucleotide mismatch.
#' @param gap A negative number associated with the gap penalty.
#' @return A list with the filled Scoring Matrix and the filled Traceback Matrix.

fill_matrices <- function(scoringMat, tracebackMat, seqA, seqB, match, mismatch, gap) {

  #sequence length
  n <- length(seqA)
  m <- length(seqB)

  for (i in 2:n) {
    for (j in 2:m) {
      # compute the scores for each possible case
      paired <- scoringMat[i-1, j-1] + pairwise_score(seqA[i], seqB[j], match, mismatch)
      deletion <- scoringMat[i-1, j] + gap
      insertion <- scoringMat[i, j-1] + gap
      scoringMat[i,j] <- max(paired, deletion, insertion, 0)

      # save directions
      if (scoringMat[i,j] == paired) {
        tracebackMat[i,j] <- "\\"
      } else if (scoringMat[i,j] == deletion) {
        tracebackMat[i,j] <- "||" # vertical - gap in seq 2
      } else if (scoringMat[i,j] == insertion) {
        tracebackMat[i,j] <- "=" # horizontal - gap in seq 1
      } else {
        tracebackMat[i,j] <- 0
      }

    }
  }
  return(list(scoring_matrix = scoringMat, directions_matrix = tracebackMat))
  # the two filled matrices are needed for the traceback process
}


#' Find the best alignment
#'
#' Starting from the element with the highest score, this function performs the traceback
#' based on the source of each score recursively, until 0 is encountered.
#' The function returns the optimal local alignment between the two sequences and its score.
#'
#' @param scoringMat filled Scoring Matrix obtained from \code{\link{fill_matrices}} applying the Smith-Waterman algorithm.
#' @param tracebackMat filled Traceback Matrix obtained from \code{\link{fill_matrices}} applying the Smith-Waterman algorithm.
#' @return A list of the two aligned sequences, the score of the optimal alignment and a matrix with starting ending position.

best_alignment <- function(scoringMat, tracebackMat) {

  # the traceback starts from the cell with maximum score
  # find the highest value of alignment
  score <- max(scoringMat)
  # retrieve the indices
  coordinates <- which(scoringMat == score, arr.ind=TRUE)
  # when there is more than one cell with the maximum score, it takes the nearest to the bottom right corner
  i <- coordinates[nrow(coordinates),1]
  j <- coordinates[nrow(coordinates),2]
  endPos <- c(i,j)

  # initialize the strings containing the alignments
  align1 <- ""
  align2 <- ""

  # the traceback stops when the matrix cell has 0 score
  while (tracebackMat[i,j] != "0") {
    # two nucleotides aligned - diagonal move
    if (tracebackMat[i,j] == "\\") {
      align1 <- paste0(rownames(scoringMat)[i], align1)
      align2 <- paste0(colnames(scoringMat)[j], align2)
      # update indices
      i <- i-1
      j <- j-1
    }
    # gap in seq 2 - vertical move
    else if (tracebackMat[i,j] == "||") {
      align1 <- paste0(rownames(scoringMat)[i], align1)
      align2 <- paste0("-", align2)
      i <- i-1 # update index of seqA
    }
    # gap in seq 1 - horizontal move
    else if (tracebackMat[i,j] == "=") {
      align1 <- paste0("-", align1)
      align2 <- paste0(colnames(scoringMat)[j], align2)
      j <- j-1 # update index of seqB
    }
  }

  startPos <- c(i,j)
  # create a matrix with the starting and ending coordinates of the alignment
  positions <- matrix(c(startPos,endPos),nrow=2,ncol=2)
  rownames(positions) <- c("Sequence A", "Sequence B")
  colnames(positions) <-c("start", "end")

  # return the two aligned sequences with the relative score and the positions
  return(list(align1, align2, score, positions))

}





#' Local Sequence Alignment
#'
#' This function performs pairwise local sequence alignment of two nucleotide sequences
#' given as input. It uses a scoring system and the Smith Waterman greedy algorithm to
#' compute the scoring and traceback matrices and returns one of the possible optimal alignments.
#'
#' @param sequenceA First nucleotide sequence. Can be a string, DNAString or DNAStringSet object.
#' @param sequenceB Second nucleotide sequence. Can be a string, DNAString or DNAStringSet object.
#' @param match A positive number associated with a nucleotide match. Default: 2
#' @param mismatch A negative number associated with a nucleotide mismatch. Default: -1
#' @param gap A negative number associated with the gap penalty. Default: -1
#' @return A list containing the optimal alignment (DNAString object), the resulting score and a matrix with the start and end positions.
#' If the two sequences are shorter than 20 nucleotides, the scoring matrix and the traceback matrix are added to the list.
#' @author Sara Ometto\cr Politecnico di Milano\cr Maintainer: Sara
#' Ometto\cr E-Mail: <sara.ometto@@mail.polimi.it> or <sara.ometto@@studenti.unimi.it>
#' @references \url{https://en.wikipedia.org/wiki/Smith–Waterman_algorithm}\cr
#' @seealso \code{\link{pairwise_score}}\cr
#' \code{\link{fill_matrices}}\cr
#' \code{\link{best_alignment}}\cr
#' @import Biostrings
#' @importFrom "methods" "is"
#' @examples
#' library(Biostrings)
#' smith_waterman('GCATGCG', DNAString('GATTACA'))
#' @export

smith_waterman <- function(sequenceA, sequenceB, match = 2, mismatch = -1, gap = -1) {
  # the default values for gap, mismatch and match are respectively -1, -1 and 2

  ## stop conditions
  # check if the sequences are strings or DNAString objects
  stopifnot("Sequence A must be a string or a DNAString" = (is.character(sequenceA) |
                                                              is(sequenceA,"DNAString") |
                                                              is(sequenceA,"DNAStringSet")))
  stopifnot("Sequence B must be a string or a DNAString" = (is.character(sequenceB) |
                                                              is(sequenceB,"DNAString") |
                                                              is(sequenceB,"DNAStringSet")))
  # make sure that all the scores are numbers
  stopifnot("All the scores must be numeric" = (is.numeric(match) &
                                                  is.numeric(mismatch) &
                                                  is.numeric(gap)))
  # check that mismatch and gap are negative
  if(mismatch > 0) mismatch <- -mismatch
  if(gap > 0) gap <- -gap



  ## prepare the sequences
  # in case of DNAString objects, convert them in strings
  if(is(sequenceA,"DNAString") | is(sequenceA,"DNAStringSet"))
    sequenceA <- toString(sequenceA)
  if(is(sequenceB,"DNAString") | is(sequenceB,"DNAStringSet"))
    sequenceB <- toString(sequenceB)

  # convert the strings in vectors, then add "-" as first character to build the matrix
  seqA <- c("-", unlist(strsplit(sequenceA, '')))
  seqB <- c("-", unlist(strsplit(sequenceB, '')))

  # compute the length of the sequences
  n <- length(seqA)
  m <- length(seqB)



  ## initialize the scoring matrix
  scoringMat <- matrix(nrow = n, ncol = m)
  dimnames(scoringMat) <- list(seqA, seqB)  # set seq nucleotides as matrix names
  # first row and column are set to 0
  scoringMat[,1] <- 0
  scoringMat[1,] <- 0

  ## initialize the traceback matrix to store the directions of the alignment
  tracebackMat <- matrix("0", nrow = n, ncol = m)
  dimnames(tracebackMat) <- list(seqA, seqB)  # set seq nucleotides as matrix names



  ## compute the scores and fill the two matrices
  matrices <- fill_matrices(scoringMat, tracebackMat, seqA, seqB,
                            match, mismatch, gap)



  ## find the best alignment
  results <- best_alignment(scoringMat = matrices[[1]],
                            tracebackMat = matrices[[2]])



  ## prepare output
  # store the final aligned sequences in a DNAString object
  alignment <- DNAStringSet(c(results[[1]], results[[2]]))
  names(alignment) <- c("Sequence A", "Sequence B")

  final_alignment <- list(optimal_alignment = alignment,
                          score = paste("Alignment score: ", results[3]),
                          positions = results[[4]])

  # if the two sequences are shorter than 20 nucleotides, add the matrices to the output
  if(nrow(scoringMat) <= 20 & ncol(scoringMat) <= 20) {
    final_alignment <- append(final_alignment, matrices)}

  return(final_alignment)

}
