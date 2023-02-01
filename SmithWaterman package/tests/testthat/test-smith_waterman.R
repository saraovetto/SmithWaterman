# Input sequences
test_that("Input error: sequence 1 not a string", {
  seq_a <- 3572
  seq_b <- 'ACGCGTT'
  expect_error(smith_waterman(seq_a, seq_b))
  })


test_that("Input error: sequence 2 not a string", {
  seq_a <- DNAString('GCATGCG')
  seq_b <- list(1, matrix(2:4, nrow=1, ncol=3))
  expect_error(smith_waterman(seq_a, seq_b))
  })



# Scoring parameters
test_that("Input error: match score not a number", {
  seq_a <- DNAString('GCATGCG')
  seq_b <- DNAString('GATTACA')
  match_score <- 'five'
  mismatch_score <- 3
  gap_penalty <- 0
  expect_error(smith_waterman(seq_a, seq_b, match_score, mismatch_score, gap_penalty))
})


test_that("Input error: mismatch score not a number", {
  seq_a <- DNAString('GCATGCG')
  seq_b <- DNAString('GATTACA')
  match_score <- 5
  mismatch_score <- FALSE
  gap_penalty <- 0
  expect_error(smith_waterman(seq_a, seq_b, match_score, mismatch_score, gap_penalty))
})


test_that("Input error: gap penalty not a number", {
  seq_a <- DNAString('GCATGCG')
  seq_b <- DNAString('GATTACA')
  match_score <- 5
  mismatch_score <- 3
  gap_penalty <- 'none'
  expect_error(smith_waterman(seq_a, seq_b, match_score, mismatch_score, gap_penalty))
})



# Function output
test_that("Function results: aligned sequences", {
  seq_a <- 'GCATGCG'
  seq_b <- DNAString('GATTACA')
  # for match, mismatch and gap use the default options

  expected_alignment <- DNAStringSet(c('GCA-TGC', 'G-ATTAC'))
  names(expected_alignment) <- c("Sequence A", "Sequence B")

  expect_equal(smith_waterman(seq_a, seq_b)$optimal_alignment, expected_alignment)
})


test_that("Function results: alignment score", {
  seq_a <- 'GCATGCG'
  seq_b <- DNAString('GATTACA')

  expected_score = paste("Alignment score: ", 5)
  expect_equal(smith_waterman(seq_a, seq_b)$score, expected_score)
})


test_that("Function results: alignment start and end", {
  seq_a <- 'GCATGCG'
  seq_b <- DNAString('GATTACA')

  expected_positions <- matrix(c(1,1,7,7), nrow=2, ncol=2)
  rownames(expected_positions) <- c("Sequence A", "Sequence B")
  colnames(expected_positions) <-c("start", "end")

  expect_equal(smith_waterman(seq_a, seq_b)$positions, expected_positions)
})


test_that("Function results: scoring matrix", {
  seq_a <- 'GCATG'
  seq_b <- DNAString('GATTA')

  expected_scores <- matrix(c(0,0,0,0,0,0,
                              0,2,1,0,0,0,
                              0,1,1,0,0,0,
                              0,0,3,2,1,2,
                              0,0,2,5,4,3,
                              0,2,1,4,4,3), nrow=6, ncol=6, byrow=TRUE)

  dimnames(expected_scores)<- list(c("-", "G", "C", "A", "T", "G"),
                                   c("-", "G", "A", "T", "T", "A"))

  expect_equal(smith_waterman(seq_a, seq_b)$scoring_matrix,
               expected_scores)
})


test_that("Function results: directions matrix", {
  seq_a <- 'GCATG'
  seq_b <- DNAString('GATTA')

  expected_directions <- matrix(c("0", "0", "0", "0", "0", "0",
                                  "0", "\\", "=", "=", "0", "0",
                                  "0", "||", "\\", "\\", "0", "0",
                                  "0", "||", "\\", "=", "=", "\\",
                                  "0", "0", "||", "\\", "\\", "=",
                                  "0", "\\", "||", "||", "\\", "\\"),
                                nrow=6, ncol=6, byrow=TRUE)

  dimnames(expected_directions)<- list(c("-", "G", "C", "A", "T", "G"),
                                   c("-", "G", "A", "T", "T", "A"))

  expect_equal(smith_waterman(seq_a, seq_b)$directions_matrix,
               expected_directions)
})
