test_that("align_sw: C++ documented example", {
    r <- align_sw("TTTACGTGTT", "ACGTGT", scoring = c(2, -3, 5, 2))
    expect_equal(r$score, 12)
    expect_equal(r$start_i, 3)
    expect_equal(r$start_j, 0)
    expect_equal(r$end_i, 8)
    expect_equal(r$end_j, 5)
    expect_equal(as.character(r$cigar), "MMMMMM")
    expect_true(r$rescore_ok)
    expect_true(r$wellformed_ok)
    expect_false(r$too_large)
})

test_that("align_sw: no positive alignment -> score 0, coord -1", {
    r <- align_sw("AAAA", "TTTT", scoring = c(2, -5, 3, 1))
    expect_equal(r$score, 0)
    expect_equal(r$start_i, -1)
    expect_equal(r$end_i, -1)
    expect_equal(nchar(as.character(r$cigar)), 0L)
})

test_that("align_sw: batch + recycling + NA", {
    r <- align_sw(c("ACGTACGT", "GGGG", NA), "ACGTTCGT", c(2, -3, 5, 2))
    expect_s3_class(r, "data.frame")
    expect_equal(nrow(r), 3L)
    # whole word aligns: 7 matches (14) + 1 mismatch (-3) = 11
    expect_equal(r$score[1], 11)
    expect_equal(r$start_i[1], 0)
    expect_equal(r$end_i[1], 7)
    expect_true(is.na(r$score[3]))
    expect_true(is.na(r$cigar[3]))
})

test_that("align_sw: with_cigar = FALSE omits cigar and keeps score + ends", {
    r <- align_sw("ACGTACGT", "ACGTTCGT", c(2, -3, 5, 2), with_cigar = FALSE)
    expect_true(r$score >= 0)
    expect_true(is.na(r$cigar))
})

test_that("align_sw: scoring must length 4", {
    expect_error(align_sw("ACGT", "ACGT", scoring = c(2, -3)))
})

test_that("align_sw: self-check rescore/wellformed on random pairs", {
    set.seed(42)
    alpha <- c("A", "C", "G", "T")
    for (i in seq_len(100)) {
        L <- sample(2:14, 1)
        a <- paste0(sample(alpha, L, replace = TRUE), collapse = "")
        b <- paste0(sample(alpha, L, replace = TRUE), collapse = "")
        s <- align_sw(a, b, scoring = c(3, -2, 4, 1))
        expect_true(s$score >= 0)
        if (s$score > 0) {
            expect_true(s$rescore_ok)
            expect_true(s$wellformed_ok)
        }
    }
})