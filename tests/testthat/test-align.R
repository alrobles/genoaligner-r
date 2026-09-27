test_that("align: known edit distance + CIGAR", {
    r <- align_edit("ACGTACGT", "ACGTTCGT", smax = 8)
    expect_equal(r$score, 1)
    expect_true(r$resolved)
    expect_equal(as.character(r$cigar), "MMMMXMMM")
    expect_true(r$rescore_ok)
    expect_true(r$wellformed_ok)
})

test_that("align: exact match", {
    r <- align_edit("ACGT", "ACGT", smax = 8)
    expect_equal(r$score, 0)
    expect_true(r$resolved)
    expect_equal(as.character(r$cigar), "MMMM")
})

test_that("align: smax bound leaves a pair unresolved", {
    r <- align_edit("AAAAAAA", "TTTTTTT", smax = 2)   # distance 7 > 2
    expect_false(r$resolved)
    expect_true(is.na(r$score))
    expect_true(is.na(r$cigar))
})

test_that("align: batch + recycling + NA", {
    r <- align_edit(c("AAAACCC", "ACGT", NA), "AAAATCC", smax = 8)
    expect_s3_class(r, "data.frame")
    expect_equal(nrow(r), 3L)
    expect_true(r$resolved[1])
    expect_true(r$resolved[2])
    expect_false(r$resolved[3])
    expect_true(is.na(r$score[3]))
    expect_true(is.na(r$cigar[3]))
})

test_that("align: with_cigar = FALSE omits cigar but keeps score", {
    r <- align_edit("ACGT", "ACCT", smax = 8, with_cigar = FALSE)
    expect_equal(r$score, 1)
    expect_true(is.na(r$cigar))
})

test_that("align: smax outside [0,511] errors", {
    expect_error(align_edit("ACGT", "ACGT", smax = -1))
    expect_error(align_edit("ACGT", "ACGT", smax = 512))
})

test_that("align: oracle parity on random pairs", {
    set.seed(123)
    alpha <- c("A", "C", "G", "T")
    for (i in seq_len(100)) {
        L <- sample(1:15, 1)
        a <- paste0(sample(alpha, L, replace = TRUE), collapse = "")
        b <- paste0(sample(alpha, L, replace = TRUE), collapse = "")
        r <- align_edit(a, b, smax = 30)
        ed <- as.numeric(adist(a, b, costs = c(1, 1, 1)))
        expect_true(r$resolved)
        expect_equal(r$score, ed)
        expect_true(r$rescore_ok)
        expect_true(r$wellformed_ok)
    }
})

test_that("align: empty inputs", {
    r <- align_edit("", "", smax = 8)
    expect_equal(r$score, 0)
    expect_true(r$resolved)
    r2 <- align_edit("ACGT", "", smax = 8)
    expect_equal(r2$score, 4)
    expect_true(r2$resolved)
})