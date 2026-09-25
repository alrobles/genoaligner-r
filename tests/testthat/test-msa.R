# genoaligner — MSA contract tests.
#
# The engine's algorithmic truth is verified upstream (bit-exact CPU/GPU
# parity gates in the C++ repo). These tests gate the R-side CONTRACT:
# shape, codon-mode invariants, QC plumbing, refusal behaviour.

test_that("msa_align: rectangular output, one row per input, names applied", {
    a <- msa_align(c(aa = "ACGTACGTAA", bb = "ACGTTCGT", cc = "ACGTACG"))
    expect_s3_class(a, "genoaligner_msa")
    expect_equal(nrow(a$aligned), 3L)
    expect_equal(ncol(a$aligned), a$params$width)
    expect_equal(rownames(a$aligned), c("aa", "bb", "cc"))
    expect_true(all(nchar(apply(a$aligned, 1, paste0, collapse = "")) ==
                    a$params$width))
})

test_that("msa_align: deterministic — identical input, identical output", {
    s <- c("ACGTACGTAA", "ACGTTCGT", "ACGTACG", "ACGTACGTACGT")
    expect_identical(as.matrix(msa_align(s)), as.matrix(msa_align(s)))
})

test_that("msa_align: a single sequence comes back unchanged", {
    a <- msa_align("ACGTACGT")
    expect_equal(paste0(a$aligned[1, ], collapse = ""), "ACGTACGT")
})

test_that("msa_align: protein mode runs and reports no codon QC", {
    a <- msa_align(c("MALWMRLLPLL", "MALWTRLLPLL"), mode = "protein")
    expect_s3_class(a, "genoaligner_msa")
    expect_true(all(is.na(a$qc$frame)))
    expect_true(all(is.na(a$qc$n_stops)))
})

test_that("msa_codon: whole-codon output keeps every row a multiple of 3", {
    s <- c("ATGATAATCACC", "ATGATTATCACCTGA", "ATGATAATCACCTAA")
    a <- msa_codon(s, gc = 2)
    expect_equal(a$params$width %% 3L, 0L)
    for (i in seq_len(nrow(a$aligned))) {
        row <- paste0(a$aligned[i, ], collapse = "")
        expect_equal(nchar(row) %% 3L, 0L)
        # indels are whole codons: every run of '-' is a multiple of 3
        runs <- rle(strsplit(row, "", fixed = TRUE)[[1]])
        gaplens <- runs$lengths[runs$values == "-"]
        expect_true(all(gaplens %% 3L == 0L),
                    info = paste0("row ", i, " gap run: ",
                                  paste(gaplens, collapse = ",")))
    }
})

test_that("msa_codon: QC reports clean coding sequences as stop-free", {
    # gc=2: TGA codes Trp (mitochondrial), so these are stop-free in frame 0
    s <- c("ATGATAATCACC", "ATGATTATCACCTGA")
    a <- msa_codon(s, gc = 2)
    expect_equal(nrow(a$qc), 2L)
    expect_true(all(a$qc$frame %in% 0:2))
    expect_equal(a$qc$n_stops, c(0L, 0L))
    expect_equal(a$qc$partial, c(0L, 0L))
})

test_that("msa_codon: an unavoidable stop under the chosen code is counted", {
    # "TAAGTAAGTAA" (gc=1) carries a stop codon in ALL THREE forward frames
    # (TAA at 0-2 for f0, at 4-6 for f1, at 8-10 for f2), so the fewest-stops
    # frame still reports n_stops >= 1 -- the QC cannot dodge it.
    a <- msa_codon("TAAGTAAGTAA", gc = 1)
    expect_gte(a$qc$n_stops[1], 1L)
})

test_that("msa_codon: a trailing partial codon is reported in qc$partial", {
    a <- msa_codon("ATGAAACCCG")   # 10 nt: frame 0 leaves a 1-nt partial
    expect_gte(a$qc$partial[1], 1L)
})

test_that("msa_codon: ambiguous bases (N) are accepted input", {
    a <- msa_codon(c("ATGNNNATCACC", "ATGATTATCACCTGA"), gc = 2)
    expect_s3_class(a, "genoaligner_msa")
})

test_that("msa_codon: refine=TRUE emits a frameshift-aware NT alignment", {
    s <- c("ATGATAATCACCTGA", "ATGATTACACCTGA")
    a <- msa_codon(s, gc = 2, refine = TRUE)
    expect_s3_class(a, "genoaligner_msa")
    expect_equal(nrow(a$aligned), 2L)
    # rows still consume exactly their input nts
    for (i in 1:2) {
        nts <- paste0(a$aligned[i, ][a$aligned[i, ] != "-"], collapse = "")
        expect_equal(nts, s[i])
    }
})

test_that("msa_*: bad input is refused, never silently served", {
    expect_error(msa_align(character(0)), "non-empty")
    expect_error(msa_align(c("ACGT", NA)), "NA")
    expect_error(msa_align(c("ACGT", "ACG-T")), "alphabet|sequence")  # pre-gapped
    expect_error(msa_align(c("ACGT", "TTTT"), names = "one"), "names")
    expect_error(msa_align("ACGT", mode = "codon"), "should be")
    expect_error(msa_codon("ACGT", gc = 3), "gc must")
    expect_error(msa_align("ACGT", guide = "mafft"), "nj")
})

test_that("msa S3 methods: print/as.matrix/as.character/as_dnabin", {
    a <- msa_align(c(x = "ACGTACGT", y = "ACGTTCGT"))
    expect_output(print(a), "genoaligner MSA: 2 sequences")
    expect_equal(as.matrix(a), a$aligned)
    fa <- as.character(a)
    expect_match(fa[1], "^>x\nACGT", perl = TRUE)
    skip_if_not_installed("ape")
    expect_s3_class(as_dnabin(a), "DNAbin")
})

test_that("align() is deprecated but answers identically to align_edit()", {
    expect_warning(r_old <- align("ACGTACGT", "ACGTTCGT", smax = 8),
                   "deprecated")
    r_new <- align_edit("ACGTACGT", "ACGTTCGT", smax = 8)
    expect_identical(r_old, r_new)
})
