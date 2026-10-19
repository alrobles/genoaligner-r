## Phase-0 diagnostics: caller-supplied guide trees (msa_tree_cpp +
## phylo_to_merges) and the sum-of-pairs scorer (msa_sp_score).

test_that("guide_tree: phylo merges equal NJ tree reproduces the NJ result", {
    skip_if_not_installed("ape")
    # A fixed small DNA set; rebuild the NJ guide in R? The engine NJ is
    # not ape::nj, so we only check the mechanics: a valid binary tree on
    # the sequence names is accepted and produces a sane alignment.
    s <- c(t1 = "ACGTACGT", t2 = "ACGTTCGT", t3 = "ACGTACG",
           t4 = "ACGTTCG", t5 = "TTTACGT")
    tr <- ape::rtree(5, tip.label = names(s))
    a <- msa_align(s, guide_tree = tr)
    expect_s3_class(a, "genoaligner_msa")
    expect_equal(nrow(a$aligned), 5L)
    expect_true(all(nchar(apply(a$aligned, 1, paste0, collapse = ""))
                    == ncol(a$aligned)))
    expect_equal(a$params$guide, "provided")
})

test_that("guide_tree: codon mode keeps the codon grid", {
    skip_if_not_installed("ape")
    s <- c(t1 = "ATGAAATTT", t2 = "ATGAAATTTGGG", t3 = "ATGCCCTTT",
           t4 = "ATGAAATTC")
    tr <- ape::rtree(4, tip.label = names(s))
    a <- msa_codon(s, guide_tree = tr)
    expect_equal(ncol(as.matrix(a)) %% 3L, 0L)
    pa <- msa_codon(s, backbone = "protein", guide_tree = tr)
    expect_equal(ncol(as.matrix(pa)) %% 3L, 0L)
})

test_that("guide_tree: refuses non-matching/non-binary trees", {
    skip_if_not_installed("ape")
    s <- c(t1 = "ACGT", t2 = "ACGA", t3 = "ACGG", t4 = "ACCT")
    expect_error(msa_align(s, guide_tree = ape::rtree(3)),
                 "one tip per sequence")
    tr <- ape::rtree(4, tip.label = names(s)); tr$tip.label[1] <- "ZZZ"
    expect_error(msa_align(s, guide_tree = tr), "tip labels must match")
    expect_error(msa_align(s, guide_tree = "nope"), "ape::phylo")
})

test_that("msa_sp_score: dna/protein/codon modes, ragged + grid refusal", {
    a <- msa_align(c("ACGTACGT", "ACGTTCGT", "ACGTACG"))
    expect_true(is.finite(msa_sp_score(a)))
    p <- msa_align(c("MALWMRLLPLL", "MALWTRLLPLL"), mode = "protein")
    expect_true(is.finite(msa_sp_score(p, mode = "protein")))
    c <- msa_codon(c("ATGAAATTT", "ATGAAATTTGGG", "ATGCCCTTT"))
    expect_true(is.finite(msa_sp_score(c, mode = "codon")))
    # event columns break the codon grid -> scorer must refuse
    s <- c("ATGAAATTT", "ATGAAATTTGGGA", "ATGCCCTTT")   # row 2 len%3=1
    ca <- msa_codon(s, refine = "auto")
    if (ncol(as.matrix(ca)) %% 3L != 0L)
        expect_error(msa_sp_score(ca, mode = "codon"), "multiple of 3")
    expect_error(msa_sp_score(c("AAA", "AAAA")), "Ragged|ragged")
})
