# genoaligner — synthetic alignment smoke test (runs on a Q6000 GPU node).
# The v1.0.0 CPU core does the alignment; --nv makes the Q6000 visible in the
# container. This validates the controlled environment end to end.
suppressMessages(library(genoaligner))
cat("=== genoaligner synthetic alignment smoke test ===\n")
cat("host      :", system("hostname", intern = TRUE), "\n")
cat("R         :", R.version.string, "\n")
cat("genoaligner:", as.character(packageVersion("genoaligner")), "\n")

# GPU visibility (informational at this stage; becomes load-bearing in stage 2)
cat("\nGPU (via --nv):\n")
cat(tryCatch(paste(system("nvidia-smi --query-gpu=name --format=csv,noheader",
                          intern = TRUE), collapse = "; "),
             error = function(e) "nvidia-smi unavailable"),
    "\n")

# --- Synthetic pairs with known answers -------------------------------------
set.seed(7)
alpha <- c("A", "C", "G", "T")
mk <- function(l) paste0(sample(alpha, l, replace = TRUE), collapse = "")

cat("\n-- edit distance (align) --\n")
ok <- 0L; n <- 20L
for (i in seq_len(n)) {
    l   <- sample(10:40, 1)
    ref <- mk(l)
    # introduce ~10% random substitutions
    s   <- strsplit(ref, "")[[1]]
    k   <- max(1L, floor(l * 0.1))
    s[sample(seq_len(l), k)] <- sample(alpha, k, replace = TRUE)
    qry <- paste0(s, collapse = "")
    r   <- align(qry, ref, smax = 100)
    ed  <- as.numeric(adist(qry, ref, costs = c(1, 1, 1)))
    if (isTRUE(r$resolved) && r$score == ed &&
        isTRUE(r$rescore_ok) && isTRUE(r$wellformed_ok)) ok <- ok + 1L
    if (i <= 3) cat(sprintf("  len=%3d score=%d (expected %d) resolved=%s\n",
                            l, r$score, ed, r$resolved))
}
cat(sprintf("  parity: %d/%d synthetic pairs match the edit-distance oracle\n",
            ok, n))

cat("\n-- Smith-Waterman local (align_sw) --\n")
r <- align_sw("TTTACGTGTT", "ACGTGT", scoring = c(2, -3, 5, 2))
cat(sprintf("  score=%d start=(%d,%d) end=(%d,%d) cigar=%s\n",
            r$score, r$start_i, r$start_j, r$end_i, r$end_j, r$cigar))
stopifnot(r$score == 12L)

cat("\nSMOKE_TEST_OK\n")