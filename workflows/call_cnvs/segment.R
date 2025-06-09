library(readr)
library(dplyr)
library(HMMcopy)

set.seed(0)

rfile <- commandArgs(TRUE)[1] # coverage
gfile <- commandArgs(TRUE)[2] # gc_wig
mfile <- commandArgs(TRUE)[3] # map_wig
sampleid <- commandArgs(TRUE)[4]
lower_thresh <- 0
mappability_thresh <- 0.9

if (is.na(rfile) | !file.exists(rfile)) {
    stop("Read coverage wig formatted should be provided.\nRscript segment.R <coverage wig> <gc wig> <mappability wig> <sample ID> <coverage min(default:100)>")
}
if (is.na(gfile) | !file.exists(gfile)) {
    stop("GC Percent wig file should be formatted.\nRscript segment.R <coverage wig> <gc wig> <mappability wig> <sample ID> <coverage min(default:100)>")
}
if (is.na(mfile) | !file.exists(mfile)) {
    stop("Mappability wig file should be provided.\nRscript segment.R <coverage wig> <gc wig> <mappability wig> <sample ID> <coverage min(default:100)>")
}
if (is.na(sampleid)) {
    stop("Sample ID should be provided.\nRscript segment.R <coverage wig> <gc wig> <mappability wig> <sample ID> <coverage min(default:100)>")
}

# Modified from HMMcopy correctReadCounts function, since it was
# overcorrecting for mappability Now the function only corrects for
# GC
mycorrect <- function(x, mappability = 0.9, samplesize = 50000, verbose = TRUE) {
    if (length(x$reads) == 0 | length(x$gc) == 0 | length(x$map) == 0) {
        stop("Missing one of required columns: reads, gc, map")
    }

    if (verbose) {
        message("Applying filter on data...")
    }

    x$valid <- TRUE
    x$valid[x$reads <= 0 | x$gc < 0] <- FALSE
    x$ideal <- TRUE
    routlier <- 0.01
    range <- quantile(x$reads[x$valid], prob = c(0, 1 - routlier), na.rm = TRUE)
    dlower <- 1e-04
    dhigher <- 0.99999
    domain <- quantile(x$gc[x$valid], prob = c(dlower, dhigher), na.rm = TRUE)
    x$ideal[!x$valid | x$map < mappability | x$reads <= range[1] | x$reads >
        range[2] | x$gc < domain[1] | x$gc > domain[2]] <- FALSE

    if (verbose) {
        message("Correcting for GC bias...")
    }

    set <- which(x$ideal)
    select <- sample(set, min(length(set), samplesize))
    rough <- loess(x$reads[select] ~ x$gc[select], span = 0.03)
    i <- seq(0, 1, by = 0.01)
    final <- loess(predict(rough, i) ~ i, span = 0.3)
    x$cor.gc <- x$reads/predict(final, x$gc)

    if (verbose) {
        message("Correcting for mappability bias...")
    }

    coutlier <- 0.01
    range <- quantile(x$cor.gc[which(x$valid)], prob = c(0, 1 - coutlier),
        na.rm = TRUE)
    set <- which(x$cor.gc < range[2])
    select <- sample(set, min(length(set), samplesize))
    final <- approxfun(lowess(x$map[select], x$cor.gc[select]))
    x$cor.map <- x$cor.gc/final(x$map)
    x$copy <- x$cor.map
    x$copy[x$copy <= 0] <- 1e-06
    x$copy <- log(x$copy, 2)
    return(x)
}


mappability <- read_tsv(mfile, col_names = c("chr","start","end","map"))
gc_content <- read_tsv(gfile, col_names = c("chr","start","end","gc"))
read_cov <- read_tsv(rfile, col_names = c("chr", "start", "end","reads"))

uncorrected_reads <- read_cov %>%
  inner_join(gc_content, by = c("chr", "start","end")) %>%
  inner_join(mappability, by = c("chr", "start","end"))
corrected_copy <- mycorrect(uncorrected_reads)

param <- data.frame(strength = 1e+30, e = 0.99999999, mu = c(-2, -1, 0,
    1, 2), lambda = 2, nu = 2.1, kappa = 75, m = c(-2, -1, 0, 1, 2), eta = 0.3,
    gamma = 3, S = 0)

corrected_copy[corrected_copy$reads < lower_thresh |
                                   corrected_copy$map < mappability_thresh, "copy"] <- NA

segmented_copy <- HMMsegment(corrected_copy, maxiter = 500, param = param,
    autosomes = rep(TRUE, nrow(corrected_copy)))

write_tsv(segmented_copy$segs, paste0(dirname(rfile), "/", sampleid, ".cn_segments.tsv"))
write_tsv(corrected_copy, paste0(dirname(rfile), "/", sampleid, ".read_cov_bin.tsv"))
write_tsv(param, paste0(dirname(rfile), "/", sampleid, ".input_params.tsv"))
