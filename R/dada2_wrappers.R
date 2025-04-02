# SPDX-FileCopyrightText: 2024, Brendan Furneaux
# SPDX-License-Identifier: MIT

# Wrappers for dada2 functions

#' Class to store dada2 filter options
#' @param maxEE_R1 Maximum expected errors for forward reads
#' @param maxEE_R2 Maximum expected errors for reverse reads
#' @return A list with class "dada2_filter_options"
#' @export
dada2_filter_options <- function(maxEE_R1, maxEE_R2) {
  checkmate::assert_number(maxEE_R1, lower = 0)
  checkmate::assert_number(maxEE_R2, lower = 0)
  structure(
    c(maxEE_R1 = maxEE_R1, maxEE_R2 = maxEE_R2),
    class = "dada2_filter_options"
  )
}

#' Names of options in dada2_filter_options
#' @export
filter_option_names = c("maxEE_R1", "maxEE_R2")

#' Update method for dada2 filter options
#' @param object ([`dada2_filter_options`][dada2_filter_options()]) existing
#' dada2 filter options object to modify
#' @param new_options (named `list`, single-row `data.frame`, named
#'  `character` vector, or named `numeric` vector) new values for the options.
#'  If a `data.frame`, then the values in each column should be all the same.
#' @param ... Additional arguments (ignored)
#' @exportS3Method stats::update
update.dada2_filter_options <- function(object, new_options, ...) {
  checkmate::assert(
    checkmate::check_list(new_options, null.ok = TRUE),
    checkmate::check_data_frame(new_options, null.ok = TRUE),
    checkmate::check_character(new_options, null.ok = TRUE),
    checkmate::check_numeric(new_options, null.ok = TRUE)
  )
  if ("maxEE_R1" %in% names(new_options)) {
    checkmate::assert_number(new_options[["maxEE_R1"]], lower = 0, finite = TRUE)
    object[["maxEE_R1"]] <- new_options[["maxEE_R1"]]
  }
  if ("maxEE_R2" %in% names(new_options)) {
    checkmate::assert_number(new_options[["maxEE_R2"]], lower = 0, finite = TRUE)
    object[["maxEE_R2"]] <- new_options[["maxEE_R2"]]
  }
  object
}


#' Wrapper for dada2::filterAndTrim
#' @param fwd Forward reads
#' @param filt Forward filtered reads
#' @param rev Reverse reads
#' @param filt.rev Reverse filtered reads
#' @param ... Additional arguments (passed to [dada2::filterAndTrim()])
#' @return `character` vector of files names for output files
#'
#' @details This function is a wrapper for `dada2::filterAndTrim` that ensures
#' that the output files are created even if they end up being empty. It also
#' returns an empty character vector if no input files are given.#'
#' @export
filterAndTrim <- function(fwd, filt, rev, filt.rev, ...) {
  if (!requireNamespace("dada2", quietly = TRUE)) {
    stop("dada2 package must be installed to run filterAndTrim()")
  }
  if (length(fwd) > 0L) {
    # ensure files are created, even if they end up being empty
    file.create(c(filt, filt.rev))
    mycall <- match.call()
    mycall[[1]] <- dada2::filterAndTrim
    eval.parent(mycall)
    # return file names for samples where at least some reads passed
    purrr::keep(c(filt, filt.rev), file.exists)
  } else {
    character()
  }
}


#' Wrapper for dada2::derepFastq
#' @param fls File names
#' @param n Number of reads to process
#' @param verbose Print progress messages
#' @param qualityType Quality score type
#' @param names Sample names for the output object
#' @return A named list of 0 or more derep objects.
#'
#' @details This function is a wrapper for `dada2::derepFastq` that ensures
#' that the output is always a list, even if it contains only one element or
#' none. It also allows the user to specify sample names for the output object.
#' @export
derepFastq <- function(fls, n = 1e+06, verbose = FALSE, qualityType = "Auto",
                       names = fls) {
  if (length(fls) == 0) {
    out <- list()
  } else {
    out <- dada2::derepFastq(fls, n = n, verbose = verbose, qualityType = qualityType)
    if (methods::is(out, "derep")) {
      out <- list(out)
    }
    names(out) <- names
  }
  out
}

#' Error function for binned quality scores
#' @param binnedQ Vector of binned quality scores
#' @return A function that estimates error rates based on binned quality scores
#' @details From dada2 commit <https://github.com/benjjneb/dada2/commit/7714487b153ca133cb6f03cb01d09fc05be60159>
#' @export
makeBinnedQualErrfun <- function(binnedQ=c(2, 11, 25, 37)) {
  function(trans, binnedQuals=binnedQ) {
    qq <- as.numeric(colnames(trans))
    # Get min and max observed quality scores
    qmax <- max(qq[colSums(trans)>0])
    qmin <- min(qq[colSums(trans)>0])
    # Check for data consistency with provided binned qualities
    if(qmax > max(binnedQuals)) stop("Input data contains a higher quality score than the provided binned values.")
    if(qmin < min(binnedQuals)) stop("Input data contains a lower quality score than the provided binned values.")
    if(!qmax %in% binnedQuals) warning("Maximum observed quality score is not in the provided binned values.")
    if(!qmin %in% binnedQuals) warning("Minimum observed quality score is not in the provided binned values.")

    est <- matrix(0, nrow=0, ncol=length(qq))
    for(nti in c("A","C","G","T")) {
      for(ntj in c("A","C","G","T")) {
        if(nti != ntj) {
          errs <- trans[paste0(nti,"2",ntj),]
          tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
          p <- errs/tot
          df <- data.frame(q=qq, errs=errs, tot=tot, p=p)
          # Check and enforce that this q scores start at zero
          if(!all(df$q == seq(nrow(df))-1)) stop("Unexpected Q score series.") ###!
          pred <- rep(NA, nrow(df))
          for(i in seq(length(binnedQuals)-1)) {
            loQ <- binnedQuals[i]
            hiQ <- binnedQuals[i+1]
            loP <- df$p[loQ+1]
            hiP <- df$p[hiQ+1]
            # Linear interpolation between the binned Q scores observed in the data
            if(!is.na(loP) && !is.na(hiP)) {
              pred[(loQ+1):(hiQ+1)] <- seq(loP, hiP, length.out=(hiQ-loQ+1))
            }
          }

          maxrli <- max(which(!is.na(pred)))
          minrli <- min(which(!is.na(pred)))
          pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
          pred[seq_along(pred)<minrli] <- pred[[minrli]]
          est <- rbind(est, pred)
        } # if(nti != ntj)
      } # for(ntj in c("A","C","G","T"))
    } # for(nti in c("A","C","G","T"))

    # HACKY
    MAX_ERROR_RATE <- 0.25
    MIN_ERROR_RATE <- 1e-7
    est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
    est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE

    # Expand the err matrix with the self-transition probs
    err <- rbind(1-colSums(est[1:3,]), est[1:3,],
                 est[4,], 1-colSums(est[4:6,]), est[5:6,],
                 est[7:8,], 1-colSums(est[7:9,]), est[9,],
                 est[10:12,], 1-colSums(est[10:12,]))
    rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
    colnames(err) <- colnames(trans)
    # Return
    return(err)
  }
}

#' Choose error function for dada2
#' @param fls File names
#' @return A function that estimates error rates based on binned quality scores
#' if the input files are binned, otherwise the default dada2 error function.
#' @export
choose_dada_error_function <- function(fls, ...) {
  bins <- fastq_qual_bins(fls, ...)
  if (length(bins) < 10) {
    makeBinnedQualErrfun(bins)
  } else {
    dada2::loessErrfun
  }
}

#' Wrapper for dada2::learnErrors
#' @param fls File names
#' @param nbases Number of bases to use for error estimation
#' @param nreads Number of reads to use for error estimation
#' @param errorEstimationFunction Function to estimate error rates
#' @param multithread Use multiple threads
#' @param randomize Randomize the order of reads
#' @param MAX_CONSIST Maximum number of consistent bases
#' @param OMEGA_C Omega constant
#' @param qualityType Quality score type
#' @param verbose Print progress messages
#' @param ... Additional arguments (passed to [dada2::learnErrors()])
#' @return An error object or `NULL` if no input files are given
#' @export
learnErrors <- function(fls, nbases = 1e+09, nreads = NULL,
                        errorEstimationFunction = choose_dada_error_function,
                        multithread = FALSE, randomize = FALSE, MAX_CONSIST = 10,
                        OMEGA_C = 0, qualityType = "Auto", verbose = FALSE, ...) {
  if (length(fls) == 0) {
    NULL
  } else {
    mycall = match.call()
    mycall[[1]] <- dada2::learnErrors
    eval.parent(mycall)
  }
}

#' Wrapper for dada2::dada
#' @param derep Dereplicated reads
#' @param err Error object
#' @param errorEstimationFunction Function to estimate error rates
#' @param selfConsist Use self-consistency
#' @param pool Pool reads
#' @param priors Sequences with prior information
#' @param multithread Use multiple threads
#' @param verbose Print progress messages
#' @param ... Additional arguments (passed to [dada2::dada()])
#' @return A `list` of [`dada`][dada2::dada-class] objects, or `NULL` if no
#' `errorEstimationFunction` is given and `selfConsist` is `FALSE`.  If only 0
#' or 1 [`derep`][dada2::derep-class] objects are given, the output will still
#' be a `list`.
#' @export
dada <- function(
    derep,
    err,
    errorEstimationFunction = choose_dada_error_function,
    selfConsist = FALSE,
    pool = FALSE,
    priors = character(0),
    multithread = FALSE,
    verbose = TRUE,
    ...
) {
  if (is.null(err) && isFALSE(selfConsist)) {
    NULL
  } else if (length(derep) == 0) {
    list()
  } else {
    mycall <- match.call()
    mycall[[1]] <- dada2::dada
    out <- eval.parent(mycall)
    if (methods::is(out, "dada") && !methods::is(derep, "derep")) {
      out <- list(out)
      names(out) <- names(derep)
    }
    out
  }
}

#' Wrapper for dada2::mergePairs
#' @param dadaF Forward denoised reads
#' @param derepF Forward dereplicated reads
#' @param dadaR Reverse denoised reads
#' @param derepR Reverse dereplicated reads
#' @param minOverlap Minimum overlap
#' @param maxMismatch Maximum mismatch
#' @param returnrejects Return rejected reads
#' @param propagateCol Columns to propagate
#' @param justConcatenate Just concatenate reads
#' @param trimOverhang Trim overhang
#' @param verbose Print progress messages
#' @param ... Additional arguments (passed to [dada2::mergePairs()])
#' @return A `data.frame` as returned by [dada2::mergePairs()], or a `list` of
#' such `data.frame`s.  In partticular, if the inputs `dadaF`, `derepF`, `dadaR`,
#' and `derepR` are all lists, even of length 0 or 1, the output will also be a
#' list of the same length.
#' @export
mergePairs <- function(dadaF, derepF, dadaR, derepR,
                       minOverlap = 12,
                       maxMismatch = 0,
                       returnrejects = FALSE,
                       propagateCol = character(0),
                       justConcatenate = FALSE,
                       trimOverhang = FALSE,
                       verbose = FALSE,
                       ...) {
  if (length(dadaF) == 0 || length(dadaR) == 0) {
    list()
  } else {
    mycall <- match.call()
    mycall[[1]] <- dada2::mergePairs
    out <- eval.parent(mycall)
    if (is.data.frame(out) && !methods::is(dadaF, "dada")) {
      out <- list(out)
      names(out) <- names(dadaF)
    }
    out
  }
}
