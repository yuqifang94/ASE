library(Gviz)
library(motifStack)
plotMB<-function (results, rsid, reverseMotif = TRUE, effect = c("strong", 
    "weak")) 
{
    g <- genome(results)[[1]]
    result <- results[names(results) %in% rsid]
    result <- result[order(sapply(result$motifPos, min), sapply(result$motifPos, 
        max)), ]
    result <- result[result$effect %in% effect]
    chromosome <- as.character(seqnames(result))[[1]]
    genome.package <- attributes(result)$genome.package
    genome.bsgenome <- eval(parse(text = genome.package))
    seq.len <- max(length(result$REF[[1]]), length(result$ALT[[1]]))
    distance.to.edge <- max(abs(c(sapply(result$motifPos, min), 
        sapply(result$motifPos, max)))) + 4
    from <- start(result)[[1]] - distance.to.edge + 1
    to <- end(result)[[1]] + distance.to.edge
    pwmList <- attributes(result)$motifs
    pwm.names <- result$providerId
    getmotifs <- mcols(pwmList)$providerId %in% result$providerId & 
        mcols(pwmList)$providerName %in% result$providerName
    pwms <- pwmList[getmotifs, ]
    pwms <- pwms[order(match(paste0(mcols(pwms)$providerId, mcols(pwms)$providerName), 
        paste0(result$providerId, result$providerName)))]
    if (reverseMotif) {
        for (pwm.i in seq_along(pwms)) {
            pwm.name <- names(pwms[pwm.i])
            pwm.id <- mcols(pwms[pwm.name, ])$providerId
            pwm.name.f <- mcols(pwms[pwm.name, ])$providerName
            doRev <- as.logical(strand(result[result$providerId == 
                pwm.id & result$providerName == pwm.name.f, ]) == 
                "-")
            if (doRev) {
                #Changed
                pwm <- pwms@listData[[pwm.i]]
                
                pwm <- pwm[, rev(1:ncol(pwm))]
                rownames(pwm) <- c("T", "G", "C", 
                  "A")
                pwm <- pwm[c("A", "C", "G", 
                  "T"), ]
                #Changed
                pwms@listData[[pwm.i]] <- pwm
                names(pwms@listData)[pwm.i] <- paste0(names(pwms)[pwm.i], 
                  "-:rc")
            }
        }
    }
    else {
        for (pwm.i in seq_along(pwms)) {
            pwm.name <- names(pwms[pwm.i])
            pwm.id <- mcols(pwms[pwm.name, ])$providerId
            pwm.name.f <- mcols(pwms[pwm.name, ])$providerName
            doRev <- as.logical(strand(result[result$providerId == 
                pwm.id & result$providerName == pwm.name.f, ]) == 
                "-")
            if (doRev) {
                #Changed
                pwm <- pwmspwms@listData[[pwm.i]]
                pwm <- pwm[, rev(1:ncol(pwm))]
                #Changed
                pwms@listData[[pwm.i]] <- pwm
                names(pwms)[pwm.i] <- paste0(names(pwms)[pwm.i], 
                  "-:r")
            }
        }
    }
    pwms <- lapply(names(pwms), function(x, pwms = pwms) {
        new("pfm", mat = pwms[[x]], name = x)
    }, pwms)
    pwms <- motifbreakR:::DNAmotifAlignment.2snp(pwms, result)
    pwmwide <- max(sapply(pwms, function(x) {
        ncol(x@mat)
    }))
    markerStart <- result$motifPos[[1]][1]
    if (markerStart > 0) {
        markerEnd <- length(result$altPos[[1]]) + 1
        markerEnd <- markerEnd - markerStart
        markerStart <- 1
    }
    else {
        markerStart <- -1 * markerStart
        markerEnd <- markerStart + length(result$altPos[[1]])
    }
    varType <- result$varType[[1]]
    varType <- switch(varType, Deletion = "firebrick", 
        Insertion = "springgreen4", Other = "gray13")
    markerRect <- new("marker", type = "rect", start = markerStart, 
        stop = markerEnd, gp = gpar(lty = 2, fill = NA, lwd = 3, 
            col = varType))
    for (pwm.i in seq_along(pwms)) {
        pwms[[pwm.i]]@markers <- list(markerRect)
    }
    ideoT <- try( IdeogramTrack(genome = g, chromosome = chromosome), 
        silent = TRUE)
    if (inherits(ideoT, "try-error")) {
        backup.band <- data.frame(chrom = chromosome, chromStart = 0, 
            chromEnd = length(genome.bsgenome[[chromosome]]), 
            name = chromosome, gieStain = "gneg")
        ideoT <- IdeogramTrack(genome = g, chromosome = chromosome, 
            bands = backup.band)
    }
    altseq <- genome.bsgenome[[chromosome]]
    at <- IRanges(start = start(result[1]), width = width(result[1]))
    if (result$varType[[1]] == "Deletion") {
        reflen <- length(result$REF[[1]])
        addedN <- DNAString(paste0(rep.int(".", reflen), 
            collapse = ""))
        addedN <- replaceLetterAt(addedN, at = (1:reflen)[-result$altPos[[1]]], 
            result$ALT[[1]])
        axisT <- GenomeAxisTrack(exponent = 0)
        seqT <- SequenceTrack(genome.bsgenome, fontcolor = colorset("DNA", 
            "auto"))
        altseq <- replaceAt(x = altseq, at = at, addedN)
    }
    else if (result$varType[[1]] == "Insertion") {
        altlen <- length(result$ALT[[1]])
        addedN <- DNAString(paste0(rep.int(".", altlen), 
            collapse = ""))
        addedN <- replaceLetterAt(addedN, at = (1:altlen)[-result$altPos[[1]]], 
            result$REF[[1]])
        refseq <- genome.bsgenome[[chromosome]]
        refseq <- DNAStringSet(replaceAt(x = refseq, at = at, 
            addedN))
        altseq <- replaceAt(x = altseq, at = at, result$ALT[[1]])
        names(refseq) <- chromosome
        seqT <- SequenceTrack(refseq, fontcolor = c(colorset("DNA", 
            "auto"), N = "#FFFFFF", . = "#FFE3E6"), 
            chromosome = chromosome)
    }
    else {
        axisT <- GenomeAxisTrack(exponent = 0)
        altseq <- replaceAt(x = altseq, at = at, result$ALT[[1]])
        seqT <- SequenceTrack(genome.bsgenome, fontcolor = colorset("DNA", 
            "auto"))
    }
    altseq <- DNAStringSet(altseq)
    names(altseq) <- chromosome
    seqAltT <- SequenceTrack(altseq, fontcolor = c(colorset("DNA", 
        "auto"), N = "#FFFFFF", . = "#FFE3E6"), 
        chromosome = chromosome)
    hiT <- HighlightTrack(trackList = list(seqT, seqAltT), start = start(result[1]) + 
        min(result[1]$altPos[[1]]) - 2, end = start(result[1]) + 
        min(result[1]$altPos[[1]]) - 2 + length(result[1]$altPos[[1]]), 
        chromosome = chromosome)
    selectingfun <- motifbreakR:::selcor
    detailfun <-  motifbreakR:::addPWM.stack
    getmotifs <- mcols(pwmList)$providerId %in% result$providerId & 
        mcols(pwmList)$providerName %in% result$providerName
    motif_ids <- names(pwmList)[getmotifs]
    names(motif_ids) <- mcols(pwmList)$providerName[getmotifs]
    for (mymotif_i in seq_along(result)) {
        mymotif <- result[mymotif_i]
        start(mymotif) <- start(mymotif) + min(mymotif$altPos[[1]]) - 
            1
        width(mymotif) <- length(mymotif$altPos[[1]])
        variant.start <- start(mymotif)
        variant.end <- end(mymotif)
        if (mymotif$motifPos[[1]][1] < 0) {
            start(mymotif) <- start(mymotif) + (mymotif$motifPos[[1]][1])
        }
        else {
            start(mymotif) <- start(mymotif) + (mymotif$motifPos[[1]][1] - 
                1)
        }
        if (mymotif$motifPos[[1]][2] < 0) {
            end(mymotif) <- end(mymotif) + (mymotif$motifPos[[1]][2] + 
                1)
        }
        else {
            end(mymotif) <- end(mymotif) + (mymotif$motifPos[[1]][2])
        }
        if ((result[mymotif_i]$varType == "Deletion" & 
            result[mymotif_i]$alleleDiff > 0) | (result[mymotif_i]$varType == 
            "Insertion" & result[mymotif_i]$alleleDiff < 
            0)) {
            mymotif <- c(mymotif, mymotif)
            end(mymotif)[1] <- variant.start - 1
            start(mymotif)[2] <- variant.end + 1
            mymotif[which.min(width(mymotif))]$motifPos <- NA
        }
        if (exists("mres")) {
            mres <- c(mres, mymotif)
        }
        else {
            mres <- mymotif
        }
    }
    result <- mres
    rm(mres)
    motif_ids <- motif_ids[result$providerName]
    motifT <- AnnotationTrack(result, id = motif_ids, fun = detailfun, 
        group = result$providerName, detailsFunArgs = list(pwm_stack = pwms), 
        feature = ifelse(!is.na(result$motifPos), paste0(result$geneSymbol, 
            "_motif"), ""), name = names(result)[[1]], 
        selectFun = selectingfun)
    if (exists("axisT")) {
        track_list <- list(ideoT, motifT, hiT, axisT)
    }
    else {
        track_list <- list(ideoT, motifT, hiT)
    }
    plotTracks(track_list, from = from, to = to, showBandId = TRUE, 
        cex.main = 0.8, col.main = "darkgrey", add53 = TRUE, 
        labelpos = "below", chromosome = chromosome, groupAnnotation = "group", 
        collapse = FALSE, min.width = 1, featureAnnotation = "feature", 
        cex.feature = 0.8, details.size = 0.85, detailsConnector.pch = NA, 
        detailsConnector.lty = 0, shape = "box", cex.group = 0.8, 
        fonts = c("sans", "Helvetica"))
    return(invisible(NULL))
}