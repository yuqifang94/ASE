preparePWM<-function (pwmList = pwmList, filterp = filterp, bkg = bkg, scoreThresh = threshold,
    method = "default",ncores=ncores)
{
    bkg <- bkg[c("A", "C", "G", "T")]
    scounts <- as.integer(mcols(pwmList)$sequenceCount)
    scounts[is.na(scounts)] <- 20L
    pwmList.pc <- Map(function(pwm, scount) {
        pwm <- (pwm * scount + 0.25)/(scount + 1)
    }, pwmList, scounts)
    if (method == "ic") {
        pwmOmegas <- mclapply(pwmList.pc, function(pwm, b = bkg) {
            omegaic <- colSums(pwm * log2(pwm/b))
        },mc.cores=ncores)
    }
    if (method == "default") {
        pwmOmegas <- mclapply(pwmList.pc, function(pwm) {
            omegadefault <- colMaxs(pwm) - colMins(pwm)
        },mc.cores=ncores)
    }
    if (method == "log") {
        pwmList.pc <- mclapply(pwmList.pc, function(pwm, b) {
            pwm <- log(pwm) - log(b)
        }, b = bkg,mc.cores=ncores)
        pwmOmegas <- 1
    }
    if (method == "notrans") {
        pwmOmegas <- 1
    }
    pwmList.pc <- Map(function(pwm, omega) {
        if (length(omega) == 1 && omega == 1) {
            return(pwm)
        }
        else {
            omegamatrix <- matrix(rep(omega, 4), nrow = 4, byrow = TRUE)
            pwm <- pwm * omegamatrix
        }
    }, pwmList.pc, pwmOmegas)
    if (filterp) {
        pwmRanges <- Map(function(pwm, omega) {
            x <- colSums(colRanges(pwm))
            return(x)
        }, pwmList.pc, pwmOmegas)
        pwmList.pc2 <- mclapply(pwmList.pc, round, digits = 2,mc.cores=ncores)
        #speed limiting step
        pwmThresh <- mclapply(pwmList.pc2, TFMpv2sc, pvalue = scoreThresh,
            bg = bkg, type = "PWM",mc.cores=ncores)
        pwmThresh <- Map("+", pwmThresh, -0.02)
    }
    else {
        pwmRanges <- Map(function(pwm, omega) {
            x <- colSums(colRanges(pwm))
            return(x)
        }, pwmList.pc, pwmOmegas)
        pwmThresh <- rep.int(scoreThresh, times = length(pwmRanges))
    }
    pwmList@listData <- mclapply(pwmList, function(pwm) {
        pwm <- rbind(pwm, N = 0)
        colnames(pwm) <- as.character(1:ncol(pwm))
        return(pwm)
    },mc.cores=ncores)
    pwmList.pc <- mclapply(pwmList.pc, function(pwm) {
        pwm <- rbind(pwm, N = 0)
        colnames(pwm) <- as.character(1:ncol(pwm))
        return(pwm)
    },mc.cores=ncores)
    return(list(pwmList = pwmList, pwmListPseudoCount = pwmList.pc,
        pwmRange = pwmRanges, pwmThreshold = pwmThresh))
}

