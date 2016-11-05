ptw.my <- function (ref, samp, selected.traces,
                 init.coef = c(0, 1, 0), try = FALSE,
		 warp.type = c("individual", "global"),
		 optim.crit = c("WCC", "RMS"),
                 mode = c("forward", "backward"),
		 smooth.param = ifelse(try, 0, 1e05),
		 trwdth = 20, trwdth.res = trwdth,
                 verbose = FALSE, dt.max=dt.max, 
                 ... )
{
  optim.crit <- match.arg(optim.crit)
  warp.type <- match.arg(warp.type)
  mode <- match.arg(mode)
  
  if (is.vector(ref)) ref <- matrix(ref, nrow = 1)
  if (is.vector(samp)) samp <- matrix(samp, nrow = 1)
  if (nrow(ref) > 1 && nrow(ref) != nrow(samp))
      stop("The number of references does not equal the number of samples")
  if (length(dim(ref)) > 2)
      stop("Reference cannot be an array")
  if (length(dim(samp)) > 2)
      stop("Sample cannot be an array")
  
  if (nrow(samp) == 1) warp.type <- "individual"
  
  r <- nrow(samp)
  
  if (!missing(selected.traces)) {
    samp <- samp[selected.traces,, drop=FALSE]  

    if (nrow(ref) > 1)
      ref <- ref[selected.traces,, drop=FALSE]
          
  }

  if (is.vector(init.coef)) init.coef <- matrix(init.coef, nrow = 1)
  if (warp.type == "global") {
    if (nrow(init.coef) != 1)
      stop("Only one warping function is allowed with global alignment.")
  } else {
    if (nrow(init.coef) != nrow(samp))
      if (nrow(init.coef) == 1) {
        init.coef <- matrix(init.coef, byrow = TRUE, nrow = nrow(samp),
                            ncol = length(init.coef))
      } else {
        stop("The number of warping functions does not match the number of samples")
      }
  }
  
  if (warp.type == "individual") {
    w <- matrix(0, nrow(samp), ncol(ref))   
    a <- matrix(0, nrow(samp), ncol(init.coef))
    v <- rep(0, nrow(samp))
    warped.sample <- matrix(NA, nrow=nrow(samp), ncol=ncol(samp))
    
    for (i in 1:nrow(samp)) {
      if (verbose & nrow(samp) > 1)
        cat(ifelse(nrow(ref) == 1,
                   paste("Warping sample", i, "with the reference \n"),
                   paste("Warping sample", i, "with reference \n", i)))

      if (nrow(ref) == 1) {
        rfrnc <- ref
      } else {
        rfrnc <- ref[i, , drop = FALSE]
      }
      quad.res <- pmwarp2(rfrnc, samp[i, , drop = FALSE],
                         optim.crit, init.coef[i,], try = try,
                         mode = mode,
                         smooth.param = smooth.param,
                         trwdth = trwdth, trwdth.res = trwdth.res,dt.max=dt.max,
                         ...)
      
      w[i, ] <- quad.res$w
      a[i, ] <- quad.res$a
      v[i] <- quad.res$v
      warped.sample[i, ] <- c(warp.sample(samp[i,,drop=FALSE],
                                          w[i,], mode = mode))
    }
  } else {
    if (nrow(ref)==1) 
        ref <- matrix(ref, nrow = nrow(samp), ncol = ncol(ref), byrow = TRUE)
    
    if (verbose) {
      if (nrow(ref) == 1) {
        cat("Simultaneous warping of samples with reference... \n")
      } else {
        cat("Simultaneous warping of samples with references... \n")
      }
    }

    quad.res <- pmwarp2(ref, samp, optim.crit, c(init.coef), try = try,
                       mode = mode, smooth.param = smooth.param,
                       trwdth = trwdth, trwdth.res = trwdth.res,dt.max=dt.max,
                       ...)
    
    w <- t(as.matrix(quad.res$w))
    a <- t(as.matrix(quad.res$a))
    v <- quad.res$v

    warped.sample <- t(warp.sample(samp, w, mode))
  }
  
  if (verbose) cat("\nFinished.\n")  

  result <-list(reference = ref, sample = samp,
                warped.sample = warped.sample,
                warp.coef = a, warp.fun = w,
                crit.value = v, optim.crit = optim.crit,
                mode = mode,
                warp.type = warp.type)
  class(result) <- "ptw"
  result
}

pmwarp2 <- function (ref, samp, optim.crit, init.coef, try = FALSE,
                    mode = c("forward", "backward"),
                    trwdth, trwdth.res, smooth.param, dt.max=dt.max, ...)
{
  mode <- match.arg(mode)
  
  ## Multiply coefficients to prevent them from becoming too 
  ## small for numeric precision.
  n <- length(init.coef)
  ncr <- ncol(ref)
  time <- (1:ncr) / ncr
  B <- matrix(time, nrow = ncr, ncol = n)
  B <- t(apply(B, 1, cumprod))/B
  a <- init.coef * ncr^(0:(n-1))

  if (optim.crit == "RMS" & smooth.param > 0) {
    samp.sm <- t(apply(samp, 1, difsm, smooth.param))
    ref.sm <- t(apply(ref, 1, difsm, smooth.param))
  }
  
  if (!try) { # perform optimization
    switch(optim.crit,
           RMS = {
             if (smooth.param > 0) {
               Opt <- optim(a, RMS, NULL, ref.sm, samp.sm, B, mode = mode, ...)
             } else {
               Opt <- optim(a, RMS, NULL, ref, samp, B, mode = mode, ...)
             }},
           WCC = {
             wghts <- 1 - (0:trwdth)/trwdth
             ref.acors <- apply(ref, 1, wac, trwdth = trwdth, wghts = wghts)
             Opt <- optim(a, WCC2, NULL, ref, samp, B,
                          trwdth = trwdth, wghts = wghts,
                          ref.acors = ref.acors, mode = mode, dt.max=dt.max, ...)
           })
    
    a <- c(Opt$par)
    v <- Opt$value

    ## if the optimization is done with a different smoothing or a
    ## different triangle, the final value for the optimization
    ## criterion is recalculated using the "original" data
    if ((optim.crit == "RMS" && smooth.param > 0) ||
        (optim.crit == "WCC" && trwdth != trwdth.res)) {
      v <- switch(optim.crit,
                  RMS = RMS(a, ref, samp, B, mode = mode),
                  WCC = WCC(a, ref, samp, B, trwdth.res, mode = mode))
    }
  }

  ## calculate, or possibly re-calculate, quality of current solution
  if (try) {
    if (optim.crit == "WCC") {
      v <- WCC(a, ref, samp, B, trwdth.res, mode = mode)
    } else {
      if (smooth.param > 0) {
        v <- RMS(a, ref.sm, samp.sm, B, mode = mode)
      } else {      
        v <- RMS(a, ref, samp, B, mode = mode)
      }
    }
  }

  ## back-transform coefficients
  w <- B %*% a
  a <- a/ncr^(0:(n-1))
  
  list(w = w, a = a, v = v)
} 

WCC2<-function(warp.coef, ref, samp, B, trwdth = 20, wghts, mode, ref.acors = NULL,dt.max=dt.max)
{
	WCC(warp.coef, ref, samp, B, trwdth, wghts, mode, ref.acors) * restriction(warp.coef,B,dt.max=dt.max)
}

restriction<-function(warp.coef,B,dt.max)
{
	if(length(warp.coef)>1) warp.coef[2]<-warp.coef[2]-1*dim(B)[1]
	w <- B %*% warp.coef
	foo<-function(x,dt.max){
		y<-abs(x)/dt.max
		y[y<1]<-1
		return(y)
	}
	out<-mean(foo(w,dt.max))
	return(out)
}



WCC <- function(warp.coef, ref, samp, B, trwdth = 20, wghts, mode,
                ref.acors = NULL)
{
  if (missing(wghts))
    wghts <- 1 - (0:trwdth)/trwdth

  if (is.null(ref.acors))
    ref.acors <- apply(ref, 1, wac, trwdth = trwdth, wghts = wghts)

  w <- B %*% warp.coef
  interp <- warp.sample(samp, w, mode)
  
  wccs <- sapply(1:ncol(interp),
                 function(i) {
                   wcc(ref[i, !is.na(interp[,i])],
                       interp[!is.na(interp[,i]), i],
                       trwdth = trwdth,
                       wghts = wghts,
                       acors1 = ref.acors[i])
                 })
  
  1 - mean(wccs) # so that an optimal value is zero
}

wac <- function(pattern1, trwdth, wghts = NULL)
{
  if (is.null(wghts)) 
    wghts <- 1 - (0:trwdth)/trwdth

  .C("wacdist",
     as.double(pattern1),
     as.integer(length(pattern1)),
     as.double(wghts),
     as.integer(trwdth),
     wacval = double(1),
     PACKAGE = "ptw")$wacval
}


wcc <- function(pattern1, pattern2, trwdth,
                wghts = NULL, acors1 = NULL, acors2 = NULL)
{
  if (is.null(wghts))
    wghts <- 1 - (0:trwdth)/trwdth
  
  if (is.null(acors1))
    acors1 <- wac(pattern1, trwdth, wghts)
  if (is.null(acors2))
    acors2 <- wac(pattern2, trwdth, wghts)

  .C("wccdist",
     as.double(pattern1),
     as.double(pattern2),
     as.integer(length(pattern1)),
     as.double(wghts),
     as.integer(trwdth),
     crossterm = double(1),
     PACKAGE = "ptw")$crossterm / (acors1*acors2)
}
warp.sample <- function(samp, w, mode) {
  if (mode == "backward") {
    apply(samp, 1, function(x) ##interpol(w, x))
      approx(x, NULL, w)$y)
  } else {
    apply(samp, 1, function(x) {
      approx(w, x, xout = 1:length(x))$y
    })
  }
}
edge<-function(ref,dt.max,slope.length=100){
	e1<-seq(1,dt.max,length.out=slope.length)
	e2<-seq(dt.max,1,length.out=slope.length)
	c(e1,rep(dt.max,length(ref)-length(e1)-length(e2)),e2)
}
