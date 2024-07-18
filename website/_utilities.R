library(rhdf5) # note: dev version to support complex
library(purrr) # mapping 
library(dplyr) # wrangling 
library(ggplot2) # plotting
library(grid) # plotting
library(patchwork) # plotting
library(gt) # tables display
library(glue) # templates

# library(lobstr) # viewing objects

## ---- tmatrix_combinedindex
p_index <- function(l, m){
  p <- l*(l+1) + m
  return(p)
}

q_index <- function(p, s, pmax){
  
  q <- (s-1)*pmax + p
  return(q)
  
}

tmatrix_combinedindex <- function(d, lmax=max(d$l)){
  
  pmax <- lmax*(lmax + 1) + lmax
  mutate(d, p = p_index(l, m), 
         pp = p_index(lp, mp),
         q = q_index(p,  s,  pmax),
         qp = q_index(pp,  sp,  pmax),
         u = 2*(p-1) + (3 - s),
         up = 2*(pp-1) + (3 - sp))
  
}


## ---- tmatrix_breaks
tmatrix_breaks <- function(lmax){
  
  l <- seq.int(lmax)
  qmax <- 2*(lmax*(lmax + 1)+lmax)
  list(breaks = cumsum(rep((2*l+1), 2)),
       labels = rep(cumsum((2*l+1)), 2),
       minor_breaks = seq.int(qmax))
}

## ---- read_tmat

read_tmat <- function(f){
  h5closeAll()
  
  tmatrix <- rhdf5::h5read(f, 'tmatrix', compoundAsDataFrame = FALSE, native = TRUE)
  modes <- rhdf5::h5read(f, 'modes', compoundAsDataFrame = FALSE, native = TRUE)
  wavelengths <- rhdf5::h5read(f, 'vacuum_wavelength', compoundAsDataFrame = FALSE, native = TRUE)
  tmat_dims <- dim(tmatrix)
  # message("tmatrix dimensions: ",paste(tmat_dims, collapse="x"))
  # tmatrix can be either a matrix (single wavelength), or an array with 1 dimension wavelengths
  # so we need to handle both cases separately (it's awkward to slice by first index)
  if(length(tmat_dims) == 2){
    n <- tmat_dims[1]
  } else if(length(tmat_dims) == 3){
    stopifnot(tmat_dims[1] == length(wavelengths))
    n <- tmat_dims[2] # size of matrix
    
  } else {
    error("dimensions of t-matrix should be 2 or 3")
  }
  
  qmax <- n^2
  
  modes$s <- ifelse(modes$polarization == 'magnetic', 1, 2)
  l <- matrix(modes$l, nrow=n, ncol=n, byrow=FALSE); lp=t(l)
  m <- matrix(modes$m, nrow=n, ncol=n, byrow=FALSE); mp=t(m)
  s <- matrix(modes$s, nrow=n, ncol=n, byrow=FALSE); sp=t(s)
  
  # process a single wavelength
  single_wavelength <- function(tmat){
    long_tmat<- data.frame(s = as.vector(s), sp = as.vector(sp), 
                           l = as.vector(l), lp = as.vector(lp),  
                           m = as.vector(m), mp = as.vector(mp), 
                           value = as.vector(tmat)) |> 
      mutate(Tr = Re(value), Ti = Im(value), mod = Mod(Tr + 1i* Ti)) |> 
      mutate(p = p_index(l, m), 
             pp = p_index(lp, mp),
             q = q_index(p,  s,  max(p)),
             qp = q_index(pp,  sp,  max(p)),
             u = 2*(p-1) + (3 - s), 
             up = 2*(pp-1) + (3 - sp)) |> 
      arrange(s,sp,l,lp,m,mp)
    
    # long_tmat$p <- p_index(long_tmat$l, long_tmat$m)
    # long_tmat$pp <- p_index(long_tmat$lp, long_tmat$mp)
    # long_tmat$q <- q_index(long_tmat$p, long_tmat$s, max(long_tmat$p))
    # long_tmat$qp <- q_index(long_tmat$pp, long_tmat$sp, max(long_tmat$pp))
    
    return(long_tmat)
  }
  
  if(length(tmat_dims) == 2){
    results <- single_wavelength(tmatrix)
  } else {
    # this will slice by first index, and return a list with Nl entries
    l <- apply(tmatrix, 1, single_wavelength, simplify = FALSE)
    results <- do.call(rbind, l)
    results$wavelength <- rep(wavelengths, each=qmax)
  } 
  return(results)
}

# display
## ---- tmatGrob
library(grid)
library(scales)
tmatGrob <- function(m){ 
  rc <- dim(m)
  # m <- oob_censor(m, range=c(min_value, max_value), only.finite = FALSE)
  d <- scales::rescale(m, from = range(m, na.rm = TRUE, finite = T), to=c(1,0))
  d[is.infinite(d)] <- NA
  dim(d) <- rc
  g1 <- rectGrob(width=unit(1,'snpc'),height=unit(1,'snpc'),gp=gpar(fill='cornsilk'))
  g2 <- rasterGrob(d, interpolate = FALSE)
  
  grobTree(g1,g2)
}


grid.tmat <- function(...) grid.draw(tmatGrob(...))




