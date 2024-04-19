#' Optimal granularity for quadrat counting
#' 
#' Given a point set, finds an optimal granularity for quadrat counting that
#' balances uniformity and robustness, as described in Ramos et al. (2021).
#'
#' @param point_set: the point set for which you want to find the optimal quadrat
#' size. Must be provided as a dataframe with the first column being the
#' 'easting' (e.g. x, longitude) and the second column the northing (e.g y, latitude)
#' @param random_sample: whether uniformity and robustness is estimated from a
#' random sample (T) or by generating a regular grid at the granularity being
#' tested (F). Using random samples generally takes less time and yields similar
#' results.
#' @param nsamples: number of samples taken if random_sample == T. In practice,
#' the default value seems to work fine.
#' @param signif: significance level for the Complete Spatial Randomness (CSR) test that is 
#' applied for the samples at each different granulairity considered.
#' @param tradeoff_crit: dictactes how the balance between uniformity and robustness
#' is determined in order to choose an optimal quadrat size. 'sum' m,eans that the
#' granularity with the greates sum of uniformity and robustness gets picked, 'product'
#' means that the granularity yielding the greatest product is picked.
#' @param uniformity_method: whether CSR is tested via the quadratcount method or
#' the nearest neighbor method. In practice, the result should not differ much
#' @param robustness_method: how the robustness of each granularity is estimated 
#' (via a Poisson model, a Binomail model, or by resampling the original point set)
#' In practice, any of the options yields similar results.
#' @param robustness_k: robustness of a cell is calculated by taking the estimated coefficient
#' of variation for a cell - let's call it x - and applying the function exp(k*x). This
#' parameter specificies which k is used. In practice, the final result is not very
#' sensitive to the specific value of k
#' @param verbose: whether to print messages while running the function or not
#' @param my_scales: which granularities to be tested. If not provided, a set is
#' automatically generated. The granularities, if provided, should be given as the
#' dimension (side) of a cell, in the unit used for the coordinates of point_set.
#' @param W: the window of interest to be considered when doing the analyzis of
#' point_set. If not provided, W is calculated as the minimum bounding rectangle
#' for point_set.
#' @param grid_crs: coordinate reference system of the points, will be ascribed to the
#' resulting grid
#' @keywords quadrat granularity adequate optimal
#' @author Rafael G. Ramos (main proponent and coder), Marcos Prates (contributor)
#' @references Ramos, R. G., Silva, B. F., Clarke, K. C., & Prates, M. (2021). 
#' Too Fine to be Good? Issues of Granularity, Uniformity and Error in Spatial 
#' Crime Analysis. Journal of Quantitative Criminology, 1-25.
#' robust.quadcount()
#' @import terra
#' @importFrom magrittr %>%
#' @importFrom stats chisq.test complete.cases filter rbinom runif
#' @importFrom stats smooth.spline var
#' @export
#' 
#' @examples
#' library(robustmap)
#' # Loading point data
#' burglary <- sf::read_sf(dsn = "inst/extdata", layer = "burglary")
#' 
#' burglary <- data.frame(x=burglary$lon_m,
#'                        y=burglary$lat_m)
#' 
#' # Estimating optimal granularity using robust.quadcount
#' burglary_map <- robust.quadcount(burglary,verbose = T)
#' 
#' # Retriving estimated granularity
#' burglary_map$opt_granularity
#' 
#' # Plotting resulting map
#' terra::plot(burglary_map$counts)

robust.quadcount<-function(point_set,
                           random_samples=T,
                           nsamples=500,
                           signif=0.99,
                           tradeoff_crit=c("product","sum"),#,"derivative"), #NOT doing derivative right now
                           uniformity_method=c("Quadratcount","Nearest-neighbor"),
                           robustness_method=c("Poisson","Binomial","Resampling"),
                           robustness_k = -3,
                           verbose=F,
                           my_scales=NULL,
                           W=NULL,
                           grid_crs=NULL){
  
  # ADD CRS INFO TO THE OUTPUTS
  tradeoff_crit <- match.arg(tradeoff_crit)
  uniformity_method <- match.arg(uniformity_method)
  robustness_method <- match.arg(robustness_method)
  
  if(is.null(W)) {
    W <- c(min(point_set[,1]),max(point_set[,1]),min(point_set[,2]),max(point_set[,2]))
  }
  
  # checking for duplicates. If there are, add jittering
  if(any(duplicated(point_set))) {
    point_set <- spatstat.geom::as.ppp(point_set,W,check=F) # no need to check, we know and will fix it!
    point_set <- spatstat.geom::rjitter(point_set, retry=TRUE, nsim=1, drop=TRUE)
    point_set <- data.frame(x = point_set$x,y = point_set$y)
  }
  
  # initializing a set of scales to test if that hasn't been provided
  if(is.null(my_scales)) {
    maxscale <- 0.10*min(c(W[2] - W[1],W[4] - W[3]))
    minscale <- maxscale*0.02
    my_scales <-  seq(minscale,maxscale,minscale)
  }
  # calculating robustness and uniformity for each granularity in 'my_scales'
  my_spatialstats <- get_spatialstats_all(point_set,
                                          my_scales,
                                          my_scales,
                                          nsamples=nsamples,
                                          random_samples=random_samples,
                                          signif=signif,
                                          uniformity_method=uniformity_method,
                                          robustness_method=robustness_method,
                                          robustness_k = robustness_k,
                                          W = W,
                                          verbose=verbose)
  
  if(verbose) {
    print("Estimating optimal balance between uniformity and robustness")
  }
  robust <- my_spatialstats$robustness
  unif <- my_spatialstats$uniformity
  plot(my_scales,robust)
  plot(my_scales,unif)
  
  # tradeoff analysis
  if(tradeoff_crit == "sum") {
    metric <- unif+robust
    splinefit <- smooth.spline(x=my_scales[!is.na(metric)],y=metric[!is.na(metric)],df=10)
    my_spline <- predict(splinefit,x=my_scales)
    opt_i <- which.max(my_spline$y)
    opt_granularity <- my_spline$x[opt_i]
    opt_ur <- my_spline$y[opt_i]
    
  } else if(tradeoff_crit == "product") {
    metric <- unif*robust
    splinefit <- smooth.spline(x=my_scales[!is.na(metric)],y=metric[!is.na(metric)],df=10)
    my_spline <- predict(splinefit,x=my_scales)
    opt_i <- which.max(my_spline$y)
    opt_granularity <- my_spline$x[opt_i]
    opt_ur <- my_spline$y[opt_i]
    
    #} else if(tradeoff_crit == "derivative") {
    
    #  splinefit = smooth.spline(x=robust,y=unif,df=6)
    #  my_spline = predict(splinefit,x=seq(0,1,0.01))
    #  my_derivs = predict(splinefit,x=seq(0,1,0.01),deriv=1)
    #  opt_i = which.min((my_derivs$y+1)^2)
    #  opt_robust = my_spline$x[opt_i]
    #  opt_uniformity = my_spline$y[opt_i]
    
    #  splinefit = smooth.spline(x=my_scales,y=robust,df=10)
    #  my_spline = predict(splinefit,x=seq(25,1000,5))
    #  opt_i = which.min((my_spline$y-opt_robust)^2)
    #  opt_granularity_1 = my_spline$x[opt_i]
    #  splinefit = smooth.spline(x=my_scales,y=unif,df=10)
    #  my_spline = predict(splinefit,x=seq(25,1000,5))
    #  opt_i = which.min((my_spline$y-opt_uniformity)^2)
    #  opt_granularity_2 = my_spline$x[opt_i]
    #  opt_granularity = 0.5*(opt_granularity_1 + opt_granularity_2)
    
  } else {
    print("Tradeoff criteria not recognized. Options allowed: sum, product.")#, derivative.")
    return(NULL)
  }
  map <- list()
  
  map$counts <- point_set %>% points2raster(opt_granularity,W=W,grid_crs=grid_crs)
  
  terra::ext(map$counts) <- W
  map$opt_granularity <- opt_granularity
  final_sample <- create_samples(point_set=point_set,
                                 random_samples = F,
                                 nsamples = nsamples,
                                 window_w = opt_granularity,
                                 window_h = opt_granularity,
                                 W = W)
  
  stats_final_sample <- get_spatialstats_sample(sample_set = final_sample,
                                                uniformity_method = uniformity_method,
                                                signif=signif,
                                                robustness_method = robustness_method)
  tmp <- matrix(data=stats_final_sample$samples_covar,
                nrow=nrow(map$counts),
                ncol=ncol(map$counts))
  map$covar <- tmp %>% terra::rast() %>% terra::flip(direction="horizontal")
  terra::ext(map$covar) <- terra::ext(map$counts)
  tmp <- matrix(data=stats_final_sample$samples_csr_pass,
                nrow=nrow(map$counts),
                ncol=ncol(map$counts))
  map$is.csr <- !(terra::flip(terra::rast(tmp),direction="horizontal"))
  terra::ext(map$is.csr) <- terra::ext(map$counts)
  robust <- my_spatialstats$robustness
  unif <- my_spatialstats$uniformity
  map$uniformity.curve <- unif
  map$granularities.tested <- my_scales
  map$robustness.curve <- robust
  
  if(!is.null(grid_crs)) {
    terra::crs(map$counts) <- grid_crs
    terra::crs(map$covar) <- grid_crs
    terra::crs(map$is.csr) <- grid_crs
  }
  
  return(map)
}


#
#  Function that estimates the internal uniformity and robustness to error for 'point_set',
#  according to the granularity dimensions listed in 'scales_x' and 'scales_y' (these should be the
#  same size; if one is larger than the other, the extra elements in the largest are discarded)
#
#  - random_sample' determinines whether the metric should be estimated from random samples
#    taken from point_set, or from a contiguous regular partion of point_set.
# 
#  - 'n_samples' determines how many samples are taken, if 'random_sample' == T.
# 
#  - 'signif' determines the threshold for considering a sample to be spatially random enough
#    for the sake of calculating internal uniformity.
#

get_spatialstats_all<-function(point_set,
                               scales_x,
                               scales_y,
                               random_samples=T,
                               nsamples=500,
                               signif=0.90,
                               robustness_method=c("Poisson","Binomial","Resampling"),
                               robustness_k = -3,
                               uniformity_method=c("Quadratcount","Nearest-neighbor"),
                               W = W,
                               verbose=verbose){
  
  robustness_method <- match.arg(robustness_method)
  uniformity_method <- match.arg(uniformity_method)
  
  if(robustness_k >= 0) {
    print("robust_k parameter must be less than zero. See documentation for details.")
    return(NULL)
  }
  
  len_scales <- min(c(length(scales_x),length(scales_y)))
  total_points <- nrow(point_set)
  
  # allocating vectors for our outputs
  
  # avg number of points per sample (and var)
  scales_count <- vector(length=len_scales)
  scales_var <- vector(length=len_scales)
  # mean and median for the nearest neighbor p-value test
  scales_csr_p <- vector(length=len_scales)
  scales_csr_median_p <- vector(length=len_scales)
  # proportion of samples that passed the threshold for CSR using the nearest neighbor approach
  scales_csr_pass <- vector(length=len_scales)
  # mean values for the coef_of var (later used to calculate robustness)
  scales_covar <- vector(length=len_scales)
  # number of zeros
  scales_zeros <- vector(length=len_scales)
  
  for(i in 1:len_scales) {
    if(verbose) {
      print_progress(i,len_scales,message="Estimating uniformity and robustness at various scales:")
    }
    
    # generating set of quadrats of the current granularity for taking the samples
    my_samples <- create_samples(point_set=point_set,
                                 random_samples = random_samples,
                                 nsamples = nsamples,
                                 window_w = scales_x[i],
                                 window_h = scales_y[i],
                                 W = W)
    
    
    my_stats_sample <- get_spatialstats_sample(sample_set = my_samples,
                                               signif = signif,
                                               uniformity_method = uniformity_method,
                                               robustness_method = robustness_method)
    
    # avg number of points per sample (and var)
    scales_count[i] <- mean(my_stats_sample$samples_count[!is.na(my_stats_sample$samples_count)])
    scales_var[i] <- var(my_stats_sample$samples_count[!is.na(my_stats_sample$samples_count)])
    
    # mean and median of the p-value for the CSR tests
    scales_csr_p[i] <- mean(my_stats_sample$samples_csr_p[!is.na(my_stats_sample$samples_csr_p)])
    scales_csr_median_p[i] <- median(my_stats_sample$samples_csr_p[!is.na(my_stats_sample$samples_csr_p)])
    # proportion of samples that passed the threshold for CSR
    scales_csr_pass[i] <- mean(my_stats_sample$samples_csr_pass[!is.na(my_stats_sample$samples_csr_pass)])
    
    # mean values for the coef_of var
    scales_covar[i] <- mean(my_stats_sample$samples_covar[!is.na(my_stats_sample$samples_covar)])
    # number of zeros
    scales_zeros[i] <- my_stats_sample$samples_zeros/my_samples$nsample
  }
  scales_robustness <- exp(-3*scales_covar)
  scales_uniformity <- 1-scales_csr_pass
  out = list(count = scales_count,
             var = scales_var,
             #csr_p = scales_csr_p,
             #csr_median_p = scales_csr_median_p,
             uniformity = scales_uniformity,
             robustness = scales_robustness,
             zerors = scales_zeros
  )
  return(out)
}

create_samples<-function(point_set,random_samples,nsamples,window_w,window_h,W=NULL) {
  if(is.null(W)) {
    max_b_x <- max(point_set$x)
    max_b_y <- max(point_set$y)
    min_b_x <- min(point_set$x)
    min_b_y <- min(point_set$y)
  } else {
    min_b_x <- W[1]
    max_b_x <- W[2]
    min_b_y <- W[3]
    max_b_y <- W[4]
  }
  w_b <- max_b_x-min_b_x
  h_b <- max_b_y-min_b_y
  if(random_samples == T) {
    # quadrats taken randomly
    offset_x <- ((w_b-window_w)*runif(nsamples))+min_b_x
    offset_y <- ((h_b-window_h)*runif(nsamples))+min_b_y
  }
  else {
    # contiguous quadrats
    nx <- floor(w_b/window_w)
    ny <- floor(h_b/window_h)
    nsamples <- nx*ny
    offset_x <- matrix(nrow=ny,ncol=nx)
    offset_y <- matrix(nrow=ny,ncol=nx)
    for(k in 1:ny) offset_x[k,] <- (0:(nx-1))*window_w+min_b_x
    for(k in 1:nx) offset_y[,k] <- (0:(ny-1))*window_h+min_b_y
    offset_x <- as.vector(offset_x)
    offset_y <- as.vector(offset_y)
  }
  out <- list()
  out$point_set <- point_set
  out$offset_x <- offset_x
  out$offset_y <- offset_y
  out$window_w <- window_w
  out$window_h <- window_h
  out$random_samples <- random_samples
  out$npopulation <- nrow(point_set)
  out$nsamples <- nsamples
  
  return(out)
}

get_spatialstats_sample<-function(sample_set,
                                  signif,
                                  robustness_method=c("Poisson","Binomial","Resampling"),
                                  uniformity_method=c("Quadratcount","Nearest-neighbor")){
  
  robustness_method <- match.arg(robustness_method)
  uniformity_method <- match.arg(uniformity_method)
  
  nsamples <- length(sample_set$offset_x)
  
  # allocating temporary structures
  samples_count <- vector(length=nsamples)
  samples_csr_p <- vector(length=nsamples)
  samples_csr_pass <- vector(length=nsamples)
  samples_covar <- vector(length=nsamples)
  samples_zeros <- 0
  
  # for each quadrat, take a the points inside and test for robustness and internal uniformity
  for(j in 1:nsamples) {
    # extracting the points inside the quadrat
    W_ext <- c(sample_set$offset_x[j],
               sample_set$offset_x[j]+sample_set$window_w,
               sample_set$offset_y[j],
               sample_set$offset_y[j]+sample_set$window_h)
    
    is_inside <- spatstat.geom::inside.owin(x=sample_set$point_set$x,
                                            y=sample_set$point_set$y,
                                            w=W_ext)
    
    sub_points <- sample_set$point_set[is_inside,]
    
    # estimate robustness and uniformity
    if(nrow(sub_points)==0) {
      # if there are no points inside, the metrics are NA
      samples_zeros <- samples_zeros + 1
      samples_count[j] <- NA
      samples_csr_p[j] <- NA
      samples_csr_pass[j] <- NA
      samples_covar[j] <- NA
    }
    else {
      # in case we have points, calculate the metric
      
      W_disp <- c(sample_set$offset_x[j]-0.000001,#FIX the margin to be relative
                  (sample_set$offset_x[j]+sample_set$window_w)+0.000001,
                  sample_set$offset_y[j]-0.000001,
                  (sample_set$offset_y[j]+sample_set$window_h)+0.000001)
      
      # testing Complete Spatial Randomness (CSR) with Clark-Evans nearest neighbor test (part of estimating uniformity)
      samples_count[j] <- nrow(sub_points)
      if(uniformity_method == "Nearest-neighbor") {
        if(samples_count[j] > 20) {
          my_clarkevans <- sub_points %>%
            spatstat.geom::as.ppp(W_disp) %>%
            spatstat.explore::clarkevans.test(alternative="two.sided")#,correction = "Donnelly")
        } else {
          my_clarkevans <- sub_points %>%
            spatstat.geom::as.ppp(W_disp) %>%
            spatstat.explore::clarkevans.test(alternative="two.sided",nsim=100)#,correction = "Donnelly")
        }
        samples_csr_p[j] <-  my_clarkevans$p.value
        # assign T to quadrats that pass the threshold 'signif' for CSR
        samples_csr_pass[j] <-  my_clarkevans$p.value < (1-signif)
      }
      else if(uniformity_method == "Quadratcount") {
        # test CST with the quadrat test
        my_quadrattest <- sub_points %>% 
          spatstat.geom::as.ppp(W_disp) %>%
          spatstat.explore::quadrat.test(nx=5,method="MonteCarlo",conditional=T)
        samples_csr_p[j] <- my_quadrattest$p.value
        samples_csr_pass[j] <- my_quadrattest$p.value < (1-signif)
      }
      else {
        samples_csr_p[j] <- NA
        samples_csr_pass[j] <- NA
      }
      samples_covar[j] <- calc_covar(nrow(sub_points),sample_set$npopulation,robustness_method)
    }
  }
  out <- list()
  out$samples_count <- samples_count
  out$samples_csr_p <- samples_csr_p
  out$samples_csr_pass <- samples_csr_pass
  out$samples_covar <- samples_covar
  out$samples_zeros <- samples_zeros
  Sys.sleep(1)
  return(out)
}


calc_covar<-function(nsub,ntotal,robustness_method){
  # Calculating robustness.
  # First we calculate the expected coefficient of variation for the samples, using different methods
  prob_event <- nsub/ntotal
  if(robustness_method == "Binomial") {
    # calculating coef of var using the Binomial estimation method (see paper)
    samples_covar <- sqrt(prob_event*(1-prob_event)*ntotal)/nsub#sd(sim_rates)/mean(sim_rates)
  }
  else if(robustness_method == "Poisson") {
    # calculating coef of var using the Poisson estimation method (see paper)
    tmp <- MASS::fitdistr(nsub,"Poisson")
    samples_covar <- tmp$sd/tmp$estimate
  }
  else if(robustness_method == "Resampling") {
    # calculating coef of var using the resampling estimation method (see paper)
    sim_rates <- rbinom(1000,ntotal,prob_event)
    samples_covar <- mean(sqrt(1/sim_rates[sim_rates!=0]))
  }
  return(samples_covar)
}

#' Robust gridded data(stack) from point data
#' 
#' Estimate robust grid using either counts or density (NOT IMPLEMENTED YET)
#' and with a time dimension (optional) as separate layers.
#' @export
#'
#' @param point_set: the point set for which you want to find the optimal quadrat
#' size. Must be provided as a dataframe with the first column being the
#' 'easting' (e.g. x, longitude) and the second column the northing (e.g y, latitude)
#' @param aggregation_method: how should the points be aggregated per cell?
#' 'count' will count the number of points per cell using robust.quadcount
#' (and proceed from there), while 'density' (NOT IMPLEMENTED YET) will estimate
#' the density from the point set using robust.density 
#' @param temporal: should time slices be generated from the point set? if so
#' a time stamp should be provided in a third column for point_set.
#' @param robust_timeslices: if temporal is TRUE, robust_timeslices being TRUE
#' further processes the time slices for improved accuraracy (NEED TO FORMALIZE VALIDATION!)
#' @param verbose: whether to print messages while running the function or not
#' @param W: the window of interest to be considered when doing the analyzis of
#' point_set. If not provided, W is calculated as the minimum bounding rectangle
#' for point_set.
#' 
#' The following are parameters for the 'count' aggregation method,
#' to be passed on to the robust.count 
#' 
#' @param random_sample: whether uniformity and robustness is estimated from a
#' random sample (T) or by generating a regular grid at the granularity being
#' tested (F). Using random samples generally takes less time and yields similar
#' results.
#' @param nsamples: number of samples taken if random_sample == T. In practice,
#' the default value seems to work fine.
#' @param signif: significance level for the Complete Spatial Randomness (CSR) test that is 
#' applied for the samples at each different granulairity considered.
#' @param tradeoff_crit: dictactes how the balance between uniformity and robustness
#' is determined in order to choose an optimal quadrat size. 'sum' m,eans that the
#' granularity with the greates sum of uniformity and robustness gets picked, 'product'
#' means that the granularity yielding the greatest product is picked.
#' @param uniformity_method: whether CSR is tested via the quadratcount method or
#' the nearest neighbor method. In practice, the result should not differ much
#' @param robustness_method: how the robustness of each granularity is estimated 
#' (via a Poisson model, a Binomail model, or by resampling the original point set)
#' In practice, any of the options yeilds similar results.
#' @param robustness_k: robustness of a cell is calculated by taking the estimated coefficient
#' of variation for a cell - let's call it x - and applying the function exp(k*x). This
#' parameter specificies which k is used. In practice, the final result is not very
#' sensitive to the specific value of k
#' @param my_scales: which granularities to be tested. If not provided, a set is
#' automatically generated. The granularities, if provided, should be given as the
#' dimension (side) of a cell, in the unit used for the coordinates of point_set.
#' @param grid_crs: coordinates reference system of the points, will be ascribed to the
#' resulting grid
#' 
#' @keywords quadrat granularity adequate optimal time
#' @references Ramos, R. G., Silva, B. F., Clarke, K. C., & Prates, M. (2021). 
#' Too Fine to be Good? Issues of Granularity, Uniformity and Error in Spatial 
#' Crime Analysis. Journal of Quantitative Criminology, 1-25.
#' robust.quadcount()

robust.grid<-function(point_set,
                      aggregation_method=c("count","density"),
                      temporal=T,
                      robust_timeslices = T,
                      opt_granularity=NULL,
                      random_samples=T,
                      nsamples=500,
                      signif=0.99,
                      tradeoff_crit=c("product","sum"),#,"derivative"), #NOT doing derivative right now
                      uniformity_method=c("Quadratcount","Nearest-neighbor"),
                      robustness_method=c("Poisson","Binomial","Resampling"),
                      robustness_k = -3,
                      verbose=F,
                      my_scales=NULL,
                      W=NULL,
                      grid_crs=NULL) {
  
  # if extent has not been provided, it is defined from the point set
  if(is.null(W)) {
    W <- c(min(point_set$x),
           max(point_set$x),
           min(point_set$y),
           max(point_set$y))
  }
  
  # checking for duplicates. If there are, add jittering
  if(any(duplicated(point_set))) {
    times <- point_set$time
    point_set <- spatstat.geom::as.ppp(point_set,W,check=F) # no need to check, we know and will fix it!
    point_set <- spatstat.geom::rjitter(point_set, retry=TRUE, nsim=1, drop=TRUE)
    if(nrow(point_set) == length(times)) {
      point_set <- data.frame(x = point_set$x,y = point_set$y,time <- times)
    } else {
      print("adding jitter dropped some points. quitting")
      return(NULL)
    }
  }
  
  # estimating optimal granularity if it has not been provided
  if(is.null(opt_granularity)) {
    # in this case the full report of the robust.quadcount test is provided
    # including the tradeoff curves of uniformity and robustness.
    fullmap <- point_set %>% robust.quadcount(random_samples=random_samples,
                                              nsamples=nsamples,
                                              signif=signif,
                                              tradeoff_crit=tradeoff_crit,
                                              uniformity_method=uniformity_method,
                                              robustness_method=robustness_method,
                                              robustness_k = robustness_k,
                                              verbose=verbose,
                                              my_scales=my_scales,
                                              W=W,
                                              grid_crs=grid_crs)
    opt_granularity <- fullmap$opt_granularity
  } else{
    # if an optimal granularity is provided, only the quadrat count at that
    # granularity is generated.
    fullmap <- point_set %>% points2raster(opt_granularity,W=W,grid_crs)
  }
  
  count_per_timestamp <- NULL
  if(temporal) {
    time_stamps <- unique(point_set$time)
    for(i in 1:length(time_stamps)){
      print(time_stamps[i])
      print(point_set$time)
      print("HA")
      print(point_set$time[(point_set$time == time_stamps[i])])
      newlayer <- point_set %>%
        dplyr::filter(time == time_stamps[i]) %>%
        points2raster(opt_granularity,W=W)
      if(is.null(count_per_timestamp)) {
        count_per_timestamp <- newlayer
      } else {
        count_per_timestamp <- c(count_per_timestamp,newlayer)
      }
    }
    terra::ext(count_per_timestamp) <- W
    terra::crs(count_per_timestamp) <- grid_crs
    if(robust_timeslices) {
      estimated_slice <- robust.timeslices(count_per_timestamp,verbose)
      terra::ext(estimated_slice) <- W
      terra::crs(estimated_slice) <- grid_crs
    }
  }
  
  out <- list(
    estimated_rates_timestamp = estimated_slice,
    observed_counts_timestamp = count_per_timestamp,
    observed_total_counts = fullmap
  )
  return(out)
}

#' Robust time slices from gridded data stack
#' 
#' Improves the rates per cells for each timeslice by considering information
#' from the other slices. For each location, model (auto.arima only for now)
#' is estimated for the corresponding time series. The predicted value at that 
#' cell-time is assigned as the improved, more robust rate. More formal
#' validation is still required, but simulation tests I've ran show that it
#' does yield better estimates.
#' @export

robust.timeslices<-function(observed_counts_timestamp,verbose=F) {
  estimated_rates_timestamp <- observed_counts_timestamp
  estimated_rates_timestamp[] <- NA
  
  for(i in 1:nrow(estimated_rates_timestamp)) {
    for(j in 1:ncol(estimated_rates_timestamp)) {
      if(verbose) {
        print_progress(j+(i-1)*ncol(estimated_rates_timestamp[[1]]),
                       length(estimated_rates_timestamp[[1]][]),
                       "Processing time slices")
      }
      count_series <- observed_counts_timestamp[i,j,] %>% t() %>% as.numeric()
      
      if(length(unique(count_series)) == 1) {
        estimated_rates_timestamp[i,j,] = count_series
      }
      else {
        #WHAT MODEL TO USE? MAYBE HAVE OPTIONS?
        #my_arima <- arima(myrobb)
        my_arima <- forecast::auto.arima(count_series)
        estimated_rates_timestamp[i,j,] <- (count_series - my_arima$residuals)
        for(k in 1:length(estimated_rates_timestamp[i,j,])){
          if(estimated_rates_timestamp[i,j,k] < 0) {
            estimated_rates_timestamp[i,j,k] <- 0
          }
        }
      }
    }
  }
  return(estimated_rates_timestamp)
}

#' Robust per capita rates from gridded data
#'
#' Estimate robust per capita rates from a event count map (numerator)
#' and a reference population map (denominator); both need to be
#' SpatRaster grids with the same dimensions and covering the same extent.
#' Event count map can be generated using robust.grid or robust.quadcount
#' from a point set.
#' @export
#'
#' @param numerator: an event count grid (SpatRaster)
#' @param denominator: a reference population map grid (SpatRaster), of the same
#' dimensions, extent and resolution as numerator.
#' @param bd: bandwidth to be used by the GWR model. If NULL (default),
#' bandwidth is estimated by cross-validation (using spgwr::gwr.sel())
#' @param weighted: if true, denominator is used also as weights for the 
#' residuals in a Weighted Least Squares fashion. The rational is that areas 
#' with a greater denominator (e.g. more individuals experiencing a per capita
#' rate) should have greater weight in the estimation of the rate.
#' @keywords quadrat granularity adequate optimal
#' @author Rafael G. Ramos (main proponent and coder)
#' @references Ramos, R. G. (2021). Improving victimization risk estimation: 
#' A geographically weighted regression approach. ISPRS International Journal of
#'  Geo-Information, 10(6), 364.
#' @export

robust.percapita<-function(numerator,denominator,bd=NULL,weighted=F){
  
  percapita_map <- numerator
  
  my_data <- data.frame(numerator = numerator[[1]][],
                        denominator = denominator[[1]][])
  
  colnames(my_data) <- c("numerator","denominator")
  
  is_complete <- complete.cases(my_data)
  
  my_coords <- terra::crds(numerator[[1]])
  if(weighted) {
    weights <- as.numeric(denominator[[1]][])[is_complete]
  } else {
    weights <- NULL
  }
  
  if(is.null(bd)) {
    print("Estimating bandwidth...")
    bd <- spgwr::gwr.sel(formula = numerator~denominator+0,
                         data = my_data[is_complete,],
                         coords = my_coords[is_complete,],
                         weights = weights)
    print(paste("Bandwidth is ",bd,sep=""))
  } 
  print("Estimating per capita ratio...")
  percapita_model <- spgwr::gwr(formula = numerator~denominator+0,
                                data = my_data[is_complete,],
                                coords = my_coords[is_complete,],
                                bandwidth = bd,
                                weights = weights)
  percapita_map[] <- NA
  percapita_map[is_complete] <- percapita_model$SDF$denominator
  
  return(percapita_map)
}

################################################################################
#
#   Utility functions, for plotting quadrat counts etc.
#
################################################################################

quad2mat<-function(quadset) {
  tmp <- matrix(data=quadset[],nrow=nrow(quadset),ncol=ncol(quadset))
  return(tmp)
}

print_progress<-function(cur_iter,n_iters,message="Processing") {
  
  cur_percent <- round(100*(cur_iter/n_iters))
  
  cat('\014',append=T)
  cat(paste(message,":\n",sep=""))
  cat("[")
  
  for(k in 1:cur_percent){
    cat("=")
  }
  if(cur_percent < 100) {
    for(k in (cur_percent+1):100){
      cat(" ")
    }
  }
  cat("] | ")
  cat(paste(cur_percent,"%\n",sep=""),append=T)
}

points2raster<-function(point_set,my_scale,W=NULL,mask=NULL,grid_crs=NULL){
  if(is.null(W)) {
    W <- c(min(point_set$x),
           max(point_set$x),
           min(point_set$y),
           max(point_set$y))
    max_b_x <- max(point_set$x)
    max_b_y <- max(point_set$y)
    min_b_x <- min(point_set$x)
    min_b_y <- min(point_set$y)
  } else {
    min_b_x <- W[1]
    max_b_x <- W[2]
    min_b_y <- W[3]
    max_b_y <- W[4]
  }
  nx <- floor((max_b_x - min_b_x)/my_scale)
  ny <- floor((max_b_y - min_b_y)/my_scale)
  my_ret <- point_set %>%
    spatstat.geom::as.ppp(W) %>%
    spatstat.geom::quadratcount(nx,ny)
  
  if(!is.null(mask)){
    Wmask <- c(min(mask$x),
               max(mask$x),
               min(mask$y),
               max(mask$y))
    my_mask <- mask %>%
      spatstat.geom::as.ppp(Wmask) %>%
      spatstat.geom::quadratcount(nx,ny) 
    my_ret[my_mask <= 0] <- NA
  }
  
  my_ret <- my_ret %>% quad2mat() %>% terra::rast()
  terra::ext(my_ret) <- W
  if(!is.null(grid_crs)) {
    terra::crs(my_ret) <- grid_crs
  }
  return(my_ret)
}

################################################################################
#
# WORK IN PROGRESS!
#
################################################################################

#
# Density estimation using max robustness while preserving uniformity
#

robust.density.unit<-function(count_mat,mask_mat=NULL,i,j,max_mask){
  k = klast <- 0
  min_points <- 40
  max_points <- 400
  while(T){
    dim_count <- dim(count_mat)
    
    i_min <- i-k
    i_max <- i+k
    j_min <- j-k
    j_max <- j+k
    if(i_min < 1) i_min <- 1
    if(i_max > dim_count[1]) i_max <- dim_count[1]
    if(j_min < 1) j_min <- 1
    if(j_max > dim_count[2]) j_max <- dim_count[2]
    
    my_mask <- max_mask[(1001-(i-i_min)):(1001 + (i_max-i)),
                        (1001-(j-j_min)):(1001 + (j_max-j))] <= k
    my_mask[!my_mask] <- NA
    
    sample_count <- count_mat[i_min:i_max,j_min:j_max]*my_mask[]
    sample_count_noNAs <- sample_count[!is.na(sample_count)]
    sum_notNAs <- length(sample_count_noNAs)
    
    if((sum_notNAs>=2)&(sum(sample_count_noNAs)>1)) {
      csr_test <- chisq.test(sample_count_noNAs,simulate.p.value=F)
      pval <- csr_test$p.value
    }
    else {
      pval <- 1 
    }
    if(pval <= 0.001) break
    if((sum(sample_count_noNAs) > max_points)| (k>=100)) {
      break
    }
    klast <- k
    k <- k + 1
  }
  
  if(sum_notNAs<=0) {
    dens <- NA
    coef_var <- NA
  }
  else {
    dens <- mean(sample_count_noNAs)
    coef_var <- 1/sqrt(sum(sample_count_noNAs))
  }
  return(c(dens,coef_var,klast))
}

robust.density<-function(point_pattern,pointsAsMask=NULL,res_y,res_x){
  
  max_mask <- matrix(nrow=2*1000+1,ncol=2*1000+1)
  for(i in 1:(2*1000+1)) {
    for(j in 1:(2*1000+1)){
      max_mask[i,j] <- sqrt((i-1001)^2 + (j-1001)^2)
    }
  }
  
  if(is.null(pointsAsMask)) {
    max_x <- max(point_pattern$x)
    min_x <- min(point_pattern$x)
    max_y <- max(point_pattern$y)
    min_y <- min(point_pattern$y)
  } else {
    max_x <- max(pointsAsMask$x)
    min_x <- min(pointsAsMask$x)
    max_y <- max(pointsAsMask$y)
    min_y <- min(pointsAsMask$y)
  }
  
  n_col <- ceiling((max_x-min_x)/res_x)
  n_row <- ceiling((max_y-min_y)/res_y)
  
  max_x <- min_x+(n_col)*res_x
  max_y <- min_y+(n_row)*res_y
  
  Wall <- c(min_x,max_x,min_y,max_y)
  count_mat <- mask_mat <- dens_mat <- coefvar_mat <- k_mat <- matrix(nrow=n_row,ncol=n_col)
  count_mat[] <- as.matrix(spatstat.geom::quadratcount(spatstat.geom::as.ppp(cbind(point_pattern$x,point_pattern$y),Wall),n_col,n_row))
  if(is.null(pointsAsMask)) {
    mask_mat <- NULL
  } else {
    mask_mat[] <- as.matrix(spatstat.geom::quadratcount(spatstat.geom::as.ppp(cbind(pointsAsMask$x,pointsAsMask$y),Wall),n_col,n_row))
    count_mat[mask_mat <= 0] <- NA
  }
  
  for(i in 1:n_row){
    print(100*i/n_row)
    
    for(j in 1:n_col){
      x <- min_x+res_x*(j-0.5)
      y <- min_y+res_y*(i-0.5)
      est <- robust.density.unit(count_mat,mask_mat,i,j,max_mask)
      dens_mat[i,j] <- est[1]
      coefvar_mat[i,j] <- est[2]
      k_mat[i,j] <- est[3]
    }
  }
  if(!is.null(pointsAsMask)) {
    dens_mat[mask_mat <= 0] <- NA
    coefvar_mat[mask_mat <= 0] <- NA
    k_mat[mask_mat <= 0] <- NA
  }
  return(list(dens_mat,coefvar_mat,k_mat))
}

