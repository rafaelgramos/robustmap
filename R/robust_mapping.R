#' Optimal granularity for quadrat counting
#' 
#' Given a point set, finds an optimal granularity for quadrat counting that
#' balances uniformity and robustness, as described in Ramos et al. (2021).
#'
#' @param point_set: the point set for which you want to find the optimal quadrat
#' size
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
#' @param verbose: whether to print messages while running the function or not
#' @param my_scales: which granularities to be tested. If not provided, a set is
#' automatically generated. The granularities, if provided, should be given as the
#' dimension (side) of a cell, in the unit used for the coordinates of point_set.
#' @param W: the window of interest to be considered when doing the analyzis of
#' point_set. If not provided, W is calculated as the minimum bounding rectangle
#' for point_set.
#' @keywords quadrat granularity adequate optimal
#' @author Rafael G. Ramos (main proponent and coder), Marcos Prates (contributor)
#' @references Ramos, R. G., Silva, B. F., Clarke, K. C., & Prates, M. (2021). 
#' Too Fine to be Good? Issues of Granularity, Uniformity and Error in Spatial 
#' Crime Analysis. Journal of Quantitative Criminology, 1-25.
#' robust.quadcount()
#' @import terra
#' @export
#' 
#' @examples
#' library(robustmap)
#' # Loading point data
#' burglary <- sf::read_sf(dsn = "data", layer = "burglary")
#' burglary <- data.frame(lon=burglary$lon_m,
#'                        lat=burglary$lat_m)
#' # Estimating optimal granularity using robust.quadcount
#' burglary_map <- robust.quadcount(burglary,verbose = T)
#' 
#' # Retriving estimated granularity
#' burglary_map$opt_granularity
#' 
#' # Plotting resulting map
#' terra::plot(burglary_map$counts)
#' 

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
                           W=NULL){
  
  tradeoff_crit <- match.arg(tradeoff_crit)
  uniformity_method <- match.arg(uniformity_method)
  robustness_method <- match.arg(robustness_method)
  
  if(is.null(W)) {
    W <- c(min(point_set[,1]),max(point_set[,1]),min(point_set[,2]),max(point_set[,2]))
  }
  
  #HOW is this initialized?
  if(is.null(my_scales)) {
    #maxscale <- 0.10*min(c(W[2] - W[1],W[4]-W[3]))
    #minscale <- maxscale*0.01
    #my_scales <- seq(minscale,maxscale,minscale)
    my_scales <- seq(25,1650,25)
  }
  
  # calculating robustness and uniformity for each granularity in 'my_scales'
  my_spatialstats <- get_spatialstats_all(point_set,
                                          my_scales,
                                          my_scales,
                                          random_samples=random_samples,
                                          signif=signif,
                                          uniformity_method=uniformity_method,
                                          robustness_method=robustness_method,
                                          robustness_k = robustness_k,
                                          verbose=verbose)
  
  if(verbose) {
    print("Estimating optimal balance between uniformity and robustness")
  }
  robust <- my_spatialstats$robustness
  unif <- my_spatialstats$uniformity
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
  nlon <- floor((W[2] - W[1])/opt_granularity)
  nlat <- floor((W[4] - W[3])/opt_granularity)
  q <- spatstat.geom::quadratcount(spatstat.geom::as.ppp(point_set,W),nlon,nlat)
  map$counts <- terra::rast(matrix(data=q[],nrow=nrow(q),ncol=ncol(q)))
  terra::ext(map$counts) <- terra::ext(W)
  map$opt_granularity <- opt_granularity
  final_sample <- create_samples(point_set=point_set,
                                 random_samples = F,
                                 window_w = opt_granularity,
                                 window_h = opt_granularity)
  stats_final_sample <- get_spatialstats_sample(sample_set = final_sample,
                                                uniformity_method = uniformity_method,
                                                signif=signif,
                                                robustness_method = robustness_method)
  tmp <- matrix(data=stats_final_sample$samples_covar,nrow=nrow(map$counts),ncol=ncol(map$counts))
  map$covar <- terra::flip(terra::rast(tmp),direction="horizontal")
  terra::ext(map$covar) <- terra::ext(map$counts)
  tmp <- matrix(data=stats_final_sample$samples_csr_pass,nrow=nrow(map$counts),ncol=ncol(map$counts))
  map$is.csr <- !(terra::flip(terra::rast(tmp),direction="horizontal"))
  terra::ext(map$is.csr) <- terra::ext(map$counts)
  robust <- my_spatialstats$robustness
  unif <- my_spatialstats$uniformity
  map$uniformity.curve <- unif
  map$granularities.tested <- my_scales
  map$robustness.curve <- robust
  return(map)
}


#
#  Function that estimates the internal uniformity and robustness to error for 'point_set',
#  according to the granularity dimensions listed in 'scales_lon' and 'scales_lat' (these should be the
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
                               scales_lon,
                               scales_lat,
                               random_samples=T,
                               nsamples=500,
                               signif=0.90,
                               robustness_method=c("Poisson","Binomial","Resampling"),
                               robustness_k = -3,
                               uniformity_method=c("Quadratcount","Nearest-neighbor"),
                               verbose=verbose){
  
  robustness_method <- match.arg(robustness_method)
  uniformity_method <- match.arg(uniformity_method)
  
  if(robustness_k >= 0) {
    print("robust_k parameter must be less than zero. See documentation for details.")
    return(NULL)
  }
  
  len_scales <- min(c(length(scales_lon),length(scales_lat)))
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
      cat('\014',append=T)
      cat("Estimating uniformity and robustness at various scales:\n")
      cat("[")
      for(k in 1:i){
        cat("=")
      }
      if(i < len_scales) {
        for(k in (i+1):len_scales){
          cat(" ")
          }
      }
      cat("] | ")
      cat(paste(round(100*(i/len_scales)),"%\n",sep=""),append=T)
      Sys.sleep(.05)
    }
    
    # generating set of quadrats of the current granularity for taking the samples
    my_samples <- create_samples(point_set=point_set,
                                 random_samples = random_samples,
                                 nsamples = nsamples,
                                 window_w = scales_lon[i],
                                 window_h = scales_lat[i])
    
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

create_samples<-function(point_set,random_samples,nsamples,window_w,window_h) {
  max_b_lon <- max(point_set$lon)
  max_b_lat <- max(point_set$lat)
  min_b_lon <- min(point_set$lon)
  min_b_lat <- min(point_set$lat)
  w_b <- max_b_lon-min_b_lon
  h_b <- max_b_lat-min_b_lat
  
  if(random_samples == T) {
    # quadrats taken randomly
    offset_lon <- ((w_b-window_w)*runif(nsamples))+min_b_lon
    offset_lat <- ((h_b-window_h)*runif(nsamples))+min_b_lat
  }
  else {
    # contiguous quadrats
    nlon <- floor(w_b/window_w)
    nlat <- floor(h_b/window_h)
    nsamples <- nlon*nlat
    offset_lon <- matrix(nrow=nlat,ncol=nlon)
    offset_lat <- matrix(nrow=nlat,ncol=nlon)
    for(k in 1:nlat) offset_lon[k,] <- (0:(nlon-1))*window_w+min_b_lon
    for(k in 1:nlon) offset_lat[,k] <- (0:(nlat-1))*window_h+min_b_lat
    offset_lon <- as.vector(offset_lon)
    offset_lat <- as.vector(offset_lat)
  }
  out <- list()
  out$point_set <- point_set
  out$offset_lon <- offset_lon
  out$offset_lat <- offset_lat
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
  
  nsamples <- length(sample_set$offset_lon)
  
  # allocating temporary structures
  samples_count <- vector(length=nsamples)
  samples_csr_p <- vector(length=nsamples)
  samples_csr_pass <- vector(length=nsamples)
  samples_covar <- vector(length=nsamples)
  samples_zeros <- 0
  
  # for each quadrat, take a the points inside and test for robustness and internal uniformity
  for(j in 1:nsamples) {
    
    # extracting the points inside the quadrat
    W_ext = c(sample_set$offset_lon[j],
              sample_set$offset_lon[j]+sample_set$window_w,
              sample_set$offset_lat[j],
              sample_set$offset_lat[j]+sample_set$window_h)
    
    sub_points <- sample_set$point_set[spatstat.geom::inside.owin(x=sample_set$point_set$lon,
                                                                  y=sample_set$point_set$lat,
                                                                  w=W_ext),]
    
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
      
      # allocating a slightly larger quadrat just to avoid that our sampled points to fall
      # exactly in the boundary of the quadrat.
      W_disp <- c(sample_set$offset_lon[j]-0.0001,
                  sample_set$offset_lon[j]+sample_set$window_w+0.0001,
                  sample_set$offset_lat[j]-0.0001,
                  sample_set$offset_lat[j]+sample_set$window_h+0.0001)
      
      # testing Complete Spatial Randomness (CSR) with Clark-Evans nearest neighbor test (part of estimating uniformity)
      samples_count[j] <- nrow(sub_points)
      if(uniformity_method == "Nearest-neighbor") {
        if(samples_count[j] > 20) {
          my_clarkevans <- spatstat.explore::clarkevans.test(spatstat.geom::as.ppp(sub_points,W_disp),alternative="two.sided")#,correction = "Donnelly")
        } else {
          my_clarkevans <- spatstat.explore::clarkevans.test(spatstat.geom::as.ppp(sub_points,W_disp),alternative="two.sided",nsim=100)#,correction = "Donnelly")
        }
        samples_csr_p[j] <-  my_clarkevans$p.value
        # assign T to quadrats that pass the threshold 'signif' for CSR
        samples_csr_pass[j] <-  my_clarkevans$p.value < (1-signif)
      }
      else if(uniformity_method == "Quadratcount") {
        # test CST with the quadrat test
        my_quadrattest <- spatstat.explore::quadrat.test(spatstat.geom::as.ppp(sub_points,W_disp),nx=5)
        samples_csr_p[j] <- my_quadrattest$p.value
        samples_csr_pass[j] <- my_quadrattest$p.value < (1-signif)
      }
      else {
        samples_csr_p[j] <- NA
        samples_csr_pass[j] <- NA
      }
      samples_covar[j] <- calc_covar(nrow(sub_points),sample_set$npopulation,robustness_method)
      # Calculating robustness.
      # First we calculate the expected coefficient of variation for the samples, using different methods
      #prob_event = nrow(sub_points)/sample_set$npopulation
      #if(robustness == "Binomial") {
      #  # calculating coef of var using the Binomial estimation method (see paper)
      #  samples_covar[j] = sqrt(prob_event*(1-prob_event)*sample_set$npopulation)/nrow(sub_points)#sd(sim_rates)/mean(sim_rates)
      #}
      #else if(robustness == "Poisson") {
      #  # calculating coef of var using the Poisson estimation method (see paper)
      #  tmp = fitdistr(nrow(sub_points),"Poisson")
      #  samples_covar[j] = tmp$sd/tmp$estimate
      #}
      #else if(robustness == "Resampling") {
      #  # calculating coef of var using the resampling estimation method (see paper)
      #  sim_rates = rbinom(1000,sample_set$npopulation,prob_event)
      #  samples_covar[j] = mean(sqrt(1/sim_rates[sim_rates!=0]))
      #}
      
    }
  }
  out <- list()
  out$samples_count <- samples_count
  out$samples_csr_p <- samples_csr_p
  out$samples_csr_pass <- samples_csr_pass
  out$samples_covar <- samples_covar
  out$samples_zeros <- samples_zeros
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





########################




#
#   Utility functions, for ploting quadrat counts etc.
#

quad2mat<-function(quadset) {
  tmp <- matrix(data=quadset[],nrow=nrow(quadset),ncol=ncol(quadset))
  #tmp2 = apply(tmp, 2, rev)
  #tmp3 = t(tmp2)
  return(tmp)
}

points2quad<-function(point_set,my_scale,mask=NULL){
  max_b_lon <- max(point_set$lon)
  max_b_lat <- max(point_set$lat)
  min_b_lon <- min(point_set$lon)
  min_b_lat <- min(point_set$lat)
  W <- c(min(point_set$lon),
         max(point_set$lon),
         min(point_set$lat),
         max(point_set$lat))
  nlon <- floor((max_b_lon - min_b_lon)/my_scale)
  nlat <- floor((max_b_lat - min_b_lat)/my_scale)
  my_ret <- spatstat.geom::quadratcount(spatstat.geom::as.ppp(point_set,W),nlon,nlat)
  
  if(!is.null(mask)){
    Wmask = c(min(mask$lon),
              max(mask$lon),
              min(mask$lat),
              max(mask$lat))
    my_mask <- spatstat.geom::quadratcount(spatstat.geom::as.ppp(mask,Wmask),nlon,nlat) 
    my_ret[my_mask <= 0] = NA
  }
  
  return(my_ret)
}


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
    
    my_mask <- max_mask[(1001-(i-i_min)):(1001 + (i_max-i)),(1001-(j-j_min)):(1001 + (j_max-j))] <= k
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

robust.density<-function(point_pattern,pointsAsMask=NULL,res_lat,res_lon){
  
  max_mask <- matrix(nrow=2*1000+1,ncol=2*1000+1)
  for(i in 1:(2*1000+1)) {
    for(j in 1:(2*1000+1)){
      max_mask[i,j] <- sqrt((i-1001)^2 + (j-1001)^2)
    }
  }
  
  if(is.null(pointsAsMask)) {
    max_lon <- max(point_pattern$lon)
    min_lon <- min(point_pattern$lon)
    max_lat <- max(point_pattern$lat)
    min_lat <- min(point_pattern$lat)
  } else {
    max_lon <- max(pointsAsMask$lon)
    min_lon <- min(pointsAsMask$lon)
    max_lat <- max(pointsAsMask$lat)
    min_lat <- min(pointsAsMask$lat)
  }
  
  n_col <- ceiling((max_lon-min_lon)/res_lon)
  n_row <- ceiling((max_lat-min_lat)/res_lat)
  
  max_lon <- min_lon+(n_col)*res_lon
  max_lat <- min_lat+(n_row)*res_lat
  
  Wall <- c(min_lon,max_lon,min_lat,max_lat)
  count_mat <- mask_mat <- dens_mat <- coefvar_mat <- k_mat <- matrix(nrow=n_row,ncol=n_col)
  count_mat[] <- as.matrix(spatstat.geom::quadratcount(spatstat.geom::as.ppp(cbind(point_pattern$lon,point_pattern$lat),Wall),n_col,n_row))
  if(is.null(pointsAsMask)) {
    mask_mat <- NULL
  } else {
    mask_mat[] <- as.matrix(spatstat.geom::quadratcount(spatstat.geom::as.ppp(cbind(pointsAsMask$lon,pointsAsMask$lat),Wall),n_col,n_row))
    count_mat[mask_mat <= 0] <- NA
  }
  
  for(i in 1:n_row){
    print(100*i/n_row)
    
    for(j in 1:n_col){
      #print(i)
      #print(n_row)
      #print(j)
      #print(n_col)
      x <- min_lon+res_lon*(j-0.5)
      y <- min_lat+res_lat*(i-0.5)
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