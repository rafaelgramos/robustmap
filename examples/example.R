library(robustmap)
# Loading point data
burglary <- sf::read_sf(dsn = "data", layer = "burglary")
burglary <- data.frame(lon=burglary$lon_m,
                       lat=burglary$lat_m)
# Estimating optimal granularity using robust.quadcount
burglary_map <- robust.quadcount(burglary,verbose = T)

# Retriving estimated granularity
burglary_map$opt_granularity

# Plotting resulting map
terra::plot(burglary_map$counts)
