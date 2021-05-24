#Robustmap package#

##Overview##

This package implement routines for mapping point data set via quadrat counting and kernel density estimation while taking into account the robustness and internal (spatial) uniformity of the aggregated counts. Too coarse of a quadrat size (or bandwidth) may lead to lack of detail and to the Ecological Fallacy; too fine of a quadrat size (or bandwidth), however, will generate unreliable counts and densities.

The **robust.quadcount** routine estimates a quadrat size that balances both robustness and uniformity given a point set, returning the 'optimal' quadrat count along with other informations such as: which granularies (quadrat sizes) were tested, the robustness and uniformities at each granularity, etc. The methodology is based on the work of Ramos, et al. (2021) [https://link.springer.com/article/10.1007/s10940-020-09474-6]

The **robust.density routine** is still in development and should estimate the density of a point pattern via a variable bandwidth method, in which the bandwidth is the largest possible such that the point subset contained by it displays Complete Spatial Randomness. This way, for any given pixel, the number of points considered in the density estimation can be maximized (improving robustness) while not incurring into the Ecological Fallacy of averaging over an area with varying point density.

##Installing##

On Mac (has not been tested on Linux nor Windows, but should workt too), run:

```R

```