# Image visual features description
Here is the description of each covariate related to image feature:
- Entropy:                    Minimal number of bits required to encode the image                             
- Contrast:                   Difference in luminance of the image                                            
- Correlation:                How correlated neighboring pixels are                                           
- Homogeneity:                How close pixel values are from the mean pixel value                            
- Energy:                     Measure of the localized change of the image                                     
- Compactness:                How closely packed the pixels of the item are                                    
- Ratio:                      Length width ratio of the item                                                   
- \# Spectral clusters:       The variety of frequencies in the image                                          
- High frequency energy:      Energy of spectral cluster with the highest frequency                            
- Highest frequency:          Centroid of the spectral cluster of highest frequency                            
- Maximum spectral distance:  Distance between spectral clusters of lowest and highest frequency

All the computations are discribed in the [extract_cov.mat](https://github.com/numediart/Covariates_Analysis/blob/main/utils/extract_cov.m) function.
