# Image visual features
Here is a brief description of each covariate related to image feature:
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

## Detailed description
<ins>Entropy</ins>
The 'entropy' built-in Matlab function is defined as -sum(p.*log2(p)) where p contains the histogram counts returned from IMHIST function.

<ins>Contrast</ins>
Contrast is computed through the 'graycoprops' Matlab function. It is defined as the intensity contrast between a pixel and its neighbor over the whole image. Range = [0 (size(GLCM,1)-1)^2].  Contrast is 0 for a constant image.

<ins>Correlation</ins>
Correlation is computed through the 'graycoprops' Matlab function. It is the statistical measure of how correlated a pixel is to its neighbor over the whole image. 
Range = [-1 1]. Correlation is 1 or -1 for a perfectly positively or negatively correlated image. Correlation is NaN for a constant image.
