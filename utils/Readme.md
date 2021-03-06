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
Note: Contrast, correlation, homogeneity and energy are computed through the 'graycoprops' Matlab function applied on the grayscale image. This function normalizes the gray-level co-occurrence matrix (GLCM) so that the sum of its elements is one.

<ins>Entropy</ins>  
The 'entropy' built-in Matlab function is defined as -sum(p.*log2(p)) where p contains the histogram counts returned from IMHIST function.

<ins>Contrast</ins>  
The contrast is defined as the intensity contrast between a pixel and its neighbor over the whole image. Range = [0 (size(GLCM,1)-1)^2].  Contrast is 0 for a constant image.

<ins>Correlation</ins>  
The correlation is the statistical measure of how correlated a pixel is to its neighbor over the whole image. 
Range = [-1 1]. Correlation is 1 or -1 for a perfectly positively or negatively correlated image. Correlation is NaN for a constant image.

<ins>Homogeneity</ins>  
The homogeneity is the closeness of the distribution of elements in the GLCM to the GLCM diagonal. Range = [0 1]. Homogeneity is 1 for a diagonal GLCM.

<ins>Energy</ins>  
The energy is the summation of squared elements in the GLCM. Range = [0 1]. Energy is 1 for a constant image.

<ins>Compactness</ins>  
The compactness of a specific shape is defined as the ratio of the area of the shape to the area of a circle (the most compact shape) having the same perimeter. It is computed as follow:  
```
compactness = 4*pi*A/P^2;
```
where A is the area of the shape and P is its perimeter (here, the area of the circle is computed from its known perimeter.

