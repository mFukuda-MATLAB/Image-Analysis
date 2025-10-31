# Image-Analysis
MATLAB code for analyzing microscopic images using multiple channels with ResNet (U-Net & PCA are additional).

## Overview
This repository contains MATLAB scripts for:
1. Cell Segmentation
2. ResNet Training / Application
3. Mean Fluorescence Intensity Measurement

The code integrates convolution-based preprocessing, watershed segmentation, and ResNet classification to quantify fluorescence signals and particle uptake in macrophage images. U-net and PCA are additional code for image analysis.

## Features
- Automated image segmentation using convolution + watershed
- ResNet 3-classification (Single, Multiple, and None) training and application
- Particle detection via circle fitting ('imfindcircles')
- Quantitative mean fluorescence intensity measurement
- Export selected images (processed) and summary data to Excel

## Workflow
1. ** Image Separation **
- Microscopic image (.tiff): slices (red, green, blue, and phase / bright-field) are separated by ImageJ / Fiji
- Store separated images in each folder
2. ** Image Prepeocessing **
- Contrast enhancement and noise filtering (bright-field)
- Edge detection (convolution)
- Create masks for each channel (red, green, blue, bright-field)
3. ** Image Selection / segmentation **
- Isolated single cell image: selection by blue channel
- Aggregated cell images: segmentation by watershed (blue & bright-filed masks)
- 
