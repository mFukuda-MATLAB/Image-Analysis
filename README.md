## Image-Analysis
MATLAB code for analyzing microscopic images using multiple channels with **ResNet** (U-Net & PCA are additional modules).

## Overview
This repository contains MATLAB scripts for:
  1. Cell Segmentation
  2. ResNet Training / Application
  3. Cell image selection by red channel (if there is red fluorescence)
  4. Mean Fluorescence Intensity Measurement

The code integrates convolution-based preprocessing, watershed segmentation, and ResNet classification to quantify fluorescence signals and particle uptake in macrophage images.  
U-Net and PCA modules are included as additional tools for image analysis.

---

## Features
  - Automated image segmentation using convolution + watershed  
  - ResNet 3-class classification (**Single**, **Multiple**, **None**)  
  - Particle detection using `imfindcircles`  
  - Quantitative mean fluorescence intensity measurement  
  - Export of processed images and summary data to Excel  

## Workflow

### 1. **Image Separation**
  - Multi-channel `.tif` images are decomposed into red, green, blue, and bright-field channels using **ImageJ/Fiji**.  
  - Each separated image is stored in its corresponding folder (e.g., `Red/`, `Green/`, `Blue/`, `BF/`).

---

### 2. **Image Preprocessing**
  - Bright-field images: contrast enhancement and noise filtering  
  - Edge detection via convolution filters (Sobel, Gaussian, Prewitt)  
  - Mask generation for each channel (red, green, blue, bright-field)

---

### 3. **Image Selection and Segmentation**
  - **Isolated single-cell images:** selected by blue-channel mask  
  - **Aggregated cell images:** segmented using watershed on combined blue and bright-field masks  

#### a. Without Red Channel
  - ResNet (ResNet-18) 3-class model trained using PNIPAM(DNA-DiYO-1) microparticle images from **RAW264.7 macrophages**  
  - The trained model is applied to different experimental conditions and particle types  

#### b. With Red Channel
  - Images containing red fluorescence are first filtered by red-channel intensity or mask presence  
  - The filtered dataset is processed using the same segmentation and classification pipeline  
  - Used primarily to analyze **co-localization** or **dual-labeling** conditions  
