# 🧠 Image-Analysis
MATLAB code for analyzing microscopic images using multiple channels with **ResNet** (U-Net & PCA are additional modules).

---

## 📘 Overview
This repository contains MATLAB scripts for:
  1. **Cell Segmentation**
  2. **ResNet Training / Application**
  3. **Cell Image Selection by Red Channel** (based on red fluorescence)
  4. **Mean Fluorescence Intensity Measurement**

The code integrates convolution-based preprocessing, watershed segmentation, and ResNet classification to quantify fluorescence signals and particle uptake in macrophage images.  
U-Net and PCA modules are included as additional tools for image analysis.

---

## 🧬 Fluorescence Channels

  ### 🟢 Green Channel
  - **DNA** within particles labeled with *DiYO-1*  
  - **Histone** within particles labeled with *AF488*
  
  ### 🔴 Red Channel
  - **PNIPAM(nps):** red fluorescent nanoparticles incorporated into PNIPAM layer  
  - **(PNIPAM&PAH-594):** PAH conjugated with AF594 and incorporated into PNIPAM layer
  
  ### 🔵 Blue Channel
  - **Hoechst** dye used for labeling nuclear DNA

---

## ⚙️ Features
  - Automated image segmentation using convolution + watershed  
  - ResNet 3-class classification (**Single**, **Multiple**, **None**)  
  - Particle detection using `imfindcircles` *(for images without red channel)*  
  - Quantitative mean fluorescence intensity measurement  
  - Export of processed images and summary data to Excel  

---

## 🔬 Workflow

### 1. **Image Separation**
  - Multi-channel `.tif` images are decomposed into red, green, blue, and bright-field channels using **ImageJ/Fiji**.  
  - Each separated image is stored in its corresponding folder (e.g., `Red/`, `Green/`, `Blue/`, `Bright-field/`).

---

### 2. **Image Preprocessing**
  - Bright-field images: contrast enhancement and noise filtering  
  - Edge detection via convolution filters (Sobel, Gaussian, Prewitt)  
  - Mask generation for each channel (red, green, blue, bright-field)

---

### 3. **Image Selection and Segmentation**
  - **Isolated single-cell images:** selected using blue-channel mask  
  - **Aggregated cell images:** segmented using watershed on combined blue and bright-field masks  

  #### a. Without Red Channel
  - ResNet-18 3-class model trained using PNIPAM(DNA-DiYO-1) microparticle images from **RAW264.7 macrophages**  
  - The trained model is applied to different experimental conditions and particle types  

  #### b. With Red Channel
  - Images containing red fluorescence are filtered by red-channel intensity or mask presence  
  - The filtered dataset is processed using the same segmentation pipeline  
  - Used primarily to select cell images containing phagocytosed single particles
  - Used to analyze **co-localization** and **mean fluorescence intensity** within **phagosome**

---

### 4. **Cell Image Selection**
  #### a. Without Red Channel
  - `imfindcircles` is used to identify cell images containing phagocytosed single particles from those selected by the ResNet model  
  - Detected circles are filtered by their area, intensity, and solidity
  
  #### b. With Red Channel
  - Images containing red fluorescence are filtered by red-channel intensity or mask presence  
  - The filtered dataset is analyzed to evaluate **co-localization** and **mean fluorescence intensity** within **phagosome**

---

### 5. **Mean Fluorescence Intensity Measurement**
  - Green channel mean fluorescence intensities were measured from selected ROI: phagosome or cytosol (depending on assay or analysis types)

## 🧪 Requirements
- MATLAB R2022a or newer  
- Image Processing Toolbox  
- Deep Learning Toolbox  
- (Optional) Statistics and Machine Learning Toolbox
- ImafeJ / Fiji software
