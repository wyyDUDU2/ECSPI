CSP-ECSPI Cotton Early Mapping - GEE JavaScript
Project Overview
Google Earth Engine interactive UI code for **Enhanced Cumulative Spectral Phenology Index (ECSPI)** early cotton mapping in Xinjiang (Aksu, Kashgar, Tacheng-Changji). This implementation fully reproduces the paper algorithm, including rainfall noise correction, pixel-wise adaptive phenology extraction, multi-source Sentinel-2 & MODIS fusion, 4 comparison classification pipelines and automatic accuracy evaluation.

Core Functions
1. Rainfall Adaptive Noise Correction (R1-R4)
   Use CHIRPS precipitation data to screen cloud/rain contaminated observations, reconstruct Sentinel-2 SI time series guided by MODIS NDVI to eliminate spectral distortion.
2. 4 Classification Methods (M1~M4)
   M1: Single Sentinel-2 CSP (unsupervised Otsu-Sauvola segmentation)
   M2: Multi-source fused ECSPI (proposed core algorithm)
   M3: Random Forest with full multi-index temporal & phenology features
   M4: Random Forest combined ECSPI + NDVI phenology features
3. Per-pixel Adaptive Phenology Detection
   Savitzky-Golay smoothing + second derivative calculation to auto extract growth start T1 & peak T2 for every pixel, adapting to inter-regional sowing differences.
4. Visualization & Evaluation Tools
   One-click temporal curve plot by map click, CSP distribution histogram, class separability table (M/JM index), binary cotton map, false-color composite, OA/Kappa/F1 accuracy statistics.
5. Asset Import & Export
   Support tiled export of feature images and classification results to GEE Assets to avoid computation overflow.

Study Regions & Input Data
Research zones: Aksu, Kashgar, Tacheng-Changji (switch via UI dropdown)
Input datasets: Sentinel-2 L2A, MOD13Q1 NDVI, CHIRPS precipitation, ESA WorldCover, SRTM DEM, field sample point assets
Sample configuration: 70% training / 30% test split, cotton & 5 non-crop categories (corn, wheat, orchard, other crop, non-agricultural land)

Usage Guide
1. Open Google Earth Engine Code Editor, paste full script
2. Pre-requirement: Upload cotton/non-cotton sample point FeatureCollections to your GEE Asset path (`projects/wangyiyao/assets/`)
3. Select study area & target year from left UI panel
4. Optional: Run R1~R4 to execute rainfall denoising module
5. Run corresponding method step-by-step following UI button order
6. Use evaluation panel to generate charts and accuracy metrics
7. Export classification results or intermediate CSP/feature layers to Assets

Algorithm Advantage
The proposed ECSPI (M2) significantly improves separability between cotton and other crops. It realizes cotton identification at squaring stage, 60–70 days earlier than traditional boll-opening period mapping. Unsupervised classification OA reaches 83.03~92.32%,and random forest accuracy up to 96.00~97.77% across 2019–2024 multi-year validation.

File Structure
Global parameter configuration: study bounds, scale, phenology DOY defaults
 Mask & sample initialization: cropland + elevation filter, train/test split
Core utility functions: cloud masking, index calculation, SG smoothing, T1/T2 detection, CSP & ECSPI computation, Otsu/Sauvola threshold, RF training, accuracy assessment
Rainfall correction pipeline (R series buttons)
Four classification UI workflows (M1/M2/M3/M4)
Visualization & chart tools
Map click time series plot callback
Side UI panel & map legend

Citation
Zhang Q, Wang Y, Dong T, et al. Early-Season Cotton Identification in Xinjiang Based on Multi-Source Time-Series Remote Sensing Data Reconstruction and Spectral Phenology Index Enhancement.
