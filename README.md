# Variability-of-POC-Stock-in-the-South-China-Sea-with-reference-to-Profile-Structures
# POC Analysis in South China Sea

MATLAB-based analysis of Particulate Organic Carbon (POC) stocks and trends in the South China Sea using satellite-derived cp660 data.

## 📋 Overview

This repository contains two main MATLAB scripts for analyzing POC distribution and long-term trends in the South China Sea:

1. **Profile_POC_monthly_and_climatological_and_rate_estimate.m**: Calculates monthly and climatological POC profiles in the upper 125m layer, and estimate the changing rates (giving the p value at the same time) during 2003-2023 period.
2. **POC_stock_Monthly_and_Climatological_and_Rate_estimate**: Calculates monthly and climatological POC stocks in the upper 100m and euphotic layer, and estimate the changing rates (giving the p value at the same time) during 2003-2023 period.

## 🎯 Features

- **Multi-depth Analysis**: 8 depth levels (5, 10, 20, 30, 50, 75, 95, 115m)
- **Long-term Trends**: 21-year time series analysis (2002-2023)
- **Statistical Robustness**: Outlier removal and significance testing
- **Climatology Calculation**: Monthly climatological averages
- **Trend Analysis**: Linear regression with p-value significance

## 📊 Data Requirements

### Input Data
- **3D cp660 profiles**: NetCDF files containing cp660 data (`cp_profile_optcp4_NN`) could be found at (https://doi.org/10.5281/zenodo.16957344)
- **Bathymetry data**: `SCSdepth9km.mat` for South China Sea depth information    afforded in this repository
- **Region coordinates**: `SCSLatLon.mat` for study area definition     afforded in this repository
- **Chlorophyll data**: MODIS satellite chlorophyll-a products   could be found at (https://oceancolor.gsfc.nasa.gov.)

### Study Area
- Northern South China Sea (15°N-25.125°N, 105°E-121.5°E)
- Grid resolution: 122×199 points
