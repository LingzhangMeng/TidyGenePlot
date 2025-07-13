# TidyGenePlot

## Description
Single-cell RNA sequencing (scRNA-seq) analysis using the R package Seurat often involves visualizing gene expression patterns through violin plots and feature plots. However, the default `VlnPlot()` and `FeaturePlot()` functions in Seurat include elements such as the "Identity" label and legends in violin plots, as well as coordinate axes in feature plots, which can clutter the output. The TidyGenePlot package addresses these issues by generating clean and streamlined violin and feature plots, removing unnecessary elements to provide researchers with clear and professional visualizations.

Additionally, TidyGenePlot resolves a bug in Seurat’s `VlnPlot()` and `FeaturePlot()` functions, which fail to generate plots for certain genes in scaled data, making it a valuable tool for scRNA-seq analysis.

## Installation
To install TidyGenePlot from GitHub, use the following commands:

```R
library(devtools)
devtools::install_github("LingzhangMeng/TidyGenePlot")
library(TidyGenePlot)
```

## Dependencies
TidyGenePlot requires the following R packages:
- Seurat (version 3.1.2 or higher)
- patchwork (version 1.0.1 or higher)
- ggplot2 (version 3.5.1 or higher)

## User Guide
To use TidyGenePlot, load the required packages:

```R
library(TidyGenePlot)
library(Seurat)
library(ggplot2)
library(patchwork)
```

### Function 1: `tidy.VlnPlot`
Generates clean violin plots without the "Identity" label or legend.

```R
feature.list <- c("Krt14", "Krt5", "Il1b", "Cd34", "Fxyd3", "Sfn")
tidy.VlnPlot(seu_obj, features = feature.list, pt.size = 0, ncol = 3)
```
**Example output** 
<img width="923" height="554" alt="Weixin Image_20250713114713" src="https://github.com/user-attachments/assets/2deeda09-1f32-47d1-8f4e-58818df44ae6" />



#### Parameters
- `seu_obj`: A Seurat object created using the Seurat package.
- `feature.list`: A character vector specifying the genes to plot.
- `pt.size`: Point size for the plot (numeric, ranging from 0 to 2, e.g., 0, 0.5, 1, or 2).
- `cols`: Color palette for clusters. Set to `NULL` to use Seurat’s default colors, or provide a custom vector of colors (length must be equal to or greater than the number of clusters).
- `ncol`: Number of plots per row (e.g., 1, 2, 3, 4, 5, or 6).

### Function 2: `tidy.FeaturePlot`
Generates clean feature plots without coordinate axes, with an option to include or exclude a legend.

```R
tidy.FeaturePlot(seu_obj, features = feature.list, pt.size = pointsize, cols = colors, Legend = FALSE, ncol = n)
```
**Example output** 

<img width="770" height="453" alt="4" src="https://github.com/user-attachments/assets/10e470ac-6f92-46bb-95c8-66c27aedcb7d" />


#### Parameters
- `seu_obj`: A Seurat object created using the Seurat package.
- `feature.list`: A character vector specifying the genes to plot.
- `pt.size`: Point size for the plot (numeric, ranging from 0 to 2, e.g., 0.5, 1, or 2).
- `cols`: A vector of two colors. The first color represents the background, and the second represents gene expression levels. Recommended settings include `c("grey", "blue")` or `c("grey", "red")`. Do not set to `NULL`.
- `Legend`: Logical (`TRUE` or `FALSE`) to determine whether a legend is displayed.
- `ncol`: Number of plots per row (e.g., 1, 2, 3, 4, 5, or 6).

## Examples and Comparison
Below are examples comparing Seurat’s default plotting functions with TidyGenePlot’s enhanced visualizations:

```R
# Define a color palette
cb_palette <- c("#ed1299", "#09f9f5", "#246b93", "#cc8e12", "#d561dd", "#c93f00", 
                "#ddd53e", "#4aef7b", "#e86502", "#9ed84e", "#AB3282", "#CCC9E6", 
                "#8249aa", "#99db27", "#DCC1DD", "#ff523f", "#ce2523", "#f7aa5d", 
                "#cebb10", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", 
                "#d66551", "#1a918f", "#ff66fc", "#2927c4", "#7149af", "#57e559", 
                "#8e3af4")

# Seurat default violin plot
VlnPlot(seu_obj, features = c("Krt14", "Krt5", "Il1b", "Cd34", "Fxyd3", "Sfn"), 
        pt.size = 0, ncol = 3)
```
**Example output** 

<img width="928" height="557" alt="Weixin Image_20250713115004" src="https://github.com/user-attachments/assets/39f65a7c-bf6d-4961-8cce-50b7961675c5" />

```R
# TidyGenePlot violin plot
tidy.VlnPlot(Control, features = c("Krt14", "Krt5", "Il1b", "Cd34", "Fxyd3", "Sfn"), 
             pt.size = 0, cols = cb_palette, ncol = 3)
```
**Example output** 

<img width="928" height="385" alt="Weixin Image_20250713115151" src="https://github.com/user-attachments/assets/9684e8e1-a5a7-42ff-9b5a-6968ca208bcf" />


```R
# Seurat default feature plot
FeaturePlot(Control, features = c("Krt14", "Krt5", "Il1b", "Cd34", "Fxyd3", "Sfn"), 
            ncol = 3)
```
**Example output** 
<img width="751" height="487" alt="Weixin Image_20250713115316" src="https://github.com/user-attachments/assets/0ea040fd-9930-4a8a-ac2a-8473abe3f343" />




```R
# TidyGenePlot feature plot with legend
tidy.FeaturePlot(Control, features = c("Krt14", "Krt5", "Il1b", "Cd34", "Fxyd3", "Sfn"), 
                pt.size = 1, cols = c("grey", "blue"), Legend = TRUE, ncol = 3)
```
**Example output** 

![Uploading Weixin Image_20250713115431.png…]()


```R
# TidyGenePlot feature plot without legend
tidy.FeaturePlot(Control, features = c("Krt14", "Krt5", "Il1b", "Cd34", "Fxyd3", "Sfn"), 
                pt.size = 1, cols = c("grey", "blue"), Legend = FALSE, ncol = 3)
```
**Example output** 

<img width="676" height="482" alt="Weixin Image_20250713115553" src="https://github.com/user-attachments/assets/7ed79626-0c1a-421d-a2e9-3bdb2373fe18" />



