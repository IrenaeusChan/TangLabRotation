# TangLabRotation
Collection of code snippets used during the rotation in Dr. Rui Tang's Lab

# Overview
1. MOBASeq Pipeline
2. Visualizations
3. Code Snippets

# Moba-Seq
The original MobaSeq was written by Andy Xu and Dr. Rui Tang and can be found [here](https://github.com/tanglab-2024/moba-seq/).

Leveraging off their backbone, [this repository](https://github.com/IrenaeusChan/mobaseq_dev) was developed with the sole purpose of improving overall processing efficiency, optimizing internal algorithms, and create a more streamlined robust CLI that allows various users the ability to process data acquired from MOBASeq.

# Analysis
Throughout this rotation three datasets were analyzed: Moba500, MobaV, and MobaV Immune

# Moba500
- Validation --> validation.R
- Scatter Plots --> validation.R
- Power Calculation --> power_calculation.R
- Dormancy Exploration --> dormancy.R
- Outlier Removal --> outlier.R

# MobaV
- Crebbp Specific Analysis --> crebbp.R
- Heatmaps --> crebbp.R
- SpikeIn Regression Line --> validation.R
- Genotype Specific Dormancy Prevention --> crebbp.R

# MobaV Immune
- Comparing BL6, NSG, and RAG2KO Mice --> immune.R
