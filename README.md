# LPiModulons
Independent component analysis (ICA) of Lactobacillus plantarum transcriptomes

## What is iModulon ?
Each iModulon (IM) is a group of genes that represents an independently modulated signal. To learn more about iModulons, how they are computed, and what they can tell you, please visit https://imodulondb.org/about.html.

## Workflow
1. Preprocess: log10 transformation and center to reference condition(wt_pH6.2)
2. ICA decomposition and IM thershold computation via D'Agostino's K-squared test
3. Regulatory and functional annotations of iModulons (IMs)
4. Case studies of acid-active and carbon source-active IMs
5. Analysis of IM activities

## Dependencies
1. PyModulon: https://github.com/SBRG/pymodulon
2. MEME suite: https://meme-suite.org/meme/
3. eggNOG-mapper: http://eggnog-mapper.embl.de
4. Biopython: https://biopython.org/
5. DNA Features Viewer: https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer
6. Networkx: https://networkx.org/
7. Seaborn: https://seaborn.pydata.org/

## Citation
Qiu, S., Huang, Y., Liang, S., Zeng, H., Yang, A. (2023, January 1). Systematic elucidation of independently modulated genes in Lactobacillus plantarum reveals a trade-off between secondary and primary metabolism. bioRxiv. https://doi.org/10.1101/2023.11.03.565434 
