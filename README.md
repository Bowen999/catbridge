# CAT Bridge: A Versatile Tool for Multi-Omics Analysis

## Overview
CAT Bridge (Compounds And Transcripts Bridge) is a robust tool built with the goal of **uncovering biosynthetic mechanisms in multi-omics data**, such as identifying genes potentially involved in compound synthesis by incorporating metabolomics and transcriptomics data.

One of the key features of CAT Bridge is its custom-built design for **handling multi-mics time series data** - a niche that has seen a dearth of dedicated tools till now. 

And CAT Bridge is **integrated with GPT 3.5 Turbo** integration to help users dive deeper into the complex biological mechanisms at play.

**Visualisation** is key to understanding complex data. That's why CAT Bridge also includes a comprehensive suite of visualization capabilities to help you understand, explore, and present your findings more intuitively and effectively.

At the heart of CAT Bridge lies a data-driven algorithm, theoretically enabling its application across a wide spectrum of multi-omics data showcasing causal relationships. However, so far, we've tested this only between metabolome-transcriptome and metabolome-metabolome interactions.

![Pipeline](https://github.com/Bowen999/catbridge/blob/main/figures/pipeline.png)

## Installation
### dependencies
Dependencies can be installed using pip (from the Unix terminal)   

```
# Python dependencies
pip install matplotlib catbridge psutil numpy pandas statsmodels sklearn scipy tslearn seaborn bioinfokit networkx pillow adjustText textwrap datashader openai


# R dependencies  
R -e 'install.packages("BiocManager"); BiocManager::install("DESeq2"); install.packages("readr")'
```

### CAT Bridge
```pip install catbridge```


### Know More
don't konw what is **Granger Causality Test**? look at this video: [STATISTICS I Time Series I Granger Causality Test I Intuition and Example](https://www.youtube.com/watch?v=6dOnNNxRJuY)

## Contact
The package is developed and maintained by Bowen Yang (by8@ualberta.ca). Please, reach us for problems, comments or suggestions.


