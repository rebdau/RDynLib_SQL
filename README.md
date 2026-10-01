# RDynLib

**RDynLib** is an R package for building, curating, and exploiting experiment-based LC–MS spectral databases, referred to as **DynLib databases**. The RDynLib/DynLib framework enables cross-experiment and cross-platform data integration, spectral analysis, and progressive metabolite annotation.

## Why RDynLib/DynLib?

Structural knowledge generated through in-house LC–MS/MS and MSⁿ experiments accumulates over years but is often difficult to preserve, connect, and reuse across studies and analytical platforms. This challenge is particularly important for unknown or partially characterized metabolites, whose annotation may depend on combining experimental evidence acquired at different times, under different experimental conditions, or using complementary analytical platforms.

**DynLib** and **RDynLib** address this challenge by organizing experimental LC–MS data and associated structural information into dynamic, progressively curated spectral databases. Rather than treating individual experiments as isolated datasets, the framework enables experimental evidence and metabolite annotations to be accumulated, connected, and reused as new data become available.

In this framework:

- **DynLib** refers to the experiment-based spectral databases used to organize and preserve LC–MS data, fragmentation spectra, metadata, and metabolite annotations.
- **RDynLib** provides the R tools to build, query, curate, align, compare, and progressively enrich these databases.

This approach facilitates the reuse of previously acquired experimental evidence and supports the progressive characterization of metabolites across experiments and analytical platforms.

> **Status:** RDynLib is under active development. The associated manuscript is currently in preparation.