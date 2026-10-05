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

## Workflows

RDynLib provides a collection of workflows illustrating how DynLib databases can be constructed, accessed, integrated, analyzed, and progressively enriched.

The complete workflows are available in the [`vignettes/`](https://github.com/rebdau/RDynLib_SQL/tree/Ahlam/vignettes) directory.

### Getting started

#### `rdynlib_overview_workflow.qmd`

This workflow provides a global overview of the different applications of **RDynLib** for working with **DynLib data**. It illustrates how the main functionalities of RDynLib can be combined to explore DynLib databases, connect information across experiments and analytical platforms, compare fragmentation spectra, and exploit accumulated structural information for metabolite characterization and annotation.

This workflow is the recommended starting point for obtaining an overview of the RDynLib/DynLib framework.

#### `dynlib_database_download_query_workflow.qmd`

This workflow demonstrates how existing DynLib databases can be downloaded and accessed in R. It provides examples of querying the underlying SQL database and retrieving experimental, compound, and spectral information stored in DynLib.

#### `dynlib_database_access_workflow.qmd`

This workflow demonstrates how DynLib databases can be accessed and manipulated using RDynLib and Bioconductor data structures. It includes examples of database subsetting, conversion to `CompDb` and `Spectra` objects, and retrieval of compounds and their associated fragmentation spectra.

### Integrating new experimental data into DynLib

#### `dynlib_ftms_database_integration_workflow.qmd`

This workflow demonstrates how a newly generated SQL database containing processed **FT-MS/MS and MSⁿ data** can be integrated into an existing DynLib database. It enables new FT-MS experimental data, including compounds, fragmentation spectra, spectral peaks, and associated experimental information, to be incorporated into DynLib while preserving the existing database content.

#### `dynlib_qtof_database_integration_workflow.qmd`

This workflow demonstrates how a newly generated SQL database containing processed **QTOF-MS/MS data** can be integrated into an existing DynLib database. It enables new QTOF experimental data and their associated spectral information to be incorporated into DynLib and subsequently connected with information accumulated from previous experiments and analytical platforms.

### Integrating and characterizing DynLib data

#### `dynlib_alignment_workflow.qmd`

This workflow demonstrates the alignment of compounds across different LC–MS experiments and analytical platforms. Candidate associations are established using precursor *m/z* and retention-time information, followed by retention-time modeling and the exploitation of spectral information to support the association of corresponding compounds.

This workflow enables structural information acquired in different experiments to be connected and reused within the DynLib framework.

#### `dynlib_cid_spectral_matching_workflow.qmd`

This workflow demonstrates the comparison of fragmentation spectra to identify compounds sharing common CID fragmentation characteristics. It illustrates how spectral information stored in DynLib can be exploited to identify related spectra and support metabolite characterization across experiments.

#### `dynlib_characterization_visualization_workflow.qmd`

This workflow illustrates tools for exploring and visualizing structural information stored in DynLib databases. It includes examples based on CID spectral relationships, CSPP and GNPS information, MS/MS and MSⁿ fragmentation data, metabolite annotations, and the export of spectral information for downstream analyses.

### Database curation and utilities

#### `dynlib_annotation_workflow.qmd`

This workflow demonstrates the management of metabolite annotations within DynLib databases. It illustrates how compound names and structural information can be added or updated as new experimental evidence becomes available.

DynLib can also preserve successive annotation events, allowing the evolution of metabolite characterization to be tracked over time.

#### `dynlib_database_filtering_workflow.qmd`

This workflow demonstrates how DynLib databases can be filtered to generate subsets containing selected experiments, compounds, spectra, or other relevant information for downstream analyses.

#### `dynlib_spectrum_type_update_workflow.qmd`

This workflow demonstrates how spectrum-type information stored in DynLib can be managed and updated, allowing different types of experimental or processed spectra to be distinguished within the database.

#### `dynlib_xcms_grouping_workflow.qmd`

This workflow demonstrates the integration of XCMS-derived feature grouping information with the DynLib SQL framework, supporting the connection between LC–MS features and the corresponding information stored in DynLib.

## Installation

The development version of **RDynLib** can be installed from the `Ahlam` branch on GitHub:

```r
remotes::install_github(
    "rebdau/RDynLib_SQL",
    ref = "Ahlam",
    build_vignettes = TRUE
)