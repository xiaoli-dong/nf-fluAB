# nf-core/influenza: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.6 - [2026-07-14]
Bug fix release

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`
Replaced csvtk concat with concat_tables.py to preserve all Nextclade output columns. This fixes an issue where H5 Nextclade reports lacked several clade-related columns (proposedSubclade, short-clade, subclade, legacy-clade, RBD), causing csvtk concat to drop these columns when merging H5 and H1/H3 reports and resulting in downstream R script failures.
## v1.0.5 - [2026-06-16]
Bug fix release

### `Added`

### `Fixed`
Fixed a Nanopore mapping statistics bug. Corrected BAM filtering to exclude secondary alignments, supplementary alignments, and PCR duplicates before generating mapping counts and depth profiles. The issue was limited to mapping statistics and did not affect variant calling, consensus sequence generation, clade assignment, or any other analysis results.

### `Dependencies`

### `Deprecated`


## v1.0.3 - [2025-05-27]

Initial release of nf-core/influenza, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`
Updated the FluAB reference database filtering criteria by relaxing the sequence length thresholds for segments 7 and 8 from 90–110% to 75–125% of the expected segment length, allowing retention of Flu B segment 7 and 8 reference sequences; Fixed an issue where some software versions were missing from the reports and corrected a bug that prevented generation of the consensus summary report for mixed infections; Added the pipeline version as a new column in the consensus summary report to improve traceability and reproducibility.

### `Dependencies`

### `Deprecated`


## v1.0dev - [date]

Initial release of nf-core/influenza, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`
