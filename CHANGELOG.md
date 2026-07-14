# nf-core/influenza: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.6 - [07/14/2026]
Bug fix release

### `Added`

### `Fixed`
csvtk concat only output common header columns. H5 database is downloaded from coummunity website, it is missing "proposedSubclade        short-clade     subclade        legacy-clade    RBD", when concating H5 nextclade report with H1 nextclade report, the output will only output the common header columns. I added concat_tables.py script to handle the problem.
### `Dependencies`

### `Deprecated`

## v1.0.5 - [06/16/2026]
Bug fix release

### `Added`

### `Fixed`
When generating the depth, count profile for nanopore data, the BAM file was inadvertently skipping the filtering step that removes secondary alignments, supplementary alignments, and PCR duplicates.
that supplementary alignments were inflating the mapping counts. This is fixed in this version

### `Dependencies`

### `Deprecated`


## v1.0dev - [date]

Initial release of nf-core/influenza, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`
