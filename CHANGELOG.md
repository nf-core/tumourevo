# nf-core/tumourevo: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.0dev - to be released

### Added

- Support for Hartwig WiGiTS / oncoanalyser somatic outputs: SAGE/PURPLE somatic VCFs in `vcf2cnaqc` and PURPLE copy number (`purple.cnv.somatic.tsv`, `purple.purity.tsv`, `cna_caller` `PURPLE`) in `cna2cnaqc`.

- Test profile: set `sparsesignatures_k` (it was spelled `sparsesignatures_K`, which the pipeline does not read), so `SPARSE_SIGNATURES` runs with the intended `K = 2:3` instead of the default `3:10` ([#92](https://github.com/nf-core/tumourevo/issues/92))

### `Dependencies`

### `Deprecated`

## v1.0.0 - 2026-09-18

Initial release of nf-core/tumourevo, created with the [nf-core](https://nf-co.re/) template.
