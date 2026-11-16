# Provenance and rights

## Code

Both detection engines are the author's own code. The dual-repeat engine is the detection core
of the author's earlier package `tas-tigr-engine` 0.3.0, with unused parts removed.

Some parameter values of the stem-loop engine were chosen with reference to the published
description of TIGRFinder (Xu, Long et al., bioRxiv 2026, doi:10.64898/2026.09.09.750476): the
box-width range, the spacer length range and the `|box C| = 3`, `|box D| = 2` reread frame. No
TIGRFinder source code is included in or derived from this package.

## Profile HMMs

`tigrdisco/data/MITTas.hmm3`, `TasH_NTD_seed.hmm` and `TasR_NTD_seed.hmm` were built by the
author with `hmmbuild` from public database sequences listed in table S1 of Faure et al.,
TIGR-Tas: A family of modular RNA-guided DNA-targeting systems in prokaryotes and their viruses,
*Science* 388, eadv9789 (2025), doi:10.1126/science.adv9789. Please cite that paper when using
the profiles. No file from the publisher is included in this package.

## Example data

`examples/example_windows.fa` is supplied only to show that the command line runs.
`examples/example_central.gb` holds eight of the same records in GenBank form for
`tigrdisco-db`.

## Licence

`LICENSE` covers the source code only. It does not assign a licence to the profile HMMs or the
example data.
