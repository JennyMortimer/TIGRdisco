# TIGRdisco

TIGRdisco finds Tas protein-coding genes in nucleotide sequences and searches their
neighbourhood for TIGR arrays. Two independent engines scan every locus: `dr-tigr_engine` for
the dual-repeat architecture and `sl-tigr_engine` for the stem-loop architecture.

Preprint: Gu, Gong and Mortimer, *bioRxiv* (2026),
[doi:10.64898/2026.11.04.764119](https://doi.org/10.64898/2026.11.04.764119).

## Install

Python 3.10 or later. No external binary is needed.

```bash
pip install git+https://github.com/JennyMortimer/TIGRdisco.git
```

or, from a clone or a release archive, `pip install .`

## Quick start

```bash
tigrdisco-fasta examples/example_windows.fa --out demo --workers 2
tigrdisco-db examples/example_central.gb --out demo_db
```

## Usage

**FASTA mode** takes contigs, genomes or metagenome assemblies. Genes are called with
pyrodigal, screened against the bundled Tas profile HMM, and every passing gene becomes an
anchor.

```bash
tigrdisco-fasta genome.fa --out results --workers 8
```

**GenBank mode** takes records in which the Tas gene is already marked. Each record needs one
`CENTRAL` feature (`start..end` or `complement(start..end)`), a unique `LOCUS` name and an
`ORIGIN` section of the declared length. No gene calling or HMM screen is run.

```bash
tigrdisco-db input.gb --out results --workers 8
```

By default both modes search 600 bp on each side of the anchor for dual-repeat arrays
(`--dr-flank`), and the anchor plus 1000 bp on each side for stem-loop arrays (`--sl-flank`). The output directory must be new or empty.
Run any command with `--help` for its options.

## Output

| file | content |
|---|---|
| `loci.csv` | one row per anchor: `SL`, `DR` and the resulting `architecture` label |
| `loci.jsonl` | one record per anchor with every array and unit |
| `manifest.json` | settings, input SHA-256, counts and runtime |

Add `--xlsx` to also write three Excel workbooks to `results/xlsx/`:
`dr_tas_with_tigr.xlsx`, `sl_tas_with_tigr.xlsx` and `tas_without_tigr.xlsx`.
`--xlsx-style tables` writes smaller workbooks without images. A finished run can be exported
later with `tigrdisco-export --run results --out reports`.

Fields and coordinate conventions are described in `docs/OUTPUT_SCHEMA.md`.
