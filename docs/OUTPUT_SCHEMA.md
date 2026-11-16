# Output files and coordinates

| file | one row or record per | content |
|---|---|---|
| `loci.csv` | anchor | `locus`, `SL`, `DR`, `architecture` |
| `loci.jsonl` | anchor | every array and unit found by both engines |
| `manifest.json` | run | settings, input SHA-256, counts, runtime |

`manifest.json` is written after `loci.csv` and `loci.jsonl`, so its presence marks a complete
detection run. Records are written in input order and do not depend on `--workers`.

In FASTA mode `loci.csv` also carries `contig`, `anchor_start1`, `anchor_end1`,
`anchor_strand`, `tier` and `hmm_score`. The anchor coordinates are 1-based and inclusive on
the input sequence.

### `loci.jsonl`

Every record is oriented so that the anchor gene is on the forward strand.

| key | meaning |
|---|---|
| `locus` | anchor identifier |
| `contig`, `anchor` | FASTA mode only: source sequence and HMM evidence of the anchor |
| `protein_sequence` | translation of the anchor gene |
| `record_length` | length of the input sequence |
| `reverse_complemented` | true when the input sequence was flipped to orient the anchor |
| `CENTRAL_start0`, `CENTRAL_end0` | anchor gene on the oriented sequence, 0-based half-open |
| `SL_window_offset0` | 0-based start of the stem-loop search window on the oriented sequence |
| `SL` | stem-loop arrays |
| `DR` | dual-repeat arrays |

To convert an oriented 0-based half-open interval `[s, e)` to 1-based inclusive input
coordinates: `(s + 1, e)` if `reverse_complemented` is false, otherwise
`(record_length - e + 1, record_length - s)`.

**Dual-repeat arrays.** Each unit has `contig_index`, the 0-based start on the oriented
sequence, and four segments: `R1`, `spacer_a`, `R2`, `spacer_b`. The unit ends at
`contig_index` plus the total length of the four segments. `classification` is `Core`,
`Prefix` or `Suffix`. `bits` is the information content of the repeat columns.

**Stem-loop arrays.** Unit `start` and `end` are 1-based and inclusive inside the search
window; add `SL_window_offset0` to reach the oriented sequence. `extended_start` and
`extended_end` include the stems, whose length is `stem_length`. The seven segments, in order,
are `box_C`, `spacer_A`, `box_D`, `connector`, `box_C_prime`, `spacer_B`, `box_D_prime`.
`exact_units` counts units identical to the array consensus and `support_fraction` is
`exact_units / unit_count`. `array_class` is `core` for arrays of at least three units and
`fragment` otherwise.

Stem-loop arrays removed by arbitration are absent from `SL`. The arbitration setting of the
run is recorded in `manifest.json`.

### Excel workbooks (`--xlsx`)

`dr_tas_with_tigr.xlsx`, `sl_tas_with_tigr.xlsx` and `tas_without_tigr.xlsx`. An anchor
positive for both engines appears in both positive workbooks. The first worksheet is a
per-anchor summary, followed by `Anchors`, `Arrays`, `Units` and `Run_info` tables;
`tas_without_tigr.xlsx` has only `Anchors` and `Run_info`. In `Units`, `input_start1` and
`input_end1` are 1-based inclusive coordinates on the input sequence.
