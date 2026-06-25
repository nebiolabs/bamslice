# fastp JSON merge fidelity

`fastp-merge` rebuilds a single fastp report from the per-chunk reports produced
when a BAM is sliced and each slice is run through fastp. Ideally the result
matches the report fastp would have produced on the unsplit dataset. This note
records where it does and where it can't.

To check, run fastp on the whole dataset and on each slice, merge the slices, and
diff:

```sh
fastp-merge slice_1.fastp.json slice_2.fastp.json -o merged.json
diff <(python3 -m json.tool full.fastp.json) <(python3 -m json.tool merged.json)
```

## Reconstructed exactly

Totals (`total_reads`, `total_bases`, `q20_bases`, `q30_bases`,
`filtering_result`, `kmer_count`) are summed. Rates (`q20_rate`, `q30_rate`,
`gc_content`) and mean lengths are recomputed from the merged totals. Floats are
written at 6 significant figures to match fastp's C++ output, and top-level keys
follow fastp's order.

## Quality and content curves

Each chunk reports a per-cycle mean over its own reads, so the whole-dataset mean
is the read-count-weighted average of the per-chunk means (weighted by the
chunk's `total_reads`), not a flat average. Two small residuals remain:

- fastp rounds each chunk's curve values to 6 significant figures before we see
  them, so the weighted average can land ±1 in the last figure. The raw per-base
  sums needed to avoid this aren't in the JSON.
- For `after_filtering`, reads vary in length, so the exact per-cycle weight is
  "reads still long enough to reach this cycle", which fastp doesn't expose. We
  use `total_reads`, which is exact for full-length `before_filtering` reads and a
  close approximation afterward.

These account for almost all of the remaining diff.

## Can't be reproduced from the chunk JSONs

- **`duplication.rate`** — duplicates split across chunks are invisible to any
  per-chunk count, so the weighted average undercounts them. To get the true rate
  you'd have to recompute duplication on the merged FASTQ.
- **`insert_size`** — read pairs split across a chunk boundary lose their mate and
  fastp counts them as `unknown`. The merge (summing histograms and `unknown`) is
  correct; the inputs are degraded. Fix belongs in the slicing step: keep both
  mates in the same chunk.
- **`adapter_cutting.*_adapter_counts`** — fastp prunes each chunk's adapter tree
  to a top-N set and rolls the rest into `others`. Each chunk prunes differently,
  so the counts aren't additive and the discarded detail can't be recovered.
  Summing and re-pruning would improve key overlap but won't be exact.
- **`command`** — there's no single command line for a merged result; it's left
  out.
