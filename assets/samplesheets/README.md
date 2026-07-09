# Example samplesheets

One samplesheet format serves all four start points. The columns you fill in decide which one you get; unused columns may be omitted or left empty.

| File                        | Run with                          |
| --------------------------- | --------------------------------- |
| `samplesheet_bulk_fastq.csv` | `--mode bulk`                     |
| `samplesheet_bulk_bam.csv`   | `--mode bulk`                     |
| `samplesheet_sc_fastq.csv`   | `--mode sc`                       |
| `samplesheet_sc_bam.csv`     | `--mode sc`                       |

RBP-STAMP uses the bulk samplesheets with `--run_flare` added; it needs no samplesheet of its own.

Paths may be absolute or relative to the samplesheet's own directory. See [`docs/usage.md`](../../docs/usage.md#samplesheet-input) for the full column reference.
