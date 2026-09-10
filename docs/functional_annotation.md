# Functional annotation (`-entry functional`)

Protein FASTA in, functional descriptions out. No genome, no samplesheet, no
structural annotation. Built for exactly the case where a good gene model
already exists — the reference target for this pipeline, *Pinus densiflora*,
ships `HA.PEP.fa` / `HB.PEP.fa` with 44,233 predicted proteins per haplotype
and zero functional annotation.

```bash
scripts/symbiomics run . -entry functional -profile singularity,local64 \
    --proteome HA.PEP.fa --genome_id Pinde_HA --taxon plant \
    --project Pinde_HA
```

## Branching

| `--taxon` | route |
|---|---|
| `fungi` | funannotate2 + `funannotate2-addons` (eggNOG-mapper, InterProScan, antiSMASH, SignalP6 are already wrapped by funannotate2 — **not called separately**). Not yet implemented — errors clearly rather than silently running the wrong thing; see the roadmap. |
| `plant` / `animal` / `other` | the shared core below |

No AHRD anywhere in this pipeline. Product naming is a direct best-hit ladder
(`bin/assign_products.py`), not a weighted-blast synthesis step.

## Phase A — runs today, zero new downloads

Every tool here reads a database this lab already has on disk.

| module | tool | database | what it adds |
|---|---|---|---|
| `swissprot` | DIAMOND blastp | `/data/genome/db/Swiss-prot/uniprot_sprot.fasta` | best-hit description, parsed straight from the UniProt FASTA header |
| `eggnog` | eggNOG-mapper 2.1.15, diamond mode | `/data/genome/db/eggnog` (v5) | GO terms, EC numbers, KEGG KO, COG category, orthology-based description |
| `dbcan` | run_dbcan (v3 line), all three callers: HMMER + DIAMOND + eCAMI | `/data/genome/db/CAZy` (2021 vintage) | CAZyme family assignment, with a `#ofTools` consensus column across all three callers in `overview.txt` |

`--functional_modules swissprot,eggnog,dbcan` (the default). Drop any of them,
or add Phase B modules once provisioned (below).

**Why DIAMOND-vs-Swiss-Prot instead of full UPIMAPI for Phase A**: UPIMAPI's
richest columns (GO/EC/Pathway/KEGG/Pfam cross-references) come from its
ID-mapping step, which needs a downloaded ID-mapping resource this lab does
not have yet. Every UniProt FASTA header already carries the product
description in a fixed format
(`>sp|ACCESSION|NAME DESCRIPTION OS=... OX=... GN=... PE=... SV=...`), so a
DIAMOND best hit plus a header parse (`bin/parse_diamond_swissprot.py`)
delivers the single most load-bearing field today, with zero extra downloads.
Full UPIMAPI (with `--local-id-mapping` resources provisioned) is a Phase B
upgrade to the same slot in the ladder — the column set grows, the ladder
position doesn't change.

**Why an old dbCAN line for Phase A**: dbCAN's CLI and database layout changed
materially between the v3 line and v4/v5 — different positional arguments,
different HMM/DIAMOND file names. This lab's `/data/genome/db/CAZy` is a
2021-vintage v3-compatible database (`CAZyDB.09242021.fa`, `dbCAN.txt` +
pressed HMM files, `stp.hmm`) — the same one `genome_pipeline.sh` already used
successfully. Pointing a v5 container at it would silently look for files that
aren't there. `RUN_DBCAN` is pinned to `dbcan:3.0.7` specifically to match.
dbCAN v5 (with the free bundled TCDB/peptidase databases) is a Phase B item —
see the roadmap — not a drop-in replacement for this container without also
fetching the matching v5 database.

**All three callers, not HMMER alone**: dbCAN3's `overview.txt` is designed
around a `#ofTools` consensus column across HMMER, DIAMOND, and eCAMI — an
earlier pass that ran `--tools hmmer` only was throwing away exactly that
consensus mechanism. Verified directly against the pinned container
(`dbcan:3.0.7--pyh5e36f6f_0`, `run_dbcan --help`): `--tools hmmer diamond
eCAMI` is the correct invocation, and eCAMI's kmer database
(`--eCAMI_kmer_db`, default `"CAZyme"`) ships bundled inside the eCAMI
package itself — no separate download or `/data/genome/db/CAZy` layout
change was needed. A real run against this lab's existing CAZy DB directory
completed in ~35s with all three columns (`HMMER`, `eCAMI`, `DIAMOND`,
`#ofTools`) populated.

## Phase B — needs new database provisioning (~80 GB)

Not yet wired into `functional_core`. Planned:

| module | tool | database | size |
|---|---|---|---|
| `interproscan` | InterProScan 6 (native Nextflow, `--datadir`) | incremental | ~60 GB |
| `kofam` | KofamScan | `ftp.genome.jp/pub/db/kofam` | ~6.5 GB |
| `dbcan5` | run_dbcan v5 | `db_v5-2_9-13-2025` | ~10 GB |
| fungal branch | funannotate2 | `funannotate2 install -d all` | ~15 GB |

## The product-name ladder

`bin/assign_products.py`, first non-empty source wins:

```
1. funannotate2 product        (fungi branch, Phase B)
2. Swiss-Prot best hit         (Phase A -- pident/evalue filtered, "uncharacterized
                                 protein"-style placeholders in the Swiss-Prot
                                 entry itself are treated as no hit)
3. eggNOG-mapper Description / Preferred_name
4. InterProScan signature description  (Phase B)
5. "hypothetical protein"
```

`bin/qc_products.py` reports the fraction that fell through to
`hypothetical protein` (`--max_hypothetical_frac`, default 0.55) and warns —
not fails — above the threshold. A high fraction usually means: wrong
`--taxon`, the reference proteome is too distant, or the Phase B modules
(InterProScan/KofamScan) would materially help.

## Consensus tags and the HTML report

`products.tsv`'s `source` column is the *ladder winner* — the one source that
supplied the product name. Independently of that, every row also carries:

- `sources_with_evidence` / `n_sources` — every source that found *any* hit
  for that protein, not just the one that won the ladder (e.g. a protein can
  show `swissprot,eggnog,dbcan` even though `source` says `swissprot`).
- `dbcan_family` / `dbcan_tools` — CAZy family call and dbCAN's own
  `#ofTools` agreement count, carried through from `overview.txt`.

This is the free, in-pipeline alternative to Blast2GO's per-gene tag column
(Blast2GO itself is paid) — same idea (which independent programs agree on
this gene), no external tool. `bin/render_annotation_report.py` turns it into
`annotation_report.html`: a single self-contained file (embedded data, vanilla
JS, no server or CDN — this lab's nodes are offline) with a sortable,
filterable table and colored tag chips per source.

## Output

```
results_functional_HA/
├── 80_functional/
│   ├── swissprot/    swissprot.diamond.tsv, swissprot_hits.tsv
│   ├── eggnog/       symbiomics.emapper.annotations
│   ├── dbcan/         dbcan_out/overview.txt
│   ├── products.tsv          <- protein_id, product, source, GO, EC, KO, Swiss-Prot hit,
│   │                             dbcan_family, dbcan_tools, sources_with_evidence, n_sources
│   ├── products_qc.json      <- %hypothetical and the per-source breakdown
│   └── annotation_report.html <- sortable/filterable consensus-tag report (open in a browser)
└── 99_report/multiqc_report.html
```

## Running it against both haplotypes

Each haplotype is a separate protein FASTA, so run it twice:

```bash
for hap in HA HB; do
    scripts/symbiomics run . -entry functional -profile singularity,local64 \
        --proteome ${hap}.PEP.fa --genome_id Pinde_${hap} --taxon plant \
        --project Pinde_${hap}
done
```

There is no haplotype-specific logic in this entry point — it is exactly the
`FUNCTIONAL_CORE` subworkflow both times. What differs between HA and HB
happens downstream, in the allele-pairing table
(`docs/allele_specific_expression.md`), which is a separate concern from
functional annotation.
