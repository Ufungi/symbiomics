# Clade notes

Some annotation tools look applicable to a genome but have no trained model or
declared support for its clade. Running them anyway produces output that looks
plausible and is not.

This file exists so those exclusions are recorded once rather than
rediscovered — and so that a future contributor who proposes adding one of these
tools finds the reason before writing the module. The policy engine enforces
the table below (`CLADE_UNSUPPORTED` in `bin/decide_strategy.py`), and
`nextflow_schema.json` deliberately has **no parameter** for the excluded tools.

## Excluded

| Tool | Excluded for | Reason |
|---|---|---|
| **Tiberius** | gymnosperms | Ships models for Mammalia, Vertebrata, Insecta, **Angiosperms (Mesangiospermae)**, Fungi, Chlorophyta and Bacillariophyta. Gymnosperms are in none of them. |
| **NCBI EGAPx** | gymnosperms, fungi | Supported taxa are Arthropoda, Vertebrata, **Magnoliopsida** (flowering plants), Cnidaria and Echinodermata. The documentation states fungi, protists and nematodes are out of scope; gymnosperms are not flowering plants. |

*Pinus densiflora*, the reference plant target for this pipeline, is a conifer —
so both are unavailable for it, and Helixer's `land_plant` model is the only
deep-learning track that generalises.

## Clade-conditional behaviour

| Tool | Behaviour |
|---|---|
| **Helixer** | lineage is one of `fungi`, `land_plant`, `vertebrate`, `invertebrate`. `land_plant` covers gymnosperms; there is no finer model. |
| **BRAKER** | `--fungus` switches on the GeneMark branch-point model. Set for fungi only; never for plants. |
| **funannotate2** | Documented as fungal-first: default `--max-intron` is 3000 bp and the bundled databases (MEROPS, dbCAN, MIBiG) are fungal/microbial. The policy engine blocks it above 2 Gb regardless of clade. |
| **antiSMASH** | Fungal branch only. Plant BGC rules are weak; plantiSMASH would be a separate integration. |
| **EDTA** | Designed and benchmarked on plants (LTR-dominated landscapes). Runs on fungi but EarlGrey gives cleaner libraries there. |

## Clade detection

`--taxon` is the coarse switch (`fungi` / `plant` / `animal` / `other`).
`--clade` is the finer hint that gates the table above. When `--species` begins
with a recognised genus the clade is inferred — `Pinus` → `gymnosperm` — but the
inference is deliberately shallow. **Set `--clade` explicitly** when it matters;
the resolved value is recorded in `strategy.yml` so it is auditable.

## Adding a tool

Before adding a predictor or annotation module, check what clades its released
models actually cover, not what the abstract implies. If coverage is partial,
add it to `CLADE_UNSUPPORTED` with the reason and a citation, and the policy
engine will exclude it with an explanation instead of producing quiet nonsense.
