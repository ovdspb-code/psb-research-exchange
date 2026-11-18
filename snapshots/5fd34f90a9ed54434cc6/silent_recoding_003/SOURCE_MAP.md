# 003: primary evidence and source map

Date retrieved: 2026-10-09. Evidence extraction is new; accepted 001/002 R1 science is preserved.

## Sequence and quantitative evidence

- GC-replay source snapshot: `matsengrp/gcreplay@2d1d6592c718f3b2997e214b21f92cc2aeb62c2e`. Selected 002 tables copied byte-for-byte from the Advisor 003 package, independently matched to retained R1 tables. The focal key is `D20_M23_GC102__seq7__seq6`; protein `P_abadea0d98857de8`. Role assignment follows the frozen DNA columns; the `name` field of planned starts is not used to infer genealogy.
- DMS source: [jbloomlab/Ab-CGGnaive_DMS](https://github.com/jbloomlab/Ab-CGGnaive_DMS/tree/5ffc88950f760a297e2196ac449111e3db023e80), pinned `5ffc88950f760a297e2196ac449111e3db023e80`. Seventeen selected files fetched; Git blob SHA1 verified before resolving LFS pointers; actual LFS size and SHA256 verified. Full manifest in results/SOURCE_MANIFEST.json. No moving-main data used.
- WT assay CDS is 705 nt / 235 AA: HC1–112; linker113–127 `GGGGSGGGGSGGGGS`; shared LC128–234 (107 AA); additional terminal R235. Tags/display/constant/native-cell components are not determined by this CDS. Focal variable CDS translates to naive HC and LC with only Y36N, equal to scFv Y163N / IMGT Y42(L)N. A matched variable protein does not establish identical complete assay carrier.
- Current per-library tables: `results/binding_Kd/variant_binding.csv` and `results/expression_meanF/variant_expression.csv`, joined on library+target+AA substitution. WT reference is within each library. Unweighted mean of two library contrasts reproduces the released +.25938/−1.26992 after rounding.
- Current `Snakefile` sends binding to `Titeseq-modeling.ipynb` and expression to `compute_expression_meanF.Rmd`; fitting is after within-library aggregation by AA variant. Individual binding CSV uses columns `TiteSeq_avgcount/log10Ka`, corresponding to the earlier per-barcode binding workflow. Its values are not the current pooled fit: Y163N legacy means ~10.60/10.49 versus current pooled 8.927/8.838. Both are retained with explicit identities. A biological discrepancy or an author error is not inferred from different estimator outputs.
- Individual expression fits and barcode mapping are available. Finite Y163N fits: binding31/44, expression40/44. Pooled support20+24 is not a claim that all44 individual fits succeeded. WT has23285 mapped per-barcode rows across the two libraries. Missing fits retained, not imputed.
- `codon_variant_table.csv` reconstructs44 Y163N barcode DNA variants from the WT CDS. None matches any complete focal HC+LC coding combination. Some carry extra synonymous substitutions. Aggregation by AA can conceal encoding-specific effects; the existing data cannot certify exact a/x/R equivalence.

## Exact limits of raw-data availability

The pinned README explicitly says initial sequencing is not included, while downstream output tables are available. Fitted barcode scores, library scores and variant-to-barcode maps were actually fetched. Raw FCS/intensities and initial PacBio/Illumina read processing were not rerun. `prepped_barcode_counts.csv` is a released420089335-byte LFS downstream count table; it was NOT fetched in this focused extraction. Its public existence is distinguished from absent raw inputs. Final scores are not converted into fabricated titration curves or confidence intervals. Questions about laboratory replicate provenance, raw curves and uncertainty go to the authors.

## Read primary methodological evidence

1. [DeWitt et al., Cell2026](https://doi.org/10.1016/j.cell.2026.05.013): accessible primary Results/Methods and data repositories; same yeast-scFv assay, barcoded library analysis, and mutation/selection distinction. No new experimental result inferred.
2. [Cui et al., J Immunol2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC5161250/): primary abstract, Methods and Results available through indexed primary page. Nonfunctional mouse κ in NP-binding λ cells; clone grouping and unique events; species/chain-specific targeting. It does not provide exact first-two-event trajectories. PMID27707999, DOI10.4049/jimmunol.1502263.
3. [Baughn et al., mBio2011](https://journals.asm.org/doi/10.1128/mbio.00186-11): primary full-text Results/Discussion. RMCE replaces endogenous human Ramos IgH; endpoint sequencing after accumulation and founder-mutation handling. Heavy-chain/human context, G/C bias and low mutation/error separation constrain transfer to this mouse LC. No proven ordered readout for our three spellings.
4. [McGrath et al., Mol Cell2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12063209/): primary abstract/Summary and Discussion; synonymous deoptimization can affect antibody yield. This is motivation to measure exact encodings, not evidence that our three encodings differ.

Browser access to some PMC full pages returned challenge pages; sections above were read from indexed primary-page results. No claim to have inspected unavailable supplements or every raw dataset. Betz1993 was in the Advisor brief; it is not needed as an additional inspected platform here.
