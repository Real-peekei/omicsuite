# omicsuite roadmap

This is the working backlog for expanding omicsuite beyond its current
modules, organized by the domains in the R & Python Packages/Libraries by
Pipeline reference. Every entry follows [CONTRIBUTING.md](CONTRIBUTING.md)'s
pipeline contract (fit + diagnose + plot + verdict) when built -- this file
tracks status, not implementation.

**Status key:** `built` | `planned` | `candidate` (worth doing, not yet
scheduled) | `out of scope` (see note).

A pipeline lands `out of scope` when the underlying tool is a CLI binary
(BWA, GATK, MACS2, AlphaFold, etc.) rather than an R library -- omicsuite
wraps statistical/analytical R packages, not command-line bioinformatics
tools. Where the reference PDF itself notes "CLI-based" or "not R-native,"
that's carried over here rather than re-litigated.

## 1. Biostatistics

| Pipeline | R package(s) | Status | Notes |
|---|---|---|---|
| Survival analysis: Cox PH | `survival` | **built** | `fit_coxph_pipeline()`, v0.1.0 |
| Survival analysis: Kaplan-Meier | `survival`, `survminer`, `flexsurv` | **built** | `fit_km_pipeline()`, v0.5.0 |
| Survival analysis: competing risks | `cmprsk` | planned | Next up -- cumulative incidence + Fine-Gray, same "Survival analysis" row as Cox/KM |
| RNA-seq DE (Bayesian mixed model) | `brms` | **built** | `fit_rnaseq_nb_pipeline()`, v0.2.0 |
| Epidemic simulation | base R | **built** | `simulate_gillespie_epidemic()`, v0.3.0 |
| Multi-omics integration | `RGCCA` | **built** | `integrate_multiomics()`, v0.4.0 |
| Mixed models / longitudinal (LMM, GEE) | `lme4`, `geepack` | candidate | Natural adjacent module to survival + RNA-seq |
| Meta-analysis | `metafor` | candidate | Well-scoped, CRAN-only, no heavy deps |
| Propensity score / causal inference | `MatchIt`, `WeightIt` | candidate | Pairs naturally with the survival modules for adjusted causal estimates |
| Bayesian hierarchical modeling (general) | `brms` | overlaps with RNA-seq module | Could generalize `fit_rnaseq_nb_pipeline()`'s Bayesian machinery into a general-purpose `fit_bayesian_mixed_pipeline()` rather than duplicating |
| Penalized regression (LASSO/elastic net) | `glmnet` | candidate | Feature-selection companion to survival/RNA-seq modules |
| ROC/AUC, biomarker diagnostics | `pROC` | candidate | Small, focused, good "quick win" module |
| Group sequential / adaptive trial design | `gsDesign`, `rpact` | candidate | Different flavor (design-time, not analysis-time) -- would need contract adaptation |
| Joint longitudinal-survival modeling | `JMbayes2` | candidate | Natural extension once mixed-models module exists |
| Instrumental variables / Mendelian randomization | `TwoSampleMR` | candidate | Fits the causal-inference cluster above |
| Multiple testing correction / FDR | `qvalue` | out of scope as standalone | Better as a shared utility other modules call, not its own pipeline |

## 2. Genomics / DNA-sequencing

| Pipeline | R package(s) | Status | Notes |
|---|---|---|---|
| Variant annotation & filtering | `VariantAnnotation` (Bioconductor) | candidate | Bioconductor install adds friction (see Dependency Policy in CONTRIBUTING.md) |
| GWAS QC & association | `SNPRelate`, `GWASTools` | candidate | |
| Polygenic risk scores | `bigsnpr`, `lassosum` | candidate | |
| Population genetics (Fst, structure, PCA) | `adegenet`, `pegas` | candidate | |
| Phylogenetics | `ape`, `phangorn`, `ggtree` | candidate | Matches the user's cosmology/evolutionary-biology interests |
| CNV / structural variants | `DNAcopy`, `ExomeDepth` | candidate | |
| ChIP-seq / ATAC-seq peak analysis | `ChIPseeker`, `DiffBind` | candidate | Bioconductor |
| DNA methylation | `minfi`, `methylKit` | candidate | Bioconductor |
| Alignment, variant calling, assembly QC | BWA, GATK, SPAdes | **out of scope** | CLI tools, not R libraries |

## 3. Transcriptomics / RNA

| Pipeline | R package(s) | Status | Notes |
|---|---|---|---|
| Bulk RNA-seq DE (frequentist) | `DESeq2`, `edgeR`, `limma` | candidate | Complements the existing Bayesian NB module rather than replacing it -- different inferential philosophy, both legitimate |
| Gene set enrichment / pathway analysis | `clusterProfiler`, `fgsea` | candidate | Natural next step after any DE module |
| Single-cell RNA-seq core workflow | `Seurat`, `scran` | candidate, large scope | Seurat alone could be several modules (QC, clustering, DE) |
| Batch correction | `sva` (ComBat), `RUVSeq` | candidate | |
| Trajectory / pseudotime | `monocle3`, `slingshot` | candidate | Depends on single-cell module existing first |
| Cell-cell communication | `CellChat`, `NicheNet` | candidate | |

## 4. Proteomics

| Pipeline | R package(s) | Status | Notes |
|---|---|---|---|
| MS data processing / stats | `MSstats`, `DEP` | candidate | |
| Protein network / PPI | `STRINGdb` | candidate | |

## 5. Metagenomics / microbiome

| Pipeline | R package(s) | Status | Notes |
|---|---|---|---|
| 16S amplicon analysis | `dada2`, `phyloseq`, `vegan` | candidate | |
| Differential abundance | `ANCOMBC` | candidate | |

## 6. Structural biology / molecular modeling

| Pipeline | R package(s) | Status | Notes |
|---|---|---|---|
| Molecular dynamics analysis | `bio3d` | candidate | |
| Protein structure prediction, docking | AlphaFold, rdkit | **out of scope** | Not R-native |

## 7. Systems biology / networks

| Pipeline | R package(s) | Status | Notes |
|---|---|---|---|
| Co-expression networks | `WGCNA` | candidate | Strong fit given the user's systems-biology thesis direction |
| Gene regulatory network inference | `GENIE3`, `minet` | candidate | |
| Network visualization/analysis | `igraph` | candidate | General-purpose enough to underpin several other modules |
| Metabolic modeling / flux balance analysis | `sybil` | candidate | |
| Multi-omics integration (alternative methods) | `mixOmics` (DIABLO), `MOFA2` | candidate | Alternative/complementary to the existing RGCCA-based module -- worth a second integration method rather than replacing RGCCA |

## 8. Cancer genomics

| Pipeline | R package(s) | Status | Notes |
|---|---|---|---|
| Mutational signatures | `MutationalPatterns` | candidate | |
| Clonal evolution / subclonal architecture | `sciClone` | candidate | |
| Copy number / purity-ploidy | `ASCAT`, `sequenza` | candidate | |
| Immune repertoire (TCR/BCR) | `immunarch` | candidate | |

## 9. Imaging / digital pathology

| Pipeline | R package(s) | Status | Notes |
|---|---|---|---|
| High-content image analysis | `EBImage` (Bioconductor) | candidate | |

## 10. Workflow management

Not applicable -- Nextflow/Snakemake/WDL/Galaxy/CWL are pipeline orchestration
engines in their own right, not R packages omicsuite would wrap.

## How this list gets worked

This isn't a queue -- entries get promoted from `candidate` to `planned` to
`built` one at a time, based on what's actually useful next (a thesis
deadline, a portfolio article, a real dataset in hand), not the order they
appear here. See NEWS.md for what's shipped so far.
