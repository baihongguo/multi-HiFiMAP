# multi-HiFiMAP
Multi-phenotype high-resolution fast identity-by-descent (IBD) mapping test

Current version: 1.1.0 ([version history](#version-history))

## Overview
**multi-HiFiMAP** tests whether local IBD sharing at a SNP is associated with many phenotypes jointly, at single-SNP testing resolution. It builds on [HiFiMAP](#citation), our single-phenotype IBD mapping test.

The test statistic is the RV coefficient between the local IBD matrix and a phenotype kernel built from the scaled residuals of all phenotypes (the double-kernel association test, DKAT). P-values come from a Pearson type III approximation to the permutation distribution, matched on its first three moments, so no permutations are needed.

Three algorithms make this feasible genome-wide in large samples:
- **Sparse algebra:** every moment is computed from the non-zero elements of the local IBD matrix, so no dense *n* × *n* matrix is ever formed.
- **Randomized traces:** the third-order trace is estimated with Hutchinson's estimator, using sparse matrix–vector products.
- **Dynamic updating:** adjacent SNPs differ by a sparse change in the IBD matrix, so the moments are updated rather than recomputed.

A stateful C++ backend (`RcppArmadillo`) keeps the sparse IBD matrix in memory as the scan moves along the chromosome.

## Installation
```bash
git clone https://github.com/baihongguo/multi-HiFiMAP.git
cd multi-HiFiMAP
```
There is no separate build step: the R script compiles the C++ backend with `Rcpp::sourceCpp()` when it starts. This needs a C++ compiler (e.g. `g++`).

## Requirements
The software is developed and tested on Linux.

**Python (>= 3.6):** `numpy`

**R (>= 4.0.0):** `data.table`, `Matrix`, `Rcpp`, `RcppArmadillo`, `PearsonDS`, `MASS`
```r
install.packages(c("data.table", "Matrix", "Rcpp", "RcppArmadillo", "PearsonDS", "MASS"))
```

## Repository layout
```text
example/                 toy data: 5,000 individuals, 128 phenotypes, chromosome 21
src/hapIBD_parsing/      parsing_hapIBD.py: hap-ibd segments -> sparse checkpoints and per-SNP differences
src/multi-HiFiMAP/       run_multi_HiFiMAP.R (scan of one chunk) and multi-HiFiMAP_helper.cpp (C++ backend)
src/Shell script/        Run_multi-HiFiMAP_parallel.sh (runs all chunks in parallel and merges the results)
```

## Usage
Run all commands from the repository root. The toy example uses the files in `example/`.

### Step 0: Parse the IBD segments (Python)
```bash
python3 src/hapIBD_parsing/parsing_hapIBD.py \
    --ibd example/chr21_toy.ibd.gz \
    --vcf example/chr21_toy.vcf.gz \
    --output example/ibd_prep/chr21 \
    --n-checkpoints 20
```
- **`--ibd`** is the output of [hap-ibd](https://github.com/browning-lab/hap-ibd): plain or gzipped, with columns sample 1, haplotype 1, sample 2, haplotype 2, chromosome, start, end, length in cM. Run hap-ibd with the IBD length threshold you want to test (we use 1 cM).
- **`--vcf`** gives the SNP positions and the sample order. The toy VCF is too large for GitHub: download [chr21_toy.vcf.gz](https://drive.google.com/file/d/1TcwBfpEx_jzbr3zR4XcmNEBJfPWfiXO9/view?usp=drive_link) and place it in `example/`.
- **`--output`** should be one folder per chromosome. It receives the checkpoint matrices (`ibd_mat_<site>.mtx`), the per-SNP differences (`delta_<site>.diff`), `sites.txt` and `samples.txt`.
- **`--n-checkpoints`** sets how many checkpoint matrices are written. Use the same number as `NUM_CHUNKS` in Step 2, so that each chunk starts close to a checkpoint.
- **Optional:** `--mask` excludes regions (a BED file or `chr:start-end`), and `--resume` continues an interrupted run.

### Step 1: Prepare the phenotype files
multi-HiFiMAP works directly on pre-computed phenotype residuals and their correlation matrix; no fitted null-model object is needed. Prepare three plain-text files without headers:

1. **Residuals (`toy_lmscaledresiduals.txt`):** an *N* × *K* matrix of scaled residuals for *N* individuals and *K* phenotypes, e.g. from a linear or linear mixed model of each phenotype on its covariates. Residuals are mean-centered by construction.
2. **Correlation matrix (`toy_cor_matrix.txt`):** the *K* × *K* correlation matrix of the residuals. It must be invertible, because the script uses its inverse.
3. **Subject IDs (`toy_subject_ids.txt`):** one ID per line, in the row order of the residuals file. The IDs must match the sample IDs in the VCF (`samples.txt`).

The individuals in the residuals file can be a subset of the samples in the IBD data: the R script maps the IBD matrices to the IDs, in the given order. If the residual rows already follow `samples.txt` exactly, you can set the ID file to `NONE`.

### Step 2: Run the parallel scan
Edit the configuration block at the top of `src/Shell script/Run_multi-HiFiMAP_parallel.sh`:
- input and output paths (`IBD_BASE_DIR`, `RES_FILE`, `COR_FILE`, `ID_FILE`, `RESULT_DIR`, `LOG_DIR`);
- `chromosomes` and `NUM_CHUNKS`, the number of chunks per chromosome;
- `MAX_JOBS` and `MIN_JOBS`, which control how many chunks run at once. Set them from your cores and memory; see [Computational resources](#computational-resources).

Then run:
```bash
bash "src/Shell script/Run_multi-HiFiMAP_parallel.sh"
```
Each chunk runs `src/multi-HiFiMAP/run_multi_HiFiMAP.R` in the background. When all chunks finish, the script merges them into one file per chromosome. Logs go to `results/logs/`, and each chunk's main-loop run time is written to a `.time` file in `results/`.

To run one chunk directly:
```bash
Rscript src/multi-HiFiMAP/run_multi_HiFiMAP.R <chr> <chunk_start> <chunk_end> <chunk_id> \
    <ibd_dir> <residuals> <cor_matrix> <subject_ids or NONE> <output_file>
```

### Output
`results/Multi-HiFiMAP_chr21.txt` has one row per tested SNP:

```text
chr     pos       n.ibd.segs    p.value
21      14240779  1052          0.0047846
21      14242245  1054          0.0066947
21      14250602  1056          0.0067795
21      14256336  1057          0.0065234
```
- **pos:** base-pair position (from `sites.txt`).
- **n.ibd.segs:** number of IBD segments covering the SNP among the analyzed individuals.
- **p.value:** multi-phenotype p-value, written with 7 significant digits.

## Computational resources
The cost at each SNP grows with the number of IBD-sharing pairs, not with *n*², because no dense *n* × *n* matrix is formed. It therefore rises with sample size and with shorter IBD segments.

In the UK Biobank UDIP cohort (*N* ≈ 36,000; 128 phenotypes; 20 chunks per chromosome; one thread per job):

| IBD threshold | CPU time per SNP | Genome-wide |
|---|---|---|
| 1 cM | about 16 s | 649,793 SNPs in 2,938 CPU hours |
| 2 cM | about 2.6 s | 640,900 SNPs in 460 CPU hours |

These times were measured with version 1.0.0; version 1.1.0 gives identical p-values and is slightly faster.

Peak memory per chunk is typically 5–8 GiB. It reaches about 31 GiB in a few regions where IBD sharing is unusually dense. Plan `MAX_JOBS` around the largest chunks.

We recommend one thread per job (the `OMP_NUM_THREADS=1` settings in the shell script) and parallelizing across chunks.

## Version history
**1.1.0** (October 2026)
- **Lower memory.**
  - The phenotype term *U*<sub>s</sub> = Σ<sub>ij</sub> *w*<sub>ij</sub><sup>3</sup> is now accumulated in blocks of 2,000 rows, so the dense *n* × *n* phenotype kernel *W* = *r*Σ<sup>−1</sup>*r*<sup>T</sup> (about 10 GB at *n* = 36,000, plus a temporary copy) is no longer formed.
  - The C++ update no longer copies the IBD matrix at each SNP.
  - In the UK Biobank UDIP analysis, typical peak memory per chunk fell from about 24 GiB to 5–8 GiB, and the maximum from about 47 GiB to about 31 GiB. The p-values are unchanged.
- **Accurate very small p-values.** When the normal approximation is used (skewness below 10<sup>−6</sup>), the upper tail is now computed directly. The previous `1 - pnorm(Z)` returned exactly 0 for p-values below about 10<sup>−16</sup>.
- **Fixed file paths.**
  - The shell script now calls `src/multi-HiFiMAP/run_multi_HiFiMAP.R`.
  - The R script loads `multi-HiFiMAP_helper.cpp` from its own folder, so it runs from any working directory.
  - The example files are now named `toy_cor_matrix.txt` and `toy_subject_ids.txt`.

**1.0.0**: initial release.

## Citation
If you use multi-HiFiMAP, please cite:
- Guo B, Naseri A, Xie Z, Zhang M, Salehi Nasab S, Sarnowski C, Zhi D, Chen H. Multi-phenotype high-resolution fast identity-by-descent mapping test. Manuscript in preparation; presented at the ASHG 2026 Annual Meeting.
- Guo B, et al. HiFiMAP: high-resolution fast identity-by-descent mapping test. *Am J Hum Genet*, in revision.

## License
GPL-2.0; see [LICENSE](LICENSE).

## Questions and issues
Please open an issue on this repository.
