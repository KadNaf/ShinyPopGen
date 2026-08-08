<p align="center">
  <img src="inst/app/www/shinypopgen_logo%20and%20name.svg" alt="ShinyPopGen" width="420"/>
</p>

<p align="center">
  Interactive R/Shiny application for population genetics analyses<br/>
  <sub>IRD · UCAD · CIRAD · INTERTRYP</sub>
</p>

ShinyPopGen covers the full workflow for multilocus, individual-based
genotype datasets (microsatellites and other codominant markers): data
import and formatting, allele frequencies, general statistics, panmixia
and subdivision testing (Weir & Cockerham 1984), genetic diversities,
linkage disequilibrium, null allele screening (FreeNA), and isolation by
distance — all computed locally, with no data leaving your machine.

## Installation

All dependencies install automatically:

```r
install.packages("remotes")

remotes::install_github("KadNaf/SPG1")

shinypopgen::run_app()
```

### macOS prerequisites (compile from source)

ShinyPopGen contains C++ code compiled with OpenMP. Apple clang does **not**
include OpenMP or gfortran by default. Install these **before** running
`remotes::install_github()`:

**1. gfortran** — download from <https://mac.r-project.org/tools/> and
install `gfortran-14.2-universal.pkg` (or the current version).

**2. libomp** — OpenMP runtime:

```bash
brew install libomp
```

**3. Apple Silicon only** — add to `~/.R/Makevars` (create the file if
absent):

```
LDFLAGS += -L/opt/homebrew/opt/libomp/lib -lomp
CPPFLAGS += -I/opt/homebrew/opt/libomp/include -Xclang -fopenmp
```

## Docker

```bash
git clone https://github.com/KadNaf/SPG1.git
cd SPG1
docker compose up
# open http://localhost:3838
```

## Features

| Module | What it does |
|---|---|
| **Import Data** | CSV/TXT import (auto-detected comma/semicolon/tab separators), single- or two-column allele encoding, auto- or manual column assignment, sampling-locality map. 500 MB upload limit. |
| **Allele Freq** | Per-locus allele frequency tables/plots and missing-data overview. |
| **General Stats** | N, Na, Ne, Ho, He per locus and per population. |
| **Local Panmixia** | Within-population FIS (Weir & Cockerham 1984), bootstrap CI and permutation p-value. |
| **Global Panmixia** | Multilocus FIT across all populations, bootstrap CI and permutation p-value. |
| **Subdivision** | FST (Weir & Cockerham 1984) per locus and overall, bootstrap CI (loci or population blocks), permutation and G-test p-values. |
| **Diversities** | HS and HT (Nei 1987) per locus and multilocus, with resampling over individuals, populations, or loci. |
| **LD** | Pairwise linkage disequilibrium between all locus pairs, G-test with permutation p-values. |
| **Null Alleles** | Null allele frequency estimation per locus × population (FreeNA EM algorithm, Chapuis & Estoup 2007), with raw and ENA-corrected FST. |
| **Isolation by Distance** | Pairwise FST/(1−FST) vs. geographic distance (Rousset 1997) and a Mantel permutation test. |

Every results table can be exported (.csv / .txt).

## Getting started

See the package vignette for a full worked example on the bundled
*Boophilus* tick dataset:

```r
vignette("shinypopgen", package = "shinypopgen")
```

## Statistical methods

F-statistics follow the unbiased moment estimators of **Weir & Cockerham
(1984)**. Confidence intervals are obtained by non-parametric bootstrap and
p-values by Monte Carlo permutation (5,000 replicates by default; 10,000
for linkage disequilibrium), parallelised in C++ via Rcpp and OpenMP. See
the in-app **Help** tab for full references.

## Citation

> ShinyPopGen: an interactive Shiny application for population genetics
> data import, exploration, and descriptive analyses. IRD / CIRAD /
> INTERTRYP.

## Credits

**Programming:** Naffiou Kadiri and Vincent Manzanilla
**Conception:** Thierry de Meeûs

## License

MIT — see [LICENSE.md](LICENSE.md).

## Bugs & support

<https://github.com/KadNaf/SPG1/issues>
