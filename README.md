<p align="center">
  <img src="inst/app/www/shinypopgen_logo and name.svg" alt="ShinyPopGen" width="420"/>
</p>

Interactive R/Shiny application for population genetics analyses (IRD · CIRAD · INTERTRYP).

## Installation

All dependencies install automatically:

```r
install.packages("remotes")

remotes::install_github("KadNaf/SPG1")

shinypopgen::run_app()
```

### macOS prerequisites (compile from source)

ShinyPopGen contains C++ code compiled with OpenMP. Apple clang does **not** include OpenMP or gfortran by default. Install these **before** running `remotes::install_git()`:

**1. gfortran** — download from <https://mac.r-project.org/tools/> and install `gfortran-14.2-universal.pkg`.

**2. libomp** — OpenMP runtime:

```bash
brew install libomp
```

**3. Apple Silicon only** — add to `~/.R/Makevars` (create the file if absent):

```
LDFLAGS += -L/opt/homebrew/opt/libomp/lib -lomp
CPPFLAGS += -I/opt/homebrew/opt/libomp/include -Xclang -fopenmp
```

## Docker

```bash
cd shinypopgen
docker compose up
# open http://localhost:3838
```

## Features

- Data import (CSV/TXT, single- or two-column allele encoding)
- Allele frequencies and missing data overview
- General statistics (Na, Ne, Ho, He, F-statistics WC84)
- Local and global panmixia (FIS, FIT with bootstrap CI and permutation p-values)
- Population subdivision (FST WC84)
- Genetic diversities (HS, HT — Nei 1987)
- Linkage disequilibrium (pairwise permutation tests)
- Null allele screening (FreeNA EM algorithm)
- Interactive map when GPS metadata is available

## Citation

> IRD / CIRAD / INTERTRYP.

## Bugs & support

<https://forge.ird.fr/intertryp/shiny_pop_gen/-/issues>
