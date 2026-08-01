# ── 03_deploy.R ──────────────────────────────────────────────────────────
# Deployment helpers for Posit Connect / shinyapps.io / Docker.

# ── Option A: rsconnect (shinyapps.io or Posit Connect) ──────────────────
# rsconnect::deployApp(
#   appDir      = ".",
#   appName     = "SPG1",
#   appTitle    = "SPG1",
#   forceUpdate = TRUE
# )

# ── Option B: build source package for server installation ────────────────
# devtools::build(".")                   # creates SPG1_x.y.z.tar.gz
# install.packages("SPG1_x.y.z.tar.gz", repos = NULL, type = "source")

# ── Option C: Docker (see Dockerfile in project root) ─────────────────────
# From the SPG1/ directory:
#   docker build -t SPG1:latest .
#   docker run --rm -p 3838:3838 SPG1:latest
