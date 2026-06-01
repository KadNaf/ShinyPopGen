# ── 03_deploy.R ──────────────────────────────────────────────────────────
# Deployment helpers for Posit Connect / shinyapps.io / Docker.

# ── Option A: rsconnect (shinyapps.io or Posit Connect) ──────────────────
# rsconnect::deployApp(
#   appDir      = ".",
#   appName     = "shinypopgen",
#   appTitle    = "ShinyPopGen",
#   forceUpdate = TRUE
# )

# ── Option B: build source package for server installation ────────────────
# devtools::build(".")                   # creates shinypopgen_x.y.z.tar.gz
# install.packages("shinypopgen_x.y.z.tar.gz", repos = NULL, type = "source")

# ── Option C: Docker (see Dockerfile in project root) ─────────────────────
# From the shinypopgen/ directory:
#   docker build -t shinypopgen:latest .
#   docker run --rm -p 3838:3838 shinypopgen:latest
