# This file is part of the standard testthat setup for R packages.
# It is run by R CMD check / devtools::test() and by the CI "test" stage
# added in .gitlab-ci.yml (see docker-test job).

library(testthat)
library(shinypopgen)

test_check("shinypopgen")
