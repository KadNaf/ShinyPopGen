# Smoke tests
#
# These do NOT validate scientific results (see test-numeric-validation-*.R
# for that). They only check that the app assembles without erroring, which
# is the minimum a CI "test" stage should verify before an image is built and
# deployed (previously the pipeline had no test stage at all - see
# .gitlab-ci.yml).

test_that("the package namespace loads and exports the expected entry points", {
  expect_true(is.function(shinypopgen::run_app))
  expect_true(is.function(shinypopgen::app_ui))
  expect_true(is.function(shinypopgen::app_server))
})

test_that("app_ui() builds without error and returns renderable HTML", {
  ui <- shinypopgen::app_ui()
  expect_false(is.null(ui))
  # bslib::page_navbar() returns tag-like content; htmltools must be able to
  # render it to HTML without raising, which is what actually matters here.
  expect_type(htmltools::doRenderTags(ui), "character")
})

test_that("app_server() initializes its reactiveValues without error", {
  testServer(shinypopgen::app_server, {
    # rv is the reactiveValues() created at the top of app_server(); testServer
    # evaluates assertions in that function's local environment, so it is
    # directly in scope here.
    expect_null(rv$raw)
    expect_null(rv$data)
    expect_null(rv$formatted_data)
  })
})
