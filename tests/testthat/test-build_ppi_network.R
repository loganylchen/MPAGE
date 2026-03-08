# Test build_ppi_network function
# These tests will need implementation once the function exists

test_that("build_ppi_network validates min_confidence parameter", {
  skip("Function not yet implemented")
  expect_error(
    build_ppi_network(proteins = "TP53", min_confidence = 1.5),
    "min_confidence must be between 0 and 1"
  )
  expect_error(
    build_ppi_network(proteins = "TP53", min_confidence = -0.1),
    "min_confidence must be between 0 and 1"
  )
})

test_that("build_ppi_network requires at least one protein", {
  skip("Function not yet implemented")
  expect_error(
    build_ppi_network(proteins = character(0)),
    "must provide at least one protein"
  )
})

test_that("build_ppi_network returns igraph object", {
  skip("Function not yet implemented")
  result <- build_ppi_network(
    proteins = c("TP53", "MDM2"),
    data_sources = "STRING",
    min_confidence = 0.4
  )
  expect_s3_class(result, "igraph")
})
