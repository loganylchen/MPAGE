test_that("get_rna_mod_proteins returns correct structure", {
  proteins <- get_rna_mod_proteins()
  
  expect_true(is.data.frame(proteins))
  expect_true(ncol(proteins) >= 4)
  expect_true(all(c("gene_symbol", "modification_type", "functional_role", "evidence_source") %in% colnames(proteins)))
})

test_that("get_rna_mod_proteins filters work correctly", {
  proteins <- get_rna_mod_proteins(
    modification_types = c("m6A", "m5C"),
    include_writers = TRUE,
    include_readers = FALSE,
    include_erasers = FALSE
  )
  
  expect_true(all(proteins$functional_role == "writer"))
})

test_that("add_custom_proteins adds custom proteins", {
  base_proteins <- get_rna_mod_proteins()
  custom_proteins <- data.frame(
    gene_symbol = "CUSTOM1",
    modification_type = "m6A",
    functional_role = "writer",
    evidence_source = "Custom"
  )
  
  enhanced_proteins <- add_custom_proteins(base_proteins, custom_proteins)
  
  expect_true("CUSTOM1" %in% enhanced_proteins$gene_symbol)
  expect_true(nrow(enhanced_proteins) > nrow(base_proteins))
})

test_that("save_rna_mod_proteins creates file", {
  proteins <- get_rna_mod_proteins()
  temp_file <- tempfile(fileext = ".csv")
  
  expect_silent(save_rna_mod_proteins(proteins, temp_file))
  expect_true(file.exists(temp_file))
  
  unlink(temp_file)
})