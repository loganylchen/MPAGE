test_that("RNA modification proteins data is loaded correctly", {
  proteins <- get_rna_mod_proteins()
  
  expect_true(is.data.frame(proteins))
  expect_true(ncol(proteins) >= 7)
  expect_true(all(c("gene_symbol", "modification_type", "functional_role", 
                   "evidence_source", "pmid", "reference", "organism") %in% colnames(proteins)))
})

test_that("RNA modification proteins has expected columns", {
  proteins <- get_rna_mod_proteins()
  
  expect_type(proteins$gene_symbol, "character")
  expect_type(proteins$modification_type, "character")
  expect_type(proteins$functional_role, "character")
  expect_type(proteins$evidence_source, "character")
  expect_type(proteins$pmid, "integer")
  expect_type(proteins$reference, "character")
  expect_type(proteins$organism, "character")
})

test_that("modification types are correct", {
  mod_types <- list_modification_types()
  expected_types <- c("m6A", "m5C", "Psi", "m1A", "m7G", "A-to-I")
  
  expect_true(all(expected_types %in% mod_types))
})

test_that("references can be retrieved", {
  refs <- get_rna_mod_references()
  
  expect_true(is.data.frame(refs))
  expect_true(ncol(refs) == 2)
  expect_true(all(c("reference", "pmid") %in% colnames(refs)))
})

test_that("organism filtering works", {
  human_proteins <- get_rna_mod_proteins(organism = "Homo sapiens")
  
  expect_true(all(human_proteins$organism == "Homo sapiens"))
})

test_that("functional role filtering works", {
  writers_only <- get_rna_mod_proteins(
    include_writers = TRUE,
    include_readers = FALSE,
    include_erasers = FALSE
  )
  
  expect_true(all(writers_only$functional_role == "writer"))
})

test_that("modification type filtering works", {
  m6a_only <- get_rna_mod_proteins(modification_types = "m6A")
  
  expect_true(all(m6a_only$modification_type == "m6A"))
})