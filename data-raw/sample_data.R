# Sample data generation for MPAGE package
set.seed(123)

# Sample expression data
sample_expression_data <- matrix(rnorm(500 * 20), nrow = 500)
rownames(sample_expression_data) <- paste0("GENE", 1:500)
colnames(sample_expression_data) <- paste0("SAMPLE", 1:20)

# Add some structure to make it more realistic
# Create 5 modules with correlated genes
module_genes <- list(
  module1 = paste0("GENE", 1:50),
  module2 = paste0("GENE", 51:100),
  module3 = paste0("GENE", 101:150),
  module4 = paste0("GENE", 151:200),
  module5 = paste0("GENE", 201:250)
)

# Make genes within modules more correlated
for (module in module_genes) {
  base_expr <- rnorm(20, mean = 5, sd = 1)
  for (gene in module) {
    if (gene %in% rownames(sample_expression_data)) {
      sample_expression_data[gene, ] <- base_expr + rnorm(20, mean = 0, sd = 0.5)
    }
  }
}

# Sample modules
sample_modules <- list(
  list(
    module_id = "m6A_writer",
    genes = c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5"),
    size = 5,
    classification = "m6A"
  ),
  list(
    module_id = "m6A_reader", 
    genes = c("GENE6", "GENE7", "GENE8", "GENE9", "GENE10"),
    size = 5,
    classification = "m6A"
  ),
  list(
    module_id = "m5C_writer",
    genes = c("GENE11", "GENE12", "GENE13", "GENE14", "GENE15"),
    size = 5,
    classification = "m5C"
  )
)

# Create sample metadata
sample_metadata <- data.frame(
  sample_id = paste0("SAMPLE", 1:20),
  group = rep(c("Control", "Treatment"), each = 10),
  stringsAsFactors = FALSE
)

# Save the data to the package
dir.create("data", showWarnings = FALSE)
save(sample_expression_data, file = "data/sample_expression_data.rda")
save(sample_modules, file = "data/sample_modules.rda")
save(sample_metadata, file = "data/sample_metadata.rda")

message("Sample data created successfully")