#!/bin/bash

# Post-create script for MPAGE development container

echo "Setting up MPAGE development environment..."

# Set git safe directory
git config --global --add safe.directory /home/vscode/workspace

# Install package dependencies
echo "Installing package dependencies..."
R -e "
if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')
BiocManager::install(c('STRINGdb', 'clusterProfiler', 'GSEABase', 'GSVA', 'Biobase', 'org.Hs.eg.db', 'AnnotationDbi'), ask=FALSE, update=FALSE)
"

# Install additional development packages
echo "Installing development tools..."
R -e "
install.packages(c('covr', 'pkgdown', 'styler', 'lintr', 'spelling'), repos='https://cran.rstudio.com/')
"

# Build and install the package locally for development
echo "Building MPAGE package..."
cd /home/vscode/workspace
R CMD build .

# Install the package for development
echo "Installing MPAGE package for development..."
R CMD INSTALL MPAGE_*.tar.gz

# Install git pre-commit hooks (optional)
if [ -f .pre-commit-config.yaml ]; then
    echo "Setting up pre-commit hooks..."
    pip install pre-commit
    pre-commit install
fi

# Set up R environment for development
echo "Configuring R development environment..."
mkdir -p ~/.R

cat > ~/.Rprofile << 'EOF'
# Custom R profile for MPAGE development
options(repos = c(CRAN = "https://cran.rstudio.com/"))
options(devtools.desc_author = "Logan Chen <loganylchen@gmail.com>")
options(devtools.desc_license = "MIT")

# Load devtools on startup for development
if (interactive()) {
  suppressMessages(require(devtools))
  suppressMessages(require(testthat))
  cat("MPAGE development environment ready!\n")
}
EOF

# Create useful aliases for development
echo "Creating development aliases..."
cat > ~/.bash_aliases << 'EOF'
alias rcheck='R CMD check --as-cran'
alias rbuild='R CMD build .'
alias rinstall='R CMD install'
alias rtest='R -e "devtools::test()"'
alias rcoverage='R -e "covr::package_coverage()"'
alias rdemo='R -e "devtools::run_examples()"'
alias rdoc='R -e "devtools::document()"'
EOF

# Set proper permissions
chmod +x ~/.bash_aliases

# Create a welcome message
echo "Creating welcome message..."
cat > ~/.welcome.txt << 'EOF'
╔══════════════════════════════════════════════════════════════════════════════╗
║                           MPAGE Development Container                        ║
╠══════════════════════════════════════════════════════════════════════════════╣
║  Welcome to the MPAGE R package development environment!                     ║
║                                                                              ║
║  Available commands:                                                         ║
║  • rcheck - Run R CMD check                                                  ║
║  • rbuild - Build package                                                    ║
║  • rtest - Run tests                                                         ║
║  • rcoverage - Run coverage analysis                                         ║
║  • rdoc - Generate documentation                                             ║
║                                                                              ║
║  The package is already installed and ready for development.                 ║
║  Start editing files and use devtools::load_all() to reload changes.         ║
╚══════════════════════════════════════════════════════════════════════════════╝
EOF

# Add welcome message to bash profile
echo "cat ~/.welcome.txt" >> ~/.bashrc

echo "✅ MPAGE development environment setup complete!"
echo "Run 'devcontainer open' or reopen in container to start developing."