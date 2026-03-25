# ===============================================================
# Build a custom BSgenome package directly from NCBI assembly
# ===============================================================

# -------------------------
# 0. Load required libraries
# -------------------------
library(BSgenomeForge)
library(devtools)

# -------------------------
# 1. Set working directory for genome package creation
# -------------------------
setwd("~/Escritorio/Genomes/")  # output directory for the generated BSgenome package

# -------------------------
# 2. Create BSgenome data package from NCBI assembly
# -------------------------
# This will automatically download the FASTA sequences and metadata
# from NCBI for the specified assembly and prepare a BSgenome package structure.
forgeBSgenomeDataPkgFromNCBI(
  assembly_accession = "GCF_016699485.2",   # NCBI assembly accession for Gallus gallus
  pkg_maintainer = "Alejandro Gonzalez <...@umh.es>",  # package maintainer info
  organism = "Gallus gallus"                # organism name
)

# -------------------------
# 3. Define the generated package directory
# -------------------------
# The directory name is automatically created based on the assembly and organism.
pkg_dir <- "./BSgenome.Ggallus.NCBI.bGalGal1.mat.broiler.GRCg7b"

# -------------------------
# 4. Build the package tarball (equivalent to 'R CMD build')
# -------------------------
tarball <- devtools::build(pkg_dir)

# -------------------------
# 5. (Optional) Run checks on the built package
# -------------------------
# This step validates that the package structure and metadata are correct.
devtools::check_built(tarball)

# -------------------------
# 6. Install the package locally from the generated tar.gz file
# -------------------------
devtools::install_local(tarball)
