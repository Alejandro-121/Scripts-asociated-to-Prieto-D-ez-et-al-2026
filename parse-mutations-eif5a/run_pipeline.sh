#!/bin/bash
# VCF Variant Analysis Pipeline
# Pipeline for variant extraction, genotype-based filtering, and visualization

set -e  # Exit on error

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored messages
print_message() {
    echo -e "${GREEN}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

# Check if required files exist
check_dependencies() {
    print_message "Checking dependencies..."
    
    if [ ! -f "merged.vcf" ]; then
        print_error "merged.vcf not found in current directory"
        exit 1
    fi
    
    # Check Python
    if ! command -v python3 &> /dev/null; then
        print_error "Python3 is not installed"
        exit 1
    fi
    
    # Check R
    if ! command -v Rscript &> /dev/null; then
        print_error "R is not installed"
        exit 1
    fi
    
    print_message "All dependencies found!"
}

# Step 1: Extract variants from VCF
run_variant_extraction() {
    print_message "============================================"
    print_message "STEP 1: Extracting variants from VCF file"
    print_message "============================================"
    
    if [ ! -f "01_variant_extraction.py" ]; then
        print_error "01_variant_extraction.py not found"
        exit 1
    fi
    
    python3 01_variant_extraction.py
    
    if [ -f "mutation_wide_table_with_alleles.csv" ]; then
        print_message "✓ Variant extraction completed successfully"
        print_message "  Output: mutation_wide_table_with_alleles.csv"
    else
        print_error "Variant extraction failed - output file not created"
        exit 1
    fi
}

# Step 2: Genotype-based filtering
run_genotype_filtering() {
    print_message "============================================"
    print_message "STEP 2: Genotype-based filtering"
    print_message "============================================"
    
    if [ ! -f "02_genotype_based_filtering.R" ]; then
        print_error "02_genotype_based_filtering.R not found"
        exit 1
    fi
    
    if [ ! -f "mutation_wide_table_with_alleles.csv" ]; then
        print_error "mutation_wide_table_with_alleles.csv not found (run step 1 first)"
        exit 1
    fi
    
    Rscript 02_genotype_based_filtering.R
    
    if [ -f "common_mutations_2.1.csv" ]; then
        print_message "✓ Genotype filtering completed successfully"
        print_message "  Generated filtered tables and plots"
    else
        print_error "Genotype filtering failed"
        exit 1
    fi
}

# Main pipeline execution
main() {
    print_message "======================================================"
    print_message "VCF Variant Analysis Pipeline"
    print_message "======================================================"
    
    # Check dependencies
    check_dependencies
    
    # Run pipeline steps
    run_variant_extraction
    run_genotype_filtering
    
    print_message "======================================================"
    print_message "Pipeline completed successfully!"
    print_message "======================================================"
    print_message ""
    print_message "Generated files:"
    print_message "  1. mutation_wide_table_with_alleles.csv"
    print_message "  2. common_mutations_2.1.csv"
    print_message "  3. common_mutations_2.3.csv"
    print_message "  4. common_reversions_2.1.csv"
    print_message "  5. common_reversions_2.3.csv"
    print_message "  6. common_mutations_101_2.1.csv"
    print_message "  7. common_mutations_010_2.1.csv"
    print_message "  8. common_mutations_101_2.3.csv"
    print_message "  9. common_mutations_010_2.3.csv"
    print_message " 10. mutation_frequency_010_2.1.pdf"
    print_message " 11. mutation_heatmap_101_2.1.pdf"
}

# Run main function
main
