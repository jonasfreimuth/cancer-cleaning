# To be run in a subshell from project root.
set -euo pipefail

mkdir -p datasets

cd datasets

DL_URL="ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/suppl/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz"

wget -qO- "$DL_URL" | tar -xz

echo "Downloaded and extracted archive from
$DL_URL
on $(date)."
