test_dir="/project/shared/DSSR/s229294/hoshida_lab_rnaseq/hbv_fusion/FastViFi"
repo_dir="/project/shared/DSSR/s229294/hoshida_lab_rnaseq/hbv_fusion/FastViFi/FastViFi"
out_dir="/project/shared/DSSR/s229294/hoshida_lab_rnaseq/hbv_fusion/FastViFi/testing"

module load gcc/12.2.0 # 5 or greater
# cd ${test_dir}/kraken2
# wget "https://ftp.ncbi.nih.gov/toolbox/ncbi_tools++/CURRENT/ncbi_cxx--25_2_0.tar.gz"
# gunzip -c ncbi_cxx--25_2_0.tar.gz | tar xvf -
# cd ncbi_cxx--25_2_0
# ./configure --without-debug --with-optimization --with-projects=src/app/dustmasker/dustmasker.lst
# cd GCC610-ReleaseMT64/build
# make all_p

cd ${test_dir}
wget "https://github.com/ncbi/ncbi-cxx-toolkit-public/archive/refs/tags/release-28.0.3.tar.gz"
gunzip -c release-28.0.3.tar.gz | tar xvf -
cd ncbi-cxx-toolkit-public-release-28.0.3
./configure --without-debug --with-optimization --with-projects=src/app/dustmasker/dustmasker.lst
cd GCC1220-ReleaseMT64/build
make all_p


# singularity shell ${out_dir}/fastvifi_v1.1.sif /bin/bash
# cd ${test_dir}/kraken2_build/kraken2/
conda activate run_fastvifi
cd ${test_dir}/kraken2
bash build_custom_kraken_database.sh \
    "hbv" \
    ${test_dir}/ViFi/viral_data/hbv/hbv.unaligned.fas \
    9000000

# check if 3 required files are present: hash.k2d, opts.k2d, and taxo.k2d
./kraken2-inspect --db Kraken2StandardDB_k_25_hbv_hg
./kraken2-inspect --db Kraken2StandardDB_k_18_hbv
./kraken2-inspect --db Kraken2StandardDB_k_22_hbv

# fix hbv taxid issue
conda activate gs
gdown "https://drive.google.com/uc?id=1VUjxYXxNgjzwla8BQqF9Kt9hrlvbrrdo"
conda deactivate
cd ${test_dir}/kraken2

bash build_custom_kraken_database.sh \
    "hbv" \
    ./hbv.unaligned_with_tax_id.fasta \
    9000000
