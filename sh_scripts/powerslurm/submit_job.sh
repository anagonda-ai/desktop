#!/bin/bash
#SBATCH --job-name=statistical_scanner
#SBATCH --output=/groups/itay_mayrose/alongonda/desktop/example_jobs/out.OU
#SBATCH --error=/groups/itay_mayrose/alongonda/desktop/example_jobs/error.ER
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --partition=itaym
#SBATCH --time=24:00:00

conda init
conda deactivate
conda activate base

cd /groups/itay_mayrose/alongonda/tools/plant_mgc_existing_tools/PlantClusterFinder/PlantClusterFinder

./run_PlantClusterFinder.sh /usr/local.cc/lib64/matlabr2024b/runtime/glnxa64/   -pgdb "/groups/itay_mayrose/alongonda/tools/plant_mgc_existing_tools/PlantClusterFinder/PlantClusterFinder/csubellipsoidea/pgdb/csubellipsoideacyc/1.0/data"   -rmdf "/groups/itay_mayrose/alongonda/tools/plant_mgc_existing_tools/PlantClusterFinder/PlantClusterFinder/Inputs/ReactionMetabolicDomainClassification.txt"   -md "{'Amines and Polyamines Metabolism'; 'Amino Acids Metabolism'; 'Carbohydrates Metabolism'; 'Cofactors Metabolism'; 'Detoxification Metabolism'; 'Energy Metabolism'; 'Fatty Acids and Lipids Metabolism'; 'Hormones Metabolism'; 'Inorganic Nutrients Metabolism'; 'Nitrogen-Containing Compounds'; 'Nucleotides Metabolism'; 'Phenylpropanoid Derivatives'; 'Polyketides'; 'Primary-Specialized Interface Metabolism'; 'Redox Metabolism'; 'Specialized Metabolism'; 'Sugar Derivatives'; 'Terpenoids'}"   -psf "/groups/itay_mayrose/alongonda/tools/plant_mgc_existing_tools/PlantClusterFinder/PlantClusterFinder/csubellipsoidea/CsubellipsoideaC_169_227_v2.0.protein.pcf13.fa"   -gtpf "/groups/itay_mayrose/alongonda/tools/plant_mgc_existing_tools/PlantClusterFinder/PlantClusterFinder/csubellipsoidea/gtpf_CsubellipsoideaC_169_227_v2.0.annotation_info.txt.txt"   -glof "/groups/itay_mayrose/alongonda/tools/plant_mgc_existing_tools/PlantClusterFinder/PlantClusterFinder/csubellipsoidea/glof_CsubellipsoideaC_169_227_v2.0.gene.gff3.txt"   -dnaf "/groups/itay_mayrose/alongonda/tools/plant_mgc_existing_tools/PlantClusterFinder/PlantClusterFinder/csubellipsoidea/CsubellipsoideaC_169_227_v2.0.hardmasked.fa"   -sitf "/groups/itay_mayrose/alongonda/tools/plant_mgc_existing_tools/PlantClusterFinder/PlantClusterFinder/Inputs/scaffold-tailoring-reactions-05082016.tab"   -gout "/groups/itay_mayrose/alongonda/tools/plant_mgc_existing_tools/PlantClusterFinder/PlantClusterFinder/results/csubellipsoidea_gene_output.tsv"   -cout "/groups/itay_mayrose/alongonda/tools/plant_mgc_existing_tools/PlantClusterFinder/PlantClusterFinder/results/csubellipsoidea_cluster_output.tsv"   SeqGapSizesChromBreak "[10000]" PGDBIdsToMap GTP