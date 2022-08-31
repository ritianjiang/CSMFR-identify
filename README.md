# CSMFR-identify
The CSM-FR identification script

The CSMFR.R contain the pipeline. Users need to read the mutation location file into R session as GRanges object first, then run the commands in CSMFR.R.

Mutect2-PoN.sh is the script for tumor-only mutation calling. The preprocess pipeline can be found at [https://github.com/Konglab404/Ritian_GATK_preprocess]

Filter.sh is the script for mutation filtering.
