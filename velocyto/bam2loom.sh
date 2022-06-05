

## following https://github.com/basilkhuder/Seurat-to-RNA-Velocity
## run the code directly not using the velocyto run10X


## Chow1
velocyto run -o loom/output/chow1 -m reference/mm10_rmsk.gtf bam/Chow_1/possorted_genome_bam.bam  Tools/CellRanger/refdata-gex-mm10-2020-A/genes/genes.gtf

## Chow2
velocyto run -o loom/output/chow2/ -m reference/mm10_rmsk.gtf bam/chow_2/possorted_genome_bam.bam  Tools/CellRanger/refdata-gex-mm10-2020-A/genes/genes.gtf

## Chow3
velocyto run -o loom/output/chow3/ -m reference/mm10_rmsk.gtf bam/chow_3/possorted_genome_bam.bam  Tools/CellRanger/refdata-gex-mm10-2020-A/genes/genes.gtf


## NASH1
velocyto run -o loom/output/nash1/ -m reference/mm10_rmsk.gtf bam/NASH_1/possorted_genome_bam.bam  Tools/CellRanger/refdata-gex-mm10-2020-A/genes/genes.gtf


## NASH2
velocyto run -o loom/output/nash2/ -m reference/mm10_rmsk.gtf bam/NASH_2/possorted_genome_bam.bam  Tools/CellRanger/refdata-gex-mm10-2020-A/genes/genes.gtf

## NASH3
velocyto run -o loom/output/nash3/ -m reference/mm10_rmsk.gtf bam/NASH_3/possorted_genome_bam.bam  Tools/CellRanger/refdata-gex-mm10-2020-A/genes/genes.gtf
