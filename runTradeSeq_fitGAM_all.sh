#!/bin/bash

for rds in linA_sce.rds linB_sce.rds precc1_cc1A_sce.rds precc1_cc1B1_sce.rds cc1B1_cc1B2_sce.rds precc2_cc2_sce.rds
do
	job2pbs.py -l "Rscript tradeSeq_fitGAM_sling.R ${rds}"
done

#for rds in mc3_opc_precc1.rds mc3_linA.rds mc3_linB.rds mc3_precc1_cc1A.rds mc3_precc1_cc1B1.rds mc3_cc1B1_cc1B2.rds mc3_precc2_cc2.rds
#do
#	job2pbs.py -l "Rscript tradeSeq_fitGAM_mc3.R ${rds}"
#done

