#!/bin/bash

npcs=2
ncores=15
subsamples=100
glm="nb"

#job2pbs.py -c ${ncores} -l "Rscript mc3.sling.tradeSeq.pseudotimeDE.R opc_precc1_sce.rds OPC 20 ${ncores} ${subsamples} 6 ${glm}"
#job2pbs.py -c ${ncores} -l "Rscript mc3.sling.tradeSeq.pseudotimeDE.R linA_sce.rds OPC 20 ${ncores} ${subsamples} 6 ${glm}"
#job2pbs.py -c ${ncores} -l "Rscript mc3.sling.tradeSeq.pseudotimeDE.R linB_sce.rds OPC 20 ${ncores} ${subsamples} 6 ${glm}"
#job2pbs.py -c ${ncores} -l "Rscript mc3.sling.tradeSeq.pseudotimeDE.R precc2_cc2_sce.rds pre-CC2 20 ${ncores} ${subsamples} 6 ${glm}"

#job2pbs.py -c ${ncores} -l "Rscript mc3.sling.tradeSeq.pseudotimeDE_on_seurat.R opc_precc1_seurat.rds OPC 20 ${ncores} ${subsamples} 6 ${glm}"
#job2pbs.py -c ${ncores} -l "Rscript mc3.sling.tradeSeq.pseudotimeDE_on_seurat.R linA_seurat.rds OPC 20 ${ncores} ${subsamples} 6 ${glm}"
#job2pbs.py -c ${ncores} -l "Rscript mc3.sling.tradeSeq.pseudotimeDE_on_seurat.R linB_seurat.rds OPC 20 ${ncores} ${subsamples} 6 ${glm}"
#job2pbs.py -c ${ncores} -l "Rscript mc3.sling.tradeSeq.pseudotimeDE_on_seurat.R precc2_cc2_seurat.rds pre-CC2 20 ${ncores} ${subsamples} 6 ${glm}"

#job2pbs.py -c 10 -l "Rscript mc3.sling.tradeSeq.pseudotimeDE.R opc_precc1_sce.rds OPC 20 10 ${subsamples} 6 ${glm}"
#job2pbs.py -c 15 -l "Rscript mc3.sling.tradeSeq.pseudotimeDE.R precc1_cc1A_sce.rds pre-CC1 20 15 ${subsamples} 6 ${glm}"
#job2pbs.py -c 20 -l "Rscript mc3.sling.tradeSeq.pseudotimeDE.R precc1_cc1B_sce.rds pre-CC1 20 20 ${subsamples} 6 ${glm}"
#job2pbs.py -c 20 -l "Rscript mc3.sling.tradeSeq.pseudotimeDE.R linB_sce.rds OPC 20 20 ${subsamples} 6 ${glm}"
#job2pbs.py -c 15 -l "Rscript mc3.sling.tradeSeq.pseudotimeDE.R linA_sce.rds OPC 20 15 ${subsamples} 6 ${glm}"
#job2pbs.py -c 10 -l "Rscript mc3.sling.tradeSeq.pseudotimeDE.R opc_precc1_sce.rds OPC 20 10 ${subsamples} 6 ${glm}"

nohup Rscript mc3.sling.tradeSeq.pseudotimeDE.R opc_precc1_sce.rds OPC 20 5 ${subsamples} 6 ${glm} > opc_precc1.out &
nohup Rscript mc3.sling.tradeSeq.pseudotimeDE.R precc1_cc1A_sce.rds pre-CC1 20 10 ${subsamples} 6 ${glm} > precc1_cc1A.out &
nohup Rscript mc3.sling.tradeSeq.pseudotimeDE.R precc1_cc1B_sce.rds pre-CC1 20 15 ${subsamples} 6 ${glm} > precc1_cc1B.out &

#job2pbs.py -c ${ncores} -l "Rscript mc3.sling.tradeSeq.pseudotimeDE.R linA_sce.rds OPC 2 ${ncores} ${subsamples} 6 ${glm}"
#job2pbs.py -c ${ncores} -l "Rscript mc3.sling.tradeSeq.pseudotimeDE.R linB_sce.rds OPC 2 ${ncores} ${subsamples} 6 ${glm} CC1-B2"
#job2pbs.py -c ${ncores} -l "Rscript mc3.sling.tradeSeq.pseudotimeDE.R linB_sce.rds OPC 2 ${ncores} ${subsamples} 6 ${glm}"
#job2pbs.py -c ${ncores} -l "Rscript mc3.sling.tradeSeq.pseudotimeDE.R precc2_cc2_sce.rds pre-CC2 2 ${ncores} ${subsamples} 6 ${glm}"

#job2pbs.py -c ${ncores} -l "Rscript mc3.sling.tradeSeq.pseudotimeDE.R opc_precc1_sce.rds OPC 2 ${ncores} ${subsamples} 6 ${glm}"
#job2pbs.py -c ${ncores} -l "Rscript mc3.sling.tradeSeq.pseudotimeDE.R linA_sce.rds OPC 2 ${ncores} ${subsamples} 6 ${glm}"
#job2pbs.py -c ${ncores} -l "Rscript mc3.sling.tradeSeq.pseudotimeDE.R linB_sce.rds OPC 2 ${ncores} ${subsamples} 6 ${glm}"
#job2pbs.py -c ${ncores} -l "Rscript mc3.sling.tradeSeq.pseudotimeDE.R precc1_cc1A_sce pre-CC1 2 ${ncores} ${subsamples} 6 ${glm}"
#job2pbs.py -c ${ncores} -l "Rscript mc3.sling.tradeSeq.pseudotimeDE.R precc1_cc1B1_sce pre-CC1 2 ${ncores} ${subsamples} 6 ${glm}"
#job2pbs.py -c ${ncores} -l "Rscript mc3.sling.tradeSeq.pseudotimeDE.R cc1B1_cc1B2_sce.rds CC1-B1 2 ${ncores} ${subsamples} 6 ${glm}"
#job2pbs.py -c ${ncores} -l "Rscript mc3.sling.tradeSeq.pseudotimeDE.R precc2_cc2_sce.rds pre-CC2 2 ${ncores} ${subsamples} 6 ${glm}"

#job2pbs.py -c ${ncores} -l "Rscript mc3.sling.tradeSeq.pseudotimeDE_on_seurat.R linA_seurat.rds OPC 20 ${ncores} ${subsamples} 6 ${glm}"
#job2pbs.py -c ${ncores} -l "Rscript mc3.sling.tradeSeq.pseudotimeDE_on_seurat.R linB_seurat.rds OPC 20 ${ncores} ${subsamples} 6 ${glm}"

#job2pbs.py -c ${ncores} -l "Rscript slingshot_iter100_tradeSeq_pseudotimeDE.R opc_precc1_sce.rds OPC ${npcs} ${ncores} ${subsamples} 6 ${glm}"
#job2pbs.py -c ${ncores} -l "Rscript slingshot_iter100_tradeSeq_pseudotimeDE.R linA_sce.rds OPC ${npcs} ${ncores} ${subsamples} 6 ${glm}"
#job2pbs.py -c ${ncores} -l "Rscript slingshot_iter100_tradeSeq_pseudotimeDE.R linB_sce.rds OPC ${npcs} ${ncores} ${subsamples} 6 ${glm}"
#job2pbs.py -c ${ncores} -l "Rscript slingshot_iter100_tradeSeq_pseudotimeDE.R precc1_cc1A_sce.rds pre-CC1 ${npcs} ${ncores} ${subsamples} 6 ${glm}"
#job2pbs.py -c ${ncores} -l "Rscript slingshot_iter100_tradeSeq_pseudotimeDE.R precc1_cc1B1_sce.rds pre-CC1 ${npcs} ${ncores} ${subsamples} 6 ${glm}"
#job2pbs.py -c ${ncores} -l "Rscript slingshot_iter100_tradeSeq_pseudotimeDE.R cc1B1_cc1B2_sce.rds CC1-B1 ${npcs} ${ncores} ${subsamples} 6 ${glm}"

#job2pbs.py -c 10 -q short -l "Rscript slingshot_iter100_tradeSeq_pseudotimeDE.R opc_precc1_seurat.rds OPC ${npcs} ${ncores} ${subsamples} 6 ${glm}"
#job2pbs.py -c 10 -q short -l "Rscript slingshot_iter100_tradeSeq_pseudotimeDE.R precc1_cc1A_seurat.rds pre-CC1 ${npcs} ${ncores} ${subsamples} 6 ${glm}"
#job2pbs.py -c ${ncores} -l "Rscript slingshot_iter100_tradeSeq_pseudotimeDE.R precc1_cc1B_seurat.rds pre-CC1 ${npcs} ${ncores} ${subsamples} 6 ${glm}"

#job2pbs.py -c ${ncores} -l "Rscript slingshot_iter100_tradeSeq_pseudotimeDE_robustPT_only.R opc_precc1_sce.rds OPC ${npcs} ${ncores} ${subsamples}"
#job2pbs.py -c ${ncores} -l "Rscript slingshot_iter100_tradeSeq_pseudotimeDE_robustPT_only.R linA_sce.rds OPC ${npcs} ${ncores} ${subsamples}"
#job2pbs.py -c ${ncores} -l "Rscript slingshot_iter100_tradeSeq_pseudotimeDE_robustPT_only.R linB_sce.rds OPC ${npcs} ${ncores} ${subsamples}"
#job2pbs.py -c ${ncores} -l "Rscript slingshot_iter100_tradeSeq_pseudotimeDE_robustPT_only.R precc1_cc1A_sce.rds pre-CC1 ${npcs} ${ncores} ${subsamples}"
#job2pbs.py -c ${ncores} -l "Rscript slingshot_iter100_tradeSeq_pseudotimeDE_robustPT_only.R precc1_cc1B1_sce.rds pre-CC1 ${npcs} ${ncores} ${subsamples}"
#job2pbs.py -c ${ncores} -l "Rscript slingshot_iter100_tradeSeq_pseudotimeDE_robustPT_only.R cc1B1_cc1B2_sce.rds CC1-B1 ${npcs} ${ncores} ${subsamples}"
