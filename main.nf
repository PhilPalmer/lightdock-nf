#!/usr/bin/env nextflow

def helpMessage() {
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --input_pdbs input_pdbs.csv

    Inputs Options:
    --input_pdbs    CSV input file containing the list of pairwise receptor ligand PDB filepaths to be processed.
    --swarms        Number of swarms/simulations to run
    --glowworms     Number of glowworms/agents to run per swarm
    --steps         Number of steps of the simulation
    """.stripIndent()
}

// Show help message
if (params.help) {
  helpMessage()
  exit 0
}

Channel
  .fromPath(params.pdb_paths)
  .ifEmpty { exit 1, "Cannot find CSV file : ${params.pdb_paths}" }
  .splitCsv(skip:1)
  .map { receptor_pdb, ligand_pdb -> [file(receptor_pdb), file(ligand_pdb)] }
  .map { receptor_pdb, ligand_pdb -> ["${receptor_pdb.baseName}_${ligand_pdb.baseName}", receptor_pdb, ligand_pdb] }
  .set { pdbs }

process lightdock {
    tag "$sample_name"
    // publishDir "${params.outdir}/${sample_name}", mode: 'copy'

    input:
    set val(sample_name), file(receptor_pdb), file(ligand_pdb) from pdbs
    
    output:
    set val(sample_name), file(receptor_pdb), file(ligand_pdb), file("lightdock*"), file("swarm_*") into swarms

    script:
    """
    lightdock3_setup.py $receptor_pdb $ligand_pdb -s $params.swarms -g $params.glowworms
    lightdock3.py setup.json $params.steps
    """
}

swarms
    .map { sample_name, receptor_pdb, ligand_pdb, lightdock, swarms -> 
        def samples_swarms = []
        swarms.each { swarm ->
            samples_swarms.add([sample_name, receptor_pdb, ligand_pdb, lightdock, swarm])
        }
        samples_swarms
    }
    .flatMap()
    .set { samples_swarms}


process ldg {
    tag "$sample_name $swarm"
    publishDir "${params.outdir}/${sample_name}", mode: 'copy'

    input:
    set val(sample_name), file(receptor_pdb), file(ligand_pdb), file(lightdock), file(swarm) from samples_swarms
    
    output:
    set val(sample_name), file("${swarm}/*") into lgd_out

    script:
    """
    cd $swarm
    lgd_generate_conformations.py ../${receptor_pdb} ../${ligand_pdb} gso_${params.steps}.out $params.glowworms
    lgd_cluster_bsas.py gso_${params.steps}.out
    lgd_rank.py ${params.swarms} ${params.steps}
    """
}