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
  .map { receptor_pdb, ligand_pdb -> [receptor_pdb.baseName, ligand_pdb.baseName, receptor_pdb, ligand_pdb] }
  .set { pdbs }

process lightdock {
    tag "$receptor_name $ligand_name"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    set val(receptor_name), val(ligand_name), file(receptor_pdb), file(ligand_pdb) from pdbs
    
    output:
    set val(receptor_name), val(ligand_name), file(receptor_pdb), file(ligand_pdb), file("swarm_*") into swarms

    script:
    """
    lightdock3_setup.py $receptor_pdb $ligand_pdb -s $params.swarms -g $params.glowworms
    lightdock3.py setup.json $params.steps
    """
}