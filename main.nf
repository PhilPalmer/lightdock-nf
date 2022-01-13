#!/usr/bin/env nextflow

def helpMessage() {
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --pdb_paths input_pdbs.csv

    Inputs Options:
    --pdb_paths     CSV input file containing the list of pairwise receptor ligand PDB filepaths to be processed.
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
    set val(sample_name), file(receptor_pdb), file(ligand_pdb), file("{setup,lightdock,init}*"), file("swarm_*") into swarms

    script:
    """
    lightdock3_setup.py $receptor_pdb $ligand_pdb -s $params.swarms -g $params.glowworms --noxt --noh --now
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

process cluster_conformations {
    tag "$sample_name $swarm"
    publishDir "${params.outdir}/${sample_name}", mode: 'copy'

    input:
    set val(sample_name), file(receptor_pdb), file(ligand_pdb), file(lightdock), file(swarm) from samples_swarms
    
    output:
    set val(sample_name), file(swarm) into clustered_conformations

    script:
    """
    cd $swarm
    lgd_generate_conformations.py ../${receptor_pdb} ../${ligand_pdb} gso_${params.steps}.out $params.glowworms > /dev/null 2> /dev/null >> generate_lightdock.list
    lgd_cluster_bsas.py gso_${params.steps}.out > /dev/null 2> /dev/null >> cluster_lightdock.list
    """
}

clustered_conformations
    .groupTuple()
    .set { all_clustered_conformations }

process rank_conformations {
    tag "$sample_name"
    publishDir "${params.outdir}/${sample_name}", mode: 'copy'

    input:
    set val(sample_name), file(swarms) from all_clustered_conformations
    
    output:
    set val(sample_name), file("*") into ranked_conformations
    file("top_ranked_by_scores.csv") into scores

    script:
    // TODO: Add filtering
    // lgd_filter_restraints.py --cutoff 5.0 --fnat 0.4 rank_by_scoring.list restraints.list A B > /dev/null 2> /dev/null
    """
    cat swarm*/generate_lightdock.list > generate_lightdock.list
    cat swarm*/cluster_lightdock.list > cluster_lightdock.list
    ant_thony.py -c $task.cpus generate_lightdock.list
    ant_thony.py -c $task.cpus cluster_lightdock.list
    lgd_rank.py ${params.swarms} ${params.steps}
    get_top.py --scores_path rank_by_scoring.list --pdb_code $sample_name --out_path top_ranked_by_scores.csv
    """
}

scores
  .collectFile(name: "${params.outdir}/top_ranked_by_scores.csv", keepHeader: true, skip: 1)