#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/


def helpMessage() {
    log.info"""
    # snptree

    ## Usage

    ```bash
    ```

    ## Parameters

    ## Output

    ## Exit codes

    - 0: All ok.
    - 1: Incomplete parameter inputs.

    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

// The genomes to run buscos on.
params.genomes = false

// Precomputed buscos to run instead of predicting them.
// Should be a directory containing many directoried of busco results for each genome.
// Folder names containing busco results will be used as species names.
params.buscos = false

// Orthologous groups to perform analyses on.
// Should be a directory of ".fasta" files.
// The basename will be used as the OG name.
// The fasta ids will be used as species.
params.ogs = false

// OG codon alignments to use.
// Should be a directory of ".fasta" files.
// File basenames will be used as OG names.
// The fasta ids will be used as species.
params.alignments = false

// Tidied OG codon alignment to use.
// Should be a single fasta file containing the concatenated codon alignments.
// No additional filtering will be done on this file.
// This should be accompanied by a partitioned file or a model file.
params.joined = false

// A nexus or RaxML formatted partition file to use alongside the tidied fasta.
params.partitions = false

// A nexus file containing the partitions and models to use for iqtree.
// Should be used alongside tidied.
params.model = false


// A folder containing the BUSCO lineage files.
params.busco_lineage = false

// The augustus species to use instead of the default for the busco lineage.
params.augustus_species = false

// The augustus config directory to use instead of the one distributed with augustus.
params.augustus_config = false

// The maximum proportion of gaps in an alignment column that is allowed.
// 1 means only filter out columns that are all gaps.
// 0 means filter all columns any gaps.
params.max_missing = 0.2
params.min_entropy = 1

// The genetic code to use to translate codons.
// Should be a number corresponding to one of the NCBI translation tables.
params.gencode = 1


def run_busco = !( params.buscos || params.ogs || params.alignments || params.joined )
def run_ogs = !( params.ogs || params.alignments || params.joined )
def run_align = !( params.alignments ||  params.joined )
def run_join = ! params.joined
def run_model_finder = ! params.model


if ( params.genomes ) {
    Channel.fromPath(params.genomes, checkIfExists: true, type: "file")
        .map { g -> [g.baseName, g] }
        .set { genomes }
} else if ( run_busco ) {
    log.error "Please provide some genomes to predict from."
    exit 1
} else {
    // We don't need them
    genomes = Channel.empty()
}


if ( params.busco_lineage ) {
    Channel
        .fromPath(params.busco_lineage, checkIfExists: true, type: "dir")
        .first()
        .set { buscoLineage }
} else if ( run_busco ) {
    log.error "Please provide a busco lineage folder."
    exit 1
} else {
    buscoLineage = Channel.empty()
}


if ( params.buscos ) {
    Channel
        .fromPath(params.buscos, checkIfExists: true, type: "dir")
        .first()
        .set { userBuscos }
}


if ( params.ogs ) {
    Channel
        .fromPath(params.ogs, checkIfExists: true, type: "dir")
        .first()
        .set { userOGs }
}


if ( params.alignments ) {
    Channel
        .fromPath(params.alignments, checkIfExists: true, type: "dir")
        .first()
        .set { userAlignments }
}


if ( params.joined ) {
    Channel
        .fromPath(params.joined, checkIfExists: true, type: "file")
        .first()
        .set { userJoined }
} else if ( !run_align && !(params.partitions || params.model) ) {
    log.error "You provided a tidied/joined fasta but not a partition or model file."
    exit 1
}


if ( params.partitions ) {
    Channel
        .fromPath(params.partitions, checkIfExists: true, type: "file")
        .first()
        .set { userPartitions }
} else {
    userPartitions = Channel.empty()
}


if ( params.model ) {
    Channel
        .fromPath(params.model, checkIfExists: true, type: "file")
        .first()
        .set { userModel }
} else {
    userModel = Channel.empty()
}


// Because augustus requires the config folder to be editable, it's easier
// to keep track of it in the pipeline.
// Assumes the correct environment variable is set.
if ( params.augustus_config ) {
    Channel
        .fromPath(params.augustus_config, checkIfExists: true, type: "dir")
        .first()
        .set { augustusConfig }

} else {
    process getAugustusConfig {

        label "busco"
        label "small_task"

        when:
        run_busco

        output:
        file "config" into augustusConfig

        script:
        """
        cp -r \${AUGUSTUS_CONFIG_PATH} ./config
        """
    }
}

genomes.set { genomes4Busco }


/*
 * TODO: Busco seems to return 0 even on error, need to check stderr for error logs.
 */
process runBusco {
    label "busco"
    label "medium_task"

    tag "${name}"

    publishDir "${params.outdir}/buscos"

    when:
    run_busco

    input:
    set val(name), file(fasta) from genomes4Busco
    file "lineage" from buscoLineage
    file "augustus_config" from augustusConfig

    output:
    file "${name}" into buscoResultsComputed

    script:
    species = params.augustus_species ? "--species ${params.augustus_species} " : ""

    """
    export AUGUSTUS_CONFIG_PATH="\${PWD}/augustus_config"

    run_BUSCO.py \
      --in "${fasta}" \
      --out "${name}" \
      --cpu ${task.cpus} \
      --mode "genome" \
      ${species} \
      --lineage_path "lineage"

    mv "run_${name}" "${name}"
    """
}


if ( params.buscos ) {
    buscoResults = userBuscos
} else {
    process processBuscoDirs {

        label "posix"
        label "small_task"

        input:
        file "*" from buscoResultsComputed.collect()

        output:
        file "buscos" into buscoResults

        script:
        """
        mkdir buscos

        find \
          -L \
          . \
          -maxdepth 1 \
          -type d \
          -path ./buscos -prune \
          -o \
          -not -path . \
          -printf '%f\\0' \
        | xargs -0 -I {} -- ln -sf "\${PWD}/{}" "\${PWD}/buscos/{}"
        """
    }
}


process selectBuscos {

    label "python3"
    label "small_task"

    publishDir "${params.outdir}"

    when:
    run_ogs

    input:
    file "buscos" from buscoResults.collect()

    output:
    file "ogs" into computedOGs

    script:
    """
    get_good_buscos.py \
      --outdir ogs \
      --max-missing "${params.max_missing}" \
      buscos
    """
}


if ( params.ogs ) {
    selectedOGs = userOGs
} else {
    selectedOGs = computedOGs
}


process alignOGs {

    label "decipher"
    label "big_task"

    publishDir "${params.outdir}"

    when:
    run_align

    input:
    file "families" from selectedOGs

    output:
    file "og_alignments" into computedOGAlignments

    script:
    """
    mkdir og_alignments

    find \
      families/ \
      \\( -name "*.fasta" -or -name "*.fa" -or -name "*.fna" -or -name "*.fas" \\) \
      -printf '%f\\0' \
    | xargs \
        -0 \
        -P "${task.cpus}" \
        -I {} \
        -- \
        run_decipher.R \
          --infile "families/{}" \
          --outfile "og_alignments/{}" \
          --gencode "${params.gencode}" \
          --maxgap "${params.max_missing}" \
          --minentropy "${params.min_entropy}"
    """
}


if ( params.alignments ) {
    ogAlignments = userAlignments
} else {
    ogAlignments = computedOGAlignments
}


process tidyCodonAlignments {

    label "python3"
    label "big_task"

    publishDir "${params.outdir}"

    input:
    file "alignments" from ogAlignments

    output:
    file "tidied_alignments" into tidiedAlignments

    script:
    """
    mkdir tidied_alignments
    find \
      alignments/ \
      \\( -name "*.fasta" -or -name "*.fa" -or -name "*.fna" -or -name "*.fas" \\) \
      -printf '%f\\0' \
    | xargs \
        -0 \
        -P "${task.cpus}" \
        -I {} \
        -- \
        tidy_alignments.py \
          --max-missing "${params.max_missing}" \
          --gencode "${params.gencode}" \
          -o "tidied_alignments/{}" \
          "alignments/{}"
    """
}


process joinAlignments {

    label "python3"
    label "small_task"

    publishDir "${params.outdir}/joined_alignments"

    when:
    run_join

    input:
    file "alignments" from tidiedAlignments

    output:
    set file("joined_alignment.fasta"),
        file("partitions.txt") into computedJoinedAlignments

    script:
    """
    join_alignments.py \
      --trim \
      --outfasta "joined_alignment.fasta" \
      --outpartition "partitions.txt" \
      --type DNA \
      alignments/*
    """
 }


if ( params.joined && run_model_finder ) {
    joinedAlignments = userJoined.combine(userPartitions)
} else {
    joinedAlignments = computedJoinedAlignments
}


joinedAlignments.into {
    joinedAlignments4FindBestModelForTree;
    joinedAlignments4RunPartitionTreeBootstraps;
}

process findBestModelForTree {

    label "iqtree"
    label "big_task"

    publishDir "${params.outdir}/selected_models"

    when:
    run_model_finder

    input:
    set file("alignment.fasta"),
        file("partitions.txt") from joinedAlignments4FindBestModelForTree

    output:
    file "partitions.txt.best_scheme.nex" into computedModel

    script:
    def cmax = 10
    def rcluster = 10
    def gencode = params.gencode

    """
    iqtree \
      -nt AUTO \
      -ntmax "${task.cpus}" \
      -st "CODON${gencode}" \
      -m "MF+MERGE" \
      -mset "GY2K,MG2K,ECMK07,ECMK07_GY2K" \
      -cmax "${cmax}" \
      -rcluster "${rcluster}" \
      -safe \
      -s alignment.fasta \
      -spp partitions.txt
    """
}


if ( params.model ) {
    model = userModel
} else {
    model = computedModel
}


/*
 */
process runPartitionTreeBootstraps {

    label "iqtree"
    label "big_task"

    publishDir "${params.outdir}/buscotree"

    input:
    set val(chunk),
        file("alignment.fasta"),
        file("partitions.txt"),
        file("partitions.nex") from joinedAlignments4RunPartitionTreeBootstraps
            .combine(model)

    script:
    """
    iqtree \
      -nt "${task.cpus}" \
      -s alignment.fasta \
      -spp partitions.nex \
      -bb 10000 \
      -alrt 1000 \
      -bspec GENESITE \
      -bnni \
      -wbt \
      -wsr \
      -wspr \
      -wslr \
      -alninfo \
      -st DNA \
      -pre "iqtree_busco_partition"
    """
}
