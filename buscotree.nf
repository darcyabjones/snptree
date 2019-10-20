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

// Tidied OG codon alignments to use.
// Should be a directory of ".fasta" files.
// File basenames will be used as OG names.
// The fasta ids will be used as species.
params.alignments = false

// A nexus or RaxML formatted partition file to use alongside the tidied fastas.
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

// Compute the gene trees using raxml-ng instead of IQtree with rapid bootstrapping.
params.raxml = false

// Collapse gene tree clades with < this support value into polytomies before running astral.
// Uses 50% for IQtree UFboot values and 10% for RAxML non-parametric bootstraps.
params.gene_bs_collapse = params.raxml ? 10 : 50

// The genetic code to use to translate codons.
// Should be a number corresponding to one of the NCBI translation tables.
params.gencode = 1

params.seed = 123

def run_busco = !( params.buscos || params.ogs || params.alignments )
def run_ogs = !( params.ogs || params.alignments )
def run_align = !( params.alignments )
def run_partitions = !( params.partitions )

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
        .fromPath(params.alignments, checkIfExists: true, type: "file")
        .set { userAlignments }
}


if ( params.partitions ) {
    Channel
        .fromPath(params.partitions, checkIfExists: true, type: "file")
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
    file "og_alignments" into ogAlignments

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



process tidyCodonAlignments {

    label "python3"
    label "big_task"

    publishDir "${params.outdir}"

    input:
    file "alignments" from ogAlignments

    output:
    file "tidied_alignments/*" into computedTidiedAlignments mode flatten

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


if ( params.alignments ) {
    tidiedAlignments = userAlignments
} else {
    tidiedAlignments = computedTidiedAlignments
}

tidiedAlignments.into {
    tidiedAlignments4GetPartitions;
    tidiedAlignments4MergeWithPartitions;
    tidiedAlignments4MergeWithPartitionedModels;
}


process getPartitions {

    label "python3"
    label "big_task"

    publishDir "${params.outdir}"

    when:
    run_partitions

    input:
    file "alignments/*" from tidiedAlignments4GetPartitions.collect()

    output:
    file "partitions/*" into computedPartitions mode flatten

    script:
    """
    mkdir partitions

    find \
      alignments/ \
      \\( -name "*.fasta" -or -name "*.fa" -or -name "*.fna" -or -name "*.fas" \\) \
      -printf '%f\\0' \
    | xargs \
        -0 \
        -P "${task.cpus}" \
        -I {} \
        -- \
        get_codon_partitions.py -o "partitions/{}" "alignments/{}"

    for f in partitions/*;
    do
      mv "\${f}" "\${f%.*}.partitions"
    done
    """
}


if ( params.partitions ) {
    partitions = userPartitions
} else {
    partitions = computedPartitions
}


tidiedAlignments4MergeWithPartitions
    .map { [ it.baseName, it ] }
    .join( partitions.map { [ it.baseName, it ] }, by: 0 )
    .set { alignmentsWithPartitions }


if (params.raxml) {

    /*
     * ModelTest-NG
     * DOI:
     *
     * NOTE: the program fails with error code 114 if it can't parse the MSA for some reason.
     * NOTE: the program fails with error code 136 with a floating point exception.
     *       It seems that this is related to MSAs without divergent sites?
     */
    process findPartitionedModel {

        label "modeltest"
        label "small_task"
        validExitStatus 0,114,136
        publishDir "${params.outdir}/partitioned_models"
        tag "${name}"

        input:
        set val(name),
            file("msa.fasta"),
            file("msa.partitions") from alignmentsWithPartitions

        output:
        set val(name),
            file("${name}.part.aic") optional true into partitionedModels

        file "${name}.part.bic" optional true
        file "${name}.part.aicc" optional true

        script:
        """
        modeltest-ng \
          --datatype nt \
          --input "msa.fasta" \
          --partitions "msa.partitions" \
          --output "${name}" \
          -f ef \
          --model-het uigf \
          --template raxml
        """
    }

    /*
     */
    tidiedAlignments4MergeWithPartitionedModels
        .map { [ it.baseName, it ] }
        .join( partitionedModels, by: 0 )
        .set { alignmentsWithPartitionedModels }


    process computeRAxMLGeneTrees {

        label "raxml"
        label "medium_task"

        publishDir "${params.outdir}/gene_trees"

        tag "${name}"

        input:
        set val(name),
            file("msa.fasta"),
            file("msa.partitions") from alignmentsWithPartitionedModels

        output:
        file "${name}.raxml.*"

        //file("T15.raxml.bestModel")
        //T15.raxml.bestTree
        //T15.raxml.bootstraps
        //T15.raxml.log
        //T15.raxml.mlTrees
        //T15.raxml.rba
        //T15.raxml.startTree
        //T15.raxml.supportFBP
        //T15.raxml.supportTBE

        script:
        """
        sed '/^>/s/[[:space:]].*\$//' "msa.fasta" > tidied.fasta
        raxml-ng --parse --msa "tidied.fasta" --model msa.partitions --prefix T1
        raxml-ng \
          --all \
          --msa T1.raxml.rba \
          --model msa.partitions \
          --prefix "${name}" \
          --threads "${task.cpus}" \
          --seed "${params.seed}" \
          --bs-cutoff 0.3 \
          --bs-metric fbp,tbe
        """
    }

} else {

    process computeIQTreeGeneTrees {
        label "iqtree"
        label "small_task"

        publishDir "${params.outdir}/gene_trees"

        tag "${name}"

        input:
        set val(name),
            file("msa.fasta"),
            file("msa.partitions") from alignmentsWithPartitions

        output:
        set val(name),
            file("${name}_iqtree.treefile") into geneTrees

	file "${name}_iqtree.alninfo"
	file "${name}_iqtree.best_scheme"
	file "${name}_iqtree.best_scheme.nex"
	file "${name}_iqtree.contree"
	file "${name}_iqtree.iqtree"
	file "${name}_iqtree.mldist"
	file "${name}_iqtree.model.gz"
	file "${name}_iqtree.rate"
	file "${name}_iqtree.sitelh"
	file "${name}_iqtree.siteprob"
	file "${name}_iqtree.ufboot"

        script:
        """
        iqtree \
          -nt 1 \
          -s msa.fasta \
          -spp msa.partitions \
          -m MFP \
          -cmax 6 \
          -st DNA \
          -bb 1000 \
          -nstop 50 \
          -bnni \
          -wbt \
          -wsr \
          -wspr \
          -wslr \
          -alninfo \
          -pre "${name}_iqtree"
        """
    }
}


process collapseGeneTrees {

    label "nwutils"
    label "small_task"

    input:
    file "gene_trees/*.nwk" from geneTrees
        .map { n, f -> f }
        .collect()


    output:
    file "collapsed.nwk" into collapsedGeneTrees

    script:
    """
    cat gene_trees/*.nwk > combined.nwk
    nw_ed combined.nwk 'i & b <= ${params.gene_bs_collapse}' o > collapsed.nwk
    """
}


process runAstral {

    label "astral"
    label "small_task"

    input:
    file "gene_trees.nwk" from collapsedGeneTrees

    script:
    """
    java -jar "\${ASTRAL_JAR}" \
      -i gene_trees.nwk \
      -o species.nwk \
    2> species_log.txt
    """
}

/*
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
*/


/*
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
*/


/*
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
 */
