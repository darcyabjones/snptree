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

params.gene_trees = false

//
params.species_tree = false


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

params.max_iterations = 200
params.max_dist = 10

params.seed = 123

def run_busco = !( params.buscos || params.ogs || params.alignments )
def run_ogs = !( params.ogs || params.alignments )
def run_align = !( params.alignments )
def run_partitions = !( params.partitions )
def run_gene_trees = !( params.gene_trees )
def run_species_tree = !( params.species_tree )

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


if ( params.gene_trees ) {
    Channel
        .fromPath(params.gene_trees, checkIfExists: true, type: 'file')
        .map { f -> [f.baseName, f] }
        .set { userGeneTrees }
} else {
    userGeneTrees = Channel.empty()
}


if ( params.species_tree ) {
    Channel
        .fromPath(params.species_tree, checkIfExists: true, type: "file")
        .first()
        .set { userSpeciesTree }
} else {
    userSpeciesTree = Channel.empty()
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
    tidiedAlignments4GetConcordanceFactors;
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

        when:
        run_gene_trees

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

        when:
        run_gene_trees

        input:
        set val(name),
            file("msa.fasta"),
            file("msa.partitions") from alignmentsWithPartitionedModels

        output:
        set val(name),
            file("${name}.raxml.bestTree") into computedGeneTrees
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

        when:
        run_gene_trees

        input:
        set val(name),
            file("msa.fasta"),
            file("msa.partitions") from alignmentsWithPartitions

        output:
        set val(name),
            file("${name}_iqtree.treefile") into computedGeneTrees

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


if ( params.gene_trees ) {
    geneTrees = userGeneTrees
} else {
    geneTrees = computedGeneTrees
}

geneTrees.into {
    geneTrees4CollapseGeneTrees;
    geneTrees4ComputeGeneTreeRFDists;
    geneTrees4GetConcordanceFactors;
    geneTrees4MergeWithRFDists;
}


process collapseGeneTrees {

    label "nwutils"
    label "small_task"

    input:
    file "gene_trees/*.nwk" from geneTrees4CollapseGeneTrees
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

    publishDir "${params.outdir}/species_trees"

    when:
    run_species_tree

    input:
    file "gene_trees.nwk" from collapsedGeneTrees

    output:
    file "species.nwk" into computedSpeciesTree
    file "species_log.nwk"

    script:
    """
    java -jar "\${ASTRAL_JAR}" \
      -i gene_trees.nwk \
      -o species.nwk \
    2> species_log.txt
    """
}


if ( params.species_tree ) {
    speciesTree = userSpeciesTree
} else {
    speciesTree = computedSpeciesTree
}

speciesTree.into {
    speciesTree4ComputeGeneTreeRFDists;
    speciesTree4GetConcordanceFactors;
    speciesTree4SelectBestGeneTrees;
}


process getConcordanceFactors {

    label "iqtree"
    label "big_task"

    input:
    file "species.nwk" from speciesTree4GetConcordanceFactors
    file "loci/*" from geneTrees4GetConcordanceFactors.map {n, f -> f} .collect()
    file "alignments/*" from tidiedAlignments4GetConcordanceFactors.collect()

    script:
    """
    cat loci/* > loci.nwk
    iqtree \
      -nt "${task.cpus}" \
      -t species.nwk \
      --gcf loci.nwk \
      -p alignments \
      --scf 100 \
      --prefix concord
    """
}


/*
 */
process computeGeneTreeRFDists {

    label "iqtree"
    label "small_task"

    tag "${name}"

    input:
    set val(name),
        file("gene.nwk"),
        file("species.nwk") from geneTrees4ComputeGeneTreeRFDists
            .combine(speciesTree4ComputeGeneTreeRFDists) 

    output:
    set val(name),
        file("${name}_rfdist.tsv") into geneTreeRFDists

    script:
    """
    iqtree -t species.nwk -rf gene.nwk
    tail -n+2 species.nwk.rfdist \
    | awk -v name="${name}" 'BEGIN {OFS="\\t"} {print name, \$2}' \
    > "${name}_rfdist.tsv"
    """
}


/*
 * This next section is a bit weird.
 * What we're doing is progressively combining the best trees
 * (by rfdist) and creating a new "species" tree from them.
 * if the new species tree is within x distance of the complete species tree,
 * then we have a valid subset and we can stop adding more trees,
 * otherwise we add the next best tree.
 *
 * The approach is based on https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3985644/
 *
 * First sort by rfdistance (ascending)
 */
geneTreeRFDists
    .map { n, f -> f }
    .collectFile(name: "${params.outdir}/gene_trees.tsv")
    .splitCsv(sep: '\t')
    .toSortedList( { a, b -> a[1] <=> b[1] } )
    .flatMap()
    .combine(geneTrees4MergeWithRFDists, by: 0) 
    .map { n, d, f -> [n, f] }
    .set { combinedGeneTreeRFDists }


/*
 * Here select the best rftree as an initial tree.
 * The "choice" hack is because .head() returns a value channel and i also need
 * the remainder after that, which doesn't have an operator
 * and i don't seem to be able to use the channel again after head.
 */
def is_first = true

// We get a depreciation warning for this, but it is unclear what
// the appropriate replacement is. .empty() doesn't work.
def head = Channel.create()
def tail = Channel.create()

combinedGeneTreeRFDists
    .choice(head, tail) { a ->
        if (is_first) {
            is_first = false
            0
        } else {
            1
        }
    }


/*
 * Now we construct the initial accumulator.
 * Essentially this is a fold operation, and we want to keep track
 * of the gene ids that make up the trees and the gene trees themselves.
 */
process initFold {

    label "posix"
    label "small_task"

    input:
    set val(name), file("gene.nwk") from head 

    output:
    set file("names.tsv"), file("gene.nwk") into initAccumulator

    script:
    """
    echo "${name}" > names.tsv
    """
}


/*
 * This is necessary to complete the loop.
 * We're combining the initial accumulator with the
 * accumulator generated after each iteration.
 */
recurseAccumulator = Channel.create()
accumulator = initAccumulator.mix(recurseAccumulator)

/*
 * This takes our accumulator and the next best gene from the
 * gene tree channel, and adds the new gene to the accumulator.
 */
process combineFold {

    label "posix"
    label "small_task"

    input:
    set val(name),
        file("gene.nwk"),
        file("names.tsv"),
        file("accum.nwk") from tail
            .merge(accumulator)

    output:
    set file("new_names.tsv"),
        file("new_accum.nwk") into combinedFold

    script:
    """
    cp -L names.tsv new_names.tsv
    cp -L accum.nwk new_accum.nwk

    echo "${name}" >> new_names.tsv
    cat gene.nwk >> new_accum.nwk
    """
}


/*
 * With our set of gene trees in accumulator
 * construct a composite tree.
 */
process selectBestGeneTrees {

    label "astral"
    label "small_task"

    input:
    set file("names.tsv"), file("gene.nwk") from combinedFold

    output:
    set file("names.tsv"), file("gene.nwk"), file("species.nwk") into combinedFoldTree

    script:
    """
    java -jar "\${ASTRAL_JAR}" \
      -i gene.nwk \
      -o species.nwk \
    2> species_log.txt
    """
}


/*
 * Evaluate the composite tree from this subset of gene trees against
 * the species tree constructed from the full set.
 * If the distance is small enough, or we've reached the maximum number
 * of iterations, we stop the recursion by not feeding new files back into
 * the accumulator channel.
 */
process findGeneTreeDist {

    label "iqtree"
    label "small_task"

    input:
    set file("names.tsv"),
        file("gene.nwk"),
        file("gene_species.nwk"),
        file("species.nwk") from combinedFoldTree
            .combine(speciesTree4SelectBestGeneTrees)

    output:
    set file("new_names.tsv"), file("new_gene.nwk") into recurseAccumulator optional true
    set file("names.tsv"), file("gene.nwk"), file("gene_species.nwk"), file("distance.txt")

    script:
    """
    iqtree -t species.nwk -rf gene_species.nwk
    DIST=\$(tail -n+2 species.nwk.rfdist | awk '{print \$2}')
    NITERS=\$(wc -l < gene.nwk)

    if [ \${DIST} -gt ${params.max_dist} ] && [ \${NITERS} -le ${params.max_iterations} ]
    then
        cp -L gene.nwk new_gene.nwk
        cp -L names.tsv new_names.tsv
    fi 

    echo \${DIST} > distance.txt
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
