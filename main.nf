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

// The vcf to use for the analysis
params.vcf = false

// The reference isolate name.
// If this is not provided, the reference isolate will be called "REF"
params.ref_name = false

// A regular gene-format gff3 file containing CDS features to find SNPs in genes from.
// Will also form partitions if neither occultercut or nopartitions is used.
params.genes = false

// A bed-file with the ID in the 4th column to use as partitions instead of the genes or occultercut.
params.partitions = false

// A fasta formatted file of SNPs
params.fasta = false

// If fasta is and you want partitions this should be an IQtree formatted partitions file for the fasta.
params.fasta_partitions = false

// A fasta file of the reference genome to find partitions in, based on on nucleotide content, and use those as partitions.
params.occultercut = false

// Don't use any partitions.
params.nopartition = false

// Minimum non-major allele frequency.
params.min_af = 0.05

// The maximum amount of missing data (or gaps) that can be present per locus.
params.max_missing = 0.05

// Default cmax
params.cmax = params.nopartition ? 15 : 10

// Default rcluster.
// Reduce to ~10 if finding model takes a long time.
params.rcluster = 50

run_select_vcf = !params.fasta

use_genes_as_partitions = (
    !params.partitions &&
    !params.occultercut &&
    run_select_vcf
)

run_occultercut = (
    params.occultercut &&
    !params.nopartition &&
    !params.partitions &&
    run_select_vcf
)


if ( params.vcf ) {
    Channel
        .fromPath(params.vcf, checkIfExists: true, type: "file")
        .first()
        .set {vcf}
} else if ( run_select_vcf ) {
    log.error "Either a VCF of Fasta file of the SNPs must be provided."
    exit 1
} else {
    Channel.empty()
}


if ( params.genes ) {
    Channel
        .fromPath(params.genes, checkIfExists: true, type: "file")
        .first()
        .set {genes}
} else if ( run_select_vcf ) {
    log.error "A gene GFF3 file is required if selecting SNPs from a vcf."
    exit 1
} else {
    Channel.empty()
}


if ( params.partitions ) {
    Channel
        .fromPath(params.partitions, checkIfExists: true, type: "file")
        .first()
        .set {userPartitions}
} else {
    userPartitions = Channel.empty()
}


if ( params.fasta ) {
    Channel
        .fromPath(params.fasta, checkIfExists: true, type: "file")
        .first()
        .set {userFasta}
} else {
    userFasta = Channel.empty()
}


if ( params.fasta_partitions ) {
    Channel
        .fromPath(params.fasta_partitions, checkIfExists: true, type: "file")
        .first()
        .set {userFastaPartitions}
} else if ( params.fasta && !params.nopartition ) {
    log.error "You provided a fasta file but no fasta_partitions. " +
              "Please provide the partitions or use --nopartitions flag " +
              "to disable."
    exit 1
} else {
    userFastaPartitions = Channel.empty()
}


if ( params.occultercut ) {
    Channel
        .fromPath(params.occultercut, checkIfExists: true, type: "file")
        .first()
        .set {genome}
} else {
    genome = Channel.empty()
}



process applyFilters {
    label "htslib"
    label "small_task"

    when:
    run_select_vcf

    input:
    file input_vcf from vcf

    output:
    file "filtered.vcf" into filteredVcf

    script:
    """
    bcftools view \
        --include 'TYPE="snp" && F_PASS(GT ~ "\\.") < ${params.max_missing}' \
        --exclude-types "indels,mnps,bnd,other" \
        --apply-filters ".,PASS" \
        --min-af "${params.min_af}:nonmajor" \
        --output-type v \
        "${input_vcf}" \
    | uniq \
    > "filtered.vcf"
    """
}


process selectInformative {
    label "snpeff"
    label "small_task"

    publishDir "${params.outdir}/select_vcf"

    when:
    run_select_vcf

    input:
    file vcf from filteredVcf

    output:
    file "informative.vcf" into informativeVcf

    script:
    """
    SnpSift \
      filter \
      -s <(echo -e 'LOW\\nMODIFIER\\n') \
      "(ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].IMPACT in SET[0])" \
      "${vcf}" \
    > "informative.vcf"
    """
}


process vcfToTsv {
    label "htslib"
    label "small_task"

    publishDir "${params.outdir}/select_vcf"

    when:
    run_select_vcf

    input:
    file vcf from informativeVcf

    output:
    file "informative.tsv" into informativeTsv

    script:
    convert_ref = params.ref_name ? "| sed '1s/REF/${params.ref_name}/'" : ""

    """
    bcftools query \
      --print-header \
      --format "%CHROM\\t%POS0\\t%END\\t%REF[\\t%TGT]\\n" \
      "${vcf}" \
    | sed '1s/\\[[[:digit:]]*\\]//g' \
    | sed '1s/:GT//g' \
    ${convert_ref} \
    > "informative.tsv"
    """
}


process genesGFFToBed {

    label "bedtools"
    label "small_task"

    publishDir "${params.outdir}/select_vcf"

    when:
    use_genes_as_partitions

    input:
    file "genes.gff3" from genes

    output:
    file "genes.bed" into geneBED

    script:
    """
    gawk '
      BEGIN {OFS="\\t"}
      \$3 == "CDS" {
        n=gensub(/.*Parent\\s*=\\s*([^;\\n]+).*/, "\\\\1", "g", \$9);
        print \$1, \$4 - 1, \$5, n
      }' \
      genes.gff3 \
    | bedtools groupby -g 4 -c 1,2,3 -o first,min,max \
    | awk 'BEGIN {OFS="\\t"} {print \$2, \$3, \$4, \$1}' \
    | sort -k1,1 -k2,2n -k3,3n \
    | bedtools merge -c 4 -o first \
    > genes.bed
    """
}


process runOcculterCut {

    label "occultercut"
    label "small_task"

    publishDir "${params.outdir}/occultercut"

    when:
    run_occultercut

    input:
    file fasta from genome

    output:
    file "occultercut_regions.bed" into occulterCutBED
    file "occultercut_regions.gff3"
    file "occultercut.png"
    file "occultercut_composition_gc.txt"
    file "occultercut_my_genome.txt"
    file "occultercut_grouped_regions.gff3" optional true
    file "occultercut_nuc_frequencies.R*" optional true

    script:
    """
    OcculterCut -f "${fasta}"


    sed -i '1i set terminal pngcairo size 900,600 enhanced font "Helvetica,20"' plot.plt
    sed -i '1a set output "plot.png"' plot.plt
    gnuplot plot.plt

    mv plot.png "occultercut.png"

    mv compositionGC.txt "occultercut_composition_gc.txt"
    mv regions.gff3 "occultercut_regions.gff3"
    mv myGenome.txt "occultercut_my_genome.txt"

    if [ -e groupedRegions.gff3 ]
    then
      mv groupedRegions.gff3 "occultercut_grouped_regions.gff3"
    fi

    for f in nuc_frequencies.R*
    do
      mv "\${f}" "occultercut_\${f}"
    done

    gawk '
      BEGIN {OFS="\\t"}
      {
        n=gensub(/.*ID\\s*=\\s*([^;\\n]+).*/, "\\\\1", "g", \$9);
        print \$1, \$4 - 1, \$5, n
      }' \
      occultercut_regions.gff3 \
    > occultercut_regions.bed
    """
}


if ( use_genes_as_partitions ) {
    partitionBed = geneBED
} else if ( run_occultercut ) {
    partitionBed = occulterCutBED
} else {
    partitionBed = userPartitions
}


process getIntersectingGenes {
    label "bedtools"
    label "small_task"

    publishDir "${params.outdir}/select_vcf"

    when:
    run_select_vcf

    input:
    set file("informative.tsv"),
        file("partitions.bed") from informativeTsv
            .combine(partitionBed)

    output:
    file "informative_with_partitions.tsv" into informativePartitionsTSV

    script:
    """
    bedtools intersect \
      -a informative.tsv \
      -b partitions.bed \
      -wa \
      -wb \
      -header \
      -sortout \
    > informative_with_partitions.tsv
    """
}


process getFasta {

    label "python3"
    label "small_task"

    publishDir "${params.outdir}/select_vcf"

    when:
    run_select_vcf

    input:
    file "informative_with_partitions.tsv" from informativePartitionsTSV

    output:
    set file("snps.fasta"), file("partitions.txt") into fastaForTree

    script:
    """
    tsv2fasta.py \
      -o snps.fasta \
      -p partitions.txt \
      informative_with_partitions.tsv
    """
}


fastaForTree.into {
    fastaForTree4FindBestPartitions;
    fastaForTree4FindBestPartitionsModel;
    fastaForTree4RunPartitionTreeBootstraps;
    fastaForTree4FindBestModelForTree;
    fastaForTree4RunTreeBootstraps;
}


process findBestPartitions {

    label "iqtree"
    label "big_task"

    publishDir "${params.outdir}/tree"

    when:
    !params.nopartition

    input:
    set file("snps.fasta"),
        file("partitions.txt") from fastaForTree4FindBestPartitions

    output:
    file "best_partitions.nex" into bestPartitions

    script:
    def cmax = 5
    def rcluster = params.rcluster

    """
    iqtree \
      -nt AUTO \
      -ntmax "${task.cpus}" \
      -s snps.fasta \
      -st DNA \
      -m "TESTNEWONLY+ASC+MERGE" \
      -mset "GTR" \
      -cmax 5 \
      -rclusterf "${rcluster}" \
      -rcluster-max 5000 \
      -safe \
      -spp partitions.txt

    mv partitions.txt.best_scheme.nex best_partitions.nex
    # -mem "${task.memory.toGiga()}G"
    """
}


process findBestPartitionsModel {

    label "iqtree"
    label "big_task"

    publishDir "${params.outdir}/tree"

    when:
    !params.nopartition

    input:
    set file("snps.fasta"),
        file("partitions.txt"),
        file("best_partitions.nex") from fastaForTree4FindBestPartitionsModel
            .combine(bestPartitions)

    output:
    file "best_model.nex" into partitionModel

    script:
    def cmax = params.cmax
    def rcluster = params.rcluster

    """
    iqtree \
      -nt AUTO \
      -ntmax "${task.cpus}" \
      -s snps.fasta \
      -st DNA \
      -m "MF+ASC" \
      -mset "GTR,JC,F81,K80,HKY,K81" \
      -cmax "${cmax}" \
      -rcluster "${rcluster}" \
      -safe \
      -spp partitions.txt

    mv partitions.txt.best_scheme.nex best_model.nex
    # -mem "${task.memory.toGiga()}G"
    """
}

process runPartitionTreeUFBootstraps {

    label "iqtree"
    label "big_task"

    publishDir "${params.outdir}/tree"

    when:
    !params.nopartition

    input:
    set file("snps.fasta"),
        file("partitions.txt"),
        file("partitions.nex") from fastaForTree4RunPartitionTreeBootstraps
            .combine(partitionModel)

    script:
    """
    iqtree \
      -nt AUTO \
      -ntmax "${task.cpus}" \
      -s snps.fasta \
      -spp partitions.nex \
      -bb 10000 \
      -nstop 50 \
      -bnni \
      -alrt 1000 \
      -bspec GENESITE \
      -wbt \
      -wsr \
      -wspr \
      -wslr \
      -alninfo \
      -st DNA \
      -pre iqtree_partition
    """
}


process findBestModelForTree {

    label "iqtree"
    label "big_task"

    publishDir "${params.outdir}/tree/model_finder"

    when:
    params.nopartition

    input:
    set file("snps.fasta"),
        file("partitions.txt") from fastaForTree4FindBestModelForTree

    output:
    file "selected_model.txt" into model
    file "snps.fasta.model.gz"
    file "snps.fasta.iqtree"
    file "snps.fasta.log"
    file "snps.fasta.treefile"

    script:
    def cmax = params.cmax
    def rcluster = params.rcluster

    """
    iqtree \
      -nt AUTO \
      -ntmax "${task.cpus}" \
      -s snps.fasta \
      -st DNA \
      -m "MF+ASC" \
      -mset "GTR,JC,F81,K80,HKY,K81" \
      -cmax "${cmax}" \
      -rcluster "${rcluster}" \
      -safe

    awk '
      /^Best-fit model/ {
        n=gensub(/.*:\\s*(\\S+).*/, "\\\\1", "g", \$0);
        print n
    }' < snps.fasta.log \
    > selected_model.txt
    """
}


/*
*/
process runTreeUFBootstraps {

    label "iqtree"
    label "big_task"

    publishDir "${params.outdir}/snptree"

    when:
    params.nopartition

    input:
    set file("snps.fasta"),
        file("partitions.txt"),
        file("selected_model.txt") from fastaForTree4RunTreeBootstraps
            .combine(model)

    output:
    file "iqtree_nopartition.*"

    script:
    """
    iqtree \
      -nt AUTO \
      -ntmax "${task.cpus}" \
      -s snps.fasta \
      -m "\$(cat selected_model.txt)" \
      -bb 10000 \
      -bnni \
      -alrt 1000 \
      -wbt \
      -wsr \
      -wspr \
      -wslr \
      -alninfo \
      -st DNA \
      -pre iqtree_nopartition
    """
}


/*
process tree {
    container "quay.io/biocontainers/raxml:8.2.10--h470a237_1"
    publishDir "tree"

    input:
    file fasta from varFasta

    output:
    file "RAxML_bootstrap.${fasta.baseName}" into bootstrapFile
    """
    raxmlHPC-PTHREADS \
        -T 15 \
        -p 12345 \
        -b 12345 \
        -f d \
        -# autoMRE \
        -m ASC_GTRCAT \
        --asc-corr=lewis \
        -s ${fasta} \
        -n ${fasta.baseName}
    """
}


process consensus {
    container "quay.io/biocontainers/raxml:8.2.10--h470a237_1"
    publishDir "tree"

    input:
    file trees from bootstrapFile

    output:
    file "RAxML_MajorityRuleExtendedConsensusTree.consensus" into consensusTree

"""
raxmlHPC-PTHREADS \
    -T 15 \
    -J MRE \
    -z ${tree} \
    -m ASC_GTRCAT \
    --asc-corr=lewis \
    -n "consensus"
"""

}


process ml {
    container "quay.io/biocontainers/raxml:8.2.10--h470a237_1"
    publishDir "tree"

    input:

    output:
    file "RAxML_bestTree.likelihoods" into bestTreeFile
    file "RAxML_result.likelihoods.RUN.*" into resultFile
    file "RAxML_parsimonyTree.likelihoods.RUN.*" into parsimonyFile

    # most likely tree
    """
    raxmlHPC-PTHREADS \
        -T 16 \
        -f d \
        -m ASC_GTRCAT \
        -asc-corr=lewis \
        -n likelihoods \
        -s ${fasta} \
        -p 12345 \
        -# 50
    """
}

process annot {
    container "quay.io/biocontainers/raxml:8.2.10--h470a237_1"
    publishDir "tree"

    input:

    output:
    file "RAxML_bipartitions.best_annotated" into
    file "RAxML_bipartitionsBranchLabels.best_annotated"

    """
    raxmlHPC-PTHREADS \
        -T 15 \
        -f b \
        -m ASC_GTRCAT \
        --asc-corr=lewis \
        -z bootstrapped_trees \
        -t best_tree \
        -n bs_tree
    """
}
*/
