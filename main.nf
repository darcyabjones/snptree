#!/usr/bin/env nextflow


params.vcf = "$baseDir/data/filtered_ann.vcf.gz"
params.genes = "$baseDir/data/SN15v9_OM_ChromOnly.gff3"
params.bed = "$baseDir/data/single_97pc_regions.bedgraph"

vcfFile = file(params.vcf)
bedFile = file(params.bed)


process applyFilters {
    container "quay.io/biocontainers/bcftools:1.9--h4da6232_0"

    input:
    file vcf from vcfFile

    output:
    file "${vcf.baseName}_filtered.vcf" into filteredVcf

    """
    bcftools view --types snps ${vcf} \
    | bcftools view -f".,PASS" -O v - \
    | uniq \
    > "${vcf.baseName}_filtered.vcf"
    """
}

process intersectVcfBed {
    container "quay.io/biocontainers/bedtools:2.27.1--1"
    publishDir "intersectBedVcf"

    input:
    file vcf from filteredVcf
    file bed from bedFile

    output:
    file "${vcf.baseName}_intersect.vcf" into intersectedVcf

    """
    bedtools intersect -a ${vcf} -b ${bed} -u -header | uniq > ${vcf.baseName}_intersect.vcf
    """
}


process selectInformative {
    container "quay.io/biocontainers/snpsift:4.3.1t--1"
    publishDir "selectInformative"

    input:
    file vcf from intersectedVcf

    output:
    file "${vcf.baseName}_informative.vcf" into informativeVcf

    """
    SnpSift filter "(ANN[*].EFFECT == 'missense_variant') || (ANN[*].EFFECT == 'synonymous_variant')" ${vcf} > "${vcf.baseName}_informative.vcf"
    """
}


process snpToFasta {
    container "quay.io/biocontainers/vcfkit:0.1.6--py27h24bf2e0_2"
    publishDir "snpToFasta"

    input:
    file vcf from informativeVcf

    output:
    file "${vcf.baseName}.fasta" into informativeFasta

    """
    vk phylo fasta ${vcf} > ${vcf.baseName}.fasta
    """
}


process variantFasta {
    container "quay.io/biocontainers/biopython:1.70--np112py36_1"
    publishDir "snpToFasta"

    input:
    file vcf from informativeFasta

    output:
    file "${vcf.baseName}_variant.fasta" into varFasta

    """
    #!/usr/bin/env python3
    
    from Bio import AlignIO
    
    msa = AlignIO.read("${vcf}", format="fasta")
    
    to_kill = []
    for i, col in enumerate(zip(*msa)):
        col = list(filter(lambda x: x in ('A', 'T', 'G', 'C'), col))
        if all(col[0] == base for base in col[1:]):
            to_kill.append(i)
    
    # Later we add +1 to i so use -1 to get 0.
    to_kill.insert(0, -1)
    to_kill.append(msa.get_alignment_length())
    
    # Run through the columns and construct new msa excluding killed sites.
    sub_msas = None
    for i, j in zip(to_kill, to_kill[1:]):
        i = i + 1
        if i == j:
            continue
        
        if sub_msas is None:
            sub_msas = msa[:, i:j]
        else:
            sub_msas += msa[:, i:j]
    
    AlignIO.write(sub_msas, "${vcf.baseName}_variant.fasta", format="fasta")
    """
}


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

