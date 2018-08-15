#!/usr/bin/env nextflow


params.vcf = "$baseDir/data/filtered_ann.vcf.gz"
params.genes = "$baseDir/data/SN15v9_OM_ChromOnly.gff3"
params.bed = "$baseDir/data/single_97pc_regions.bedgraph"

vcfFile = file(params.vcf)
bedFile = file(params.vcf)


process applyFilters {
    input:
    file vcf from vcfFile

    output:
    file "${vcf.baseName}_filtered.vcf.gz" into filteredVcf

    """
    bcftools view -f".,PASS" -O z ${vcf} > "${vcf.baseName}_filtered.vcf.gz"
    """
}

process intersectVcfBed {
    //container "quay.io/biocontainers/bedtools:2.27.1--1"
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
    publishDir "selectInformative"

    input:
    file vcf from intersectedVcf

    output:
    file "${vcf.baseName}_informative.vcf" into informativeVcf

    """
    SnpSift filter "(ANN[*].EFFECT == 'missense_variant') || (ANN[*].EFFECT == 'synonymous_variant')" ${vcf} > "${vcf.baseName}_informative.vcf"
    """
}

process coordsToBED {
    publishDir "bed_alignments"

    input:
    file coord from mcoords

    output:
    file "${coord.baseName}.bed" into coordsBEDs

    """
    tail -n +5 ${coord} \
    | awk '{ printf "%s\\t%s\\t%s\\t\\.\\t+\\t%s\\tfiltered\\n", \$12, \$1-1, \$2, \$7 }' \
    | sort -k1,1 -k2,2n \
    > ${coord.baseName}.bed
    """
}


process genomeIndex {
    container "quay.io/biocontainers/samtools:1.9--h46bd0b3_0"
    input:
    file reference from reference_file

    output:
    file "${reference}.fai" into referenceIndex

    """
    samtools faidx ${reference} > ${reference}.fai
    """
}


process bedCoverage {
    container "quay.io/biocontainers/bedtools:2.27.1--1"
    publishDir "bedgraph_coverages"

    input:
    file bed from coordsBEDs
    file index from referenceIndex

    output:
    file "${bed.baseName}.bedgraph" into coverageBEDs

    """
    bedtools genomecov -bga -g ${index} -i ${bed} > ${bed.baseName}.bedgraph
    """
}



percentages = [1, 2, 3, 4, 5, 7, 10, 20, 30, 40, 50, 60, 70, 80, 90, 93, 95, 96, 97, 98, 99]

process compareCoverage {
    publishDir "compare_coverages"

    input:
    file bg from combinedCoverage
    val perc from percentages

    output:
    file "core_${perc}pc_regions.bedgraph" into core_regions
    file "duplicate_${perc}pc_regions.bedgraph" into duplicate_regions

    """
    /usr/bin/env python3

    import os
    import pandas as pd

    table = pd.read_table("${bg}")
    table.columns = [os.path.split(os.path.splitext(c)[0])[1] for c in table.columns]

    # Cheating a little bit here
    table.drop(columns=["15FG031_contigs", "15FG039_contigs", "15FG046_contigs", "203FG217_contigs"], inplace=True) 
    pc = ${perc} / 100

    table2 = table[(table.iloc[:, 3:] > 0).mean(axis=1) > pc]
    table2.to_csv("core_${perc}pc_regions.bedgraph", index=False, sep="\t")

    table3 = table2[(table2.iloc[:, 3:] > 1).mean(axis=1) > pc]
    table3.to_csv("duplicate_${perc}pc_regions.bedgraph", index=False, sep="\t")

    table4 = table2[(table2.iloc[:, 3:] == 1).mean(axis=1) > pc]
    table4.to_csv("single_${perc}pc_regions.bedgraph", index=False, sep="\t")
    """
}

