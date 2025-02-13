nextflow.enable.dsl=2

//sample name
sample = "GJM_rev"


process callVariants {
    container 'broadinstitute/gatk:4.2.3.0'
    
    publishDir "${params.outpath}/results/", mode: 'copy'

    input:
    
    file("${params.bam}")
    file("${params.bam}.bai")
    file("${params.ref}")

    output:
    file("${sample}.raw_variants.vcf")

    script:
    """
    gatk --java-options "-Xmx128g" HaplotypeCaller \
    -R ${params.ref} \
    -I ${params.bam} \
    --output-mode EMIT_ALL_ACTIVE_SITES \
    -O ${sample}.raw_variants.vcf
    """
}

process extractSNPs {
    container 'broadinstitute/gatk:4.2.3.0'

    publishDir "${params.outpath}/results/", mode: 'copy'

    input:
    file("${params.ref}")
    file("${sample}.raw_variants.vcf")

    output:
    file("${sample}.raw_snps.vcf")

    script:
    """
    gatk --java-options "-Xmx128g" SelectVariants \
    -R ${params.ref} \
    -V ${sample}.raw_variants.vcf \
    --select-type-to-include SNP \
    -O ${sample}.raw_snps.vcf
    """
}

process extractIndels {
    container 'broadinstitute/gatk:4.2.3.0'

    publishDir "${params.outpath}/results/", mode: 'copy'

    input:
    file("${params.ref}")
    file("${sample}.raw_variants.vcf")

    output:
    file("${sample}.raw_indels.vcf")

    script:
    """
    gatk --java-options "-Xmx128g" SelectVariants \
    -R ${params.ref} \
    -V ${sample}.raw_variants.vcf \
    --select-type-to-include INDEL \
    -O ${sample}.raw_indels.vcf
    """
}

process filterSNPs {
    container 'broadinstitute/gatk:4.2.3.0'

    publishDir "${params.outpath}/results/", mode: 'copy'

    input:
    file("${params.ref}")
    file("${sample}.raw_snps.vcf")

    output:
    file("${sample}.filtered_snps.vcf")

    script:
    """
    gatk --java-options "-Xmx128g" VariantFiltration \
    -R ${params.ref} \
    -V ${sample}.raw_snps.vcf \
    -O ${sample}.filtered_snps.vcf \
    -filter-name "QD_filter" -filter "QD < 2.0" \
    -filter-name "FS_filter" -filter "FS > 60.0" \
    -filter-name "MQ_filter" -filter "MQ < 40.0" \
    -filter-name "QUAL30" -filter "QUAL < 30.0" \
    -filter-name "SOR_filter" -filter "SOR > 4.0" \
    -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
    -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"
    """
}

process filterIndels {
    container 'broadinstitute/gatk:4.2.3.0'

    publishDir "${params.outpath}/results/", mode: 'copy'

    input:
    file("${params.ref}")
    file("${sample}.raw_indels.vcf")

    output:
    file("${sample}.filtered_indels.vcf")

    script:
    """
    gatk --java-options "-Xmx128g" VariantFiltration \
    -R ${params.ref} \
    -V ${sample}.raw_indels.vcf \
    -O ${sample}.filtered_indels.vcf \
    -filter-name "QD_filter" -filter "QD < 2.0" \
    -filter-name "FS_filter" -filter "FS > 200.0" \
    -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -20.0"
    """
}

process compressAndIndexSNPs {
    container 'biocontainers/bcftools:v1.9-1-deb_cv1'

    publishDir "${params.outpath}/results/", mode: 'copy'

    input:
    file("${sample}.filtered_snps.vcf")

    output:
    file("${sample}.filtered_snps_final.vcf.gz")

    script:
    """
    bcftools view ${sample}.filtered_snps.vcf -Oz -o ${sample}.filtered_snps_final.vcf.gz
    bcftools index ${sample}.filtered_snps_final.vcf.gz
    """
}

process compressAndIndexIndels {
    container 'biocontainers/bcftools:v1.9-1-deb_cv1'

    publishDir "${params.outpath}/results/", mode: 'copy'

    input:
    file("${sample}.filtered_indels.vcf")

    output:
    file("${sample}.filtered_indels_final.vcf.gz")

    script:
    """
    bcftools view ${sample}.filtered_indels.vcf -Oz -o ${sample}.filtered_indels_final.vcf.gz
    bcftools index ${sample}.filtered_indels_final.vcf.gz
    """
}


workflow {
    foo = callVariants()

    snps = extractSNPs(foo.out.file("${sample}.raw_variants.vcf"))
    indels = extractIndels(foo.out.file("${sample}.raw_variants.vcf"))

    filtered_snps = filterSNPs(snps.out.file("${sample}.raw_snps.vcf"))
    filtered_indels = filterIndels(indels.out.file("${sample}.raw_indels.vcf"))

    compressAndIndexSNPs(filtered_snps.out.file("${sample}.filtered_snps.vcf"))
    compressAndIndexIndels(filtered_indels.out.file("${sample}.filtered_indels.vcf"))
}