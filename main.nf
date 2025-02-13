nextflow.enable.dsl=2

process callVariants {
    container 'broadinstitute/gatk:4.2.3.0'
    
    publishDir "${params.outpath}/results/", mode: 'copy'

    input:
    tuple val(sample), val(bam), val(bai)

    output:
    path '*'

    script:
    """
    gatk --java-options "-Xmx128g" HaplotypeCaller \
    -R ${params.ref} \
    -I ${bam} \
    --output-mode EMIT_ALL_ACTIVE_SITES \
    -O ${sample}.raw_variants.vcf
    """
}

process extractSNPs {
    container 'broadinstitute/gatk:4.2.3.0'

    publishDir "${params.outpath}/results/", mode: 'copy'

    input:
    file raw_vcf

    output:
    path '*'

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
    file raw_vcf2

    output:
    path '*'

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
    file raw_snps_vcf

    output:
    path '*'

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
    file raw_indels_vcf

    output:
    path '*'

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
    file filtered_snps_vcf

    output:
    path '*'

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
    file filtered_indels_vcf

    output:
    path '*'

    script:
    """
    bcftools view ${sample}.filtered_indels.vcf -Oz -o ${sample}.filtered_indels_final.vcf.gz
    bcftools index ${sample}.filtered_indels_final.vcf.gz
    """
}


workflow {
    sample_ch = Channel.empty()
    sample_tsv = file(params.inputs)    
    sample_ch = extractInput(sample_tsv)

    foo = callVariants(sample_ch)

    snps = extractSNPs(foo)
    indels = extractIndels(foo)

    filtered_snps = filterSNPs(snps)
    filtered_indels = filterIndels(indels)

    compressAndIndexSNPs(filtered_snps)
    compressAndIndexIndels(filtered_indels)
}

def returnFile(it) {
    if (!file(it).exists()) exit 1, "Missing file in TSV file: ${it}, see --help for more information"
    return file(it)
}

def extractInput(tsvFile) {
    Channel.from(tsvFile)
        .splitCsv(sep: '\t')
        .map { row ->           
            def sample    = row[0]
            def bam       = returnFile(row[1])
            def bai       = returnFile(row[2])

            [sample, bam, bai]
        }
}