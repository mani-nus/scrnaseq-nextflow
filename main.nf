#!/usr/bin/env nextflow
////////////////////////////////////////////////////
/*
 * Configuration parameters
 */
 ////////////////////////////////////////////////////

// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(', ')}"
}

//Check if one of the available aligners is used (alevin, kb, star)
if (params.aligner != 'star' && params.aligner != 'alevin' && params.aligner != 'kb'){
    exit 1, "Invalid aligner option: ${params.aligner}. Valid options: 'star', 'alevin', 'kb'"
}

////////////////////////////////////////////////////
/*
 * Create channels
 */
 ////////////////////////////////////////////////////

// Input read files
Channel
    .fromPath( params.input )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.input}\nNB: Path needs to be enclosed in quotes!\n\NB: Path requires at least one * wildcard!\n" }
    .set { read_files_kb }


// GTF file
if (params.gtf){
    Channel
        .fromPath(params.gtf)
        .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
        .set {gtf_kb}
}

// FASTA file
if (params.fasta){
    Channel
        .fromPath(params.fasta)
        .ifEmpty { exit 1, "FASTA file not found: ${params.fasta}" }
        .set {fasta_kb}
}

////////////////////////////////////////////////////
/*
 * Run kb-python pipeline
 */ 
 ////////////////////////////////////////////////////

params.kb_ref_files = "${params.outdir}/reference/kb"
params.tech = params.type + params.chemistry

process kb_ref {
    label params.aligner
    publishDir params.kb_ref_files, mode: "copy"
    
    input:
        file gtf from gtf_kb
        file fasta from fasta_kb
    
    output:
        file "index.idx" into kb_index
        file "t2g.txt" into kb_t2g
        // file "cDNA.fa" into kb_cDNA_fa

    when: 
        params.aligner == "kb" && !params.kb_prebuilt

    script:
    """
    kb ref --tmp ${params.tmpdir} \\
    --verbose \\
    --workflow standard \\
    -d human \\
    -i index.idx \\
    -g t2g.txt \\
    -f1 cDNA.fa \\
    ${fasta} \\
    ${gtf}
    """
}

process kb_count {
    label params.aligner
    tag "${params.tech}"
    publishDir params.outdir, mode: "copy"
    cpus 3

    input:
        file reads from read_files_kb.collect()
        file index from kb_index
        file t2g from kb_t2g

    // output:
    //   file results 

    when:
        params.aligner == "kb"

    script:
    """
    kb count --tmp ${params.tmpdir} \\
    -i ${index} \\
    -g ${t2g} \\
    -x ${params.tech} \\
    -o ${params.outdir} \\
    --filter bustools \\
    ${reads}
    """
    
}
