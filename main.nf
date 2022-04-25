#!/usr/bin/env nextflow
////////////////////////////////////////////////////
/*
 * Configuration parameters
 */
 ////////////////////////////////////////////////////
// If kb-aligner used and user specifies that indices have been prebuilt, check if index AND t2g files are supplied
if ( params.aligner == 'kb' && params.kb_prebuilt && !(params.kb_index || params.kb_t2g) ){
    exit 1, "If kb has been prebuilt, you need to supply its precomputed index and the transcript to gene .txt file! Otherwise, leave kb_prebuilt to false to run kb ref"
}

// If kb-aligner used and is not prebuilt, check if species are human or mouse. Otherwise, gtf and fasta files must BOTH be supplied.
if ( params.aligner == "kb" && !params.kb_prebuilt && !(params.gtf || params.fasta) && !(params.species == 'human' || params.species == 'mouse') ){
    exit 1, "If kb has not been prebuilt for non human/mouse species, you need to supply the gtf and genome fasta files."
}

// If genome is specified, check if genome exists in the config file
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
} else {

  gtf_kb = Channel.empty()

}

// FASTA file
if (params.fasta){
    Channel
        .fromPath(params.fasta)
        .ifEmpty { exit 1, "FASTA file not found: ${params.fasta}" }
        .set {fasta_kb}
} else {
  fasta_kb = Channel.empty()
}

// index file
if (params.kb_index){
    Channel
	.fromPath(params.kb_index)
	.ifEmpty { exit 1, "index file not found: ${params.kb_index}"}
	.set {kb_prebuilt_index}
} else {
  kb_prebuilt_index = Channel.empty()
}

// t2g file
if (params.kb_t2g){
   Channel
	.fromPath(params.kb_t2g)
	.ifEmpty { exit 1, "t2g file not found: ${params.kb_t2g}"}
	.set {kb_prebuilt_t2g}
}  else {
   kb_prebuilt_t2g = Channel.empty()
}
////////////////////////////////////////////////////
/*
 * Run kb-python pipeline
 */ 
 ////////////////////////////////////////////////////

params.kb_ref_files = "${params.outdir}/reference/kb"
params.tech = params.type + params.chemistry

process kb_build_index {
    label params.aligner
    label "mid_memory"
    publishDir params.kb_ref_files, mode: "copy"
    
    input:
        file gtf from gtf_kb
        file fasta from fasta_kb
    
    output:
        file "index.idx" into kb_built_index
        file "t2g.txt" into kb_built_t2g

    when: 
        params.aligner == "kb" && !params.kb_prebuilt

    script:
    """
    kb ref --tmp ${params.tmpdir} --verbose --workflow standard -i index.idx -g t2g.txt -f1 cDNA.fa ${fasta} ${gtf}
    """
}

process kb_download_ref {
    label params.aligner
    label "low_memory"
    publishDir params.kb_ref_files, mode: "copy"

    output:
	file "index.idx" into kb_download_index
	file "t2g.txt" into kb_download_t2g

    when:
	params.aligner == "kb" && !params.kb_prebuilt && !(params.gtf && params.fasta)

    script:
    """
    kb ref --tmp ${params.tmpdir} --verbose --workflow standard -d ${params.species} -i index.idx -g t2g.txt
    """
}

process kb_count {
    label params.aligner
    label "high_memory"
    tag "${params.tech}"
    publishDir params.outdir, mode: "copy"

    input:
        file reads from read_files_kb.collect()
        file index from kb_built_index.mix(kb_download_index,kb_prebuilt_index).collect()
        file t2g from kb_built_t2g.mix(kb_download_t2g,kb_prebuilt_t2g).collect()

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
    -t 8 \\
    --filter bustools \\
    ${reads}
    """    
}
