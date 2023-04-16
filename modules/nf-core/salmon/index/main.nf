process SALMON_INDEX {
    tag "$transcript_fasta"
    label "process_medium"

    conda "bioconda::salmon=1.9.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/salmon:1.9.0--h7e5ed60_1' :
        'quay.io/biocontainers/salmon:1.9.0--h7e5ed60_1' }"

    input:
    path genome_fasta
    tuple val(meta), path(transcript_fasta)

    output:
    tuple val(meta), path("salmon") , emit: index
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gentrome      = "gentrome.fa"
    def get_decoy_ids = ''
    def decoy_arg = ''
    def make_decoy_txt_and_gentrome = ''
    if(params.salmon_decoy_fasta){
        get_decoy_ids = "grep '^>' $genome_fasta | cut -d ' ' -f 1 > decoys.txt"
        make_decoy_txt_and_gentrome = "sed -i.bak -e 's/>//g' decoys.txt; cat $transcript_fasta $genome_fasta > $gentrome"
        decoy_arg = "-d decoys.txt"
        if (genome_fasta.endsWith('.gz')) {
            get_decoy_ids = "grep '^>' <(gunzip -c $genome_fasta) | cut -d ' ' -f 1 > decoys.txt"
            gentrome      = "gentrome.fa.gz"
        }
    } else {
        make_decoy_txt_and_gentrome = "cp $transcript_fasta $gentrome"
    }
    """
    $get_decoy_ids
    $make_decoy_txt_and_gentrome

    salmon \\
        index \\
        --threads $task.cpus \\
        -t $gentrome \\
        $decoy_arg \\
        $args \\
        -i salmon
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
    END_VERSIONS
    """
}
