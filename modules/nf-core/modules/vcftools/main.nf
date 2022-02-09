process VCFTOOLS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::vcftools=0.1.16" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vcftools:0.1.16--he513fc3_4' :
        'quay.io/biocontainers/vcftools:0.1.16--he513fc3_4' }"

    input:
    // Owing to the nature of vcftools we here provide solutions to working with optional bed files and optional
    // alternative variant files, for use with the 'diff' suite of tools.
    // Other optional input files can be utilised in a similar way to below but we do not exhaustively itterate through all
    // possible options. Instead we leave that to the user.
    tuple val(meta), path(variant_file)
    path  bed
    path  diff_variant_file

    output:
    tuple val(meta), path("*.vcf")                    , optional:true, emit: vcf
    tuple val(meta), path("*.bcf")                    , optional:true, emit: bcf
    tuple val(meta), path("*.frq")                    , optional:true, emit: frq
    tuple val(meta), path("*.frq.count")              , optional:true, emit: frq_count
    tuple val(meta), path("*.idepth")                 , optional:true, emit: idepth
    tuple val(meta), path("*.ldepth")                 , optional:true, emit: ldepth
    tuple val(meta), path("*.ldepth.mean")            , optional:true, emit: ldepth_mean
    tuple val(meta), path("*.gdepth")                 , optional:true, emit: gdepth
    tuple val(meta), path("*.hap.ld")                 , optional:true, emit: hap_ld
    tuple val(meta), path("*.geno.ld")                , optional:true, emit: geno_ld
    tuple val(meta), path("*.geno.chisq")             , optional:true, emit: geno_chisq
    tuple val(meta), path("*.list.hap.ld")            , optional:true, emit: list_hap_ld
    tuple val(meta), path("*.list.geno.ld")           , optional:true, emit: list_geno_ld
    tuple val(meta), path("*.interchrom.hap.ld")      , optional:true, emit: interchrom_hap_ld
    tuple val(meta), path("*.interchrom.geno.ld")     , optional:true, emit: interchrom_geno_ld
    tuple val(meta), path("*.TsTv")                   , optional:true, emit: tstv
    tuple val(meta), path("*.TsTv.summary")           , optional:true, emit: tstv_summary
    tuple val(meta), path("*.TsTv.count")             , optional:true, emit: tstv_count
    tuple val(meta), path("*.TsTv.qual")              , optional:true, emit: tstv_qual
    tuple val(meta), path("*.FILTER.summary")         , optional:true, emit: filter_summary
    tuple val(meta), path("*.sites.pi")               , optional:true, emit: sites_pi
    tuple val(meta), path("*.windowed.pi")            , optional:true, emit: windowed_pi
    tuple val(meta), path("*.weir.fst")               , optional:true, emit: weir_fst
    tuple val(meta), path("*.het")                    , optional:true, emit: heterozygosity
    tuple val(meta), path("*.hwe")                    , optional:true, emit: hwe
    tuple val(meta), path("*.Tajima.D")               , optional:true, emit: tajima_d
    tuple val(meta), path("*.ifreqburden")            , optional:true, emit: freq_burden
    tuple val(meta), path("*.LROH")                   , optional:true, emit: lroh
    tuple val(meta), path("*.relatedness")            , optional:true, emit: relatedness
    tuple val(meta), path("*.relatedness2")           , optional:true, emit: relatedness2
    tuple val(meta), path("*.lqual")                  , optional:true, emit: lqual
    tuple val(meta), path("*.imiss")                  , optional:true, emit: missing_individual
    tuple val(meta), path("*.lmiss")                  , optional:true, emit: missing_site
    tuple val(meta), path("*.snpden")                 , optional:true, emit: snp_density
    tuple val(meta), path("*.kept.sites")             , optional:true, emit: kept_sites
    tuple val(meta), path("*.removed.sites")          , optional:true, emit: removed_sites
    tuple val(meta), path("*.singletons")             , optional:true, emit: singeltons
    tuple val(meta), path("*.indel.hist")             , optional:true, emit: indel_hist
    tuple val(meta), path("*.hapcount")               , optional:true, emit: hapcount
    tuple val(meta), path("*.mendel")                 , optional:true, emit: mendel
    tuple val(meta), path("*.FORMAT")                 , optional:true, emit: format
    tuple val(meta), path("*.INFO")                   , optional:true, emit: info
    tuple val(meta), path("*.012")                    , optional:true, emit: genotypes_matrix
    tuple val(meta), path("*.012.indv")               , optional:true, emit: genotypes_matrix_individual
    tuple val(meta), path("*.012.pos")                , optional:true, emit: genotypes_matrix_position
    tuple val(meta), path("*.impute.hap")             , optional:true, emit: impute_hap
    tuple val(meta), path("*.impute.hap.legend")      , optional:true, emit: impute_hap_legend
    tuple val(meta), path("*.impute.hap.indv")        , optional:true, emit: impute_hap_indv
    tuple val(meta), path("*.ldhat.sites")            , optional:true, emit: ldhat_sites
    tuple val(meta), path("*.ldhat.locs")             , optional:true, emit: ldhat_locs
    tuple val(meta), path("*.BEAGLE.GL")              , optional:true, emit: beagle_gl
    tuple val(meta), path("*.BEAGLE.PL")              , optional:true, emit: beagle_pl
    tuple val(meta), path("*.ped")                    , optional:true, emit: ped
    tuple val(meta), path("*.map")                    , optional:true, emit: map_
    tuple val(meta), path("*.tped")                   , optional:true, emit: tped
    tuple val(meta), path("*.tfam")                   , optional:true, emit: tfam
    tuple val(meta), path("*.diff.sites_in_files")    , optional:true, emit: diff_sites_in_files
    tuple val(meta), path("*.diff.indv_in_files")     , optional:true, emit: diff_indv_in_files
    tuple val(meta), path("*.diff.sites")             , optional:true, emit: diff_sites
    tuple val(meta), path("*.diff.indv")              , optional:true, emit: diff_indv
    tuple val(meta), path("*.diff.discordance.matrix"), optional:true, emit: diff_discd_matrix
    tuple val(meta), path("*.diff.switch")            , optional:true, emit: diff_switch_error
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args_list = args.tokenize()

    def bed_arg  = (args.contains('--bed')) ? "--bed ${bed}" :
        (args.contains('--exclude-bed')) ? "--exclude-bed ${bed}" :
        (args.contains('--hapcount')) ? "--hapcount ${bed}" : ''
    args_list.removeIf { it.contains('--bed') }
    args_list.removeIf { it.contains('--exclude-bed') }
    args_list.removeIf { it.contains('--hapcount') }

    def diff_variant_arg = (args.contains('--diff')) ? "--diff ${diff_variant_file}" :
        (args.contains('--gzdiff')) ? "--gzdiff ${diff_variant_file}" :
        (args.contains('--diff-bcf')) ? "--diff-bcf ${diff_variant_file}" : ''
    args_list.removeIf { it.contains('--diff') }
    args_list.removeIf { it.contains('--gzdiff') }
    args_list.removeIf { it.contains('--diff-bcf') }

    def input_file = ("$variant_file".endsWith(".vcf")) ? "--vcf ${variant_file}" :
        ("$variant_file".endsWith(".vcf.gz")) ? "--gzvcf ${variant_file}" :
        ("$variant_file".endsWith(".bcf")) ? "--bcf ${variant_file}" : ''

    """
    vcftools \\
        $input_file \\
        --out $prefix \\
        ${args_list.join(' ')} \\
        $bed_arg \\
        $diff_variant_arg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcftools: \$(echo \$(vcftools --version 2>&1) | sed 's/^.*VCFtools (//;s/).*//')
    END_VERSIONS
    """
}
