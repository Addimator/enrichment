rule fgsea:
    input:
        # samples=f"{sleuth_sample}",
        input_effects=f"{input_fgsea}.tsv",
        gene_sets=config["enrichment"]["fgsea"]["gene_sets_file"],
        common_src=workflow.source_path("../scripts/common.R"),
    output:
        enrichment=report(
            f"{output_fgsea}.all-gene-sets.tsv",
            caption="../report/fgsea-table-all.rst",
            category="Gene set enrichment analysis",
            labels=lambda wildcards: get_dynamic_labels(wildcards),
        ),
        rank_ties=report(
            f"{output_fgsea}.rank-ties.tsv",
            caption="../report/fgsea-rank-ties.rst",
            category="Gene set enrichment analysis",
            labels=lambda wildcards: get_dynamic_labels(wildcards),
        ),
        significant=report(
            f"{output_fgsea}.sig-gene-sets.tsv",
            caption="../report/fgsea-table-significant.rst",
            category="Gene set enrichment analysis",
            labels=lambda wildcards: get_dynamic_labels(wildcards),
        ),
        plot=report(
            f"{plots_fgsea}.table-plot.pdf",
            caption="../report/fgsea-table-plot.rst",
            category="Gene set enrichment analysis",
            labels=lambda wildcards: get_dynamic_labels(wildcards),
        ),
        plot_collapsed=report(
            f"{plots_fgsea}.collapsed_pathways.table-plot.pdf",
            caption="../report/fgsea-collapsed-table-plot.rst",
            category="Gene set enrichment analysis",
            labels=lambda wildcards: get_dynamic_labels(wildcards),
        ),
    params:
        bioc_species_pkg=bioc_species_pkg,
        # model=get_model,
        effect_col=lambda wildcards: get_effect_col(wildcards, config),
        gene_set_fdr=config["enrichment"]["fgsea"]["fdr_gene_set"],
        eps=config["enrichment"]["fgsea"]["eps"],
    conda:
        enrichment_env
    log:
        f"logs/enrichment/fgsea/{logs_fgsea}.log",
    threads: 25
    script:
        "../scripts/fgsea.R"


rule fgsea_plot_gene_sets:
    input:
        # samples="results/sleuth/{group}.samples.tsv",
        input_effects=f"{input_fgsea}.tsv",
        gene_sets=config["enrichment"]["fgsea"]["gene_sets_file"],
        sig_gene_sets=f"{output_fgsea}.sig-gene-sets.tsv",
        common_src=workflow.source_path("../scripts/common.R"),
    output:
        report(
            directory(f"{plots_fgsea}"),
            patterns=["{*.gene-set-plot.pdf"],
            caption="../report/plot-fgsea-gene-set.rst",
            category="Gene set enrichment analysis",
            labels=lambda wildcards: get_dynamic_labels(wildcards),
        ),
    params:
        effect_col=lambda wildcards: get_effect_col(wildcards, config),
    conda:
        enrichment_env
    log:
        f"logs/enrichment/fgsea_plot_gene_sets/{logs_fgsea}.log",
    script:
        "../scripts/plot-fgsea-gene-sets.R"


# gene ontology term enrichment analysis


rule ens_gene_to_go:
    input:
        common_src=workflow.source_path("../scripts/common.R"),
    output:
        "resources/ontology/ens_gene_to_go.tsv",
    params:
        bioc_species_pkg=bioc_species_pkg,
    conda:
        enrichment_env
    log:
        "logs/enrichment/resources/ens_gene_to_go.log",
    script:
        "../scripts/ens_gene_to_go.R"


rule download_go_obo:
    output:
        "resources/ontology/gene_ontology.obo",
    params:
        download=config["resources"]["ontology"]["gene_ontology"],
    conda:
        "../envs/curl.yaml"
    log:
        "logs/enrichment/resources/curl.download_go_obo.log",
    shell:
        "( curl --silent -o {output} {params.download} ) 2> {log}"


rule goatools_go_enrichment:
    input:
        obo="resources/ontology/gene_ontology.obo",
        ens_gene_to_go="resources/ontology/ens_gene_to_go.tsv",
        input_effects=f"{input_goatools}.tsv",
    output:
        enrichment=f"{output_goatools}.tsv",
        enrichment_sig_terms=f"{output_goatools}.sig_terms.tsv",
        plot=expand(
            f"{plots_goatools}.pdf",
            ns=["BP", "CC", "MF"],
        ),
    wildcard_constraints:
        # go_term_fdr=".*(?<!_postprocessed)$"
        go_term_fdr="[^_]+",  # Einschränkung, dass 'go_term_fdr' keine Unterstriche enthält
        gene_fdr="[^_]+",
    params:
        species=get_bioc_species_name(),
        effect_col=lambda wildcards: get_effect_col(wildcards, config),
        gene_fdr=config["enrichment"]["goatools"]["fdr_genes"],
        go_term_fdr=config["enrichment"]["goatools"]["fdr_go_terms"],
    conda:
        "../envs/goatools.yaml"
    log:
        f"logs/enrichment/goatools.log",
    script:
        "../scripts/goatools-go-enrichment-analysis.py"
