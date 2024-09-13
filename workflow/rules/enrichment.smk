## gene set enrichment analysis


rule fgsea:
    input:
        samples="results/sleuth/{model}.samples.tsv",
        diffexp="results/tables/diffexp/{model}.genes-representative.diffexp.tsv",
        gene_sets=config["enrichment"]["fgsea"]["gene_sets_file"],
        common_src=workflow.source_path("../scripts/common.R"),
    output:
        enrichment=report(
            "results/tables/fgsea/{model}.all-gene-sets.tsv",
            caption="../report/fgsea-table-all.rst",
            category="Gene set enrichment analysis",
            labels={"model": "{model}"},
        ),
        rank_ties=report(
            "results/tables/fgsea/{model}.rank-ties.tsv",
            caption="../report/fgsea-rank-ties.rst",
            category="Gene set enrichment analysis",
            labels={"model": "{model}"},
        ),
        significant=report(
            "results/tables/fgsea/{model}.sig-gene-sets.tsv",
            caption="../report/fgsea-table-significant.rst",
            category="Gene set enrichment analysis",
            labels={"model": "{model}"},
        ),
        plot=report(
            "results/plots/fgsea/{model}.table-plot.pdf",
            caption="../report/fgsea-table-plot.rst",
            category="Gene set enrichment analysis",
            labels={"model": "{model}"},
        ),
        plot_collapsed=report(
            "results/plots/fgsea/{model}.collapsed_pathways.table-plot.pdf",
            caption="../report/fgsea-collapsed-table-plot.rst",
            category="Gene set enrichment analysis",
            labels={"model": "{model}"},
        ),
    params:
        bioc_species_pkg=bioc_species_pkg,
        model=get_model,
        gene_set_fdr=config["enrichment"]["fgsea"]["fdr_gene_set"],
        eps=config["enrichment"]["fgsea"]["eps"],
        covariate=lambda w: config["diffexp"]["models"][w.model]["primary_variable"],
    conda:
        enrichment_env
    log:
        "logs/tables/fgsea/{model}.gene-set-enrichment.log",
    threads: 25
    script:
        "../scripts/fgsea.R"


rule fgsea_plot_gene_sets:
    input:
        samples="results/sleuth/{model}.samples.tsv",
        diffexp="results/tables/diffexp/{model}.genes-representative.diffexp.tsv",
        gene_sets=config["enrichment"]["fgsea"]["gene_sets_file"],
        sig_gene_sets="results/tables/fgsea/{model}.sig-gene-sets.tsv",
        common_src=workflow.source_path("../scripts/common.R"),
    output:
        report(
            directory("results/plots/fgsea/{model}"),
            patterns=["{model}.{gene_set}.gene-set-plot.pdf"],
            caption="../report/plot-fgsea-gene-set.rst",
            category="Gene set enrichment analysis",
            labels={"model": "{model}"},
        ),
    params:
        model=get_model,
        covariate=lambda w: config["diffexp"]["models"][w.model]["primary_variable"],
    conda:
        enrichment_env
    log:
        "logs/plots/fgsea/{model}.plot_fgsea_gene_set.log",
    script:
        "../scripts/plot-fgsea-gene-sets.R"


## gene ontology term enrichment analysis


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
        "logs/resources/ens_gene_to_go.log",
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
        "logs/resources/curl.download_go_obo.log",
    shell:
        "( curl --silent -o {output} {params.download} ) 2> {log}"


def get_effect_col(wildcards, config):
    effect_col_config = lookup(dpath="pathvars/effect_col", within=config)

    if effect_col_config["dynamic"]:
        return eval(effect_col_config["dynamic_expression"])
    else:
        return effect_col_config["static_value"]


input_prefix = lookup(dpath="pathvars/input_prefix", within=config)
input_name = lookup(dpath="pathvars/input_name", within=config)
output_prefix = lookup(dpath="pathvars/output_prefix", within=config)
output_name = lookup(dpath="pathvars/output_name", within=config)
plot_prefix = lookup(dpath="pathvars/plot_prefix", within=config)
plot_name = lookup(dpath="pathvars/plot_name", within=config)
logfile = lookup(dpath="pathvars/log_name", within=config)


rule goatools_go_enrichment:
    input:
        obo="resources/ontology/gene_ontology.obo",
        ens_gene_to_go="resources/ontology/ens_gene_to_go.tsv",
        diffexp=f"{input_prefix}{input_name}.tsv",
    output:
        enrichment=f"{output_prefix}{output_name}.tsv",
        enrichment_sig_terms=f"{output_prefix}{output_name}.sig_terms.tsv",
        plot=expand(
            f"{plot_prefix}{plot_name}.pdf",
            ns=["BP", "CC", "MF"],
        ),
    params:
        species=get_bioc_species_name(),
        effect_col=lambda wildcards: get_effect_col(wildcards, config),
        gene_fdr=config["enrichment"]["goatools"]["fdr_genes"],
        go_term_fdr=config["enrichment"]["goatools"]["fdr_go_terms"],
    conda:
        "../envs/goatools.yaml"
    log:
        f"{logfile}.log",
    script:
        "../scripts/goatools-go-enrichment-analysis.py"
