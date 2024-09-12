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


input_prefix = lookup(dpath="pathvars/input_prefix", within=config)
output_prefix = lookup(dpath="pathvars/output_prefix", within=config)
file_name = lookup(dpath="pathvars/file_name", within=config)


rule goatools_go_enrichment:
    input:
        obo="resources/ontology/gene_ontology.obo",
        ens_gene_to_go="resources/ontology/ens_gene_to_go.tsv",
        # diffexp="results/tables/diffexp/{model}.genes-representative.diffexp.tsv",
        # diffexp="results/dmr_calls/BC02-ref-sorted/genes_transcripts/chipseeker_postprocessed.tsv",
        diffexp=f"{input_prefix}{file_name}.tsv",
    output:
        enrichment=f"{output_prefix}{file_name}.tsv",
        # enrichment="results/tables/go_terms/{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_fdr_{go_term_fdr}.tsv",
        # enrichment_sig_terms="results/tables/go_terms/{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_fdr_{go_term_fdr}.sig_terms.tsv",
        # plot=expand(
        # "results/plots/go_terms/{{model}}.go_term_enrichment_{ns}.gene_fdr_{{gene_fdr}}.go_term_fdr_{{go_term_fdr}}.pdf",
        # ns=["BP", "CC", "MF"],
        # ),
        enrichment_sig_terms=f"{output_prefix}{file_name}.sig_terms.tsv",
        plot=expand(f"{output_prefix}_{{ns}}_{file_name}.tsv", ns=["BP", "CC", "MF"]),
    params:
        species=get_bioc_species_name(),
        # model=lambda w: config["diffexp"]["models"][w.model]["primary_variable"],
        effect_col="mean_methylation_difference",
        gene_fdr=0.05,
        go_term_fdr=0.05,
        # gene_fdr=lambda wc: wc.gene_fdr.replace("-", "."),
        # go_term_fdr=lambda wc: wc.go_term_fdr.replace("-", "."),
    conda:
        "../envs/goatools.yaml"
    log:
        f"logs/goatools/tables_and_plots.{file_name}.log",
    script:
        "../scripts/goatools-go-enrichment-analysis.py"
