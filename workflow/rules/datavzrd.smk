# Postprocessing GO Enrichment Data
rule postprocess_go_enrichment:
    input:
        enrichment="results/tables/go_terms/{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_fdr_{go_term_fdr}.tsv",
        significant_terms="results/tables/go_terms/{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_fdr_{go_term_fdr}.sig_terms.tsv",
    output:
        "results/tables/go_terms/{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_sig_study_fdr_{go_term_fdr}.tsv",
    conda:
        "../envs/polars.yaml"
    log:
        "logs/yte/postprocess_go_enrichment/{model}_{gene_fdr}.go_term_fdr_{go_term_fdr}.log",
    script:
        "../scripts/postprocess_go_enrichment.py"


# Postprocessing Spia Data
rule postprocess_spia:
    input:
        spia="results/tables/pathways/{model}.pathways.tsv",
    output:
        "results/tables/pathways/{model}.pathways_postprocessed.tsv",
    conda:
        "../envs/pandas.yaml"
    params:
        model=get_model,
    log:
        "logs/yte/postprocess_spia/{model}.log",
    script:
        "../scripts/postprocess_spia.py"


# Generating SPIA Datavzrd Report
rule spia_datavzrd:
    input:
        config=workflow.source_path("../resources/datavzrd/spia-template.yaml"),
        # files required for rendering the given configs
        vega_circle=workflow.source_path(
            "../resources/custom_vega_plots/circle_diagram_genes.json"
        ),
        spia_table="results/tables/pathways/{model}.pathways_postprocessed.tsv",
        vega_waterfall=workflow.source_path(
            "../resources/custom_vega_plots/waterfall_plot_study_items.json"
        ),
    output:
        report(
            directory("results/datavzrd-reports/spia-{model}"),
            htmlindex="index.html",
            caption="../report/spia.rst",
            category="Pathway enrichment",
            patterns=["index.html"],
            labels={"model": "{model}"},
        ),
    log:
        "logs/datavzrd-report/spia-{model}/spia-{model}.log",
    params:
        offer_excel=lookup(within=config, dpath="report/offer_excel", default=False),
        pathway_db=config["enrichment"]["spia"]["pathway_database"],
    wrapper:
        "v3.13.8/utils/datavzrd"


# Generating GO Enrichment Datavzrd Report
rule go_enrichment_datavzrd:
    input:
        config=workflow.source_path("../resources/datavzrd/go-enrichment-template.yaml"),
        vega_bars=workflow.source_path(
            "../resources/custom_vega_plots/horizontal_bars.json"
        ),
        vega_waterfall=workflow.source_path(
            "../resources/custom_vega_plots/waterfall_plot_study_items.json"
        ),
        enrichment="results/tables/go_terms/{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_sig_study_fdr_{go_term_fdr}.tsv",
    output:
        report(
            directory(
                "results/datavzrd-reports/go_enrichment-{model}_{gene_fdr}.go_term_fdr_{go_term_fdr}"
            ),
            htmlindex="index.html",
            caption="../report/go-enrichment-sig_terms.rst",
            category="GO term enrichment",
            subcategory="{model}",
            patterns=["index.html"],
            labels={
                "model": "{model}",
                "gene_fdr": "{gene_fdr}",
                "go_term_fdr": "{go_term_fdr}",
            },
        ),
    log:
        "logs/datavzrd-report/go_enrichment-{model}/go_enrichment-{model}_{gene_fdr}.go_term_fdr_{go_term_fdr}.log",
    params:
        offer_excel=lookup(within=config, dpath="report/offer_excel", default=False),
        samples=get_model_samples,
    wrapper:
        "v3.13.8/utils/datavzrd"
