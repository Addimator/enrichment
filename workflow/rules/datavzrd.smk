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


# Generating Meta Comparison Datavzrd Reports
rule meta_compare_go_terms_datavzrd:
    input:
        config=lambda wildcards: workflow.source_path(
            f"../resources/datavzrd/meta_comparison-go_terms-template.yaml"
        ),
        table="results/tables/go_terms/meta_compare_{meta_comp}.tsv",
        plot="results/meta_comparison/go_terms/{meta_comp}.json",
    output:
        report(
            directory("results/datavzrd-reports/go_terms_meta_comparison_{meta_comp}"),
            htmlindex="index.html",
            caption="../report/meta_compare.rst",
            category="Comparisons",
            subcategory="{meta_comp}",
            patterns=["index.html"],
            labels=lambda wildcards: get_meta_compare_labels,
        ),
    log:
        "logs/datavzrd-report/meta_comp_go_terms.{meta_comp}.log",
    wrapper:
        "v3.13.8/utils/datavzrd"
