# Postprocessing GO Enrichment Data
rule postprocess_go_enrichment:
    input:
        enrichment=f"{output_goatools}.tsv",
        significant_terms=f"{output_goatools}.sig_terms.tsv",
    output:
        enrichment=f"postprocessed_{output_goatools}.tsv",
    conda:
        "../envs/polars.yaml"
    log:
        f"logs/enrichment/postprocess_go_enrichment/{logs_goatools}.log",
    script:
        "../scripts/postprocess_go_enrichment.py"


def get_dynamic_labels(wildcards):
    """
    Dynamically generate labels from all wildcards without hardcoding.
    The wildcards object is a dict-like structure, so we can iterate through its keys and values.
    """
    # Dynamically create labels dictionary from all wildcard keys and values
    print(f"Available wildcards: {wildcards}")
    return {key: value for key, value in wildcards.items()}

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
        enrichment=f"postprocessed_{output_goatools}.tsv",
    output:
        report(
            # directory(lambda wildcards: expand("{datavzrd_goatools}", datavzrd_goatools=wildcards.datavzrd_goatools)),
            directory(datavzrd_goatools),
            htmlindex="index.html",
            caption="../report/go-enrichment-sig_terms.rst",
            category="GO term enrichment",
            # subcategory="{model}",
            patterns=["index.html"],
            labels=lambda wildcards: get_dynamic_labels(wildcards),
            # labels={
            #     "model": "{model}",
            #     "gene_fdr": "{gene_fdr}",
            #     "go_term_fdr": "{go_term_fdr}",
            # },
        ),
    log:
        f"logs/enrichment/go_enrichment_datavzrd/{logs_goatools}.log",
    params:
        offer_excel=lookup(within=config, dpath="report/offer_excel", default=False),
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
