def all_input(wildcards):
    """
    Function defining all requested inputs for the rule all (below).
    """

    wanted_input = []

    # request goatools if 'activated' in config.yaml
    if config["enrichment"]["goatools"]["activate"]:
        wanted_input.extend(
            expand(
                [
                    "results/tables/go_terms/{model}.go_term_enrichment.gene_fdr_{gene_fdr}.go_term_fdr_{go_term_fdr}.tsv",
                    "results/plots/go_terms/{model}.go_term_enrichment_{go_ns}.gene_fdr_{gene_fdr}.go_term_fdr_{go_term_fdr}.pdf",
                    "results/datavzrd-reports/go_enrichment-{model}_{gene_fdr}.go_term_fdr_{go_term_fdr}",
                ],
                model=config["diffexp"]["models"],
                go_ns=["BP", "CC", "MF"],
                gene_fdr=str(config["enrichment"]["goatools"]["fdr_genes"]).replace(
                    ".", "-"
                ),
                go_term_fdr=str(
                    config["enrichment"]["goatools"]["fdr_go_terms"]
                ).replace(".", "-"),
            )
        )

    # request fgsea if 'activated' in config.yaml
    if config["enrichment"]["fgsea"]["activate"]:
        wanted_input.extend(
            expand(
                [
                    "results/tables/fgsea/{model}.all-gene-sets.tsv",
                    "results/tables/fgsea/{model}.sig-gene-sets.tsv",
                    "results/plots/fgsea/{model}.table-plot.pdf",
                    "results/plots/fgsea/{model}",
                ],
                model=config["diffexp"]["models"],
            )
        )

    # meta comparisons
    if config["meta_comparisons"]["activate"]:
        wanted_input.extend(
            directory(
                expand(
                    "results/datavzrd-reports/go_terms_meta_comparison_{meta_comp}",
                    meta_comp=lookup(
                        dpath="meta_comparisons/comparisons", within=config
                    ),
                )
            ),
        )
    return wanted_input
