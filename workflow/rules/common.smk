import yaml
import pandas as pd


def get_model(wildcards):
    if wildcards.model == "all":
        return {"full": None}
    return config["diffexp"]["models"][wildcards.model]


def get_model_samples(wildcards):
    samples = pd.read_csv(config["samples"], sep="\t", dtype=str, comment="#")
    units = pd.read_csv(config["units"], sep="\t", dtype=str, comment="#")
    sample_file = units.merge(samples, on="sample")
    sample_file["sample_name"] = sample_file.apply(
        lambda row: "{}-{}".format(row["sample"], row["unit"]), axis=1
    )
    gps = config["diffexp"]["models"][wildcards.model]["primary_variable"]
    sample_groups = sample_file.loc[sample_file[gps].notnull(), ["sample_name"]]
    samples = sample_groups["sample_name"].values
    return samples


def get_bioc_species_name():
    first_letter = config["resources"]["ref"]["species"][0]
    subspecies = config["resources"]["ref"]["species"].split("_")[1]
    return first_letter + subspecies


def get_bioc_species_pkg():
    """Get the package bioconductor package name for the the species in config.yaml"""
    species_letters = get_bioc_species_name()[0:2].capitalize()
    return "org.{species}.eg.db".format(species=species_letters)


def render_enrichment_env():
    species_pkg = f"bioconductor-{get_bioc_species_pkg()}"
    with open(workflow.source_path("../envs/enrichment.yaml")) as f:
        env = yaml.load(f, Loader=yaml.SafeLoader)
    env["dependencies"].append(species_pkg)
    env_path = Path("resources/envs/enrichment.yaml")
    env_path.parent.mkdir(parents=True, exist_ok=True)
    with open(env_path, "w") as f:
        yaml.dump(env, f)
    return env_path.absolute()


bioc_species_pkg = get_bioc_species_pkg()
enrichment_env = render_enrichment_env()


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
