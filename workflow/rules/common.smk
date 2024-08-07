# from snakemake.utils import validate
import pandas as pd
import yaml

# from pathlib import Path

# ##### load config and sample sheets #####


# # validate(config, schema="../schemas/config.schema.yaml")

# samples = pd.read_csv(config["samples"], sep="\t", dtype=str, comment="#").set_index(
#     "sample", drop=False
# )
# samples.index.names = ["sample_id"]


# def drop_unique_cols(df):
#     singular_cols = df.nunique().loc[(df.nunique().values <= 1)].index
#     return df.drop(singular_cols, axis=1)


# samples = drop_unique_cols(samples)
# validate(samples, schema="../schemas/samples.schema.yaml")

# units = pd.read_csv(config["units"], dtype=str, sep="\t", comment="#").set_index(
#     ["sample", "unit"], drop=False
# )
# units.index.names = ["sample_id", "unit_id"]
# units.index = units.index.set_levels(
#     [i.astype(str) for i in units.index.levels]
# )  # enforce str in index
# validate(units, schema="../schemas/units.schema.yaml")


# report: "../report/workflow.rst"


# ##### wildcard constraints #####


# wildcard_constraints:
#     sample="|".join(samples.index),
#     unit="|".join(units["unit"]),
#     model="|".join(list(config["diffexp"].get("models", [])) + ["all"]),


# ####### helpers ###########

# is_3prime_experiment = (
#     config.get("experiment", dict())
#     .get("3-prime-rna-seq", dict())
#     .get("activate", False)
# )
# three_prime_vendor = (
#     config.get("experiment", dict()).get("3-prime-rna-seq", dict()).get("vendor")
# )

# if is_3prime_experiment:
#     if three_prime_vendor != "lexogen":
#         raise ValueError(
#             f"Currently, only lexogene is supported. Please check the vendor "
#             "in the config file and try again"
#         )


# def check_config():
#     representative_transcripts_keywords = ["canonical", "mostsignificant"]
#     representative_transcripts = config["resources"]["ref"][
#         "representative_transcripts"
#     ]
#     if representative_transcripts not in representative_transcripts_keywords:
#         if not os.path.exists(representative_transcripts):
#             raise ValueError(
#                 f"Invalid value given for resources/ref/representative_transcripts in "
#                 "configuration. Must be 'canonical', 'mostsignificant' or valid path, "
#                 "but {representative_transcripts} does not exist or is not readable."
#             )


# check_config()


# def get_meta_compare_labels(method=""):
#     def _get_labels(wildcards):
#         return {
#             "comparison": method
#             + lookup(
#                 dpath=f"meta_comparisons/comparisons/{wildcards.meta_comp}/label",
#                 within=config,
#             )
#         }

#     return _get_labels


def get_model(wildcards):
    if wildcards.model == "all":
        return {"full": None}
    return config["diffexp"]["models"][wildcards.model]


# def is_single_end(sample, unit):
#     """Determine whether unit is single-end."""
#     bam_paired_not_present = pd.isnull(units.loc[(sample, unit), "bam_paired"])
#     fq2_not_present = pd.isnull(units.loc[(sample, unit), "fq2"])
#     return fq2_not_present and bam_paired_not_present


# def get_fastqs(wildcards):
#     """Get raw FASTQ files from unit sheet."""
#     if not pd.isnull(units.loc[(wildcards.sample, wildcards.unit), "bam_single"]):
#         return f"results/fastq/{wildcards.sample}-{wildcards.unit}.fq.gz"
#     elif not pd.isnull(units.loc[(wildcards.sample, wildcards.unit), "bam_paired"]):
#         fqfrombam1 = f"results/fastq/{wildcards.sample}-{wildcards.unit}.1.fq.gz"
#         fqfrombam2 = f"results/fastq/{wildcards.sample}-{wildcards.unit}.2.fq.gz"
#         return [fqfrombam1, fqfrombam2]
#     elif is_single_end(wildcards.sample, wildcards.unit):
#         return units.loc[(wildcards.sample, wildcards.unit), "fq1"]
#     else:
#         u = units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
#         return [f"{u.fq1}", f"{u.fq2}"]


# def get_all_fastqs(wildcards):
#     for item in units[["sample", "unit"]].itertuples():
#         if is_single_end(item.sample, item.unit):
#             yield f"results/trimmed/{item.sample}-{item.unit}.fastq.gz"
#         else:
#             yield f"results/trimmed/{item.sample}-{item.unit}.1.fastq.gz"
#             yield f"results/trimmed/{item.sample}-{item.unit}.2.fastq.gz"


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


# def get_trimmed(wildcards):
#     if not is_single_end(**wildcards):
#         # paired-end sample
#         return expand(
#             "results/trimmed/{sample}-{unit}.{group}.fastq.gz",
#             group=[1, 2],
#             **wildcards,
#         )
#     # single end sample
#     return expand("results/trimmed/{sample}-{unit}.fastq.gz", **wildcards)


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


# def kallisto_quant_input(wildcards):
#     if is_3prime_experiment:
#         return "results/main_transcript_3prime_reads/{sample}-{unit}.fastq"
#     elif not is_single_end(wildcards.sample, wildcards.unit):
#         return expand(
#             "results/trimmed/{{sample}}-{{unit}}.{group}.fastq.gz", group=[1, 2]
#         )
#     else:
#         return expand("results/trimmed/{sample}-{unit}.fastq.gz", **wildcards)


# def kallisto_params(wildcards, input):
#     extra = config["params"]["kallisto"]
#     if len(input.fastq) == 1 or is_3prime_experiment:
#         extra += " --single --single-overhang --pseudobam"
#         extra += (
#             " --fragment-length {unit.fragment_len_mean} " "--sd {unit.fragment_len_sd}"
#         ).format(unit=units.loc[(wildcards.sample, wildcards.unit)])
#     else:
#         extra += " --fusion"
#     return extra


# def input_genelist(predef_genelist):
#     if config["diffexp"]["genes_of_interest"]["activate"] == True:
#         predef_genelist = config["diffexp"]["genes_of_interest"]["genelist"]
#     else:
#         predef_genelist = []

#     return predef_genelist


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
    # request spia if 'activated' in config.yaml
    if config["enrichment"]["spia"]["activate"]:
        wanted_input.extend(
            expand(
                [
                    "results/tables/pathways/{model}.pathways.tsv",
                    "results/datavzrd-reports/spia-{model}/",
                ],
                model=config["diffexp"]["models"],
            )
        )

    # meta comparisons
    if config["meta_comparisons"]["activate"]:
        wanted_input.extend(
            directory(
                expand(
                    "results/datavzrd-reports/{report_type}_meta_comparison_{meta_comp}",
                    report_type=["go_terms", "pathways"],
                    meta_comp=lookup(
                        dpath="meta_comparisons/comparisons", within=config
                    ),
                )
            ),
        )

    return wanted_input
