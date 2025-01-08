import yaml
import pandas as pd
import itertools

input_goatools = lookup(dpath="enrichment/goatools/pathvars/input_file", within=config)
output_goatools = lookup(
    dpath="enrichment/goatools/pathvars/output_file", within=config
)
plots_goatools = lookup(dpath="enrichment/goatools/pathvars/plot_file", within=config)
datavzrd_goatools = lookup(
    dpath="enrichment/goatools/pathvars/datavzrd_file", within=config
)
logs_goatools = lookup(
    dpath="enrichment/goatools/pathvars/log_file_name", within=config
)

input_fgsea = lookup(dpath="enrichment/fgsea/pathvars/input_file", within=config)
output_fgsea = lookup(dpath="enrichment/fgsea/pathvars/output_file", within=config)
plots_fgsea = lookup(dpath="enrichment/fgsea/pathvars/plot_file", within=config)
logs_fgsea = lookup(dpath="enrichment/fgsea/pathvars/log_file_name", within=config)


def get_effect_col(wc, config):
    effect_col_config = lookup(dpath="enrichment/effect_col", within=config)
    if effect_col_config["dynamic"]:
        return eval(effect_col_config["dynamic_expression"])
    else:
        return effect_col_config["static_value"]


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


# Since the input and ouput files are dynamic there is a dynamic number of wildcards
def eval_wildcards(wildcard_str, local_vars):
    return eval(wildcard_str, {}, local_vars)


def get_dynamic_wildcards(method):
    wildcard_definitions = config["enrichment"][method]["pathvars"]["wildcards"]
    wildcard_values = {}

    for wildcard_name, wildcard_expr in wildcard_definitions.items():
        evaluated_value = eval(wildcard_expr)
        if not isinstance(evaluated_value, list):
            evaluated_value = [evaluated_value]

        wildcard_values[wildcard_name] = evaluated_value

    wildcard_combinations = list(itertools.product(*wildcard_values.values()))

    wildcard_combinations_named = [
        dict(zip(wildcard_values.keys(), comb)) for comb in wildcard_combinations
    ]

    return wildcard_combinations_named


# We have to create the wildcards out of the given enrichment paths and the dynamic number of wildcards
def create_paths(method):
    # All possible combinations of wildcards
    combinations = get_dynamic_wildcards(method)

    paths = set()

    # For every wildcard combination create a path
    for combo in combinations:
        if method == "goatools":
            # input_path = input_goatools.format(**combo)
            output_path = output_goatools.format(**combo)
            datavzrd_path = datavzrd_goatools.format(**combo)
            plot_path = (
                plots_goatools.replace("{{", "{").replace("}}", "}").format(**combo)
            )
            paths.update([f"{datavzrd_path}", f"{output_path}.tsv", f"{plot_path}.pdf"])
        elif method == "fgsea":
            # input_path = input_fgsea.format(**combo)
            output_path = output_fgsea.format(**combo)
            plot_path = plots_fgsea.format(**combo)
            paths.update(
                [
                    f"{output_path}.all-gene-sets.tsv",
                    f"{output_path}.sig-gene-sets.tsv",
                    f"{plot_path}.table-plot.pdf",
                    f"{plot_path}",
                ]
            )

        # paths.update([f"{input_path}.tsv", f"{output_path}.tsv", f"{plot_path}.pdf"])

    return paths


def all_input(wildcards):
    """
    Function defining all requested inputs for the rule all (below).
    """

    wanted_input = []
    # request fgsea if 'activated' in config.yaml
    if config["enrichment"]["fgsea"]["activate"]:
        wanted_input.extend(create_paths("fgsea"))

    # request goatools if 'activated' in config.yaml
    if config["enrichment"]["goatools"]["activate"]:
        wanted_input.extend(create_paths("goatools"))

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
