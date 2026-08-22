import os
import math
from dim_red_tools import *

# Scaling methods to loop over
scaling_method_list = ["standard", "robust", "minmax"]

# Number of principal components
n_components = 7
dim_ids_list = ["dim" + str(i) for i in range(n_components)]

# Paths
pca_input_path = "analysis/pca_all_af2_models/af2/"
output_path = "analysis/pca_all_af2_models_spectral_all_organisms/"

# Large palette — enough colours for the kingdom with the most organisms (bacteria, 16)
palette = sns.color_palette("tab20", 20)

# Complete organism lists per kingdom
organism_animal_list = [
    "Caenorhabditis elegans",
    "Danio rerio",
    "Drosophila melanogaster",
    "Mus musculus",
    "Rattus norvegicus",
    "Homo sapiens",
    "Brugia malayi",
    "Dracunculus medinensis",
    "Onchocerca volvulus",
    "Schistosoma mansoni",
    "Strongyloides stercoralis",
    "Trichuris trichiura",
    "Wuchereria bancrofti",
]

organism_bacteria_list = [
    "Escherichia coli",
    "Helicobacter pylori",
    "Campylobacter jejuni",
    "Enterococcus faecium",
    "Klebsiella pneumoniae",
    "Mycobacterium leprae",
    "Mycobacterium tuberculosis",
    "Mycobacterium ulcerans",
    "Neisseria gonorrhoeae",
    "Nocardia brasiliensis",
    "Pseudomonas aeruginosa",
    "Salmonella typhimurium",
    "Shigella dysenteriae",
    "Staphylococcus aureus",
    "Streptococcus pneumoniae",
    "Haemophilus influenzae",
]

organism_fungi_list = [
    "Candida albicans",
    "Saccharomyces cerevisiae",
    "Cladophialophora carrionii",
    "Fonsecaea pedrosoi",
    "Madurella mycetomatis",
    "Sporothrix schenckii",
    "Ajellomyces capsulatus",
    "Schizosaccharomyces pombe",
    "Paracoccidioides lutzii",
]

organism_plant_list = [
    "Arabidopsis thaliana",
    "Glycine max",
    "Oryza sativa",
    "Zea mays",
]

organism_protozoan_list = [
    "Plasmodium falciparum",
    "Dictyostelium discoideum",
    "Leishmania infantum",
    "Trypanosoma brucei",
    "Trypanosoma cruzi",
]

kingdom_plots = [
    ("Animal",   organism_animal_list,   "spectral_plot_animal"),
    ("Bacteria", organism_bacteria_list, "spectral_plot_bacteria"),
]

for scaling_method in scaling_method_list:

    input_path_scaled = pca_input_path + scaling_method + "/"
    output_path_scaled = output_path + scaling_method + "/"
    os.makedirs(output_path_scaled, exist_ok=True)

    pca_transformed_data = pd.read_csv(input_path_scaled + "pca_transformed_data.csv")

    # Match the same filter applied in pca_analysis.py
    pca_transformed_data = pca_transformed_data[
        ~pca_transformed_data["dssp_bin"].isin(["Hbond Turn", "Bend", "3 10 Helix"])
    ].reset_index(drop=True)

    all_organisms = [org for _, filt_list, _ in kingdom_plots for org in filt_list]

    # Compute global y limits from per-organism means — these are what the lines
    # actually show, so the axis range matches the plotted data, not the raw spread
    org_means = (
        pca_transformed_data[
            pca_transformed_data["organism_scientific_name"].isin(all_organisms)
        ]
        .groupby("organism_scientific_name")[dim_ids_list]
        .mean()
    )
    margin = 0.1 * (org_means.values.max() - org_means.values.min())
    global_ylim = (
        org_means.values.min() - margin,
        org_means.values.max() + margin,
    )

    for title, filt_list, file_name in kingdom_plots:
        ncols = 2
        spectral_plot(
            pca_data=pca_transformed_data.sort_values(
                by="organism_scientific_name", ascending=True
            ),
            group_var="organism_scientific_name",
            value_var_list=dim_ids_list,
            filt_list=filt_list,
            title=title,
            legend_title="",
            output_path=output_path_scaled,
            file_name=file_name,
            palette=palette,
            ncols=ncols,
            ylim=global_ylim,
            legend_loc="center left",
            legend_bbox_to_anchor=(1.02, 0.5),
        )

    print(f"Saved spectral plots for {scaling_method}")
