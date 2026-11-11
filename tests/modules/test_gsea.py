"""Tests for flyte_bio.modules.gsea.

The cases and expected values are upstream's gsea/gsea tests (modules commit
9a1cc9f), with their config's ``-rnd_seed 10``: the md5s of the tables and
of the correlation heat map page, and the names of the HTML reports and
plots (as upstream, their contents vary). The per-gene-set files are read
from the module's tarball of them.
"""

import hashlib
import io
import tarfile
from dataclasses import dataclass

from flyte_bio.modules.gsea import GseaResult, gsea
from tests.framework import assert_md5, env, fixture

GENE_SET_ANALYSIS = "genomics/mus_musculus/gene_set_analysis/"
PREFIX = "Condition_genotype_WT_KO"


@dataclass
class Expected:
    heat_map_corr_plot: str
    report_ref: str
    report_target: str
    ranked_gene_list: str
    gene_set_sizes: str
    gene_set_tsvs: dict[str, str]  # gene_sets_<prefix><name>.tsv md5s
    gene_set_htmls: list[str]  # gene_sets_<prefix><name>.html
    gene_set_heatmaps: list[str]  # gene_sets_<prefix><name>.png
    enplots: list[str]  # <prefix>enplot_<name>.png
    rnd_es_dists: list[str]  # <prefix>gset_rnd_es_dist_<name>.png


HALLMARK = Expected(
    heat_map_corr_plot="129c561984542746bc50c3f6db9c7c5e",
    report_ref="674b21a1646cc554d6686ddd67c55dab",
    report_target="a4e7de6cf0400befc6e52e379c89c178",
    ranked_gene_list="3ae0171ae04398f350cd3e0f183b9d67",
    gene_set_sizes="61ca77f6b13620a78170c0b9366b8425",
    gene_set_tsvs={
        "HALLMARK_ADIPOGENESIS": "0e7035e5cbaba4a35d87538da823e05c",
        "HALLMARK_ANDROGEN_RESPONSE": "46cee22eddeac12e1d8643e539c6c397",
        "HALLMARK_ANGIOGENESIS": "3b4a7f1b28282bf886c9c06d73f61e45",
        "HALLMARK_APICAL_JUNCTION": "e9d02a848a5aca73a22a6aef7473defd",
        "HALLMARK_APICAL_SURFACE": "84eed15b2aad41077b2a4c5367c05782",
        "HALLMARK_APOPTOSIS": "302c2afbf1b77c6eccd7983fb657007c",
        "HALLMARK_BILE_ACID_METABOLISM": "967e4be880742f32720419fe879b7cc7",
        "HALLMARK_COAGULATION": "6c5bece0e26c1e1b32d9d47b10651baf",
        "HALLMARK_COMPLEMENT": "0d3dfe9073ddbbfa3a6f0c7bcb4c8601",
        "HALLMARK_DNA_REPAIR": "975b4b3905685e73500ec49cd54a2c03",
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION": "122c6786499f76f5ec0c6ce604b52aa0",
        "HALLMARK_ESTROGEN_RESPONSE_EARLY": "7fea0527089b002ccf840223a60ecbf4",
        "HALLMARK_ESTROGEN_RESPONSE_LATE": "a838545de228235bc4f8830b5777e893",
        "HALLMARK_FATTY_ACID_METABOLISM": "926ddf4b74ad3a6d76bc99e82ae6aaf3",
        "HALLMARK_GLYCOLYSIS": "4dcdce8b8db83dd971800f2ff18ccbeb",
        "HALLMARK_HEDGEHOG_SIGNALING": "2d508b31d0ec00228db110a0e9e2f059",
        "HALLMARK_HEME_METABOLISM": "bb723be28feecbbc68909fdb17f8e2f4",
        "HALLMARK_HYPOXIA": "e3a7d53519ee12377377659743d9cf4c",
        "HALLMARK_IL2_STAT5_SIGNALING": "4bcc7835887602c732c85ef2db4f5584",
        "HALLMARK_IL6_JAK_STAT3_SIGNALING": "18e50d2154c81700b2c44ee848ced133",
        "HALLMARK_INTERFERON_ALPHA_RESPONSE": "fb6fcf993118271ee3340e9a12de0ed7",
        "HALLMARK_KRAS_SIGNALING_DN": "d0649d9c33c0f81bf4e8e2dcdb6a036d",
        "HALLMARK_KRAS_SIGNALING_UP": "1b4411d95e2bf08db2bae06996d6729c",
        "HALLMARK_MYOGENESIS": "564c3e924517993e6eac6819a5f53d92",
        "HALLMARK_NOTCH_SIGNALING": "e91e740dac62555b603bf44c427e671e",
        "HALLMARK_OXIDATIVE_PHOSPHORYLATION": "75095d76fc23eccb98d67e045839d4a9",
        "HALLMARK_P53_PATHWAY": "ce4a9fc08c539fc4d930f81e08a2db55",
        "HALLMARK_PANCREAS_BETA_CELLS": "9df68cbfc6e4666d18ff6d214c5f1544",
        "HALLMARK_PEROXISOME": "f0e60a2f8bc511be09ff8f04886eca61",
        "HALLMARK_PI3K_AKT_MTOR_SIGNALING": "f76d919c3e1342304bb69e50a58c81d9",
        "HALLMARK_PROTEIN_SECRETION": "f16a1fcc6e9b4fd976b3c9fefd0f8dc9",
        "HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY": "66bcbea77b774b7d0542abeb460bfef8",
        "HALLMARK_TGF_BETA_SIGNALING": "5970edfc5faad9bb73408d75a2f40350",
        "HALLMARK_UV_RESPONSE_DN": "bb751fac5d95a22e6872ef347038cd4e",
        "HALLMARK_UV_RESPONSE_UP": "8cfe52d138bc42ec8fd7d2c6fbf6761f",
        "HALLMARK_WNT_BETA_CATENIN_SIGNALING": "b45d874ea3fdd7cd4f8bd09843eaa941",
        "HALLMARK_XENOBIOTIC_METABOLISM": "ca28e2bd7b3584d1354e97b09a37f8af",
        "Symbol_to_probe_set_mapping_details": "a8bed433fe3e8081495c82226ec44917",
    },
    gene_set_htmls=[
        "HALLMARK_ADIPOGENESIS",
        "HALLMARK_ANDROGEN_RESPONSE",
        "HALLMARK_ANGIOGENESIS",
        "HALLMARK_APICAL_JUNCTION",
        "HALLMARK_APICAL_SURFACE",
        "HALLMARK_APOPTOSIS",
        "HALLMARK_BILE_ACID_METABOLISM",
        "HALLMARK_COAGULATION",
        "HALLMARK_COMPLEMENT",
        "HALLMARK_DNA_REPAIR",
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
        "HALLMARK_ESTROGEN_RESPONSE_EARLY",
        "HALLMARK_ESTROGEN_RESPONSE_LATE",
        "HALLMARK_FATTY_ACID_METABOLISM",
        "HALLMARK_GLYCOLYSIS",
        "HALLMARK_HEDGEHOG_SIGNALING",
        "HALLMARK_HEME_METABOLISM",
        "HALLMARK_HYPOXIA",
        "HALLMARK_IL2_STAT5_SIGNALING",
        "HALLMARK_IL6_JAK_STAT3_SIGNALING",
        "HALLMARK_INTERFERON_ALPHA_RESPONSE",
        "HALLMARK_KRAS_SIGNALING_DN",
        "HALLMARK_KRAS_SIGNALING_UP",
        "HALLMARK_MYOGENESIS",
        "HALLMARK_NOTCH_SIGNALING",
        "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
        "HALLMARK_P53_PATHWAY",
        "HALLMARK_PANCREAS_BETA_CELLS",
        "HALLMARK_PEROXISOME",
        "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
        "HALLMARK_PROTEIN_SECRETION",
        "HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY",
        "HALLMARK_TGF_BETA_SIGNALING",
        "HALLMARK_UV_RESPONSE_DN",
        "HALLMARK_UV_RESPONSE_UP",
        "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
        "HALLMARK_XENOBIOTIC_METABOLISM",
    ],
    gene_set_heatmaps=[
        "HALLMARK_ADIPOGENESIS_13",
        "HALLMARK_ANDROGEN_RESPONSE_100",
        "HALLMARK_ANGIOGENESIS_94",
        "HALLMARK_APICAL_JUNCTION_55",
        "HALLMARK_APICAL_SURFACE_82",
        "HALLMARK_APOPTOSIS_40",
        "HALLMARK_BILE_ACID_METABOLISM_31",
        "HALLMARK_COAGULATION_34",
        "HALLMARK_COMPLEMENT_91",
        "HALLMARK_DNA_REPAIR_43",
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_67",
        "HALLMARK_ESTROGEN_RESPONSE_EARLY_58",
        "HALLMARK_ESTROGEN_RESPONSE_LATE_70",
        "HALLMARK_FATTY_ACID_METABOLISM_16",
        "HALLMARK_GLYCOLYSIS_97",
        "HALLMARK_HEDGEHOG_SIGNALING_61",
        "HALLMARK_HEME_METABOLISM_28",
        "HALLMARK_HYPOXIA_103",
        "HALLMARK_IL2_STAT5_SIGNALING_106",
        "HALLMARK_IL6_JAK_STAT3_SIGNALING_46",
        "HALLMARK_INTERFERON_ALPHA_RESPONSE_52",
        "HALLMARK_KRAS_SIGNALING_DN_109",
        "HALLMARK_KRAS_SIGNALING_UP_85",
        "HALLMARK_MYOGENESIS_64",
        "HALLMARK_NOTCH_SIGNALING_73",
        "HALLMARK_OXIDATIVE_PHOSPHORYLATION_4",
        "HALLMARK_P53_PATHWAY_22",
        "HALLMARK_PANCREAS_BETA_CELLS_37",
        "HALLMARK_PEROXISOME_7",
        "HALLMARK_PI3K_AKT_MTOR_SIGNALING_49",
        "HALLMARK_PROTEIN_SECRETION_10",
        "HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY_25",
        "HALLMARK_TGF_BETA_SIGNALING_79",
        "HALLMARK_UV_RESPONSE_DN_76",
        "HALLMARK_UV_RESPONSE_UP_88",
        "HALLMARK_WNT_BETA_CATENIN_SIGNALING_112",
        "HALLMARK_XENOBIOTIC_METABOLISM_19",
    ],
    enplots=[
        "HALLMARK_ADIPOGENESIS_12",
        "HALLMARK_ANDROGEN_RESPONSE_99",
        "HALLMARK_ANGIOGENESIS_93",
        "HALLMARK_APICAL_JUNCTION_54",
        "HALLMARK_APICAL_SURFACE_81",
        "HALLMARK_APOPTOSIS_39",
        "HALLMARK_BILE_ACID_METABOLISM_30",
        "HALLMARK_COAGULATION_33",
        "HALLMARK_COMPLEMENT_90",
        "HALLMARK_DNA_REPAIR_42",
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_66",
        "HALLMARK_ESTROGEN_RESPONSE_EARLY_57",
        "HALLMARK_ESTROGEN_RESPONSE_LATE_69",
        "HALLMARK_FATTY_ACID_METABOLISM_15",
        "HALLMARK_GLYCOLYSIS_96",
        "HALLMARK_HEDGEHOG_SIGNALING_60",
        "HALLMARK_HEME_METABOLISM_27",
        "HALLMARK_HYPOXIA_102",
        "HALLMARK_IL2_STAT5_SIGNALING_105",
        "HALLMARK_IL6_JAK_STAT3_SIGNALING_45",
        "HALLMARK_INTERFERON_ALPHA_RESPONSE_51",
        "HALLMARK_KRAS_SIGNALING_DN_108",
        "HALLMARK_KRAS_SIGNALING_UP_84",
        "HALLMARK_MYOGENESIS_63",
        "HALLMARK_NOTCH_SIGNALING_72",
        "HALLMARK_OXIDATIVE_PHOSPHORYLATION_3",
        "HALLMARK_P53_PATHWAY_21",
        "HALLMARK_PANCREAS_BETA_CELLS_36",
        "HALLMARK_PEROXISOME_6",
        "HALLMARK_PI3K_AKT_MTOR_SIGNALING_48",
        "HALLMARK_PROTEIN_SECRETION_9",
        "HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY_24",
        "HALLMARK_TGF_BETA_SIGNALING_78",
        "HALLMARK_UV_RESPONSE_DN_75",
        "HALLMARK_UV_RESPONSE_UP_87",
        "HALLMARK_WNT_BETA_CATENIN_SIGNALING_111",
        "HALLMARK_XENOBIOTIC_METABOLISM_18",
    ],
    rnd_es_dists=[
        "101",
        "104",
        "107",
        "11",
        "110",
        "113",
        "14",
        "17",
        "20",
        "23",
        "26",
        "29",
        "32",
        "35",
        "38",
        "41",
        "44",
        "47",
        "5",
        "50",
        "53",
        "56",
        "59",
        "62",
        "65",
        "68",
        "71",
        "74",
        "77",
        "8",
        "80",
        "83",
        "86",
        "89",
        "92",
        "95",
        "98",
    ],
)
MULTIPLE = Expected(
    heat_map_corr_plot="129c561984542746bc50c3f6db9c7c5e",
    report_ref="0aadcf1751fc08ba9f3b253ea3c5d3a6",
    report_target="21297fb2b50155ed348d83f10ab1ae8b",
    ranked_gene_list="3ae0171ae04398f350cd3e0f183b9d67",
    gene_set_sizes="42fe75837b300bc9c0be6eed72135a1b",
    gene_set_tsvs={
        "HALLMARK_ADIPOGENESIS": "0e7035e5cbaba4a35d87538da823e05c",
        "HALLMARK_ANGIOGENESIS": "3b4a7f1b28282bf886c9c06d73f61e45",
        "HALLMARK_APICAL_JUNCTION": "e9d02a848a5aca73a22a6aef7473defd",
        "HALLMARK_APICAL_SURFACE": "84eed15b2aad41077b2a4c5367c05782",
        "HALLMARK_BILE_ACID_METABOLISM": "967e4be880742f32720419fe879b7cc7",
        "HALLMARK_COAGULATION": "6c5bece0e26c1e1b32d9d47b10651baf",
        "HALLMARK_COMPLEMENT": "0d3dfe9073ddbbfa3a6f0c7bcb4c8601",
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION": "122c6786499f76f5ec0c6ce604b52aa0",
        "HALLMARK_ESTROGEN_RESPONSE_EARLY": "7fea0527089b002ccf840223a60ecbf4",
        "HALLMARK_ESTROGEN_RESPONSE_LATE": "a838545de228235bc4f8830b5777e893",
        "HALLMARK_FATTY_ACID_METABOLISM": "926ddf4b74ad3a6d76bc99e82ae6aaf3",
        "HALLMARK_GLYCOLYSIS": "4dcdce8b8db83dd971800f2ff18ccbeb",
        "HALLMARK_HEDGEHOG_SIGNALING": "2d508b31d0ec00228db110a0e9e2f059",
        "HALLMARK_HEME_METABOLISM": "bb723be28feecbbc68909fdb17f8e2f4",
        "HALLMARK_KRAS_SIGNALING_UP": "1b4411d95e2bf08db2bae06996d6729c",
        "HALLMARK_MYOGENESIS": "564c3e924517993e6eac6819a5f53d92",
        "HALLMARK_NOTCH_SIGNALING": "e91e740dac62555b603bf44c427e671e",
        "HALLMARK_OXIDATIVE_PHOSPHORYLATION": "75095d76fc23eccb98d67e045839d4a9",
        "HALLMARK_P53_PATHWAY": "ce4a9fc08c539fc4d930f81e08a2db55",
        "HALLMARK_PEROXISOME": "f0e60a2f8bc511be09ff8f04886eca61",
        "HALLMARK_PROTEIN_SECRETION": "f16a1fcc6e9b4fd976b3c9fefd0f8dc9",
        "HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY": "66bcbea77b774b7d0542abeb460bfef8",
        "HALLMARK_TGF_BETA_SIGNALING": "5970edfc5faad9bb73408d75a2f40350",
        "HALLMARK_UV_RESPONSE_DN": "bb751fac5d95a22e6872ef347038cd4e",
        "HALLMARK_UV_RESPONSE_UP": "8cfe52d138bc42ec8fd7d2c6fbf6761f",
        "HALLMARK_XENOBIOTIC_METABOLISM": "ca28e2bd7b3584d1354e97b09a37f8af",
        "MP_DECREASED_INCIDENCE_OF_TUMORS_BY_CHEMICAL_INDUCTION": "7b29d2bc9cfcc1ff4bcdabd8e46a970b",
        "MP_DECREASED_METASTATIC_POTENTIAL": "c70ceae8d96717f33983442970478710",
        "MP_DECREASED_TUMOR_GROWTH_SIZE": "da6905b554d983629ef19c275fe8724a",
        "MP_DECREASED_TUMOR_NECROSIS_FACTOR_SECRETION": "f318a881241f725ad1de5febd89cd0d2",
        "MP_EMBRYO_TUMOR": "24f94c2c4456aa1f5b773eae810bffb5",
        "MP_INCREASED_ADENOMA_INCIDENCE": "88e9f1759362869157702034e0273169",
        "MP_INCREASED_CARCINOMA_INCIDENCE": "693c20644fd2d9a20f3aa5cf87e3037c",
        "MP_INCREASED_GASTROINTESTINAL_TUMOR_INCIDENCE": "af36f27bcb01c52e9e8bf4b417c26aef",
        "MP_INCREASED_HEMANGIOMA_INCIDENCE": "7741038545f1f9d7bde58f5e338dbdab",
        "MP_INCREASED_HEMANGIOSARCOMA_INCIDENCE": "2f9971d0e21ababd7480812dfda76f6d",
        "MP_INCREASED_HEPATOCELLULAR_CARCINOMA_INCIDENCE": "4ffc87c24f7403067a148d0e92ea7075",
        "MP_INCREASED_LIVER_TUMOR_INCIDENCE": "b5189d8034e7b730f9c35f85aeffe2ae",
        "MP_INCREASED_METASTATIC_POTENTIAL": "0a4a0048b3e87f093ec830aa6707dc8f",
        "MP_INCREASED_TUMOR_NECROSIS_FACTOR_SECRETION": "f6cad28a554d33db2ef536944ee5f884",
        "Symbol_to_probe_set_mapping_details": "a8bed433fe3e8081495c82226ec44917",
    },
    gene_set_htmls=[
        "HALLMARK_ADIPOGENESIS",
        "HALLMARK_ANGIOGENESIS",
        "HALLMARK_APICAL_JUNCTION",
        "HALLMARK_APICAL_SURFACE",
        "HALLMARK_BILE_ACID_METABOLISM",
        "HALLMARK_COAGULATION",
        "HALLMARK_COMPLEMENT",
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
        "HALLMARK_ESTROGEN_RESPONSE_EARLY",
        "HALLMARK_ESTROGEN_RESPONSE_LATE",
        "HALLMARK_FATTY_ACID_METABOLISM",
        "HALLMARK_GLYCOLYSIS",
        "HALLMARK_HEDGEHOG_SIGNALING",
        "HALLMARK_HEME_METABOLISM",
        "HALLMARK_KRAS_SIGNALING_UP",
        "HALLMARK_MYOGENESIS",
        "HALLMARK_NOTCH_SIGNALING",
        "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
        "HALLMARK_P53_PATHWAY",
        "HALLMARK_PEROXISOME",
        "HALLMARK_PROTEIN_SECRETION",
        "HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY",
        "HALLMARK_TGF_BETA_SIGNALING",
        "HALLMARK_UV_RESPONSE_DN",
        "HALLMARK_UV_RESPONSE_UP",
        "HALLMARK_XENOBIOTIC_METABOLISM",
        "MP_DECREASED_INCIDENCE_OF_TUMORS_BY_CHEMICAL_INDUCTION",
        "MP_DECREASED_METASTATIC_POTENTIAL",
        "MP_DECREASED_TUMOR_GROWTH_SIZE",
        "MP_DECREASED_TUMOR_NECROSIS_FACTOR_SECRETION",
        "MP_EMBRYO_TUMOR",
        "MP_INCREASED_ADENOMA_INCIDENCE",
        "MP_INCREASED_CARCINOMA_INCIDENCE",
        "MP_INCREASED_GASTROINTESTINAL_TUMOR_INCIDENCE",
        "MP_INCREASED_HEMANGIOMA_INCIDENCE",
        "MP_INCREASED_HEMANGIOSARCOMA_INCIDENCE",
        "MP_INCREASED_HEPATOCELLULAR_CARCINOMA_INCIDENCE",
        "MP_INCREASED_LIVER_TUMOR_INCIDENCE",
        "MP_INCREASED_METASTATIC_POTENTIAL",
        "MP_INCREASED_TUMOR_NECROSIS_FACTOR_SECRETION",
    ],
    gene_set_heatmaps=[
        "HALLMARK_ADIPOGENESIS_13",
        "HALLMARK_ANGIOGENESIS_115",
        "HALLMARK_APICAL_JUNCTION_64",
        "HALLMARK_APICAL_SURFACE_100",
        "HALLMARK_BILE_ACID_METABOLISM_37",
        "HALLMARK_COAGULATION_49",
        "HALLMARK_COMPLEMENT_112",
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_79",
        "HALLMARK_ESTROGEN_RESPONSE_EARLY_67",
        "HALLMARK_ESTROGEN_RESPONSE_LATE_82",
        "HALLMARK_FATTY_ACID_METABOLISM_16",
        "HALLMARK_GLYCOLYSIS_118",
        "HALLMARK_HEDGEHOG_SIGNALING_70",
        "HALLMARK_HEME_METABOLISM_31",
        "HALLMARK_KRAS_SIGNALING_UP_103",
        "HALLMARK_MYOGENESIS_76",
        "HALLMARK_NOTCH_SIGNALING_91",
        "HALLMARK_OXIDATIVE_PHOSPHORYLATION_4",
        "HALLMARK_P53_PATHWAY_25",
        "HALLMARK_PEROXISOME_7",
        "HALLMARK_PROTEIN_SECRETION_10",
        "HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY_28",
        "HALLMARK_TGF_BETA_SIGNALING_97",
        "HALLMARK_UV_RESPONSE_DN_94",
        "HALLMARK_UV_RESPONSE_UP_109",
        "HALLMARK_XENOBIOTIC_METABOLISM_22",
        "MP_DECREASED_INCIDENCE_OF_TUMORS_BY_CHEMICAL_INDUCTION_43",
        "MP_DECREASED_METASTATIC_POTENTIAL_106",
        "MP_DECREASED_TUMOR_GROWTH_SIZE_121",
        "MP_DECREASED_TUMOR_NECROSIS_FACTOR_SECRETION_55",
        "MP_EMBRYO_TUMOR_73",
        "MP_INCREASED_ADENOMA_INCIDENCE_58",
        "MP_INCREASED_CARCINOMA_INCIDENCE_52",
        "MP_INCREASED_GASTROINTESTINAL_TUMOR_INCIDENCE_85",
        "MP_INCREASED_HEMANGIOMA_INCIDENCE_19",
        "MP_INCREASED_HEMANGIOSARCOMA_INCIDENCE_61",
        "MP_INCREASED_HEPATOCELLULAR_CARCINOMA_INCIDENCE_34",
        "MP_INCREASED_LIVER_TUMOR_INCIDENCE_46",
        "MP_INCREASED_METASTATIC_POTENTIAL_40",
        "MP_INCREASED_TUMOR_NECROSIS_FACTOR_SECRETION_88",
    ],
    enplots=[
        "HALLMARK_ADIPOGENESIS_12",
        "HALLMARK_ANGIOGENESIS_114",
        "HALLMARK_APICAL_JUNCTION_63",
        "HALLMARK_APICAL_SURFACE_99",
        "HALLMARK_BILE_ACID_METABOLISM_36",
        "HALLMARK_COAGULATION_48",
        "HALLMARK_COMPLEMENT_111",
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_78",
        "HALLMARK_ESTROGEN_RESPONSE_EARLY_66",
        "HALLMARK_ESTROGEN_RESPONSE_LATE_81",
        "HALLMARK_FATTY_ACID_METABOLISM_15",
        "HALLMARK_GLYCOLYSIS_117",
        "HALLMARK_HEDGEHOG_SIGNALING_69",
        "HALLMARK_HEME_METABOLISM_30",
        "HALLMARK_KRAS_SIGNALING_UP_102",
        "HALLMARK_MYOGENESIS_75",
        "HALLMARK_NOTCH_SIGNALING_90",
        "HALLMARK_OXIDATIVE_PHOSPHORYLATION_3",
        "HALLMARK_P53_PATHWAY_24",
        "HALLMARK_PEROXISOME_6",
        "HALLMARK_PROTEIN_SECRETION_9",
        "HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY_27",
        "HALLMARK_TGF_BETA_SIGNALING_96",
        "HALLMARK_UV_RESPONSE_DN_93",
        "HALLMARK_UV_RESPONSE_UP_108",
        "HALLMARK_XENOBIOTIC_METABOLISM_21",
        "MP_DECREASED_INCIDENCE_OF_TUMORS_BY_CHEMICAL_INDUCTION_42",
        "MP_DECREASED_METASTATIC_POTENTIAL_105",
        "MP_DECREASED_TUMOR_GROWTH_SIZE_120",
        "MP_DECREASED_TUMOR_NECROSIS_FACTOR_SECRETION_54",
        "MP_EMBRYO_TUMOR_72",
        "MP_INCREASED_ADENOMA_INCIDENCE_57",
        "MP_INCREASED_CARCINOMA_INCIDENCE_51",
        "MP_INCREASED_GASTROINTESTINAL_TUMOR_INCIDENCE_84",
        "MP_INCREASED_HEMANGIOMA_INCIDENCE_18",
        "MP_INCREASED_HEMANGIOSARCOMA_INCIDENCE_60",
        "MP_INCREASED_HEPATOCELLULAR_CARCINOMA_INCIDENCE_33",
        "MP_INCREASED_LIVER_TUMOR_INCIDENCE_45",
        "MP_INCREASED_METASTATIC_POTENTIAL_39",
        "MP_INCREASED_TUMOR_NECROSIS_FACTOR_SECRETION_87",
    ],
    rnd_es_dists=[
        "101",
        "104",
        "107",
        "11",
        "110",
        "113",
        "116",
        "119",
        "122",
        "14",
        "17",
        "20",
        "23",
        "26",
        "29",
        "32",
        "35",
        "38",
        "41",
        "44",
        "47",
        "5",
        "50",
        "53",
        "56",
        "59",
        "62",
        "65",
        "68",
        "71",
        "74",
        "77",
        "8",
        "80",
        "83",
        "86",
        "89",
        "92",
        "95",
        "98",
    ],
)


async def check(result: GseaResult, expected: Expected) -> None:
    names = {
        "rpt": (result.rpt, f"{PREFIX}.Gsea.rpt"),
        "index": (result.index_html, f"{PREFIX}index.html"),
        "ref html": (result.report_html_ref, f"{PREFIX}gsea_report_for_WT.html"),
        "target html": (result.report_html_target, f"{PREFIX}gsea_report_for_KO.html"),
        "histogram": (result.histogram, f"{PREFIX}global_es_histogram.png"),
        "heatmap": (result.heatmap, f"{PREFIX}heat_map_1.png"),
        "pvalues": (result.pvalues_vs_nes_plot, f"{PREFIX}pvalues_vs_nes_plot.png"),
        "ranked corr": (result.ranked_list_corr, f"{PREFIX}ranked_list_corr_2.png"),
    }
    for label, (file, name) in names.items():
        assert file.name == name, f"{label}: {file.name}, expected {name}"
    assert result.butterfly_plot is not None and result.butterfly_plot.name == f"{PREFIX}butterfly_plot.png"
    snapshots = [f.name for f in result.snapshots]
    assert snapshots == [f"{PREFIX}neg_snapshot.html", f"{PREFIX}pos_snapshot.html"], snapshots
    assert result.archive is None, "a zip report without -zip_report"
    await assert_md5(result.heat_map_corr_plot, expected.heat_map_corr_plot, label="heat_map_corr_plot")
    await assert_md5(result.report_tsv_ref, expected.report_ref, label="report for WT")
    await assert_md5(result.report_tsv_target, expected.report_target, label="report for KO")
    assert result.ranked_gene_list.name == f"{PREFIX}ranked_gene_list_KO_versus_WT.tsv"
    await assert_md5(result.ranked_gene_list, expected.ranked_gene_list, label="ranked gene list")
    await assert_md5(result.gene_set_sizes, expected.gene_set_sizes, label="gene set sizes")

    async with result.gene_sets.open("rb") as fh:
        data = bytes(await fh.read())
    with tarfile.open(fileobj=io.BytesIO(data)) as tar:
        members = {m.name: m for m in tar.getmembers()}
        tsvs = {
            name[len(f"gene_sets_{PREFIX}") : -len(".tsv")]: hashlib.md5(tar.extractfile(m).read()).hexdigest()
            for name, m in members.items()
            if name.startswith("gene_sets_") and name.endswith(".tsv")
        }
    assert tsvs == expected.gene_set_tsvs, {k: v for k, v in tsvs.items() if expected.gene_set_tsvs.get(k) != v}

    def named(prefix: str, suffix: str) -> list[str]:
        return sorted(n[len(prefix) : -len(suffix)] for n in members if n.startswith(prefix) and n.endswith(suffix))

    assert named(f"gene_sets_{PREFIX}", ".html") == sorted(expected.gene_set_htmls), "gene set HTML pages differ"
    assert named(f"gene_sets_{PREFIX}", ".png") == sorted(expected.gene_set_heatmaps), "gene set heat maps differ"
    assert named(f"{PREFIX}enplot_", ".png") == sorted(expected.enplots), "enrichment plots differ"
    assert named(f"{PREFIX}gset_rnd_es_dist_", ".png") == sorted(expected.rnd_es_dists), "ES distributions differ"


async def run(gene_sets: list[str]) -> GseaResult:
    return await gsea(
        await fixture(GENE_SET_ANALYSIS + "Condition_treatment_Control_Treated.gct"),
        await fixture(GENE_SET_ANALYSIS + "Condition_genotype_WT_KO.cls"),
        [await fixture(GENE_SET_ANALYSIS + name) for name in gene_sets],
        reference="WT",
        target="KO",
        prefix=PREFIX,
        chip=await fixture(GENE_SET_ANALYSIS + "Mus_musculus.anno.feature_metadata.chip"),
        args="-rnd_seed 10",
    )


@env.task
async def test_gsea_hallmark() -> None:
    # upstream case: "test"
    await check(await run(["mh.all.v2022.1.Mm.symbols.gmt"]), HALLMARK)


@env.task
async def test_gsea_multiple_gmt() -> None:
    # upstream case: "test_multiple_gmt"
    await check(await run(["mh.all.v2022.1.Mm.symbols.gmt", "m5.mpt.v2022.1.Mm.symbols.gmt"]), MULTIPLE)


tests = [test_gsea_hallmark, test_gsea_multiple_gmt]
