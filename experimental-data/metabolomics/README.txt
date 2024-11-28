all_group--------------------------------------All sample information
├── ALL_sample_data_raw.xlsx-------------------Metabolite information table for all samples (The data before CV filtering)
├── ALL_sample_data.xlsx-----------------------Metabolite information table for all samples
├── sample_info.xlsx---------------------------Sample grouping information
├── ALL_sample_cor.xlsx------------------------Table of Pearson correlation coefficients for all samples
├── untargeted_metabolites.xlsx----------------List of metabolites for non-target detection
├── untargeted_metabolites_plot.*--------------Scatter plot of mixed samples with non-target detection
└── README.txt---------------------------------Document content introduction
Table description:

ALL_sample_data.xlsx
  Index：MWDB ID of metabolites (* indicates that the metabolite is an isomer of another substance with the same number.)
  Compounds：The names of metabolites
  Class I：Primary classification of metabolites
  Class II：Secondary classification of metabolites
  Q1 (Da)：The molecular weight of the precursor ions after the addition of ions by electrospray ion source
  Ionization model：The electrospray ion source (ESI) mode (either positive or negative) used to detect the metabolite
  Formula：The chemical formula of metabolites
  Level：Substance identification level (1：The substance identification is based on parameters such as Standard, Q1, Q3, RT, MS2, etc; 2：The substance identification is based on parameters such as Standard, Q1, Q3, RT or public library; 3：The substance identification is based on MetDNA or predicted library.)
  cpd_ID：ID information of metabolites in KEGG database
  CAS：A unique numerical identification number for a substance
  kegg_map：KEGG database pathway number
  Other columns：relative content of samples

untargeted_metabolites.xlsx
  Index: MetWare ID; Mass: Relative molecular mass of metabolite, mass-to-charge ratio; RT (min): Retention time; Precursor (Da): mass of precursor ion; Score: Qualitative score; Level: Qualitative level.

ALL_sample_cor.xlsx
  Rows and columns are sample names; the values in the table are Pearson correlation coefficients between samples
