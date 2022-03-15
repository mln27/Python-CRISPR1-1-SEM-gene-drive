# CRISPR1-1_gene-drive

The model package is an extended version of a previously developed model (https://github.com/jozinzapletal/Python-CRISPR1-1-SEM-gene-drive). This new version explicitly model the lifecyle of mosquitoes using an existing mosquito model (Deredec A, et al. PNAS 108(43), E874–80). We developed a deterministic and stocastic version of the model. Like the previous model, this model consists of six total components, of which five are required to run the model and output the results. The following components must be located within the same folder in order to execute the model without errors:

• CRISPR1-1_main.py

• CRISPR1_1_model.py

• CRISPR1_1_Equations.py

• CRISPR1-1 Input Parameters.xlsx

• Allele_Plots.py

Model parameters can be modified within the CRISPR1-1 Input Parameters.xlsx Excel file. The file must be saved after any changes to be read in and executed properly by Python. The model can then be executed by running the CRISPR1-1_main.py module. The primary output will consist of allele plots for all the scenarios entered in the CRISPR1-1 Input Parameters.xlsx Excel file. The quantitative results from which these plots were generated will be stored in the file CRISPR1_1 Results.xlsx. These results are overwritten each time the model is run and store the proportion of each allele and population size at each time step for all scenarios.

The CRISPR1_1_Equation_Generation.py module is not required to execute the model, but provides the framework from which the equations within the CRISPR1_1_Equations.py were generated.
