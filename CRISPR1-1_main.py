# This module will execute the model and generate the associated plots
#
# Developed by Josef Zapletal (jozinzapletal@gmail.com)

import os.path
from os import path
import CRISPR1_1_Equation_Generation


file_name = 'CRISPR1_1_Equations.py'
if path.exists(file_name):
    pass
else:
    print('CRISPR1_1_Equations.py not found. Generating equations...')
    CRISPR1_1_Equation_Generation.CRISPR1_1_Eq_gen()
    print('Equation generation complete. Running model...')

from CRISPR1_1_model import CRISPR_1_1_Progression
from Allele_Plots import generate_plots

CRISPR_1_1_Progression()
generate_plots()
