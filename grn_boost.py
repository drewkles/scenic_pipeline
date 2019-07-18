# GRNBoost
# grn_boost.py

import subprocess
import pandas as pd
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

if __name__ == '__main__':

    print("Running GRNboost...")

    in_file = 'int/1.1_exprMatrix_filtered_t.txt'
    tf_file = 'int/1.1_inputTFs.txt'
    out_file = 'int/1.2_grnoutput.txt'

    ex_matrix = pd.read_csv(in_file, sep='\t')
    tf_names = load_tf_names(tf_file)

    network = grnboost2(expression_data=ex_matrix, tf_names=tf_names)
    network.to_csv(out_file, sep='\t', index=False, header=False)

    print("...GRNboost complete")

