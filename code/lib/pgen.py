import os.path

import olga
import olga.load_model as load_model
import olga.generation_probability as generation_probability

main_folder = os.path.dirname(olga.__file__)
model_folder = os.path.join(main_folder, 'default_models', 'human_T_beta')
params_file_name = os.path.join(model_folder,'model_params.txt')
marginals_file_name = os.path.join(model_folder,'model_marginals.txt')
V_anchor_pos_file = os.path.join(model_folder,'V_gene_CDR3_anchors.csv')
J_anchor_pos_file = os.path.join(model_folder,'J_gene_CDR3_anchors.csv')
genomic_data = load_model.GenomicDataVDJ()
genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
generative_model = load_model.GenerativeModelVDJ()
generative_model.load_and_process_igor_model(marginals_file_name)
pgen_model = generation_probability.GenerationProbabilityVDJ(generative_model, genomic_data)

def calc_nt_pgen(row):
    "Calc probability of generation under the default OLGA human TCR CDRbeta model"
    ntseq = str(row['nt'])
    aaseq = str(row['aa'])
    if aaseq == 'nan':
        return np.nan
    else:
        return pgen_model.compute_nt_CDR3_pgen(ntseq[81-3*len(aaseq):81])
