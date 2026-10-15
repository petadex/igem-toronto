import numpy as np
import tensorflow as tf
import pandas as pd
import io
from tensorflow import keras
from sitetack_general_functions import convert_to_onehot, tensor_encoding, get_kmer
from functools import partial

def predict_ptm (sequence, seq_name, labeled_model, unlabeled_model, alphabet_with_labels, alphabet_without_labels, model, PR, PR2, k=53):
    folder="results"  
    PR_sites=[i + 1 for i, char in enumerate(sequence) if char == PR]
    PR_kmers= [get_kmer(sequence, s, k=k) for s in PR_sites]
    if PR2=='':
        PR2_sites=[]
        PR2_kmers=[]
    else:
        PR2_sites=[i + 1 for i, char in enumerate(sequence) if char == PR2]
        PR2_kmers= [get_kmer(sequence, s, k=k) for s in PR2_sites]

    kmers=PR_kmers+PR2_kmers
    sites=PR_sites+PR2_sites

    tensor1 = tensor_encoding(kmers, 23, 'emb', alphabet_without_labels, k=k)
    tensor2 = tensor_encoding(kmers, 23, 'emb', alphabet_with_labels, k=k)

    nl_y_pred = unlabeled_model.predict(tensor1)[:,0]
    l_y_pred = labeled_model.predict(tensor2)[:,0]

    dict={"Name": seq_name,"PTM Type": model, "Site": sites, "No labels model":nl_y_pred, "With PTM labels model":l_y_pred}
    df = pd.DataFrame(dict)
    return (df)

predict_n_glycosylation = partial(predict_ptm, model = f"N-glycosylation(N)", alphabet_with_labels = f"ARNDCEQGHILKMFPSTWYV@-UXZB", alphabet_without_labels = f"ARNDCEQGHILKMFPSTWYV-UXZB", PR = "N", PR2 = "", labeled_model = tf.keras.models.load_model(f"./sitetack/models/N-Linked glycosylation (N)/All organism/N-linked_glycosylation_N_emb_CNN_with_labels_4.h5"), unlabeled_model = tf.keras.models.load_model(f"./sitetack/models/N-Linked glycosylation (N)/All organism/N-linked_glycosylation_N_emb_CNN_no_labels_4.h5"))
predict_ubiquitination = partial(predict_ptm, model = f"Ubiquitination(K)", alphabet_with_labels = f"ARNDCEQGHILKMFPSTWYV@-UX", alphabet_without_labels = f"ARNDCEQGHILKMFPSTWYV-UX", PR = "K", PR2 = "", labeled_model = tf.keras.models.load_model(f"./sitetack/models/Ubiquitination (K)/All organism/Ubiquitin_K_emb_CNN_with_labels_7.h5"), unlabeled_model= tf.keras.models.load_model(f"./sitetack/models/Ubiquitination (K)/All organism/Ubiquitin_K_emb_CNN_no_labels_6.h5"))
predict_o_glycosylation = partial(predict_ptm, model = f"O-glycosylation(S,T)", alphabet_with_labels = f"ARNDCEQGHILKMFPSTWYV@&-UX", alphabet_without_labels = f"ARNDCEQGHILKMFPSTWYV-UX", PR = "S", PR2 = "T", labeled_model = tf.keras.models.load_model(f"./sitetack/models/O-Linked glycosylation (S,T)/All organism/O-linked glycosylation_S_T_emb_CNN_with_labels_10.h5"), unlabeled_model = tf.keras.models.load_model(f"./sitetack/models/O-Linked glycosylation (S,T)/All organism/O-linked glycosylation_S_T_emb_CNN_no_labels_6.h5"))
predict_phosphorylation_st = partial(predict_ptm, model = f"Phosphorylation(S,T)", alphabet_with_labels = f"ARNDCEQGHILKMFPSTWYV@&-UXBZ", alphabet_without_labels = f"ARNDCEQGHILKMFPSTWYV-UXBZ", PR = "S", PR2 = "T", labeled_model = tf.keras.models.load_model(f"./sitetack/models/Phosphorylation (S,T)/All organism/Phosphorylation_S_T_emb_CNN_with_labels_2.h5"), unlabeled_model = tf.keras.models.load_model(f"./sitetack/models/Phosphorylation (S,T)/All organism/Phosphorylation_S_T_emb_CNN_no_labels_6.h5"))
predict_phosphorylation_y = partial(predict_ptm, model = f"Phosphorylation(Y)", alphabet_with_labels = f"ARNDCEQGHILKMFPSTWYV@-UXBZ", alphabet_without_labels = f"ARNDCEQGHILKMFPSTWYV-UXBZ", PR = "Y", PR2 = "", labeled_model = tf.keras.models.load_model(f"./sitetack/models/Phosphorylation (Y)/All organism/Phosphorylation_Y_emb_CNN_with_labels_8.h5"), unlabeled_model = tf.keras.models.load_model(f"./sitetack/models/Phosphorylation (Y)/All organism/Phosphorylation_Y_emb_CNN_no_labels_2.h5"))
predict_sumoylation = partial(predict_ptm, model = f"SUMOylation(K)", alphabet_with_labels = f"ARNDCEQGHILKMFPSTWYV@-UX", alphabet_without_labels = f"ARNDCEQGHILKMFPSTWYV-UX", PR = "K", PR2 = "", labeled_model = tf.keras.models.load_model(f"./sitetack/models/SUMOylation (K)/All organism/SUMOylation_K_emb_CNN_with_labels_8.h5"), unlabeled_model = tf.keras.models.load_model(f"./sitetack/models/SUMOylation (K)/All organism/SUMOylation_K_emb_CNN_no_labels_8.h5"))
predict_acetylation = partial(predict_ptm, model = f"Acetylation(K)", alphabet_with_labels = f"ARNDCEQGHILKMFPSTWYV@-UX", alphabet_without_labels = f"ARNDCEQGHILKMFPSTWYV-UX", PR = "K", PR2 = "", labeled_model = tf.keras.models.load_model(f"./sitetack/models/Acetylation (K)/All organism/N6-acetyllysine_K_emb_CNN_with_labels_1.h5"), unlabeled_model = tf.keras.models.load_model(f"./sitetack/models/Acetylation (K)/All organism/N6-acetyllysine_K_emb_CNN_no_labels_10.h5"))
predict_methylation_r = partial(predict_ptm, model = f"Methylation(R)", alphabet_with_labels = f"ARNDCEQGHILKMFPSTWYV@-UX", alphabet_without_labels = f"ARNDCEQGHILKMFPSTWYV-UX", PR = "R", PR2 = "", labeled_model = tf.keras.models.load_model(f"./sitetack/models/Methylation (R)/All organism/Methylation_R_emb_CNN_with_labels_6.h5"), unlabeled_model = tf.keras.models.load_model(f"./sitetack/models/Methylation (R)/All organism/Methylation_R_emb_CNN_no_labels_6.h5"))
predict_methylation_k = partial(predict_ptm, model = f"Methylation(K)", alphabet_with_labels = f"ARNDCEQGHILKMFPSTWYV@-UX", alphabet_without_labels = f"ARNDCEQGHILKMFPSTWYV-UX", PR = "K", PR2 = "", labeled_model = tf.keras.models.load_model(f"./sitetack/models/Methylation (K)/All organism/Methylation_K_emb_CNN_with_labels_5.h5"), unlabeled_model = tf.keras.models.load_model(f"./sitetack/models/Methylation (K)/All organism/Methylation_K_emb_CNN_no_labels_1.h5"))
predict_pyroglutamylation = partial(predict_ptm, model = f"Pyroglutamylation(Q)", alphabet_with_labels = f"ARNDCEQGHILKMFPSTWYV@-UXBZ", alphabet_without_labels = f"ARNDCEQGHILKMFPSTWYV-UXBZ", PR = "Q", PR2 = "", labeled_model = tf.keras.models.load_model(f"./sitetack/models/Pyroglutamylation (Q)/All organism/Pyrrolidone-carboxylic-acid_Q_emb_CNN_with_labels_3.h5"), unlabeled_model = tf.keras.models.load_model(f"./sitetack/models/Pyroglutamylation (Q)/All organism/Pyrrolidone-carboxylic-acid_Q_emb_CNN_no_labels_1.h5"))
predict_palmitoylation = partial(predict_ptm, model = f"Palmitoylation(C)", alphabet_with_labels = f"ARNDCEQGHILKMFPSTWYV@-UX", alphabet_without_labels = f"ARNDCEQGHILKMFPSTWYV-UX", PR = "C", PR2 = "", labeled_model = tf.keras.models.load_model(f"./sitetack/models/Palmitoylation (C)/All organism/S-Palmitoylation_C_emb_CNN_with_labels_4.h5"), unlabeled_model = tf.keras.models.load_model(f"./sitetack/models/Palmitoylation (C)/All organism/S-Palmitoylation_C_emb_CNN_no_labels_9.h5"))
predict_hydroxylation_p = partial(predict_ptm, model = f"Hydroxylation(P)", alphabet_with_labels = f"ARNDCEQGHILKMFPSTWYV@-UXZB", alphabet_without_labels = f"ARNDCEQGHILKMFPSTWYV-UXZB", PR = "P", PR2 = "", labeled_model = tf.keras.models.load_model(f"./sitetack/models/Hydroxylation (P)/All organism/Hydroxyproline_P_emb_CNN_with_labels_1.h5"), unlabeled_model = tf.keras.models.load_model(f"./sitetack/models/Hydroxylation (P)/All organism/Hydroxyproline_P_emb_CNN_no_labels_7.h5"))
predict_hydroxylation_k = partial(predict_ptm, model = f"Hydroxylation(K)", alphabet_with_labels = f"ARNDCEQGHILKMFPSTWYV@-UZ", alphabet_without_labels = f"ARNDCEQGHILKMFPSTWYV-UZ", PR = "K", PR2 = "", labeled_model = tf.keras.models.load_model(f"./sitetack/models/Hydroxylation (K)/All organism/Hydroxyproline_K_emb_CNN_with_labels_5.h5"), unlabeled_model = tf.keras.models.load_model(f"sitetack/sitetack/models/Hydroxylation (K)/All organism/Hydroxyproline_K_emb_CNN_no_labels_2.h5"))