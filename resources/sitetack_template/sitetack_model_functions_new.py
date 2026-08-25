import numpy as np
import tensorflow as tf
import pandas as pd
import io
from tensorflow import keras
from sitetack_general_functions import convert_to_onehot, tensor_encoding, get_kmer
from functools import partial
import re

def predict_ptm_batch(predict_fn, records, k_default=53):
    kw = predict_fn.keywords

    labeled_model = kw['labeled_model']
    unlabeled_model = kw['unlabeled_model']
    alphabet_with_labels = kw['alphabet_with_labels']
    alphabet_without_labels = kw['alphabet_without_labels']
    model = kw['model']
    PR = kw['PR']
    PR2 = kw['PR2']
    k = kw.get('k', k_default)

    all_kmers = []
    all_sites = []
    all_names = []

    for record in records:
        sequence = str(record.seq)
        seq_name = re.search(r'\|([^|]+)\|', record.description).group(1)

        PR_sites = [i + 1 for i, char in enumerate(sequence) if char == PR]
        PR_kmers = [get_kmer(sequence, s, k=k) for s in PR_sites]

        if PR2 == '':
            PR2_sites, PR2_kmers = [], []
        else:
            PR2_sites = [i + 1 for i, char in enumerate(sequence) if char == PR2]
            PR2_kmers = [get_kmer(sequence, s, k=k) for s in PR2_sites]

        kmers = PR_kmers + PR2_kmers
        sites = PR_sites + PR2_sites

        all_kmers.extend(kmers)
        all_sites.extend(sites)
        all_names.extend([seq_name] * len(kmers))

    # one batched encode + one batched predict call for ALL sequences
    tensor1 = tensor_encoding(all_kmers, 23, 'emb', alphabet_without_labels, k=k)
    tensor2 = tensor_encoding(all_kmers, 23, 'emb', alphabet_with_labels, k=k)

    nl_y_pred = unlabeled_model.predict(tensor1)[:, 0]
    l_y_pred = labeled_model.predict(tensor2)[:, 0]

    df = pd.DataFrame({
        "Name": all_names,
        "PTM Type": model,
        "Site": all_sites,
        "No labels model": nl_y_pred,
        "With PTM labels model": l_y_pred,
    })
    return df

predict_n_glycosylation = partial(predict_ptm_batch, model = f"N-glycosylation(N)", alphabet_with_labels = f"ARNDCEQGHILKMFPSTWYV@-UXZB", alphabet_without_labels = f"ARNDCEQGHILKMFPSTWYV-UXZB", PR = "N", PR2 = "", labeled_model = tf.keras.models.load_model(f"./sitetack/models/N-Linked glycosylation (N)/All organism/N-linked_glycosylation_N_emb_CNN_with_labels_4.h5"), unlabeled_model = tf.keras.models.load_model(f"./sitetack/models/N-Linked glycosylation (N)/All organism/N-linked_glycosylation_N_emb_CNN_no_labels_4.h5"))
predict_ubiquitination = partial(predict_ptm_batch, model = f"Ubiquitination(K)", alphabet_with_labels = f"ARNDCEQGHILKMFPSTWYV@-UX", alphabet_without_labels = f"ARNDCEQGHILKMFPSTWYV-UX", PR = "K", PR2 = "", labeled_model = tf.keras.models.load_model(f"./sitetack/models/Ubiquitination (K)/All organism/Ubiquitin_K_emb_CNN_with_labels_7.h5"), unlabeled_model= tf.keras.models.load_model(f"./sitetack/models/Ubiquitination (K)/All organism/Ubiquitin_K_emb_CNN_no_labels_6.h5"))
predict_o_glycosylation = partial(predict_ptm_batch, model = f"O-glycosylation(S,T)", alphabet_with_labels = f"ARNDCEQGHILKMFPSTWYV@&-UX", alphabet_without_labels = f"ARNDCEQGHILKMFPSTWYV-UX", PR = "S", PR2 = "T", labeled_model = tf.keras.models.load_model(f"./sitetack/models/O-Linked glycosylation (S,T)/All organism/O-linked glycosylation_S_T_emb_CNN_with_labels_10.h5"), unlabeled_model = tf.keras.models.load_model(f"./sitetack/models/O-Linked glycosylation (S,T)/All organism/O-linked glycosylation_S_T_emb_CNN_no_labels_6.h5"))
predict_phosphorylation_st = partial(predict_ptm_batch, model = f"Phosphorylation(S,T)", alphabet_with_labels = f"ARNDCEQGHILKMFPSTWYV@&-UXBZ", alphabet_without_labels = f"ARNDCEQGHILKMFPSTWYV-UXBZ", PR = "S", PR2 = "T", labeled_model = tf.keras.models.load_model(f"./sitetack/models/Phosphorylation (S,T)/All organism/Phosphorylation_S_T_emb_CNN_with_labels_2.h5"), unlabeled_model = tf.keras.models.load_model(f"./sitetack/models/Phosphorylation (S,T)/All organism/Phosphorylation_S_T_emb_CNN_no_labels_6.h5"))
predict_phosphorylation_y = partial(predict_ptm_batch, model = f"Phosphorylation(Y)", alphabet_with_labels = f"ARNDCEQGHILKMFPSTWYV@-UXBZ", alphabet_without_labels = f"ARNDCEQGHILKMFPSTWYV-UXBZ", PR = "Y", PR2 = "", labeled_model = tf.keras.models.load_model(f"./sitetack/models/Phosphorylation (Y)/All organism/Phosphorylation_Y_emb_CNN_with_labels_8.h5"), unlabeled_model = tf.keras.models.load_model(f"./sitetack/models/Phosphorylation (Y)/All organism/Phosphorylation_Y_emb_CNN_no_labels_2.h5"))
predict_sumoylation = partial(predict_ptm_batch, model = f"SUMOylation(K)", alphabet_with_labels = f"ARNDCEQGHILKMFPSTWYV@-UX", alphabet_without_labels = f"ARNDCEQGHILKMFPSTWYV-UX", PR = "K", PR2 = "", labeled_model = tf.keras.models.load_model(f"./sitetack/models/SUMOylation (K)/All organism/SUMOylation_K_emb_CNN_with_labels_8.h5"), unlabeled_model = tf.keras.models.load_model(f"./sitetack/models/SUMOylation (K)/All organism/SUMOylation_K_emb_CNN_no_labels_8.h5"))
predict_acetylation = partial(predict_ptm_batch, model = f"Acetylation(K)", alphabet_with_labels = f"ARNDCEQGHILKMFPSTWYV@-UX", alphabet_without_labels = f"ARNDCEQGHILKMFPSTWYV-UX", PR = "K", PR2 = "", labeled_model = tf.keras.models.load_model(f"./sitetack/models/Acetylation (K)/All organism/N6-acetyllysine_K_emb_CNN_with_labels_1.h5"), unlabeled_model = tf.keras.models.load_model(f"./sitetack/models/Acetylation (K)/All organism/N6-acetyllysine_K_emb_CNN_no_labels_10.h5"))
predict_methylation_r = partial(predict_ptm_batch, model = f"Methylation(R)", alphabet_with_labels = f"ARNDCEQGHILKMFPSTWYV@-UX", alphabet_without_labels = f"ARNDCEQGHILKMFPSTWYV-UX", PR = "R", PR2 = "", labeled_model = tf.keras.models.load_model(f"./sitetack/models/Methylation (R)/All organism/Methylation_R_emb_CNN_with_labels_6.h5"), unlabeled_model = tf.keras.models.load_model(f"./sitetack/models/Methylation (R)/All organism/Methylation_R_emb_CNN_no_labels_6.h5"))
predict_methylation_k = partial(predict_ptm_batch, model = f"Methylation(K)", alphabet_with_labels = f"ARNDCEQGHILKMFPSTWYV@-UX", alphabet_without_labels = f"ARNDCEQGHILKMFPSTWYV-UX", PR = "K", PR2 = "", labeled_model = tf.keras.models.load_model(f"./sitetack/models/Methylation (K)/All organism/Methylation_K_emb_CNN_with_labels_5.h5"), unlabeled_model = tf.keras.models.load_model(f"./sitetack/models/Methylation (K)/All organism/Methylation_K_emb_CNN_no_labels_1.h5"))
predict_pyroglutamylation = partial(predict_ptm_batch, model = f"Pyroglutamylation(Q)", alphabet_with_labels = f"ARNDCEQGHILKMFPSTWYV@-UXBZ", alphabet_without_labels = f"ARNDCEQGHILKMFPSTWYV-UXBZ", PR = "Q", PR2 = "", labeled_model = tf.keras.models.load_model(f"./sitetack/models/Pyroglutamylation (Q)/All organism/Pyrrolidone-carboxylic-acid_Q_emb_CNN_with_labels_3.h5"), unlabeled_model = tf.keras.models.load_model(f"./sitetack/models/Pyroglutamylation (Q)/All organism/Pyrrolidone-carboxylic-acid_Q_emb_CNN_no_labels_1.h5"))
predict_palmitoylation = partial(predict_ptm_batch, model = f"Palmitoylation(C)", alphabet_with_labels = f"ARNDCEQGHILKMFPSTWYV@-UX", alphabet_without_labels = f"ARNDCEQGHILKMFPSTWYV-UX", PR = "C", PR2 = "", labeled_model = tf.keras.models.load_model(f"./sitetack/models/Palmitoylation (C)/All organism/S-Palmitoylation_C_emb_CNN_with_labels_4.h5"), unlabeled_model = tf.keras.models.load_model(f"./sitetack/models/Palmitoylation (C)/All organism/S-Palmitoylation_C_emb_CNN_no_labels_9.h5"))
predict_hydroxylation_p = partial(predict_ptm_batch, model = f"Hydroxylation(P)", alphabet_with_labels = f"ARNDCEQGHILKMFPSTWYV@-UXZB", alphabet_without_labels = f"ARNDCEQGHILKMFPSTWYV-UXZB", PR = "P", PR2 = "", labeled_model = tf.keras.models.load_model(f"./sitetack/models/Hydroxylation (P)/All organism/Hydroxyproline_P_emb_CNN_with_labels_1.h5"), unlabeled_model = tf.keras.models.load_model(f"./sitetack/models/Hydroxylation (P)/All organism/Hydroxyproline_P_emb_CNN_no_labels_7.h5"))
predict_hydroxylation_k = partial(predict_ptm_batch, model = f"Hydroxylation(K)", alphabet_with_labels = f"ARNDCEQGHILKMFPSTWYV@-UZ", alphabet_without_labels = f"ARNDCEQGHILKMFPSTWYV-UZ", PR = "K", PR2 = "", labeled_model = tf.keras.models.load_model(f"./sitetack/models/Hydroxylation (K)/All organism/Hydroxyproline_K_emb_CNN_with_labels_5.h5"), unlabeled_model = tf.keras.models.load_model(f"./sitetack/models/Hydroxylation (K)/All organism/Hydroxyproline_K_emb_CNN_no_labels_2.h5"))

predict_functions = [
    predict_n_glycosylation, predict_ubiquitination, predict_o_glycosylation,
    predict_phosphorylation_st, predict_phosphorylation_y, predict_sumoylation,
    predict_acetylation, predict_methylation_r, predict_methylation_k,
    predict_pyroglutamylation, predict_palmitoylation,
    predict_hydroxylation_p, predict_hydroxylation_k,
]

def predict_ptm_batch_chunk(predict_fn, records, k_default=53, chunk_size=1000):
    all_dfs = []
    for i in range(0, len(records), chunk_size):
        chunk = records[i:i + chunk_size]
        df = predict_ptm_batch(predict_fn, chunk, k_default=k_default)
        all_dfs.append(df)
    return pd.concat(all_dfs, ignore_index=True)