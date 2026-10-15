import numpy as np
import tensorflow as tf
import pandas as pd
from tqdm import tqdm

#define helper functions
from tensorflow import keras

def convert_to_onehot(data, alphabet):
    char_to_int = dict((c, i) for i, c in enumerate(alphabet))
    return [char_to_int[char] for char in data]

def tensor_encoding(x_data, depth, type, alphabet, k=53):
    indices = []
    for i in range(len(x_data)):
        indices.append(convert_to_onehot(x_data[i], alphabet))
        if len(convert_to_onehot(x_data[i], alphabet)) != k:
            print(x_data[i], "Length off")
    array = np.stack(indices, axis=0)
    if type == 'emb':
        return array
    t2 = []
    for i in tqdm(range(len(indices))):
        t1 = tf.one_hot(indices[i], depth)
        t2.append(t1)
    return t2

def get_kmer(seq, location, k=53):
    half = int((k - 1) / 2)
    if location > len(seq):
        print(f"Site outside of seq bounds, site: {location}, sequence length: {len(seq)}")
        return ''
    elif location <= half:
        if location > len(seq) - half:
            gap = "-" * (half - location + 1)
            gap2 = "-" * int(half - (len(seq) - location))
            kmer = seq[0:int(location + half)]
            kmer = gap + kmer + gap2
        else:
            gap = "-" * (half - location + 1)
            kmer = seq[0:int(location + half)]
            kmer = gap + kmer
    elif location > len(seq) - half:
        gap = "-" * int(half - (len(seq) - location))
        kmer = seq[int(location - half - 1): len(seq)]
        kmer = kmer + gap
    else:
        kmer = seq[int(location - half - 1): int(location + half)]
    assert len(kmer) == k
    return kmer