import numpy as np
import tensorflow as tf
import pandas as pd
import io
from tqdm import tqdm 
import zstandard as zstd
import re
from Bio import SeqIO
from itertools import islice
from tensorflow import keras
from sitetack_model_functions_new import predict_functions, predict_ptm_batch_chunk

tf.debugging.set_log_device_placement(True)

#import boto3
#session = boto3.Session(profile_name='petadex-claire')
#s3 = session.client('s3')
#obj = s3.get_object(Bucket='petadex-orf-fastaa', Key='blastnr_pazy.catalytic_orfs.fa')
#body = obj['Body']
#records = SeqIO.parse(io.TextIOWrapper(body, encoding='utf-8'), 'fasta')

#ec2 version
records = SeqIO.parse('blastnr_pazy.catalytic_orfs.fa', 'fasta')
first_10000 = list(islice(records, 10000))

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

def has_only_standard_aa(seq):
    return set(str(seq.seq)).issubset(STANDARD_AA)

first_10000 = [r for r in first_10000 if has_only_standard_aa(r)]

results_list = [predict_ptm_batch_chunk(fn, first_10000) for fn in predict_functions]
predict_results = pd.concat(results_list, ignore_index=True)
predict_results.to_csv('prediction_results_batch_gpu_10000.csv')
