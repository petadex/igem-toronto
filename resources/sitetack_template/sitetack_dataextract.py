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
from sitetack_model_functions import predict_ptm, predict_n_glycosylation, predict_ubiquitination, predict_o_glycosylation
#only for local testing, must sso to AWS first, generally sucks ass and i don't understand how it works
import boto3
session = boto3.Session(profile_name='petadex-claire')
s3 = session.client('s3')
obj = s3.get_object(Bucket='petadex-orf-fastaa', Key='blastnr_pazy.catalytic_orfs.fa')
body = obj['Body']

records = SeqIO.parse('blastnr_pazy.catalytic_orfs.fa', 'fasta')
first_100 = list(islice(records, 100))

#fix this bro
predict_results = pd.DataFrame()
results_list = []
for record in first_100:
    seq_name = re.search(r'\|([^|]+)\|', record.description).group(1)
    results_list.append(predict_n_glycosylation(record.seq, seq_name=seq_name))
    results_list.append(predict_ubiquitination(record.seq, seq_name=seq_name))
    results_list.append(predict_o_glycosylation(record.seq, seq_name=seq_name))
    results_list.append(predict_phosphorylation_st(record.seq, seq_name=seq_name))
    results_list.append(predict_phosphorylation_y(record.seq, seq_name=seq_name))
    results_list.append(predict_sumoylation(record.seq, seq_name=seq_name))
    results_list.append(predict_acetylation(record.seq, seq_name=seq_name))
    results_list.append(predict_methylation_r(record.seq, seq_name=seq_name))
    results_list.append(predict_methylation_k(record.seq, seq_name=seq_name))
    results_list.append(predict_pyroglutamylation(record.seq, seq_name=seq_name))
    results_list.append(predict_palmitoylation(record.seq, seq_name=seq_name))
    results_list.append(predict_hydroxylation_p(record.seq, seq_name=seq_name))
    results_list.append(predict_hydroxylation_k(record.seq, seq_name=seq_name))

predict_results = pd.concat(results_list, ignore_index=True)
predict_results.to_csv('prediction_results_2.csv')
