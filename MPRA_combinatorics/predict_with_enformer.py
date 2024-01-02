import tensorflow as tf
import tensorflow_hub as hub
from loguru import logger

import kipoiseq
import pandas as pd
import numpy as np

MODEL_PATH = 'https://tfhub.dev/deepmind/enformer/1'
SEQUENCE_LENGTH = 393216



def one_hot_encode(sequence):
  return kipoiseq.transforms.functional.one_hot_dna(sequence).astype(np.float32)

def crop_seq(seq):
    left_length=np.floor((SEQUENCE_LENGTH-len(seq))/2).astype(int)
    right_length=np.ceil((SEQUENCE_LENGTH-len(seq))/2).astype(int)
    return "N"*left_length+seq+"N"*right_length

def seq2x(seqs):
    seqs=list(map(crop_seq,seqs))
    X=np.array(list(map(one_hot_encode,seqs)))
    return X

def one_hot_encode(sequence):
  return kipoiseq.transforms.functional.one_hot_dna(sequence).astype(np.float32)


class Enformer:

  def __init__(self, tfhub_url):
    self._model = hub.load(tfhub_url).model

  def predict_on_batch(self, inputs):
    predictions = self._model.predict_on_batch(inputs)
    return {k: v.numpy() for k, v in predictions.items()}


if __name__ == '__main__':
  assert tf.config.list_physical_devices('GPU')
  model = Enformer(MODEL_PATH)
  i=0
  for chunk in pd.read_csv("Pd1_formatted_sequences/MPRA_seqs.csv",chunksize=4,header=0):
      i=i+1
      logger.info(f"{i*4} sequences processed")
      seqs=list(chunk.loc[:,"sequences"])
      x=seq2x(seqs)
      predictions = model.predict_on_batch(x)['human'][:,[448,449],:].mean(axis=1)
      pd.DataFrame(predictions).to_csv("Pd3_Enformer_predictions/enformer_predictions.csv", mode='a',index=False,header=False)



