import numpy as np

def generate_random_seq(length):
    return ''.join(np.random.choice(["A", "C", "G", "T"], length))

def generate_random_seqs(num_seqs, length):
    return [generate_random_seq(length) for _ in range(num_seqs)]