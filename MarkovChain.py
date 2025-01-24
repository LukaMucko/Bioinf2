from Bio import SeqIO
import numpy as np
import argparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

def load_genome(filepath):
    """
    Load a genome sequence from a FASTA file.
    """
    for record in SeqIO.parse(filepath, "fasta"):
        return str(record.seq)

def generate_kmers(sequence, k):
    """
    Generate k-mers from a given sequence.
    """
    return [sequence[i:i+k] for i in range(len(sequence) - k + 1)]

def main():
    parser = argparse.ArgumentParser(description='Process genome sequences and generate k-mers.')
    parser.add_argument('-k', type=int, default=1, required=True, help='Length of k-mers')
    parser.add_argument('-r', '--reference', required=False, help='Path to reference genome FASTA file', default="ecoli/ref.fasta")
    parser.add_argument('-o', '--output', required=False, help='Output folder name', default='samples/')
    parser.add_argument('-n', '--num_genes', required=False, help="Number of genes to sample", default=1, type=int)
    args = parser.parse_args()
    return args

def create_transition_matrix(kmers: list, n_kmers: int):
    """
    Create a Markov chain transition matrix from genome sequence.

    Args:
        kmers: List of k-m
        n_kmers: Number of unique k

    Returns:
        Tuple containing:
        - transition_matrix: Normalized transition probabilities
        - final_state: Index
    """
    # Create and fill transition matrix
    transition_matrix = np.zeros((n_kmers, 4))
    
    # Calculate transitions
    for i in range(len(kmers)-1):
        prev = kmers[i]
        curr = kmers[i+1][-1]
        i_idx = kmer_dict[prev]
        j_idx = base_encoder[curr]
        transition_matrix[i_idx, j_idx] += 1

    # Find final states
    final_state = np.where(np.sum(transition_matrix, axis=1) == 0.)[0]
    mask = np.ones(transition_matrix.shape[0], dtype=bool)
    mask[final_state] = False

    # Normalize transition probabilities
    if final_state.size != 0:
        transition_matrix[mask] = transition_matrix[mask] / np.sum(transition_matrix[mask], axis=1)[:, None]
    else:
        transition_matrix = transition_matrix / np.sum(transition_matrix, axis=1)[:, None]

    starting_state = kmer_dict[kmers[0]]

    return transition_matrix, starting_state, None if final_state.size == 0 else final_state[0]

def sample_sequence(transition_matrix, reverse_kmer_dict, starting_state, sequence_length, final_state = None):
    """
    Generate a random sequence from the Markov chain.
    
    Parameters:
    - transition_matrix: The normalized transition matrix (rows sum to 1).
    - reverse_kmer_dict: Dictionary to decode integer indices back to k-mers.
    - starting_state: The index of the starting k-mer.
    - sequence_length: The desired length of the generated sequence.
    
    Returns:
    - A string representing the generated sequence.
    """
    kmer = reverse_kmer_dict[starting_state]
    k = len(kmer)

    sequence = bytearray(sequence_length)
    sequence[:k] = kmer.encode()
    
    current_state = starting_state
    pos = k

    transition_matrix = np.cumsum(transition_matrix, axis=1)
    while pos < sequence_length:
        next_base_idx = np.searchsorted(transition_matrix[current_state], np.random.rand())
        next_base = reverse_base_decoder[int(next_base_idx)]

        sequence[pos] = ord(next_base)

        current_state = kmer_dict[sequence[pos-k+1:pos+1].decode()]

        if next_base_idx == final_state:
            print("Reached final_state")
            return sequence[:pos+1].decode()

        pos += 1
        
    return sequence.decode()

def save_to_fasta(sequence, output_folder, k, n, sequence_id=None):
    """
    Save generated sequence to a FASTA file with standardized naming.
    
    Args:
        sequence: The sequence to save
        output_folder: Folder path to save the file
        k: k-mer size used
        n: Sample number
        sequence_id: Optional sequence identifier
    """

    if sequence_id is None:
        sequence_id = f"sample_{k}_{n}_{len(sequence)}"
    
    filename = f"k{k}_n{n}_{len(sequence)}.fasta"
    output_path = os.path.join(output_folder, filename)
    
    # Create output directory if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)
    
    # Create and save SeqRecord
    record = SeqRecord(Seq(sequence),
                      id=sequence_id,
                      description=f"Generated sequence with k={k}")
    SeqIO.write(record, output_path, "fasta")

if __name__ == "__main__":
    args = main()

    # Load genome and generate k-mers
    reference = load_genome(args.reference)
    k = args.k
    kmers = generate_kmers(reference, k)
    unique_kmers = set(kmers)
    n_kmers = len(unique_kmers)
    print(f"Genome loaded: {len(reference)} bases")
    print(f"Number of {k}-mers: {len(kmers)}")
    print(f"Number of unique {k}-mers: {n_kmers}")

    # Create mappings
    kmer_dict = {kmer: idx for idx, kmer in enumerate(unique_kmers)}
    reverse_kmer_dict = {idx: kmer for kmer, idx in kmer_dict.items()}
    base_encoder = {"A": 0, "C": 1, "G": 2, "T": 3}
    reverse_base_decoder = {0: "A", 1: "C", 2: "G", 3: "T"}

    transition_matrix, starting_state, final_state = create_transition_matrix(kmers, n_kmers)
    print(f"Transition matrix shape: {transition_matrix.shape}")
    print(f"Data usage of transition matrix: {transition_matrix.nbytes / (1024 * 1024 * 1024):.2f} GB")
    del kmers
    
    for i in range(args.num_genes):
        sample = sample_sequence(transition_matrix, reverse_kmer_dict, starting_state, len(reference), final_state)
        save_to_fasta(sample, args.output, k, i)