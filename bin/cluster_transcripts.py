import numpy as np
import sys

def to_sequences(fn):
    sample_names, sequences, lengths = [], [], []
    count = 0
    sequence = ''
    if fn.split('.')[-1] == 'gz':
        with gzip.open(fn, 'r') as f:
            for line in f:
                line = line.decode().strip()
                if (line) != '':
                    if line[0] == '>':
                        sample_names.append(line[1:])
                        if count > 0:
                            lengths.append(len(sequence))
                            sequences.append(sequence)
                            sequence = ''
                    else:
                        sequence += line
                        count += 1
    else:
        with open(fn, 'r') as f:
            for line in f:
                line = line.strip()
                if (line) != '':
                    if line[0] == '>':
                        sample_names.append(line[1:])
                        if count > 0:
                            lengths.append(len(sequence))
                            sequences.append(sequence)
                            sequence = ''
                    else:
                        sequence += line
                        count += 1
    if count > 0:
        lengths.append(len(sequence))
        sequences.append(sequence)
        sequence = ''
    out = {}
    out['sample'] = sample_names
    out['sequence'] = sequences
    out['sequence_length'] = lengths
    return out


def main():
    fn_fasta = sys.argv[1]
    fn_pep = sys.argv[2]
    fn_cluster = sys.argv[3]
    output_fasta = sys.argv[4]

    seq = to_sequences(fn_fasta)
    pep = to_sequences(fn_pep)
    seq_lengths = {x:y for x, y in zip(seq['sample'], seq['sequence_length'])}
    sequences = {x:y for x, y in zip(seq['sample'], seq['sequence'])}
    pep['sample'] = [x.split(' ')[0].split('.')[0].replace('_pilon_pilon', '_pilon') for x in pep['sample']]
    pep_lengths = {x:y for x, y in zip(pep['sample'], pep['sequence_length'])}
    final_clusters = []
    final_seqs = []

    collapsed = {}
    with open(fn_cluster, 'r') as f:
        for line in f:
            cluster = np.array(line.strip().split('\t'))
            seq_length = [seq_lengths[x] for x in cluster]
            pep_length = [pep_lengths[x] if x in pep_lengths else 0 for x in cluster]
            maxpeplen = np.argwhere(np.array(pep_length) == max(pep_length)).flatten()
            order = maxpeplen[np.argsort(-np.array(seq_length)[maxpeplen])]
            ind = order[0]
            clname = cluster[ind]
            longest_transcript = sequences[clname]
            final_clusters.append(clname)
            final_seqs.append(longest_transcript)

    reorder = np.argsort([int(x.split('_')[1]) for x in final_clusters])
    final_clusters = np.array(final_clusters)[reorder]
    final_seqs = np.array(final_seqs)[reorder]

    with open(output_fasta, 'w') as w:
        for x, y in zip(final_clusters, final_seqs):
            w.write('>'+x+'\n'+y+'\n')

if __name__ == "__main__":
    main()
