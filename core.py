import os
import numpy as np
import pandas as pd
import pysam
from collections import defaultdict
import itertools

from scipy.stats import poisson, multinomial

def gather_input(bam_filepath, transcripts_gtf_filepath):
    rpf_counts = list()

    # If BAM index does not exist, create it
    if not os.path.isfile(bam_filepath + '.bai'):
        pysam.index(bam_filepath)
    bam = pysam.AlignmentFile(bam_filepath, 'rb')
    with open(transcripts_gtf_filepath) as transcripts_gtf:
        for line in transcripts_gtf:
            if line.startswith('#'):
                continue
            record = line.strip().split('\t')
            if record[2] != 'transcript':
                continue
            _counts = np.array(bam.count_coverage(contig=record[0], start=int(record[3]), stop=int(record[4])))
            rpf_counts.append(np.sum(_counts, axis=0))
    return rpf_counts


def get_triplets(rpf_counts, F_reading_frame=1, triplet_size=3):
    """
    Given a list of RPF counts for a transcript, returns the triplets (X) based on
    the given reading frame.

    Assumes the ``len(rpf_counts)`` is divisible by 3.

    :param rpf_counts: list<int> RPF counts for each base in a particular fragment, over the length of the fragment
    :param F_reading_frame: int Offset of the reading frame, must be in {1, 2, 3}
    :param triplet_size: int Size of a triplet (nearly always 3)
    :return: list<list<int>> Triplets based on the reading frame
    """
    _triplets = [rpf_counts[i:i+triplet_size] for i in range(F_reading_frame - 1, len(rpf_counts), triplet_size)]
    if F_reading_frame == 1:
        return _triplets
    elif F_reading_frame == 2:
        _triplets[-1] = [rpf_counts[0], rpf_counts[-2], rpf_counts[-1]]
        return _triplets
    elif F_reading_frame == 3:
        _triplets[-1] = [rpf_counts[0], rpf_counts[1], rpf_counts[-1]]
        return _triplets
    raise ValueError('Reading frame (F) must be in {1, 2, 3}')


def prob_T_given_theta_S_E():
    prob_F_given_theta_S_E = 1/3
    return sum([prob_triplets_given_F_theta_S_E(f) * prob_F_given_theta_S_E for f in range(1, 4)])


def prob_triplets_given_F_theta_S_E(f):
    return 0



def hmm(theta, emissions, states, expression_value):
    """
    Emissions are the vector of triplets (X)
    Hidden states are one of the 9 possible states

    :param theta:
    :param emissions:
    :return:
    """
    # Find emission probabilities
    emission_probs = defaultdict(lambda: 0)

    # Get mu for all states
    # TODO This is a point of dissonance between what I'm seeeing online and what's in the paper
    mu = theta['alpha'] / theta['beta']

    for emission, state in itertools.product(emissions, states):
        """
        emission is a vector<int> of size 3 representing RPF counts
        """
        prob_Y_given_Z = poisson.pmf(sum(emission), mu[state] * expression_value)
        prob_X_given_Y_Z = multinomial.pmf(x=emission, n=sum(emission), p=theta['pi'].loc[state])
        """
        TODO Does the emission (X) need to be uniquely identified? For example, if two separate emissions are both 
        (4,5,5) then they will come out with the same emission prob, is that desired?
        """
        emission_probs[(emission, state)] = prob_X_given_Y_Z * prob_Y_given_Z

    return emission_probs



    # Find transition probabilities


states = ['5primeUTS', '5primeUTS+', 'TIS', 'TIS+', 'TES', 'TTS-', 'TTS', '3primeUTS-', '3primeUTS']
omega_start = {'AUG', 'CUG', 'GUG', 'UUG', 'AAG', 'ACG'}
omega_stop = {'UAA', 'UAG', 'UGA'}
_pi = np.random.randint(1, 5, size=(9, 3))
_pi = _pi / np.sum(_pi, axis=1, keepdims=True)
pi = pd.DataFrame(
    _pi,
    index=states,
    columns=[1, 2, 3]
)
alpha = pd.Series(np.random.random(9), index=states)
beta = pd.Series(np.random.random(9), index=states)
xi_kozak = np.random.random()
xi_start = pd.Series(np.random.random(len(omega_start)), index=omega_start)
theta = {
    'pi': pi,
    'alpha': alpha,
    'beta': beta,
    'xi_kozak': xi_kozak,
    'xi_start': xi_start
}

# Set constant probabilities
pi.loc[['5primeUTS', '5primeUTS+', '3primeUTS-', '3primeUTS']] = 1/3


rpf_counts = list()
rna_seq = list()
expression_tpm = list()

