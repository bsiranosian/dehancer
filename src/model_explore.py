from dragonn.utils import *
from os.path import join

fastaDir='/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/fasta'
# start with one cell line for now 
cell='Dnd41'

enhancerFasta=join(fastaDir, cell+'_enhancers_windows.fa')
promoterFasta=join(fastaDir, cell+'_promoters_windows.fa')
negativeFasta=join(fastaDir, cell+'_negative_windows.fa')

# load up the sequences and one hot encode them 
eoh = encode_fasta_sequences(enhancerFasta)
poh = encode_fasta_sequences(promoterFasta)
noh = encode_fasta_sequences(negativeFasta)

