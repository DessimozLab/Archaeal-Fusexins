# importing libraries
import os
import subprocess
import glob

# listing MEs
ME_data = [(ME.replace('ME.', '').replace('.fna', ''),
            '../data/mobile_elements_sequences/{0}'.format(ME)) for ME in os.listdir('../data/mobile_elements_sequences/')]

# creating directory to allocate predicted ORFs
if not os.path.exists('../results/MEs_predicted_orfs'):
  os.mkdir('../results/MEs_predicted_orfs')
  
# running getorf (EMBOSS package) with bacterial genetic code, and between START and STOP codons
for me_data in ME_data:
  # unpacking variables
  genome_flag, me_file = me_data
  # running getorf if file doesnt exist
  if not os.path.exists('../results/MEs_predicted_orfs/ME.{0}.predicted_orfs.faa'.format(genome_flag)):
    getorf_command = 'getorf -sequence {0} -find 1 -table 11 -outseq ../results/MEs_predicted_orfs/ME.{1}.predicted_orfs.faa'.format(me_file, genome_flag).split(' ')
    subprocess.run(getorf_command)
