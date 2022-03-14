# Comparative modeling by the AutoModel class
#
# Demonstrates how to build multi-chain models, and symmetry restraints
#
from modeller import *
from modeller.automodel import *    # Load the AutoModel class

log.verbose()

# Override the 'special_restraints' and 'user_after_single_model' methods:
class MyModel(AutoModel):
    #def special_restraints(self, aln):
        # Constrain the A and B chains to be identical (but only restrain
        # the C-alpha atoms, to reduce the number of interatomic distances
        # that need to be calculated):
       # s1 = Selection(self.chains['A']).only_atom_types('CA')
       # s2 = Selection(self.chains['B']).only_atom_types('CA')
       # s3 = Selection(self.chains['C']).only_atom_types('CA')
        #self.restraints.symmetry.append(Symmetry(s1, s2, s3, 1.0))
    def user_after_single_model(self):
        # Report on symmetry violations greater than 1A after building
        # each model:
        self.restraints.symmetry.report(1.0)

env = Environ()
# directories for input atom files
env.io.atom_files_directory = ['.','../']

# Be sure to use 'MyModel' rather than 'AutoModel' here!
a = MyModel(env,
            alnfile  = ('both.ali'),     # alignment filename
            knowns   = 'fx11coot-2.pdb',              # codes of the templates
            sequence = 'HHR27186')              # code of the target

a.starting_model= 3                # index of the first model
a.ending_model  = 3                # index of the last model
                                   # (determines how many models to calculate)
a.make()                           # do comparative modeling
