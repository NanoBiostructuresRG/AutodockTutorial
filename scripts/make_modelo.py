from modeller import *
from modeller.automodel import *
from modeller.automodel import refine


env = Environ()
a = AutoModel(env, alnfile='seq_6MD4.ali',
              knowns='6MD4_receptor', 
              sequence='PPARG_P37231',
              assess_methods=(assess.DOPE,
                              assess.GA341))
a.starting_model = 1
a.ending_model = 3
a.md_level = refine.very_slow

a.make()

