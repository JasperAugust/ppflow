import pyrosetta
from pyrosetta import pose_from_sequence
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pyrosetta.rosetta.core.scoring import get_score_function

# Initialize PyRosetta
pyrosetta.init()

# Your peptide sequence
sequence = "CGNKRTRGC"

pose = pose_from_sequence(sequence, "fa_standard")  

scorefxn = get_score_function()
min_mover = MinMover()
min_mover.score_function(scorefxn)
min_mover.min_type('lbfgs_armijo_nonmonotone')
min_mover.max_iter(100)  # Maximum number of iterations

min_mover.apply(pose)

# Output PDB file path
output_pdb = "lyp1.pdb"

# Dump the pose to a PDB file
pose.dump_pdb(output_pdb)

print(f"PDB file saved as {output_pdb}")