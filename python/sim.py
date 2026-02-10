import MRzeroCore as mr0
from pathlib import Path

# Folder where sim.py is located
HERE = Path(__file__).resolve().parent

seq_path = HERE / "epi_se.seq"   # or put it in a subfolder: HERE / "seq" / "se_epi.seq"
print("seq_path:", seq_path)
print("exists:", seq_path.exists())

# That's it - automatic phantom download and simulation!
signal, ktraj_adc = mr0.util.simulate(str(seq_path))