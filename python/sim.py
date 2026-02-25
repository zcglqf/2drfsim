import MRzeroCore as mr0
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

# Folder where sim.py is located
HERE = Path(__file__).resolve().parent

seq_path = HERE.parent / "epi_gre.seq"   # or put it in a subfolder: HERE / "seq" / "se_epi.seq"
print("seq_path:", seq_path)
print("exists:", seq_path.exists())

# That's it - automatic phantom download and simulation!
signal, ktraj_adc = mr0.util.simulate(str(seq_path))

Nx = 32   # set to your ADC samples per readout
Ny = 32   # number of lines

sig = np.asarray(signal).reshape(Ny, Nx)

img = np.fft.ifftshift(np.fft.ifft2(np.fft.ifftshift(sig)))
plt.figure()
plt.imshow(np.abs(img), cmap="gray")
plt.title("Reconstructed magnitude (quick FFT)")
plt.colorbar()
plt.show()