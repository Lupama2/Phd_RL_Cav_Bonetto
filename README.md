# Cavitation Bubble Dynamics (Bonetto et al.)

Numerical model of ultrasound‑driven cavitation bubbles following Bonetto et al.  This model is derived from first principles, and no parameters are fitted when comparing with experiments. The bubble's system is divided into three regions: the bubble interior, the bubble interface, and the surrounding liquid. In the interior, the model includes transient and spatially nonuniform heat transfer, gas dissociation, chemical reactions, and vapor mass diffusion in the noncondensable gas mixture. At the interface, it includes nonequilibrium evaporation/condensation of water and a temperature jump between regions. In the liquid, it includes transient and spatially nonuniform heat transfer, liquid compressibility, and mass diffusion of gas in the liquid.

## Quick Start
- **Requirement:** Borland C++ Builder 6 (bcc32 / bcb). Project folders must **not** contain hyphens (`-`) or the toolchain will fail to link.
- **Build full project:** `bcb -B src/PulWid.bpr`
- **Build single translation unit:** `bcc32 -c src/PSULSOAR.cpp`
- **Run:** execute the generated `.exe`; output time series are written as `.dat` files in `src/` by default.

## What’s in this repository
- `src/sbsl.cpp` — entry point that advances the cavitation bubble dynamics.
- `src/ecdiff.cpp` — time‑dependent differential equation.
- `src/funcions.cpp` — utility functions used across the solver.
- `src/psulsoar.cpp` — model parameters (pressure drive, physical constants, initial conditions).
- `src/Rkqc.cpp` — adaptive RKCK integrator.
- `src/Rkutta4.cpp` — fixed‑step RK4 integrator.
- Borland project files: `PulWid.bpr`, `PulWid.bpf`, `PulWid.res`.

## Usage notes
- Adjust physical parameters in `src/psulsoar.cpp` before building to set the driving pressure, frequency, and medium properties.
- Outputs: `.dat` files contain bubble radius vs time and related quantities; visualize with your preferred plotting tool (e.g., Python/Matplotlib or gnuplot).

## Citation
If you use this code, please cite the Physica Scripta paper (BibTeX will be added once DOI is issued) and acknowledge this repository.

## License
Distributed under the terms of the LICENSE file in this repository.
