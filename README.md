![image](https://github.com/grfrederic/visavis/assets/26434160/c0ab759c-fd78-4264-9237-1f33828465bd)

Overview
--------

*VIS-A-VIS* is an agent-based simulator of viral infection spread and viral
infection self-containment in a monolayer of cell.

The simulation mimics the innate immune response to an infection with an RNA
virus. The hard-coded and externally parametrized interactions between the host
cell and virus are specific to the respiratory syncytial virus (RSV) and human
alveolar epithelial cells (A549 cell line). Infected cells attempt to produce
and secrete an interferon, which alerts the non-infected bystander cells
about the nearby threat.

*VIS-A-VIS* is a hybrid simulator as it executes alternating phases of
interferon diffusion and chemical kinetics. The interferon diffuses between
compartments above neighboring lattice nodes and its local concentrations are
expressed using real variables. The nodes are occupied by cells, in which
stochastic chemical kinetics (Gillespie algorithm) causes transitions within
a predefined set of states for molecules of several chemical species. Some of
these states are associated with higher level or higher activity of these
molecules. The lattice has periodic boundary conditions.

*VIS-A-VIS* was used in the paper "[Antagonism between viral infection and innate immunity at the single-cell level](https://doi.org/10.1371/journal.ppat.1011597)", *PLOS Pathogens* (2023). See citing section below.


Compilation
-----------

The simulator has been implemented in [Rust](https://www.rust-lang.org). If your
Rust toolchain has been set up via [rustup](https://rustup.rs), then prior to
compilation in the command line you may want to issue:

```bash
$ rustup update
```

To compile the source code, enter the `visavis-1.0.0` folder and type:
```bash
$ cargo build --release
```

The software has been tested to successfully compile under the current stable
toolchain (shipping rustc version 1.64.0). If libglib2.0-dev and libcairo2-dev
are installed, all dependencies (rust crates) should be retrieved and compiled
on the fly.


Running
-------

To run a simulation using default parameters and stimulation protocol,
you may invoke the simulator with:
```bash
$ cargo run --release parameters/default.json protocols/default.protocol --images
```
which is equivalent to:
```bash
$ target/release/vis-a-vis parameters/default.json protocols/default.protocol --images
```

At simulated-time intervals prescribed in the protocol, the simulator produces
CSV-formatted files with the current state of all molecules in the lattice.
Generating PNG images with lattice snapshots (turned on when using the `--images`
argument) is optional.


Output
------

**CSV**: Each generated CSV file fully describes the state of the simulated system
in the time point specified in the file name. Each row of the file corresponds
to a lattice node that contains a cell, whose molecular internal state is given
in comma-separated entries corresponding to molecules as given in the file header.
States of molecules are discrete (and are in the range `MIN..MAX` defined in
respective arrays in `src/cell.rs`). Last two columns give the amount of
the extracellular interferon in the cell culture medium above the cell-node
(in the lower and in the upper subcompartment separately).

**PNG**: Lattice images depict cells as circles inscribed in hexagons. The more
yellow is the hexagon fill, the higher is the amount of the extracellular
interferon in the lower medium subcompartment above the cell. Pinkish outer
cell ring indicates viral infection; reddish color of the inner circle
corresponds to a high level of p-IRF3, whereas greenish inner circle shows
STAT1/2 activity.

If in module lattice (`src/lattice.rs`) the boolean variable
`Lattice::NEIGHS_TO_FILE` is set to true, then additionally a file `neighbors.csv`
with complete information about lattice node neighborhoods is dumped.

All the output files are generated in the current working directory.


Tweaking
--------

Changes in stimulation protocols (provided as text files, see included examples
in `protocols/`) or in parameter values (provided as JSON-formatted text files,
see examples included in `parameters/`) do not require the code to be recompiled.

Modifications of the wiring of the molecular virus--host and intra-host
interactions require changes in module simulation (`src/simulation.rs`) and
code recompilation.

To change lattice size, you need not only to change `Lattice::WIDTH`,
`Lattice::HEIGHT` in module lattice (`src/lattice.rs`), but also likely tweak
the `THREAD_STACK_SIZE` parameter in module config (`src/config.rs`), and then
recompile the code.


Extra: Python wrapper
---------------------

Since running many simulations from terminal can be cumbersome, in `extra/`
we provide a simple Python wrapper and associated convenience code for working
with *Vis-A-Vis*:
  * `visavis.py`: the *Vis-A-Vis* client,
  * `simulation_result.py`: class wrapping the result of a simulation,
  * `annotate.py`: script for annotating lattice snapshots,
  * `parameters/default.py`: Python format for parameters.

Script `extra/example.py` demonstrates concisely how the convenience code may
be used. The script can be run with:
```bash
$ cd extra/
$ python example.py
```


Citing
------

The simulator, written by Marek Kochanczyk and Frederic Grabowski from the
Institute of Fundamental Technological Research of the Polish Academy of
Sciences (IPPT PAN, Warsaw, Poland), features the research article
"[Antagonism between viral infection and innate immunity at the single-cell level](https://doi.org/10.1371/journal.ppat.1011597)"
published in *PLOS Pathogens* (2023). If you use the code for research purposes,
please consider citing this work.


License
-------

The code is open source under the 3-Clause BSD license (see file LICENSE or
visit https://opensource.org/licenses/BSD-3-Clause).
