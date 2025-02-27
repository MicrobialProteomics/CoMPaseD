# Welcome to CoMPaseD
CoMPaseD - **Co**mparrison of **M**ultiple **P**rote**ase** **D**igestions

---
## Purpose
CoMPaseD is a computational tool designed to assist in the selection of optimal protease combinations for the analysis of small proteins in bacterial proteomes. 

Recent studies on *Methanosarcina mazei* {cite}`MMazei-citation` and *Bacillus subtilis* {cite}`BSU-citation` have demonstrated that the use of multiple proteases enhances both, sequence coverage and the number of detected small proteins compared to trypsin digestion alone. However, selecting the most suitable combination often relies on empirical knowledge or extensive experimental trials. CoMPaseD addresses this challenge through a Monte Carlo simulation-based approach that models a typical bottom-up proteomics workflow, including digestion, *semi*-random peptide selection in the mass spectrometer, and protein inference from identified peptides. It then predicts a protease score that objectively evaluates the effectiveness of different digestion strategies.

While particularly useful for studying small proteins, CoMPaseD can also be applied to other sub-proteomes, such as those defined by protein localisation, biochemical properties, or functional categories.

---
## Quick Start
If CoMPaseD is not installed on your system plase refer to [Setup](setup-Installation).

CoMPaseD provides both, a command-line interface (CLI), and a graphical user interface (GUI).  
```{note} The following command-line calls might slightly differ depending on your python setup.
```
For most applications the GUI is the preferred option. It can be started by running the `CoMPaseD_gui.py` command:  
```bash
python CoMPaseD_gui.py
```  
To start an analysis via the CLI, the `CoMPaseD_cli.py` command must be called together with a parameter file:  
```bash
python CoMPaseD_cli.py -p path/to/parameters.params
```

The content of the parameter file is described in the [Parameter File](reference-params) reference.  

Information on further command-line arguments can be found in the [Command-line Interface](reference-cli) reference.  


---
## Limitations  

- Do not use CoMPaseD if you would like to identify an individual protein rather than a sub-proteome.  
- Simulating real data based on actual protein abundance values needs still a lot of testing and parameter optimisation.  
- Usage of extremely large replicate numbers for random sampling (several hundreds) is not recommended due to memory issues.

---
## Citation  
The accompanying manuscript is currently under preparation.  
