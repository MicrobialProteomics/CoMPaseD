(reference-cli)=
# Command Line Interface (CLI)

The command-line interface can be used by calling the `CoMPaseD_CLI.py` command.

All parameters required for the different steps of CoMPaseD can be set in the `parameters` file provided in the GitHub repository.

```{note}
Many of these parameters can be directly overwritten through the command-line arguments.
```

A full list of command-line arguments with explanations can be displayed using the `-h` flag.

```bash
python CoMPaseD_cli.py -h
```

## Options  

| Parameter Type            | Command Line Argument       | Description                                                                                                          |
|---------------------------|-----------------------------|----------------------------------------------------------------------------------------------------------------------|
| Help                      | -h, --help                  | show the help message and exit                                                                                       |
| Parameter File            | -p, --param_file            | path to CoMPaseD parameter file                                                                                      |
| CoMPaseD Mode             | -e, --export                | export simulated protein abundance values                                                                            |
| CoMPaseD Mode             | -d, --digest                | perform in-silico digest using crux toolkit                                                                          |
| CoMPaseD Mode             | -a, --analysis              | perform analysis from simulated protein abundance and in-silico digestion                                            |
| Digestion Mode Arguments  | --use_original_proteomapper | use original perl scripts for mapping in-silico digested peptides, this might be slower but requires less memory     |
| Digestion Mode Arguments  | --differentiate_I_L         | distinguish between peptide variants containing leucine or iso-leucine (default treat as identical)                  |
| Digestion Mode Arguments  | --indexing_key_len          | length in amino acids of the indexing keys for mapping                                                               |
| Analysis Mode Arguments   | --export_result             | path to CoMPaseD export result file with simulated protein abundance values and protein group assignment             |
| Analysis Mode Arguments   | --digestion_result          | path to CoMPaseD digestion result file ('unique_peptides_table_filtered')                                            |
| Overwrite Parameter File  | --out_folder                | change output directory                                                                                              |
| Overwrite Parameter File  | --fasta                     | change fasta file                                                                                                    |
| Overwrite Parameter File  | --score_peptide             | change weight of peptide IDs for protease score calculation                                                          |
| Overwrite Parameter File  | --score_protein             | change weight of protein IDs for protease score calculation                                                          |
| Overwrite Parameter File  | --score_coverage            | change weight of protein coverage for protease score calculation                                                     |
| Overwrite Parameter File  | --enzymes                   | change enzymes                                                                                                       |
| Overwrite Parameter File  | --mc                        | change maximal missed cleavage sites                                                                                 |
| Overwrite Parameter File  | --mc_freq                   | change missed cleavage sites frequency                                                                               |
| Overwrite Parameter File  | --num_peps                  | change number of peptides to sample                                                                                  |
| Overwrite Parameter File  | --min_pep_mw                | minimal molecular weight of peptides after digestion in Da (default: 400)                                            |
| Overwrite Parameter File  | --max_pep_mw                | maximal molecular weight of peptides after digestion in Da (default: 6000)                                           |
| Overwrite Parameter File  | --min_pep_len               | minimal length of peptides after digestion in amino acids (default: 6)                                               |
| Overwrite Parameter File  | --max_pep_len               | maximal length of peptides after digestion in amino acids (default: 55)                                              |
| Overwrite Parameter File  | --frac_peps                 | change fraction of generated peptides to sample                                                                      |
| Overwrite Parameter File  | --bins                      | change protein binning                                                                                               |
| Overwrite Parameter File  | --undetectable              | change undetectable protein fraction in protein bins                                                                 |
| Overwrite Parameter File  | --DMSP_weight               | change weighting factor of Deep-MS-Peptide prediction                                                                |
| Overwrite Parameter File  | --DMSP_model                | change Deep-MS-Peptide prediction model path                                                                         |
| Overwrite Parameter File  | --samplings                 | change number of random samplings                                                                                    |
| Overwrite Parameter File  | --dynamic_range             | change dynamic range of protein abundance                                                                            |
| Overwrite Parameter File  | --use_unique_peptides_only  | change whether to use only unique peptides or assemble protein groups and consider shared peptides (default is true) |


To run the full analysis provided by CoMPaseD use the following command.

```bash
python CoMPaseD_cli.py -p path/to/parameters
```

To run specific parts of the analysis, the export (-e), digestion (-d) or result analysis (-a) flags can be used.

```bash
python CoMPaseD_cli.py -p path/to/parameter.param [-e] [-d] [-a]
```

```{note} 
The [-a] flag requires to run the rest of the analysis beforehand [-e] and [-d].
```

All result files are written to the folder/files indicated in the `parameters` file or given as command-line arguments.

## Advanced usage

For complex analyses or the comparison of different settings, CoMPaseD can be used in a command line mode by executing `CoMPaseD_CLI.py`. This also allows the automation of several analyses via shell scripts or batch files.

For example, the following commands in a Windows batch file would run an analysis with identical parameters for three different organisms and output the results to folders named `Organism_1`, `Organism_2` and `Organism_3`:

```
CoMPaseD_cli.py -p C:\CoMPaseD.param --out_folder C:\Result\Organism_1 --fasta C:\Organism_1.fasta
```

```
CoMPaseD_cli.py -p C:\CoMPaseD.param --out_folder C:\Result\Organism_2 --fasta C:\Organism_2.fasta
```

```
CoMPaseD_cli.py -p C:\CoMPaseD.param --out_folder C:\Result\Organism_3 --fasta C:\Organism_3.fasta
```

Similarly, the following command in a Linux shell script would run one analysis with each of the two Deep MS-Peptide models using an existing export and digestion file which was specified in the parameter file:

```
python3 /home/user/CoMPaseD/CoMPaseD_cli.py -p /home/user/CoMPaseD.param -a --DMSP_model /home/user/CoMPaseD/bin/DeepMSPep_Confetti_Model.h5 --out_folder /home/user/results/confetti_model/
```

```
python3 /home/user/CoMPaseD/CoMPaseD_cli.py -p /home/user/CoMPaseD.param -a --DMSP_model /home/user/CoMPaseD/bin/DeepMSPep_Original_Model.h5 --out_folder /home/user/results/original_model/
```

Command line options will generally overwrite values in the parameter file, and a parameter file with the used values will automatically be saved in the output folder together with the results.
The complete list of options available can be displayed by executing:

```bash
CoMPaseD_cli.py --help
```
