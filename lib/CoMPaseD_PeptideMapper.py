import argparse
import os

import colorama
from multiprocessing import Process
from subprocess import Popen, PIPE
from datetime import datetime
from time import perf_counter
from Bio import SeqIO
from os import path, remove, makedirs, chdir, listdir
from pandas import read_csv, DataFrame, Series, concat
from shutil import copy, rmtree
import itertools
import pickle


def main():
    parser = argparse.ArgumentParser(description="map peptide sequences to their positions in proteins in a fasta file with fasta indexing")
    parser.add_argument('--fasta', required=True)
    parser.add_argument('--out_folder', required=True)

    parser.add_argument('--crux_path', required=True)
    parser.add_argument('--enzyme_list', required=True)
    parser.add_argument('--max_mc_list', help='maximum missed cleavage sites', required=True)
    parser.add_argument('--min_mass', required=False, default=400)
    parser.add_argument('--max_mass', required=False, default=6000)
    parser.add_argument('--min_len', required=False, default=6)
    parser.add_argument('--max_len', required=False, default=55)
    parser.add_argument('--digestion_type', required=False, default='full-digest')
    parser.add_argument('--clip_n_term_met', required=False, default='T')
    parser.add_argument('--decoy_format', required=False, default='None')
    parser.add_argument('--unique_peps_file', required=False, default='')

    parser.add_argument('--indexing_key_len', required=False, default=5)
    parser.add_argument('--differentiate_I_L', action='store_false')


    # get start time
    time_0 = perf_counter()

    # output to progress tab:
    print(f"CoMPaseD - Comparison of Multiple-Protease Digestions", flush=True)
    print(f"", flush=True)
    print(f"In-silico digestion started", flush=True)
    print(f"", flush=True)
    print(f"---------------------------------------------------------------------------", flush=True)
    print(f"Using Crux mass spectrometry toolkit for digestion", flush=True)
    print(f"\tfor a detailed description please see:", flush=True)
    print(f"\tSean McIlwain, Kaipo Tamura, Attila Kertesz-Farkas, Charles E. Grant, Benjamin Diament, Barbara Frewen,")
    print(f"\tJ. Jeffry Howbert, Michael R. Hoopmann, Lukas Käll, Jimmy K. Eng, Michael, J.MacCoss and William S. Noble:",
        flush=True)
    print(f"\tCrux: rapid open source protein tandem mass spectrometry analysis.")
    print(f"\tJournal of Proteome Research. 13(10):4488-4491, 2014.", flush=True)
    print(f"\tDOI: 10.1021/pr500741y", flush=True)
    print("---------------------------------------------------------------------------", flush=True)

    # switches to remove non-existing dirs later
    tmp_exists = False
    crux_out_exists = False

    # parse arguments
    args = parser.parse_args()

    fasta = args.fasta
    out_folder = args.out_folder
    if not os.path.isabs(out_folder):
        out_folder = os.path.abspath(out_folder)

    # create output folders if necessary
    tmp_out_folder = path.join(out_folder, 'Tmp')
    if not path.isdir(tmp_out_folder):
        makedirs(tmp_out_folder, exist_ok=True)
        print(f"Created output folder", flush=True)
    else:
        tmp_exists = True
        makedirs(tmp_out_folder, exist_ok=True)
        print("WARNING: output folder exists, files may be overwritten", flush=True)

    if path.isdir(path.join(tmp_out_folder, "crux-output")):
       crux_out_exists = True

    # change directory to tmp_out_folder, all crux output will be in this folder
    chdir(tmp_out_folder)

    print("Start digest:", flush=True)
    crux_path = path.join(args.crux_path)
    # clean crux_path in case it contains escape chr's and is not found initially
    # solution from https://stackoverflow.com/questions/18682695/python-escape-character
    # by terrachild on Sep 8, 2013 at 10:31
    if not path.isfile(crux_path):
        backslash_map = {'\a': r'\a', '\b': r'\b', '\f': r'\f',
                         '\n': r'\n', '\r': r'\r', '\t': r'\t', '\v': r'\v'}

        def reconstruct_broken_string(s):
            for key, value in backslash_map.items():
                s = s.replace(key, value)
            return s

        crux_path = path.join(reconstruct_broken_string(crux_path))

    fasta = path.join(args.fasta)
    out_folder = path.join(out_folder)
    proteases_string = args.enzyme_list
    mc_string = args.max_mc_list
    protease_list = proteases_string.split(',')
    mc_list = [int(mc) for mc in mc_string.split(',')]
    min_pep_mw = args.min_mass
    max_pep_mw = args.max_mass
    min_pep_len = args.min_len
    max_pep_len = args.max_len

    # get all crux commands in a list
    crux_cmd_list, exp_protease_list, exp_mc_list, crux_out_file_list = get_crux_cmds(protease_list,
                                                                                      mc_list,
                                                                                      fasta,
                                                                                      crux_path,
                                                                                      min_pep_mw,
                                                                                      max_pep_mw,
                                                                                      min_pep_len,
                                                                                      max_pep_len)
    crux_file_list_file = path.join(tmp_out_folder, 'crux_result_files.tsv')

    with open(crux_file_list_file, 'w') as f:
        f.write('\n'.join(crux_out_file_list))

    crux_proc_list = list()

    def crux_process(crux_cmd, n, total, protease, mc):
        new_proc = Popen(crux_cmd, shell=True, stdout=PIPE, stderr=PIPE)
        print(f"Started digestion with {protease} and {mc} missed cleavages", flush=True)
        out = new_proc.stdout.read()
        err = new_proc.stderr.read()
        new_proc.wait()
        out = out.splitlines()
        err = err.splitlines()
        for line in out:
            print(line.decode('utf8'), flush=True)
        for line in err:
            print(line.decode('utf8'), flush=True)

        print(f"Finished digestion {n} of {total}.", flush=True)

    # count number of digestions
    n = 0
    total = len(crux_cmd_list)

    # init all processes
    for crux_cmd, protease, mc in zip(crux_cmd_list, exp_protease_list, exp_mc_list):
        n += 1
        proc = Process(target=crux_process(crux_cmd, n, total, protease, mc))
        crux_proc_list.append(proc)

    # start multiprocessing queue
    for proc in crux_proc_list:
        proc.start()

    # wait for finishing
    for proc in crux_proc_list:
        proc.join()

    print("", flush=True)
    print("Finished all In-silico digestions", flush=True)
    print("", flush=True)
    print("", flush=True)
    print("---------------------------------------------------------------------------", flush=True)
    print("", flush=True)
    print("", flush=True)

    ILEquivalence = args.differentiate_I_L
    splitLen = int(args.indexing_key_len)
    # generate fasta index
    idxing_result = generate_index(fasta, splitLen=splitLen, aaAlphabet = 'extended', ILEquivalence=ILEquivalence)
    # check index generation
    if not idxing_result == 0:
        print(f"{colorama.Fore.RED}ERROR: Index generation for {fasta} failed. Please check. Stopping.{colorama.Style.RESET_ALL}", flush=True)
        print("", flush=True)
        return 1

    # get crux-output file names
    crux_file_list_file = path.join(tmp_out_folder, 'crux_result_files.tsv')

    # lists for protease, crux_file, out_file and MCs
    crux_file_list = list()
    mapped_crux_file_list = list()
    protease_list = list()
    mc_list = list()

    # pre-load fasta index to avoid repeated loading, this saves some time
    fastaIndex = path.join(path.dirname(fasta), (path.basename(fasta) + ".idx.pickle"))
    with open(fastaIndex, 'rb') as handle:
        fasta_idx = pickle.load(handle)
        print(f"\t \t Loaded fasta index", flush=True)

    # open list with crux result files
    with open(crux_file_list_file, 'r') as crux_f:
        # process line-by-line
        while (line := crux_f.readline().rstrip()):
            # crux_file_list_file is tabular with file names in column [0]
            tmp_col = line.split(sep='\t')

            # generate complete file name from current line
            crux_file = path.join(tmp_out_folder, 'crux-output', tmp_col[0])
            # add protease and MC information from current line to lists
            protease_list.append(str(tmp_col[1]))
            mc_list.append(str(tmp_col[2]))
            # generate output file name
            mapped_crux_file = path.join(tmp_out_folder, str('Mapped_' + tmp_col[0]))

            # read peptide column (first col) from crux file into memory and convert to list
            # df should be sorted alphabetically
            # use pandas read_csv, header=None is required to not loose first peptide
            df = read_csv(crux_file, delimiter='\t', usecols=[0], header=None)[0].to_list()

            # output format should be:
            # header line, 6 cols: peptide \t protein \t location \t prevAA \t in_fasta \t nextAA \n
            # body:                AAAAFRVVK \t lcl|AL009126.3_prot_2464 \t 19 \t \t \t \n
            mapping_result_list = map_peptides(df, fasta, splitLen=splitLen, ILEquivalence=ILEquivalence, protease= str(tmp_col[1]), MC= str(tmp_col[2]), fasta_idx=fasta_idx)

            # check mapping result
            if mapping_result_list == 1:
                print(f"{colorama.Fore.RED}ERROR: Mapping failed for {crux_file}. Please check. Stopping.{colorama.Style.RESET_ALL}", flush=True)
                print("", flush=True)
                return 1

            # save result to file
            with open(mapped_crux_file, 'w') as f:
                for entry in mapping_result_list:
                    f.write(entry)

            # log finished crux files
            crux_file_list.append(crux_file)
            mapped_crux_file_list.append(mapped_crux_file)

    # return from open crux_file_list_file
    print("", flush=True)
    print("Finished peptide mapping", flush=True)
    print("", flush=True)
    print("", flush=True)
    print("---------------------------------------------------------------------------", flush=True)
    print("", flush=True)
    print("", flush=True)

    # generate list with mapped file names and protease / mc combination
    crux_out_folder = path.join(out_folder, 'Tmp', 'crux-output')
    mapped_file_list_file = path.join(crux_out_folder, 'File_List.tsv')
    mapped_out_file_list = list()
    for mf, pr, mc in zip(mapped_crux_file_list, protease_list, mc_list):
        ln = str(mf) + "\t" + str(pr) + "\t" + str(mc)
        mapped_out_file_list.append(ln)
    # # save this list
    with open(mapped_file_list_file, 'w') as mf:
        mf.write('\n'.join(mapped_out_file_list))
    # remove files from tmp output and keep only 'Mapped_' peptide lists and 'File_List'
    print("Removing tmp files", flush=True)
    for f in listdir(crux_out_folder):
        if not f.startswith('Mapped_'):
            if not f.startswith('File_List'):
                remove(path.join(crux_out_folder, f))
    print("Finished removing tmp files ", flush=True)

    # merge mapped files and annotate table
    print("", flush=True)
    print("", flush=True)
    print("Started pooling digests", flush=True)
    df_complete = DataFrame()
    with open(mapped_file_list_file, 'r') as mfl:
        while(line := mfl.readline().rstrip()):
            # file names with mapped peptides and annotation is tabular
            tmp_col = line.split(sep='\t')
            f = tmp_col[0]
            protease = tmp_col[1]
            mc = tmp_col[2]
            # read relevant columns to pandas df
            df = read_csv(f, delimiter='\t', usecols=["peptide", "protein", "location"])
            # keep = False to keep only unique peptides within a certain protease / mc combination
            df_reduced = df.drop_duplicates(subset="peptide", keep=False, ignore_index=True)

            df_reduced = df_reduced.assign(MC=Series([mc] * len(df_reduced.index)))
            df_reduced = df_reduced.assign(Enzyme=Series([protease] * len(df_reduced.index)))

            # df_complete = df_complete.append(df_reduced, ignore_index=True) # deprecated
            df_complete = concat([df_complete, df_reduced], axis=0, join='outer', ignore_index=True)

    complete_df_file = path.join(out_folder, "unique_peptides_table_unfiltered.tsv")
    df_complete.to_csv(complete_df_file, index=False, sep="\t")

    # remove redundant peptides from multiple missed cleavages
    df_complete_sorted = df_complete.sort_values(by='MC', ascending=False, kind='mergesort')
    # keep=last to remove duplicates with less than max mc in each file (crux generates peptides with zero to max mc's)
    df_filtered_sorted = df_complete_sorted.drop_duplicates(subset=['peptide', 'Enzyme'], keep='last',
                                                            ignore_index=True)

    if not path.join(args.unique_peps_file) == "":
        filtered_df_file = path.join(out_folder, "unique_peptides_table_filtered.tsv")
    else:
        filtered_df_file = path.join(out_folder, "unique_peptides_table_filtered.tsv")

    df_filtered_sorted.to_csv(filtered_df_file, index=False, sep="\t")

    print("Finished pooling digests", flush=True)
    chdir(out_folder)

    # try to remove folders that eventually were created
    if not crux_out_exists:
        rmtree(path.join(tmp_out_folder, "crux-output"), ignore_errors=True)
    if not tmp_exists:
        rmtree(tmp_out_folder, ignore_errors=True)

    print("", flush=True)
    print("", flush=True)
    print("---------------------------------------------------------------------------", flush=True)
    time_1 = perf_counter()
    print(f"Digestion took {time_1 - time_0:0.1f} seconds")
    print("", flush=True)
    print("", flush=True)


def clean_protease_names(protease_list: list) -> list:
    '''remove special characters from protease names in preparation of file name creation'''
    protease_list_clean = list()
    for each_protease in protease_list:
        for character in '\\/:*?"<>|,;-':
            each_protease = each_protease.replace(character, '_')
        # ensure that no duplicates are created
        while each_protease in protease_list_clean:
            each_protease = each_protease + "_1"

        protease_list_clean.append(each_protease)
    return protease_list_clean


def generate_peptides_call(out_folder,
                           fasta,
                           crux_path,
                           min_mass=400,
                           max_mass=6000,
                           min_len=6,
                           max_len=55,
                           missed_cleavages=0,
                           digestion="full-digest",
                           clip_n_term_met="T",
                           decoy_format="none",
                           enzyme="trypsin") -> str:
    """command to execute crux generate-peptides for one enzyme/MCs/fasta/... combination"""
    # the digestion specificity could be set to semi-specific (see "Optimised Proteomics Workflow for the Detection
    # of Small Proteins", Journal of Proteome Research, 2020)
    # check parameters
    # with gui-usage there is no possibility to manually select other than full-digest, kept for future compatibility
    allowed_digestion = {"full-digest", "partial-digest", "non-specific-digest"}
    if digestion not in allowed_digestion:
        raise ValueError("digestion must be one of %r." % allowed_digestion)
    # Caution: the following enzymes use the proline rule, i.e. do not cleave if cleavage site is followed by P:
    # "trypsin", "chymotrypsin", "elastase","arg-c","glu-c","pepsin-a","elastase-trypsin-chymotrypsin","lys-c"
    allowed_enzyme = {"trypsin", "trypsin/p", "chymotrypsin", "elastase", "clostripain", 'cyanogen-bromide',
                      "iodosobenzoate", "proline-endopeptidase", "staph-protease", "asp-n", "lys-c", "lys-n",
                      "arg-c",
                      "glu-c", "pepsin-a", "elastase-trypsin-chymotrypsin", "lysarginase"}
    if enzyme not in allowed_enzyme:
        raise ValueError("Protease not defined. Must be one of %r." % allowed_enzyme)

    # generate crux command
    peptide_cutter_command = str(f'"{crux_path}"') + " generate-peptides --overwrite T --min-mass " + str(
        min_mass) + " --max-mass " + str(
        max_mass) + " --min-length " + str(min_len) + " --max-length " + str(max_len) + " --digestion " + str(
        digestion) + " --missed-cleavages " + str(missed_cleavages) + " --clip-nterm-methionine " + str(
        clip_n_term_met) + " --decoy-format " + str(decoy_format) + " --enzyme " + str(
        enzyme) + " --fileroot " + str(f'"{out_folder}"') + " " + str(f'"{fasta}"')

    return peptide_cutter_command


def get_crux_cmds(protease_list, mc_list, fasta, crux_path, min_pep_mw=400, max_pep_mw=6000, min_pep_len=6, max_pep_len=55):
    # get plain list with possible protease / MCs combinations
    expanded_protease_list = list()
    expanded_mc_list = list()
    # list out files to file for later read back
    crux_out_file_list = list()
    # use protease_list instead of protease_list_clean for crux command
    for each_protease, each_max_mc in zip(protease_list, mc_list):
        for each_mc in range(each_max_mc + 1):
            expanded_protease_list.append(each_protease)
            expanded_mc_list.append(each_mc)
    crux_cmd_list = list()
    crux_opt_list = list()
    # loop through both expanded lists and generate crux call
    for protease, mc in zip(expanded_protease_list, expanded_mc_list):
        protease_ls = list()
        protease_ls.append(protease)
        protease_clean = clean_protease_names(protease_ls)
        result_name = str(str(protease_clean[0]) + "_" + str(mc) + "_MCs")
        crux_cmd = generate_peptides_call(out_folder=result_name,
                                          fasta=path.join(fasta),
                                          crux_path=path.join(crux_path),
                                          enzyme=protease,
                                          missed_cleavages=mc,
                                          min_mass=min_pep_mw,
                                          max_mass=max_pep_mw,
                                          min_len=min_pep_len,
                                          max_len=max_pep_len)
        crux_cmd_list.append(crux_cmd)
        crux_out_file_str = str(protease_clean[0]) + "_" + str(mc) + \
                            "_MCs.generate-peptides.target.txt" + "\t" + protease + \
                            "\t" + str(mc)

        crux_out_file_list.append(crux_out_file_str)
    return crux_cmd_list, expanded_protease_list, expanded_mc_list, crux_out_file_list


def generate_index(fastaFile, splitLen=5, aaAlphabet='ACDEFGHIKLMNPQRSTVWYBXZJUO', ILEquivalence = True):
    '''
    This is a simplified python implementation of ProteoMappers Clips.pl script

    Function to generate a pickled fasta index based on aa-combinations of length=splitLen
        protein identifier are encoded as six-digit hexadecimal values and
         start position of each aa-combination is assigned (1-based), e.g:

         'AAACD': ['000001, 1', '0000CD, 18', '0000E3, 21', '0000EF, 2']

         would be one entry indicating matches of the sequence AAACD to four proteins at
         position 1, 18, 21 and 2 respectively.

         As the index is saved as a python dict memory may become an issue for very large fasta files.
         An additional dict will be created to save ID mappings and log index generation settings.

         Extended alphabet from Bio.Alphabet is set as default for aaAlphabet to ensure complete mapping
         but may be modified by direct setting ('ABC...') or the following key-words: 20-aa; 22-aa, extended

    '''

    # correct aaAlphabet for key-words
    if aaAlphabet == '20-aa':
        if ILEquivalence:
            aaAlphabet = 'ACDEFGHKLMNPQRSTVWY'
        else:
            aaAlphabet = 'ACDEFGHIKLMNPQRSTVWY'

    if aaAlphabet == '22-aa':
        if ILEquivalence:
            aaAlphabet = 'ACDEFGHKLMNPQRSTVWYOU'
        else:
            aaAlphabet = 'ACDEFGHIKLMNPQRSTVWYOU'

    if aaAlphabet == 'extended':
        if ILEquivalence:
            aaAlphabet = 'ACDEFGHKLMNPQRSTVWYBXZJUO'
        else:
            aaAlphabet = 'ACDEFGHIKLMNPQRSTVWYBXZJUO'

    # inform user on indexing start and currently used aaAlphabet
    print(f"---------------------------------------------------------------------------", flush=True)
    print(f"Using a simplified python implementation of ProteoMappers clips.pl and promast.pl for peptide to protein mapping", flush=True)
    print(f"\tfor a detailed description of ProteoMapper please see:", flush=True)
    print(f"\tLuis Mendoza, Eric W. Deutsch, Zhi Sun, David S. Campbell, David D. Shteynberg and Robert L. Moritz:",
        flush=True)
    print(f"\tFlexible and fast mapping of peptides to a proteome with ProteoMapper")
    print(f"\tJournal of Proteome Research. 17(12):4337-4344, 2018.", flush=True)
    print(f"\tDOI: 10.1021/acs.jproteome.8b00544", flush=True)
    print(f"---------------------------------------------------------------------------", flush=True)

    print(f"\t Started index generation for {fastaFile}.", flush=True)
    print(f"\t \t indexed amino acids: {aaAlphabet}.", flush=True)
    if ILEquivalence:
        print(f"\t \t treating L and I amino acids as the same.", flush=True)

    # check fasta file existence
    if not path.isfile(fastaFile):
        print(f"{colorama.Fore.RED}ERROR: Index generation for {fastaFile} failed. Please check if this file exists. \n Stopping. {colorama.Style.RESET_ALL}", flush = True)
        print("", flush=True)
        return 1
    else:
        # init entry ID to hexadecimal conversion for reduced file size compared to decimal numbering
        n = 0
        id_dict = dict()
        inv_id_dict = dict()
        id_seq_dict = dict()

        # generate all possible combinations, result is list of tuples
        aa_idxs = [p for p in itertools.product(aaAlphabet, repeat=splitLen)]

        # create dict to hold all combinations
        split_dict = dict()
        # loop through aa_idxs, convert each idx tuple of individual characters to one string and add empty list
        for idx in aa_idxs:
            split_dict[''.join(idx)] = list()

        # open fasta
        with open(fastaFile) as handle:
            # process each entry separately
            for record in SeqIO.parse(handle, "fasta"):

                # get sequence
                rec_seq = record.seq
                # replace I with L
                if ILEquivalence:
                    rec_seq = rec_seq.replace('I', 'L')

                # fill ID translation dicts
                id_dict[f'{n:06x}'] = record.id
                inv_id_dict[record.id] = f'{n:06x}'
                # fill id_seq_dict
                id_seq_dict[f'{n:06x}'] = rec_seq

                # split sequence into substrings of split_len length
                splits = [rec_seq[idx:idx + splitLen] for idx in range(len(rec_seq) - splitLen + 1)]

                # add splits to split_dict
                ii = 1
                for split in splits:
                    split_dict[split].append(str(inv_id_dict[record.id]) + ', ' + str(ii))
                    ii += 1

                # increase ID counter
                n += 1

        # generate standardised idx file name
        fastaIndex = path.join(path.dirname(fastaFile), (path.basename(fastaFile) + ".idx.pickle"))
        indexAnnot = path.join(path.dirname(fastaFile), (path.basename(fastaFile) + ".annot.pickle"))

        # remove any existing index but warn the user about this
        if path.isfile(fastaIndex):
            print(f"\t\tWARNING: Existing index ({fastaIndex}) was found and will be replaced. Please check", flush=True)
            try:
                remove(fastaIndex)
            except Exception as e:
                print(f"{colorama.Fore.RED}ERROR: Existing fasta file index {fastaIndex} could not be deleted due to {e}. Please check permission. \n Stopping. {colorama.Style.RESET_ALL}", flush=True)
                print("", flush=True)
                return 1

        # remove any existing index annotation but warn the user about this
        if path.isfile(indexAnnot):
            print(f"\t\tWARNING: Existing index annotation ({indexAnnot}) was found and will be replaced. Please check", flush=True)
            try:
                remove(indexAnnot)
            except Exception as e:
                print(f"{colorama.Fore.RED}ERROR: Existing fasta file index annotation for {indexAnnot} could not be deleted due to {e}. Please check permission. \n Stopping. {colorama.Style.RESET_ALL}", flush=True)
                print("", flush=True)
                return 1

        with open(fastaIndex, 'wb') as handle:
            pickle.dump(split_dict, handle)

        with open(indexAnnot, 'wb') as handle:
            pickle.dump(fastaFile, handle) # fasta file
            pickle.dump(fastaIndex, handle) # index file
            pickle.dump(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S")), handle) # creation date and time
            pickle.dump(splitLen, handle)  # index splitLen
            pickle.dump(id_dict, handle)  # id-to-hexid dict
            pickle.dump(inv_id_dict, handle)  # hexid-to-id dict
            pickle.dump(id_seq_dict, handle) # id-to-sequence dict
            pickle.dump(aaAlphabet, handle)  # aa-alphabet
            if ILEquivalence:
                pickle.dump(True, handle) # ILEquivalence
            else:
                pickle.dump(False, handle) # ILEquivalence

        print(f"\t Finished index generation.", flush=True)
        print("", flush=True)

        return 0

def map_peptides(peptideList, fastaFile, splitLen, ILEquivalence, protease, MC, fasta_idx = None):
    '''
    This is a simplified python implementation of ProteoMappers Promast.pl script

    Function to map a list of peptides to their positions in a fasta file
        using the pickled index file as template
         start position of mappings is 1-based, e.g in:
         >Protein
         MPEPTIDE
         the peptide EPTIDE is at pos 3-8

    '''

    # check whether fasta file exists
    if not path.isfile(fastaFile):
        print(f"{colorama.Fore.RED}ERROR: Fasta file not existing. Please check: {fastaFile}. \n Stopping. {colorama.Style.RESET_ALL}",
            flush=True)
        print("", flush=True)
        return 1

    # init indexing result
    indexing_result = 0

    # retrieve standardised idx file name
    fastaIndex = path.join(path.dirname(fastaFile), (path.basename(fastaFile) + ".idx.pickle"))
    indexAnnot = path.join(path.dirname(fastaFile), (path.basename(fastaFile) + ".annot.pickle"))

    if (not path.isfile(fastaIndex)) or (not path.isfile(indexAnnot)):
        print(f"\t\tWARNING: No existing index was found. Start indexing fasta file.", flush=True)
        try:
            indexing_result = generate_index(fastaFile=fastaFile, splitLen=splitLen, aaAlphabet='extended', ILEquivalence=ILEquivalence)
        except Exception as e:
            print(
                f"{colorama.Fore.RED}ERROR: Could not create fasta file index for {fastaFile} due to {e}. Please check permission. \n Stopping. {colorama.Style.RESET_ALL}",
                flush=True)
            print("", flush=True)
            return 1

    if indexing_result == 0:
        # load annotations
        with open(indexAnnot, 'rb') as handle:
            # re-load fasta file and index file names from annotation to ensure correct index-fasta mapping
            fastaFile = pickle.load(handle)
            fastaIndex = pickle.load(handle)
            # show idx_creation time, splitLen and len of fasta file in log later
            idx_creation_date_time = pickle.load(handle)
            splitLen = pickle.load(handle)
            id_dict = pickle.load(handle)  # id-to-hexid dict
            inv_id_dict = pickle.load(handle)  # hexid-to-id dict
            id_seq_dict = pickle.load(handle)  # hexid-to-sequence dict
            aaAlphabet = pickle.load(handle)  # indexed amino acids
            ILEquivalence = pickle.load(handle) # bool ILEquivalence

        print(" ", flush=True)
        print(f"\t Started peptide mapping for {len(peptideList)} peptides.", flush=True)
        print(f"\t \t Protease: {protease}", flush=True)
        print(f"\t \t Missed cleavages: {MC}", flush=True)
        print(f"\t \t Indexed fasta file: {fastaFile}", flush=True)
        print(f"\t \t Index creation date: {idx_creation_date_time}", flush=True)
        print(f"\t \t Protein entries: {len(id_dict)}", flush=True)
        print(f"\t \t Index len: {str(splitLen)}", flush=True)
        print(f"\t \t Indexed amino acids: {str(aaAlphabet)}", flush=True)

        # load fasta index if not provided as function parameter
        if fasta_idx == None:
            with open(fastaIndex, 'rb') as handle:
                fasta_idx = pickle.load(handle)
                print(f"\t \t Loaded fasta index", flush=True)

        print(f"\t \t Mapping...", flush=True)

        result_list = list()
        result_list.append(f'peptide\tprotein\tlocation\tprevAA\tin_fasta\tnextAA\n')
        # loop through peptide list and map to all possible proteins
        for pept in peptideList:
            if len(pept) < splitLen:
                print(f"\t \t \t WARNING: Peptide '{pept}' is shorter than index len and will be removed. Please check digestion settings.", flush=True)
            else:
                # handle ILEquivalence on-the fly
                if ILEquivalence:
                    pept = pept.replace('I', 'L')

                pept_idx = pept[0:splitLen]
                idx_mapping_list = fasta_idx[pept_idx]

                for idx_mapping in idx_mapping_list:
                    # split individual mappings to key and position
                    idx_mapping_split = idx_mapping.split(', ')
                    idx_key = str(idx_mapping_split[0])
                    idx_pos = int(idx_mapping_split[1])

                    # correct for zero-position
                    idx_start = idx_pos - 1
                    idx_end = idx_start + len(pept)

                    retrieved_seq = id_seq_dict[idx_key][idx_start:idx_end]
                    if retrieved_seq == pept:

                        # output format should be:
                        # header line, 6 cols: peptide \t protein \t location \t prevAA \t in_fasta \t nextAA
                        # body:                AAAAFRVVK \t lcl|AL009126.3_prot_2464 \t 19 \t \t \t
                        result_list.append(f'{pept}\t{id_dict[idx_key]}\t{idx_pos}\t\t\t\n')
                        #print(f"Mapped: '{pept}' to {id_dict[idx_key]} starting at position {idx_pos} and ending at position {idx_end}.")

                    else:
                        #print(f"\t '{pept}' is not equal to '{retrieved_seq}'")
                        pass

        return result_list


if __name__ == "__main__":
    main()
