import argparse
import colorama
import re
from multiprocessing import Process
from sys import platform
from pandas import read_csv, DataFrame, Series, concat
from subprocess import Popen, PIPE
from time import perf_counter
from os import path, makedirs, chdir, listdir, remove
from shutil import copy, rmtree


def main():
    parser = argparse.ArgumentParser(description="run in-silico digestions using crux toolkit")
    parser.add_argument('--out_folder', required=True)
    parser.add_argument('--fasta', required=True)
    parser.add_argument('--crux_path', required=True)
    parser.add_argument('--clips_path', required=True)
    parser.add_argument('--promast_path', required=True)
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

    # get start time
    time_0 = perf_counter()

    # output to progress tab:
    print("", flush=True)
    print(f"\t\tCoMPaseD - Comparison of Multiple-Protease Digestions", flush=True)
    print("", flush=True)
    print("", flush=True)
    print("---------------------------------------------------------------------------", flush=True)
    print("", flush=True)
    print("Starting In-silico digestion", flush=True)
    print("", flush=True)
    print(f"\tusing Crux mass spectrometry toolkit for digestion", flush=True)
    print(f"\tfor a detailed description please see:", flush=True)
    print(f"\tSean McIlwain, Kaipo Tamura, Attila Kertesz-Farkas, Charles E. Grant, Benjamin Diament, Barbara Frewen,")
    print(f"\tJ. Jeffry Howbert, Michael R. Hoopmann, Lukas Käll, Jimmy K. Eng, Michael, J.MacCoss and William S. Noble:", flush=True)
    print(f"\tCrux: rapid open source protein tandem mass spectrometry analysis.")
    print(f"\tJournal of Proteome Research. 13(10):4488-4491, 2014.", flush=True)
    print(f"\tDOI: 10.1021/pr500741y", flush=True)
    print("", flush=True)

    # switches to remove non-existing dirs later
    tmp_exists = False
    crux_out_exists = False

    # parse out_folder and create if necessary
    args = parser.parse_args()
    tmp_out_folder = path.join(args.out_folder, 'Tmp')
    if not path.isdir(tmp_out_folder):
        makedirs(tmp_out_folder, exist_ok=True)
        print("created output folder", flush=True)
    else:
        tmp_exists = True
        makedirs(tmp_out_folder, exist_ok=True)
        print("WARNING: output folder exists, files may be overwritten", flush=True)

    if path.isdir(path.join(tmp_out_folder, "crux-output")):
       crux_out_exists = True

    # change directory to tmp_out_folder, all crux output will be in this folder
    chdir(tmp_out_folder)

    print("", flush=True)
    print(f"\tStart Digest:", flush=True)
    crux_path = path.join(args.crux_path)
    clips_path = path.join(args.clips_path)
    promast_path = path.join(args.promast_path)
    fasta = path.join(args.fasta)
    out_folder = path.join(args.out_folder)
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
        print("", flush=True)
        print(f"\t\tStarted digestion with {protease} and {mc} missed cleavages", flush=True)
        out = new_proc.stdout.read()
        err = new_proc.stderr.read()
        new_proc.wait()
        out = out.splitlines()
        err = err.splitlines()
        for line in out:
            print(line.decode('utf8'), flush=True)
        for line in err:
            print(line.decode('utf8'), flush=True)

        print(f"\t\tFinished digestion {n} of {total}.", flush=True)

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
    print("---------------------------------------------------------------------------", flush=True)
    print("", flush=True)
    print(f"\tusing ProteoMapper's perl scripts 'clips.pl' and 'promast.pl' for mapping", flush=True)
    print(f"\tfor a detailed description please see:", flush=True)
    print(f"\tLuis Mendoza, Eric W. Deutsch, Zhi Sun, David S. Campbell, David D. Shteynberg and Robert L. Moritz:", flush=True)
    print(f"\tFlexible and fast mapping of peptides to a proteome with ProteoMapper")
    print(f"\tJournal of Proteome Research. 17(12):4337-4344, 2018.", flush=True)
    print(f"\tDOI: 10.1021/acs.jproteome.8b00544", flush=True)
    print("", flush=True)

    # generate fasta index with clips
    def clips_process(clips_cmd):
        new_proc = Popen(clips_cmd, shell=True, stdout=PIPE, stderr=PIPE)
        print(f"\tStart fasta file index generation", flush=True)
        out = new_proc.stdout.read()
        err = new_proc.stderr.read()
        new_proc.wait()
        out = out.splitlines()
        err = err.splitlines()
        for line in out:
            print(line.decode('utf8'), flush=True)
        for line in err:
            print(line.decode('utf8'), flush=True)
        print(f"\tFinished fasta file index generation.", flush=True)
        print("", flush=True)

    # clips needs to run only once to index the fasta file but will be treated in a queue like the other cmds
    clips_proc_list = list()
    # copy fasta for clips to crux-output dir
    chdir(path.join(tmp_out_folder, 'crux-output'))
    tmp_fasta_loc = path.join(path.join(tmp_out_folder, 'crux-output'), path.basename(fasta))
    copy(fasta, tmp_fasta_loc)
    # init clips process
    clips_proc_list.append(
        Process(target=clips_process(clips_cmd=generate_clips_call(tmp_fasta_loc, clips_path))))
    # start multiprocessing queue
    for proc in clips_proc_list:
        proc.start()
    # wait for finishing
    for proc in clips_proc_list:
        proc.join()

    # map peptides with promast
    promast_proc_list = list()

    def promast_process(promast_cmd, n, total, protease, mc):
        new_proc = Popen(promast_cmd, shell=True, stdout=PIPE, stderr=PIPE)
        print(f"\tStarted peptide mapping for {protease} and {mc} missed cleavages", flush=True)
        out = new_proc.stdout.read()
        err = new_proc.stderr.read()
        new_proc.wait()
        out = out.splitlines()
        err = err.splitlines()
        # promast reports via stderr, thus filter there
        for line in out:
            print(line.decode('utf8'), flush=True)
        for line in err:
            out_ln = line.decode('utf8')
            # filter un-necessary information here
            if "Cannot provide peptide mapping context" not in out_ln:
                if not out_ln.startswith("Reading index file header"):
                    print(out_ln, flush=True)
        print(f"\tFinished peptide mapping for digest {n} of {total}.", flush=True)
        print("", flush=True)

    promast_cmd_list, protease_list, mc_list, mapped_file_list = map_peptides(out_folder=path.join(out_folder, 'Tmp', 'crux-output'),
                                                                              tmp_out_folder=path.join(out_folder, 'Tmp'),
                                                                              fasta=tmp_fasta_loc,
                                                                              promast_path=promast_path)

    n = 0
    total = len(promast_cmd_list)

    # init all processes
    for promast_cmd, protease, mc in zip(promast_cmd_list, protease_list, mc_list):
        n += 1
        proc = Process(target=promast_process(promast_cmd, n, total, protease, mc))
        promast_proc_list.append(proc)

    # start multiprocessing queue
    for proc in promast_proc_list:
        proc.start()

    # wait for finishing
    for proc in promast_proc_list:
        proc.join()

    print("", flush=True)
    print("Finished peptide mapping", flush=True)
    print("---------------------------------------------------------------------------", flush=True)
    print("", flush=True)
    # generate list with mapped file names and protease / mc combination
    crux_out_folder = path.join(out_folder, 'Tmp', 'crux-output')
    mapped_file_list_file = path.join(crux_out_folder, 'File_List.tsv')
    mapped_out_file_list = list()
    for mf, pr, mc in zip(mapped_file_list, protease_list, mc_list):
        ln = str(mf) + "\t" + str(pr) + "\t" + str(mc)
        mapped_out_file_list.append(ln)
    # # save this list
    with open(mapped_file_list_file, 'w') as mf:
        mf.write('\n'.join(mapped_out_file_list))
    # remove files from tmp output and keep only 'Mapped_' peptide lists and 'File_List'
    print(f"\tRemoving tmp files", flush=True)
    for f in listdir(crux_out_folder):
        if not f.startswith('Mapped_'):
            if not f.startswith('File_List'):
                remove(path.join(crux_out_folder, f))
    print(f"\tFinished removing tmp files ", flush=True)
    print("---------------------------------------------------------------------------", flush=True)
    print("", flush=True)
    # merge mapped files and annotate table
    print(f"\tStart pooling digests", flush=True)
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
        filtered_df_file = path.join(args.unique_peps_file)
    else:
        filtered_df_file = path.join(out_folder, "unique_peptides_table_filtered.tsv")

    df_filtered_sorted.to_csv(filtered_df_file, index=False, sep="\t")
    print(f"\tFinished pooling digests", flush=True)
    print("---------------------------------------------------------------------------", flush=True)
    print("", flush=True)
    chdir(out_folder)

    # try to remove folders that eventually were created
    if not crux_out_exists:
        rmtree(path.join(tmp_out_folder, "crux-output"), ignore_errors=True)
    if not tmp_exists:
        rmtree(tmp_out_folder, ignore_errors=True)


    time_1 = perf_counter()
    print(f"Digestion took {time_1 - time_0:0.1f} seconds")
    print("---------------------------------------------------------------------------", flush=True)

def clean_protease_names(protease_list: list) -> list:
    '''remove special characters from protease names in preparation of file name creation'''
    protease_list_clean = list()
    for each_protease in protease_list:
        for character in '\\/:*?"<>|,;-[]{}':
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

    if enzyme.lower().startswith("custom"):
        return_code, cleavage_spec = handle_custom_proteases(enzyme)

        if return_code == -1:
            str_help = "{P}"
            raise ValueError(f"Custom proteases must contain exactly two pairs of squared or curly brackets, e.g. 'custom,trypsin,[KR]|{str_help}'. '{enzyme}' is mal-formatted. Please check. Stopping.")
        elif return_code == -2:
            raise ValueError(f"Inconsistent bracket type in first bracket pair of '{enzyme}'. Please check. Stopping.")
        elif return_code == -3:
            raise ValueError(f"Inconsistent bracket type in second bracket pair of '{enzyme}'. Please check. Stopping.")
        elif return_code == 0:
            enzyme_custom_switch = " --custom-enzyme "
            # !Important: Use inner double quotation marks to ensure windows handles the cleavage specificity correctly
            tmp_enzyme = f'"{cleavage_spec}"'

    elif enzyme not in allowed_enzyme:
        raise ValueError(f"Protease not defined. Must be one of {allowed_enzyme}.")

    else:
        enzyme_custom_switch = " --enzyme "
        tmp_enzyme = enzyme

    # generate crux command
    peptide_cutter_command = str(crux_path) + " generate-peptides --overwrite T --min-mass " + str(
        min_mass) + " --max-mass " + str(
        max_mass) + " --min-length " + str(min_len) + " --max-length " + str(max_len) + " --digestion " + str(
        digestion) + " --missed-cleavages " + str(missed_cleavages) + " --clip-nterm-methionine " + str(
        clip_n_term_met) + " --decoy-format " + str(decoy_format) + str(enzyme_custom_switch) + str(
        tmp_enzyme) + " --fileroot " + str(out_folder) + " " + str(fasta)

    return peptide_cutter_command


def generate_clips_call(fasta, clips_path, segment_size=5):
    """
    clips.pl provides further command line options: \n
        -V = do not use PEFF variants (default: use them) \n
        -f = force index file overwriting (default clips.pl: do not; here set to -f in order to avoid errors) \n
	    -I = do not convert I->L (default: convert) \n
	    -A = do not generate all possible keys in index \n

    example clips call:
    C:/TPP/bin/clips2.pl -V -f -s 5 protein.fasta
    """
    if "linux" in platform:
        clips_cmd = "perl " + str(clips_path) + " -V -f -s " + str(segment_size) + " " + str(fasta)
    else:
        clips_cmd = str(clips_path) + " -V -f -s " + str(segment_size) + " " + str(fasta)
    return clips_cmd


def generate_promast_call(fasta, peptide_list, out_name, promast_path, number_cpus=1):
    """
    promast_write.pl provides further command line options: \n
        -c = provide sequence context of mapped peptide \n
        -t = number of threads to use \n
        -n = output name
        ... \n

    example promast call:
    C:/TPP/bin/promast_write.pl -c protein.fasta trimmed_generate-peptides-output.txt
    """
    if "linux" in platform:
        promast_cmd = "perl " + str(promast_path) + " -c -t " + str(number_cpus) + " -n " + str(out_name) + " " + str(
            fasta) + " " + str(peptide_list)
    else:
        promast_cmd = str(promast_path) + " -c -t " + str(number_cpus) + " -n " + str(out_name) + " " + str(
            fasta) + " " + str(peptide_list)
    return promast_cmd


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


def map_peptides(out_folder, tmp_out_folder, fasta, promast_path):
    """
    Use promast_write.pl to map peptide positions in protein.fasta \n\n

    :param out_folder: path.join(args.out_folder, 'Tmp', 'crux-output')
    :param tmp_out_folder: path.join(args.out_folder, 'Tmp')
    :param fasta: path.join(args.fasta)
    :param promast_path: path.join(args.promast_path)
    :return:
    """
    # get crux-output file names
    crux_file_list_file = path.join(tmp_out_folder, 'crux_result_files.tsv')
    crux_file_list = list()
    mapping_crux_file_list = list()
    mapped_crux_file_list = list()
    protease_list = list()
    mc_list = list()
    with open(crux_file_list_file, 'r') as f:
        while (line := f.readline().rstrip()):
            # crux_file_list_file is tabular with file names in column [0]
            tmp_col = line.split(sep='\t')

            # remove MW and protein ID from file and save under modified name
            crux_file = path.join(out_folder, tmp_col[0])
            mapping_crux_file = path.join(out_folder, str('Mapping_' + tmp_col[0]))
            mapped_crux_file = path.join(out_folder, str('Mapped_' + tmp_col[0]))
            df = read_csv(crux_file, delimiter='\t', usecols=[0])  # use pandas read_csv
            df.to_csv(mapping_crux_file, index=False)

            crux_file_list.append(crux_file)
            mapping_crux_file_list.append(mapping_crux_file)
            mapped_crux_file_list.append(mapped_crux_file)
            protease_list.append(str(tmp_col[1]))
            mc_list.append(str(tmp_col[2]))

    promast_cmd_list = list()
    for f, o in zip(mapping_crux_file_list, mapped_crux_file_list):
        promast_cmd_list.append(generate_promast_call(fasta, peptide_list=f, out_name=o, promast_path=promast_path))

    return promast_cmd_list, protease_list, mc_list, mapped_crux_file_list


def handle_custom_proteases(protease_string):
    """generate string for custom enzyme in crux"""

    backup_protease_string = protease_string

    # protease cleavage string can only be in the form []|[], []|{}, {}|[] or {}|{} filled with any letter in between
    # there is no check whether letters are part of bio-alphabet but upper-case is required
    # X indicates any amino acid
    # input would be, e.g. 'custom [X]|[RKD]' for lysarginase and asp-n
    # whitespaces are stripped

    # remove all whitespaces, including e.g. tab which might originate from pasting
    protease_string = "".join(protease_string.split())

    # custom - trypsin + chymotrypsin [FWYLKR]|{P}
    # regex that extracts the opening bracket, letters and closing bracket before and after the pipe as individual groups
    pattern = r'.*?([\[{])([A-Za-z]*)([\]}])\|([\[{])([A-Za-z]*)([\]}])'

    m = re.match(pattern, protease_string)
    if not m:
        # if pattern was not found, there is no way to handle this - directly raise error
        raise ValueError(f"Invalid custom enzyme syntax: {protease_string}")

    pre_open, pre_content, pre_close, post_open, post_content, post_close = m.groups()

    # change amino acids to upper-case and handle empty brackets as any amino acid
    if pre_content:
        pre_content = pre_content.upper()
    elif pre_open == "{":
        print(
            f"WARNING: First bracket pair in {backup_protease_string} is empty. Will be replaced with X for any amino acid.",
            flush=True)
        print("This results in undigested proteins sequences.", flush=True)
        print("Please check if this was intended.", flush=True)
        pre_content = "X"
    else:
        print(
            f"WARNING: First bracket pair in {backup_protease_string} is empty. Will be replaced with X for any amino acid.",
            flush=True)
        print("Please check if this was intended.", flush=True)
        pre_content = "X"

    if post_content:
        post_content = post_content.upper()
    elif post_open == "{":
        print(
            f"WARNING: Second bracket pair in {backup_protease_string} is empty. Will be replaced with X for any amino acid.",
            flush=True)
        print("This results in undigested proteins sequences.", flush=True)
        print("Please check if this was intended.", flush=True)
        post_content = "X"
    else:
        print(
            f"WARNING: Second bracket pair in {backup_protease_string} is empty. Will be replaced with X for any amino acid.",
            flush=True)
        print("Please check if this was intended.", flush=True)
        post_content = "X"

    # check that opening and closing brackets match, return error code otherwise
    bracket_dict = {"[": "]", "{": "}"}
    if bracket_dict.get(pre_open) != pre_close:
        return -2, "", ""
    if bracket_dict.get(post_open) != post_close:
        return -3, "", ""

    return 0, f"{pre_open}{pre_content}{pre_close}|{post_open}{post_content}{post_close}"


if __name__ == "__main__":
    main()
