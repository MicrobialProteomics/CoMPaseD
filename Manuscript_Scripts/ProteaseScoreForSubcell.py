import argparse
from numpy import mean, median, sum
from Bio import Seq
from collections import namedtuple

from multiprocessing import cpu_count, Pool
from itertools import combinations
from os import environ, rename, path, makedirs
from datetime import datetime
from time import perf_counter
from Bio import SeqIO
from pandas import read_csv, merge, DataFrame

# define classes
class ProteinClass(Seq.Seq):
    """
    Class providing extentended functionality for proteins;
    inherits from biopythons sequence class with amino acids
    """

    def __init__(self, identifier, sequence):

        self.id = identifier
        self.seq = sequence
        self.length = len(self.seq)
        self.peps = []
        self.pep_pos = []
        self.pep_length = []
        self.coverage = 0


    def add_pep(self, peptide, peptide_position, peptide_length):
        """
        Assign peptides to a protein; this variant might be minimally faster
        as peptide length is obtained from input and not calculated
        """
        self.peps.append(peptide)
        self.pep_pos.append(peptide_position)
        self.pep_length.append(peptide_length)


    def add_pep_variant_2(self, peptide, peptide_position):
        """
        Assign peptides to a protein; this variant might be minimally slower
        as peptide length is obtained from input and not calculated
        as peptide length is calculated for each peptide
        """
        self.peps.append(peptide)
        self.pep_pos.append(peptide_position)
        self.pep_length.append(len(peptide))


    def calcCoverage(self):
        """
        Calculate sequence coverage for a particular protein,
        based on all currently assigned peptides
        """
        sequence_copy = list(self.seq)
        for each_position, each_length in zip(self.pep_pos, self.pep_length):
            each_position_corrected = each_position - 1
            for ii in range(each_position_corrected, (each_position_corrected+each_length)):
                sequence_copy[ii] = "z"
        "".join(sequence_copy)
        covered_seq = sequence_copy.count("z")

        calc_coverage = covered_seq / self.length
        self.coverage = calc_coverage
        return calc_coverage


class CoMPaseD_results():
    """
    Class to hold results ;
    inherits from biopythons sequence class with amino acids
    """

    def __init__(self, protease_combination: list, replicate: str, group: str, min_peps_per_prot = 2):
        '''
        if len(protease_combination) > 1:
            self.combination = " - ".join(protease_combination)
        elif len(protease_combination) == 1:
            self.combination = str(protease_combination)
        '''

        self.combination = " - ".join(protease_combination)

        self.random_sampling = replicate
        self.group = group
        self.min_peps_per_prot = min_peps_per_prot

        self.number_proteins = int()
        self.number_peptides_total = int()
        self.number_peptides_mean = float()
        self.number_peptides_median = float()
        self.coverage_mean = float()
        self.coverage_median = float()

        self.number_proteins_filtered = int()
        self.number_peptides_total_filtered = int()
        self.number_peptides_mean_filtered = float()
        self.number_peptides_median_filtered = float()
        self.coverage_mean_filtered = float()
        self.coverage_median_filtered = float()

        self.score = float()
        self.score_filtered = float()


    def get_results(self, protein_list, update_coverage = False):
        """get results based on list of ProteinClass objects, set update_coverage to True if coverage was not calculated before"""
        self.get_number_proteins(protein_list)
        self.get_number_peptides(protein_list)
        if update_coverage:
            self.get_coverage_result(protein_list, calc_coverage=True)
        else:
            self.get_coverage_result(protein_list, calc_coverage=False)


    def get_number_proteins(self, protein_list):
        """count number of proteins with one or min_peps_per_prot peptides"""
        total = int()
        filtered = int()

        for prot in protein_list:
            if len(prot.peps) >= 1:
                total += 1
            if len(prot.peps) >= self.min_peps_per_prot:
                filtered += 1

        # update result obj
        self.number_proteins = total
        self.number_proteins_filtered = filtered


    def get_number_peptides(self, protein_list):
        """calculate number of peptides per protein, and total peptides"""
        total = list()
        filtered = list()

        for prot in protein_list:
            if len(prot.peps) >= 1:
                total.append(len(prot.peps))
            if len(prot.peps) >= self.min_peps_per_prot:
                filtered.append(len(prot.peps))

        # update result obj
        self.number_peptides_total = sum(total)
        self.number_peptides_mean = mean(total)
        self.number_peptides_median = median(total)
        self.number_peptides_total_filtered = sum(filtered)
        self.number_peptides_mean_filtered = mean(filtered)
        self.number_peptides_median_filtered = median(filtered)


    def get_coverage_result(self, protein_list, calc_coverage = False):
        """calculate average and median protein coverage"""
        total = list()
        filtered = list()

        for prot in protein_list:
            # update coverage if calc_coverage = True
            if calc_coverage:
                prot.calcCoverage()

            if len(prot.peps) >= 1:
                total.append(prot.coverage)
            if len(prot.peps) >= self.min_peps_per_prot:
                filtered.append(prot.coverage)

        # update result obj
        self.coverage_mean = mean(total)
        self.coverage_median = median(total)
        self.coverage_mean_filtered = mean(filtered)
        self.coverage_median_filtered = median(filtered)



def main():
    parser = argparse.ArgumentParser(description="run CoMPaseD score calculation on existing RandomSampling output file")
    parser.add_argument('-r','--RandomSampling', required=True, help="path to RandomSampling.tsv file")
    parser.add_argument('-o', '--Output_directory', required=True, help="path to Output directory")
    parser.add_argument('-f', '--Fasta', required=True, help="path to fasta file")
    parser.add_argument('--max_sampling_number', required=False, default=10)
    parser.add_argument('--Protein_IDs_weight', required=False, default=1.0)
    parser.add_argument('--Peptide_IDs_weight', required=False, default=1.0)
    parser.add_argument('--Coverage_weight', required=False, default=1.0)

    # get start time
    time_0 = perf_counter()

    # parse arguments
    args = parser.parse_args()

    # check existing param file
    if not path.isfile(args.RandomSampling):
        # raised errors are not displayed on progress tab when running from gui,
        # thus print and raise error afterwards to output some useful error info
        print(f"RandomSampling file not found, please save file under {args.RandomSampling}")
        raise FileNotFoundError(f"RandomSampling file not found, please save file under {args.RandomSampling}")

    # read file to df
    else:
        try:
            pep_df = read_csv(path.join(args.RandomSampling), sep='\t')
        except Exception as e:
            print(f"ERROR: {e}")
            print(f"ERROR: Could not open {path.join(args.RandomSampling)}. \n Stopped!")
            raise RuntimeError(f"ERROR: Could not open {path.join(args.RandomSampling)}. \n Stopped!")

    # find possible protease combinations to analyse
    max_proteases = len(set(pep_df.Enzyme))
    proteases = set(pep_df.Enzyme)
    combin_list = get_protease_combinations(protease_list=proteases,
                                            max_proteases=max_proteases)
    # get groups
    tmp_grouping_df = pep_df[["Group", "protein"]]
    tmp_grouping_df = tmp_grouping_df.drop_duplicates()
    groups_list = list(set(tmp_grouping_df["Group"].to_list()))
    groups_list.sort()

    if not path.isfile(args.Fasta):
        # raised errors are not displayed on progress tab when running from gui,
        # thus print and raise error afterwards to output some useful error info
        print(f"Fasta file not found, please save file under {args.Fasta}")
        raise FileNotFoundError(f"Fasta file not found, please save file under {args.Fasta}")
    else:
        fasta = args.Fasta

    max_sampling_number = args.max_sampling_number

    # analyse results and store as nested list of CoMPaseD_results objs
    result_list = list()
    for curr_group in groups_list:

        print(f"Started sampling result analysis for group {curr_group}", flush=True)

        # subset pep_df by protein group
        group_pep_df = pep_df.loc[pep_df["Group"] == curr_group]

        # counter for protease combinations
        n = 1
        tot_n = len(combin_list)

        # generate analyse_sampling argument tuple for all protease combinations
        analysis_args_list = list()
        for combin in combin_list:
            analysis_args = (group_pep_df, combin, curr_group, max_sampling_number, fasta, n, tot_n)
            analysis_args_list.append(analysis_args)
            n += 1

        # use starmap from multiprocessing.Pool to run analysis in parallel
        # pool_n = multiprocessing.cpu_count() - 1
        pool_n = cpu_count() - 1
        # with multiprocessing.Pool(pool_n) as analysis_pool:
        with Pool(pool_n) as analysis_pool:
            analysis_result_list = analysis_pool.starmap(func=analyse_sampling, iterable=analysis_args_list)

        # append list with result objs from all protease combinations
        result_list.append(analysis_result_list)
        print(f"Finished sampling result analysis for group {curr_group}", flush=True)

    # reduce dimensionality of result_list
    # a) unpack combination level
    result_list = [lst for sublist_result in result_list for lst in sublist_result]
    # b) unpack random_sampling level
    result_list = [lst for sublist_result in result_list for lst in sublist_result]

    # find trypsin results, extract and remove from result_list
    trypsin_result = list()
    trypsin_idx_list = list()
    for idx, result in enumerate(result_list):
        if str(result.combination) == "trypsin":
            trypsin_result.append(result)
            trypsin_idx_list.append(idx)
    for idx in sorted(trypsin_idx_list, reverse=True):
        del result_list[idx]

    # extract raw results
    res_df_trypsin = make_result_df(trypsin_result)
    res_df_other = make_result_df(result_list)

    # merge trypsin results by random sampling
    final_res_df = merge(left=res_df_other, right=res_df_trypsin,
                                on=['Random sampling', 'Protein group'], suffixes=['', ' trypsin'], how='left')

    # calculate ratios used for scores
    final_res_df['Protein ID ratio (unfiltered)'] = final_res_df['Total proteins identified (unfiltered)'] / \
                                                    final_res_df['Total proteins identified (unfiltered) trypsin']
    final_res_df['Protein ID ratio (filtered)'] = final_res_df['Total proteins identified (filtered)'] / final_res_df[
        'Total proteins identified (filtered) trypsin']
    final_res_df['Protein ID weight'] = float(args.Protein_IDs_weight)

    final_res_df['Peptide ID ratio (unfiltered)'] = final_res_df['Total number of peptides identified (unfiltered)'] / \
                                                    final_res_df[
                                                        'Total number of peptides identified (unfiltered) trypsin']
    final_res_df['Peptide ID ratio (filtered)'] = final_res_df['Total number of peptides identified (filtered)'] / \
                                                  final_res_df['Total number of peptides identified (filtered) trypsin']
    final_res_df['Peptide ID weight'] = float(args.Peptide_IDs_weight)

    final_res_df['Protein coverage ratio (unfiltered)'] = final_res_df['Mean protein coverage (unfiltered)'] / \
                                                          final_res_df['Mean protein coverage (unfiltered) trypsin']
    final_res_df['Protein coverage ratio (filtered)'] = final_res_df['Mean protein coverage (filtered)'] / final_res_df[
        'Mean protein coverage (filtered) trypsin']
    final_res_df['Protein coverage weight'] = float(args.Coverage_weight)

    # calculate scores and insert as fourth and fifth column
    final_res_df.insert(3, 'Protease score (unfiltered)', (
            (final_res_df['Protein ID ratio (unfiltered)'] ** final_res_df['Protein ID weight']) * (
            final_res_df['Peptide ID ratio (unfiltered)'] ** final_res_df['Peptide ID weight']) * (
                    final_res_df['Protein coverage ratio (unfiltered)'] ** final_res_df[
                'Protein coverage weight'])) ** (1 / (
            final_res_df['Protein ID weight'] + final_res_df['Peptide ID weight'] + final_res_df[
        'Protein coverage weight'])))

    final_res_df.insert(4, 'Protease score (filtered)', (
            (final_res_df['Protein ID ratio (filtered)'] ** final_res_df['Protein ID weight']) * (
            final_res_df['Peptide ID ratio (filtered)'] ** final_res_df['Peptide ID weight']) * (
                    final_res_df['Protein coverage ratio (filtered)'] ** final_res_df[
                'Protein coverage weight'])) ** (1 / (
            final_res_df['Protein ID weight'] + final_res_df['Peptide ID weight'] + final_res_df[
        'Protein coverage weight'])))

    # generate output file name and save
    final_res_df_file_name = path.join(args.Output_directory, "CoMPaseD_results.tsv")

    if not path.isdir(path.join(args.Output_directory)):
        makedirs(path.join(args.Output_directory))

    # if file exists, try to rename existing file with last modification date and time
    if path.isfile(final_res_df_file_name):
        mti = datetime.fromtimestamp(path.getmtime(final_res_df_file_name))
        rename_f_name = path.join(args.Output_directory, mti.strftime("%Y-%m-%d_%Hh%Mmin%Ssec_CoMPaseD_results.tsv"))
        try:
            rename(final_res_df_file_name, rename_f_name)
        except Exception as e:
            print(f"ERROR: Could not rename existing file {final_res_df_file_name}. File will be overwritten.")
            print(f"{e}")

    print(f"Write results to {final_res_df_file_name}")
    final_res_df.to_csv(final_res_df_file_name, sep='\t', index=False)

def make_result_df(result_list):
    """Combine CoMPaseD result objects that are stored as a list to a pandas df"""

    # init result df
    result_col_names = ['Protease combination',
                        'Protein group',
                        'Random sampling',
                        'Total proteins identified (unfiltered)',
                        'Total number of peptides identified (unfiltered)',
                        'Mean number of peptides per protein identified (unfiltered)',
                        'Median number of peptides per protein identified (unfiltered)',
                        'Mean protein coverage (unfiltered)',
                        'Median protein coverage (unfiltered)',
                        'Min peptides per protein for filtering',
                        'Total proteins identified (filtered)',
                        'Total number of peptides identified (filtered)',
                        'Mean number of peptides per protein identified (filtered)',
                        'Median number of peptides per protein identified (filtered)',
                        'Mean protein coverage (filtered)',
                        'Median protein coverage (filtered)'
                        ]
    tmp_res_df = DataFrame(columns=result_col_names)

    # generate list with results for each {group / protease combination / random sampling} combination
    for res in result_list:
        tmp_res_list = list()
        tmp_res_list.append(str(res.combination))
        tmp_res_list.append(str(res.group))
        tmp_res_list.append(str(res.random_sampling))
        tmp_res_list.append(int(res.number_proteins))
        tmp_res_list.append(int(res.number_peptides_total))
        tmp_res_list.append(float(res.number_peptides_mean))
        tmp_res_list.append(float(res.number_peptides_median))
        tmp_res_list.append(float(res.coverage_mean))
        tmp_res_list.append(float(res.coverage_median))
        tmp_res_list.append(int(res.min_peps_per_prot))
        tmp_res_list.append(int(res.number_proteins_filtered))
        tmp_res_list.append(int(res.number_peptides_total_filtered))
        tmp_res_list.append(float(res.number_peptides_mean_filtered))
        tmp_res_list.append(float(res.number_peptides_median_filtered))
        tmp_res_list.append(float(res.coverage_mean_filtered))
        tmp_res_list.append(float(res.coverage_median_filtered))

        # add list as last row to tmp_res_df
        tmp_res_df.loc[len(tmp_res_df)] = tmp_res_list
    return tmp_res_df


def analyse_sampling(pep_df, protease_combin, curr_group, max_sampling_number, fasta, n, tot_n) -> list:
    """Analyse results for one (combination of) protease(s)"""

    # convert potential list to string prior print
    combination_str = " - ".join(protease_combin)
    print(f"\t Started analysing protease combination {n} of {tot_n} ({combination_str})", flush=True)

    # result_list
    combin_result_list = list()
    # subset pep_df
    tmp_df = pep_df.loc[pep_df["Enzyme"].isin(protease_combin)]

    # number of random sampling is taken from params to allow fewer samplings than available in df
    sampling_col_list = list()
    for i in range(1, int(max_sampling_number) + 1):
        sampling_col_list.append("sampling_" + str(i))

    # generate list of "identified" proteins and fill protein objs with peptides for each random_sampling
    for sampling_col in sampling_col_list:
        protein_list = makeProteinList(SeqIO.parse(fasta, "fasta"))
        tmp_df_smp = tmp_df.loc[tmp_df[sampling_col] == 1, ["peptide", "protein", "location"]].reset_index(drop=True)
        protein_list = fillProteinList(protein_list, tmp_df_smp)

        # after each sampling generate new result obj and put to list
        tmp_result = CoMPaseD_results(protease_combin, sampling_col, curr_group, min_peps_per_prot=2)
        tmp_result.get_results(protein_list, update_coverage=True)
        combin_result_list.append(tmp_result)

    print(f"\t Finished analysing protease combination {n} of {tot_n}", flush=True)

    return combin_result_list



def makeProteinList(fasta_file):
    """
    Import proteins from fasta file; the returned list has empty containers
    for each protein to hold corresponding peptides, etc.
    """
    # make a list of all proteins in the fasta, this can be filled with the corresponding peptides later
    protein_list = list()
    for record in fasta_file:
        tmp = ProteinClass(record.id, record.seq)
        protein_list.append(tmp)
    return protein_list


def fillProteinList(protein_list_to_fill, fill_df):
    """
    Faster variant that uses a dict to fill protein_list
    """
    # sort df by protein and extract lists with peptide information, moved to Analysis_MPD, analyse_results function
    # outside the loop, thus only done once fill_df.sort_values(by="protein", inplace=True)
    protein_col = fill_df["protein"]
    location_col = fill_df["location"]
    peptide_col = fill_df["peptide"]
    # convert/combine lists to tuple and than to namedtuple
    peptide_tuples = list(zip(protein_col, location_col, peptide_col))
    tuple_list = namedtuple('entry', 'protein location pepseq')
    tuple_list = [tuple_list(*el) for el in peptide_tuples]
    # generate dict from tuple_list
    tmp_peps = {}
    tmp_loc = {}
    for protein, location, pepseq in tuple_list:
        if protein not in tmp_peps:
            tmp_peps[protein] = list()
            tmp_loc[protein] = list()
        tmp_peps[protein].append(pepseq)
        tmp_loc[protein].append(location)
    # loop through protein list and fill peptides by dict-key
    for prot in protein_list_to_fill:
        # do not try to fill proteins without any peptide to avoid 'NoneType' error
        if prot.id in tmp_peps.keys():
            for peptide_seq, pep_location in zip(tmp_peps.get(prot.id), tmp_loc.get(prot.id)):
                prot.add_pep_variant_2(peptide_seq, pep_location)
    return protein_list_to_fill


def get_protease_combinations(protease_list: list, max_proteases: int) -> list:
    """Obtain all possible protease combinations up to max_protease items"""

    # convert list to tuple for usage with itertools
    protease_list = tuple(protease_list)

    # fix cases where more proteases shall be combined than available
    if max_proteases >= len(protease_list) + 1:
        max_proteases = len(protease_list)

    # calculate all combinations with 1 to N proteases
    combination_list = list()
    for ii in range(1, max_proteases + 1):
        # tmp_combinations = list(itertools.combinations(protease_list, ii))
        tmp_combinations = list(combinations(protease_list, ii))
        # itertools.combinations returns a list of tuples
        # append each tuple to combination_list
        for protease_ecombination in tmp_combinations:
            combination_list.append(protease_ecombination)
    # convert tuple elements in combination_list to lists and return
    return [list(combin) for combin in combination_list]

if __name__ == "__main__":
    main()
