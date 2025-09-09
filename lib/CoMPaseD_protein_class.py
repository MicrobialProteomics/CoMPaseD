from numpy import mean, median, sum
from Bio import Seq
from collections import namedtuple


class ProteinClass(Seq.Seq):
    """
    Class providing extentended functionality for proteins;
    inherits from biopythons sequence class with amino acids
    """

    def __init__(self, identifier, sequence):
        super().__init__(sequence)
        self.id = identifier
        self.seq = sequence
        self.length = len(self.seq)
        self.peps = []
        self.pep_pos = []
        self.pep_length = []
        self.coverage = 0

    def __hash__(self):
        return hash(self.id)

    def __eq__(self, other):
        if isinstance(other, ProteinClass):
            return self.id == other.id
        return False


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

    def __init__(self, protease_combination: list, replicate: str, group: str, min_peps_per_prot = 2, use_unique_peps_only = "True"):
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
        self.use_unique_peps_only = use_unique_peps_only

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
        if self.use_unique_peps_only == "True":
            self.get_number_proteins(protein_list)
            self.get_number_peptides(protein_list)
            self.get_coverage_result(protein_list, calc_coverage=update_coverage)
        # whenever the parameter is not explicitly the True string, fall back to grouped version
        else:
            self.get_number_proteins_groups(protein_list)
            self.get_number_peptides_groups(protein_list)
            self.get_coverage_result_groups(protein_list, calc_coverage=update_coverage)

    def get_number_proteins_groups(self, protein_groups):
        total = int()
        filtered = int()
        for gr in protein_groups:
            if any(len(p.peps) >= 1 for p in gr):
                total += 1

            if any(len(p.peps) >= self.min_peps_per_prot for p in gr):
                filtered += 1

        self.number_proteins = total
        self.number_proteins_filtered = filtered



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


    def get_number_peptides_groups(self, protein_groups):
        total = list()
        filtered = list()

        for group in protein_groups:
            group_total = list()
            group_filtered = list()
            for p in group:
                if len(p.peps) >= 1:
                    group_total.append(len(p.peps))
                if len(p.peps) >= self.min_peps_per_prot:
                    group_filtered.append(len(p.peps))


            # treat the median number of peptides for a group as the representative value
            if group_total:
                total.append(median(group_total))
            if group_filtered:
                filtered.append(median(group_filtered))

        # aggregate on current result level (using only one representative value for each group)
        # convert totals back to int as median might have introduced float numbers

        self.number_peptides_total = int(round(sum(total))) if total else 0
        self.number_peptides_mean = mean(total)
        self.number_peptides_median = median(total)
        self.number_peptides_total_filtered = int(round(sum(filtered))) if filtered else 0
        self.number_peptides_mean_filtered = mean(filtered)
        self.number_peptides_median_filtered = median(filtered)


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


    def get_coverage_result_groups(self, protein_groups, calc_coverage=False):
        coverages_total = list()
        coverages_filtered = list()

        for group in protein_groups:
            if calc_coverage:
                for p in group:
                    p.calcCoverage()

            if group:
                group_coverages = [p.coverage for p in group]
                group_cov = median(group_coverages)
                if group_cov:
                    coverages_total.append(group_cov)

                # filter groups in the same way as get_number_proteins_groups
                if any(len(p.peps) >= self.min_peps_per_prot for p in group):
                    coverages_filtered.append(group_cov)

            self.coverage_mean = mean(coverages_total) if coverages_total else 0
            self.coverage_median = median(coverages_total) if coverages_total else 0
            self.coverage_mean_filtered = mean(coverages_filtered) if coverages_filtered else 0
            self.coverage_median_filtered = median(coverages_filtered) if coverages_filtered else 0


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
