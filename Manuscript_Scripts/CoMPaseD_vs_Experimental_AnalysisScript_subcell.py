# imports
import argparse
import os
import pandas
import warnings
import time
import colorama
from pathlib import Path

colorama.init()


# class from: https://stackoverflow.com/questions/29484443/python-argparse-insert-blank-line-between-help-entries
# which helps formating help with additional empty lines between the arguments
# "formatter_class=BlankLinesHelpFormatter" in parser = argparse.ArgumentParser (... belongs to this class.
class BlankLinesHelpFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
        return super()._split_lines(text, width) + ['']


# parse arguments from the command line
parser = argparse.ArgumentParser(description="Prepare MS2PIP input from fasta file.",
                                 formatter_class=BlankLinesHelpFormatter)

parser.add_argument('--insilico',
                    help="Input in silico file.",
                    type=str,
                    required=False)

parser.add_argument('--insilicolist',
                    help="Input in silico files from a .txt file. Provide one line per file.",
                    type=str,
                    required=False)

parser.add_argument('--experimental',
                    help="Input experimental file folder. File names should contain all proteases used in this analysis (possible options: 'chymotrypsin', 'gluc', 'lysc', 'lysarginase', 'trypsin').",
                    type=str,
                    required=True)

parser.add_argument('--out_file',
                    help="Output file name. \n",
                    type=str,
                    required=True)

parser.add_argument('--prot_weight',
                    help="Protein ID weighting factor (default = 1.0)",
                    type=float,
                    required=False,
                    default= 1.0)

parser.add_argument('--pep_weight',
                    help="Peptide ID weighting factor (default = 1.0)",
                    type=float,
                    required=False,
                    default= 1.0)

parser.add_argument('--cov_weight',
                    help="Protein coverage weighting factor (default = 1.0)",
                    type=float,
                    required=False,
                    default= 1.0)

parser.add_argument('--grouping',
                    help="Protein grouping bins (default =  0,50,100,99999), i.e. protein length in amino acids for proteins belonging to one group. Provide multiple integer values, SPACE separated. E.g.: '--grouping 0 50 100 99999'",
                    type=int,
                    nargs='+',
                    required=False,
                    default= [0,50,100,99999])

parser.add_argument('--positive_tag',
                    help="String common in all target sequences (default: None for no filtering; usage e.g. 'BACSU' for Bacillus subtils). \n",
                    type=str,
                    default=None,
                    required=False)

parser.add_argument('--decoy_tag',
                    help="Decoy prefeix string (default: 'rev_'). \n",
                    type=str,
                    default='rev_',
                    required=False)

parser.add_argument('--iRT_tag',
                    help="String common to all iRT sequences (default: 'Biognosys'). \n",
                    type=str,
                    default='Biognosys',
                    required=False)

parser.add_argument('--format',
                    help='Format output and calculate correlation between experimental and in-silico result(s)',
                    action='store_true',
                    default=False,
                    required=False)

parser.add_argument('--format_file',
                   help='File containing folder names organised as matrix. Must contain header in first row and index in first column. Correlations will be sorted in additional formatted output accordingly. Only valid when "--format" was provided.',
                   default=None,
                   required=False)

parser.add_argument('--alternative_grouping',
                    help='Use alternative grouping indicated by annotation file.',
                    action='store_true',
                    default=False,
                    required=False)

parser.add_argument('--grouping_file',
                   help='File containing a column with protein identifiers and a grouping column. Columns must be named as "Identifier" and "Grouping", respectively. Will be effective only when --alternative_grouping flag was set.',
                   default=None,
                   required=False)


args = parser.parse_args()

# import single result file...
if (not args.insilico is None) and (args.insilicolist is None):
    if os.path.isfile(os.path.join(args.insilico)):
        # test file and check format roughly?
        in_silico_ls = [os.path.join(args.insilico)]
    else:
        print(
            f"{colorama.Fore.RED}Error: in silico file ({args.insilico}) is no a valid file. Please check.{colorama.Style.RESET_ALL}")
        raise OSError
# ...or import list of result files
elif (not args.insilicolist is None) and (args.insilico is None):
    try:
        # open file and store all files as list
        in_silico_ls = list()
        with open(args.insilicolist, 'r') as f:
            tmp_lst = [line.rstrip() for line in f]
        for ln in tmp_lst:
            if os.path.isfile(os.path.join(ln)):
                in_silico_ls.append(os.path.join(ln))
            else:
                print(
                    f"{colorama.Fore.LIGHTBLUE_EX}Warning: in silico file ({ln}) is no a valid file. Will not use this file. Please check.{colorama.Style.RESET_ALL}")
        if len(in_silico_ls) < 1:
            print(f"{colorama.Fore.RED}Error: No in silico file found. Please check.{colorama.Style.RESET_ALL}")
            raise OSError
    except:
        print(f"{colorama.Fore.RED}Error: Could not open in silico file list. Please check.{colorama.Style.RESET_ALL}")
        raise OSError
# avoid ambiguity when user provides both
else:
    print(
        f"{colorama.Fore.RED}Error: Please provide either, --insilico (one file) or --insilicolist (file containing multiple in silico result files).{colorama.Style.RESET_ALL}")
    raise OSError

# get all file names in experimental folder and filter for tsv files
if os.path.exists(os.path.join(args.experimental)):
    exp_ls = os.listdir(args.experimental)
    exp_ls = [fl for fl in exp_ls if fl.endswith(".tsv")]
    if len(exp_ls) < 1:
        print(f"{colorama.Fore.RED}Error: No experimental result file (must have extension '.tsv') found. Please check.{colorama.Style.RESET_ALL}")
        raise OSError
else:
    print(f"{colorama.Fore.RED}Error: Could not find experimental result folder ({args.experimental}). Please check.{colorama.Style.RESET_ALL}")
    raise OSError

# enzyme names in file names and CoMPaseD result table, must be ordered identically
exp_valid_enzymes = ('chymotrypsin', 'gluc', 'lysc', 'lysarginase', 'trypsin')
insilico_valid_enzymes = ('chymotrypsin', 'glu-c', 'lys-c', 'lysarginase', 'trypsin')

# extract proteases from experimental file names and convert to id-string
exp_contains_ls = list()
for fl in exp_ls:
    fl_edit = fl.lower()
    # handle 'trypsin-contained-in-chymotrypsin' ambiguity
    fl_edit = fl_edit.replace('chymotrypsin', 'chymTrps')
    exp_valid_enzymes_edit = ('chymTrps', 'gluc', 'lysc', 'lysarginase', 'trypsin')

    curr_fl_contains = list()
    for idx, enz in enumerate(exp_valid_enzymes_edit):
        if enz in fl_edit:
            curr_fl_contains.append("1")
        else:
            curr_fl_contains.append("0")

    exp_contains_ls.append(''.join(curr_fl_contains))

# open experimental files and collect data for calculation of protease scores
exp_dict = dict()   # contains file-wise sub-dicts with protease_combination and group-wise values for score calc.
# result lists
grp_lst = list()
grp_lower_lst = list()
grp_upper_lst = list()
combin_lst = list()
score_unfiltered_lst = list()
score_filtered_lst = list()

# load annotation file once at the beginning if required
if args.alternative_grouping:
    # read alternative annotation
    if os.path.isfile(os.path.join(args.grouping_file)):
        annot_file = pandas.read_csv(os.path.join(args.grouping_file), sep='\t')
    else:
        print(
            f"{colorama.Fore.RED}Error: in grouping annotation file ({args.grouping_file}) is no a valid file. Please check.{colorama.Style.RESET_ALL}")
        raise OSError

    # quickly check format
    if not 'Identifier' in annot_file.columns:
        print(
            f"{colorama.Fore.RED}Error: in grouping annotation file. Column 'Identifier' not found. Please check.{colorama.Style.RESET_ALL}")
        raise OSError
    if not 'Grouping' in annot_file.columns:
        print(
            f"{colorama.Fore.RED}Error: in grouping annotation file. Column 'Identifier' not found. Please check.{colorama.Style.RESET_ALL}")
        raise OSError

    # list of unique group categories
    group_lst = list(set(annot_file.Grouping))

    # shorten annotation table to minimum and remove whitespace behind identifier
    annot_file = annot_file[["Identifier", "Grouping"]]
    annot_file["Identifier"] = annot_file.Identifier.str.rstrip()

for idx, fl in enumerate(exp_ls):
    # read experimental file
    curr_exp_res = pandas.read_csv(os.path.join(args.experimental, fl), sep = '\t')

    # filter out reverse and iRT entries
    curr_exp_res = curr_exp_res.loc[~curr_exp_res.Protein.str.startswith(args.decoy_tag)].reset_index(drop=True)
    curr_exp_res = curr_exp_res.loc[~curr_exp_res.Protein.str.contains(args.iRT_tag)].reset_index(drop=True)

    # filter everything that is not of BACSU (Bacillus subtilis) or user-set origin
    if not args.positive_tag is None:
        curr_exp_res = curr_exp_res.loc[curr_exp_res.Protein.str.contains(args.positive_tag)].reset_index(drop=True)
    else:
        print("not filtered any organism")

    key_list = list()
    val_list = list()
    key_list.append('protease_combination')
    val_list.append(exp_contains_ls[idx])

    if not args.alternative_grouping:
        for idx2, grp in enumerate(args.grouping):
            if idx2 < len(args.grouping)-1:

                # filter results for current grp interval
                curr_exp_res_grp = curr_exp_res[curr_exp_res.Length.between(args.grouping[idx2]+1, args.grouping[idx2+1]+1)].reset_index(drop=True)

                # print(f"Min len: {curr_exp_res_grp.Length.min()}; \t Max len: {curr_exp_res_grp.Length.max()}")

                # get key - value pairs for different properties and store them in lists
                tmp_var_name = 'group_' + str(eval(str(idx2))) + '_lower_bound'
                key_list.append(tmp_var_name)
                val_list.append(args.grouping[idx2]+1)

                tmp_var_name = 'group_' + str(eval(str(idx2))) + '_upper_bound'
                key_list.append(tmp_var_name)
                val_list.append(args.grouping[idx2+1])

                tmp_var_name = 'group_' + str(eval(str(idx2))) + '_coverage_unfiltered'
                key_list.append(tmp_var_name)
                val_list.append(curr_exp_res_grp['Percent Coverage'].mean())

                tmp_var_name = 'group_' + str(eval(str(idx2))) + '_protein_ids_unfiltered'
                key_list.append(tmp_var_name)
                val_list.append(curr_exp_res_grp['Percent Coverage'].count())

                tmp_var_name = 'group_' + str(eval(str(idx2))) + '_peptide_ids_unfiltered'
                key_list.append(tmp_var_name)
                val_list.append(curr_exp_res_grp['Stripped Peptides'].sum())

                curr_exp_res_grp = curr_exp_res_grp.loc[curr_exp_res_grp['Stripped Peptides'] > 1]

                tmp_var_name = 'group_' + str(eval(str(idx2))) + '_coverage_filtered'
                key_list.append(tmp_var_name)
                val_list.append(curr_exp_res_grp['Percent Coverage'].mean())

                tmp_var_name = 'group_' + str(eval(str(idx2))) + '_protein_ids_filtered'
                key_list.append(tmp_var_name)
                val_list.append(curr_exp_res_grp['Percent Coverage'].count())

                tmp_var_name = 'group_' + str(eval(str(idx2))) + '_peptide_ids_filtered'
                key_list.append(tmp_var_name)
                val_list.append(curr_exp_res_grp['Stripped Peptides'].sum())
    else:
        # annotate grouping
        curr_exp_res = curr_exp_res.merge(annot_file, how='left',
                           left_on='Protein',
                           right_on='Identifier', copy=None)

        # filter rows without grouping information
        curr_exp_res = curr_exp_res[curr_exp_res.Grouping.isin(group_lst)]

        for idx2, grp in enumerate(group_lst):
            # assuming that group list is allways present
            curr_exp_res_grp = curr_exp_res[curr_exp_res.Grouping == grp]

            tmp_var_name = grp + '_coverage_unfiltered'
            key_list.append(tmp_var_name)
            val_list.append(curr_exp_res_grp['Percent Coverage'].mean())

            tmp_var_name = grp + '_protein_ids_unfiltered'
            key_list.append(tmp_var_name)
            val_list.append(curr_exp_res_grp['Percent Coverage'].count())

            tmp_var_name = grp + '_peptide_ids_unfiltered'
            key_list.append(tmp_var_name)
            val_list.append(curr_exp_res_grp['Stripped Peptides'].sum())

            curr_exp_res_grp = curr_exp_res_grp.loc[curr_exp_res_grp['Stripped Peptides'] > 1]

            tmp_var_name = grp + '_coverage_filtered'
            key_list.append(tmp_var_name)
            val_list.append(curr_exp_res_grp['Percent Coverage'].mean())

            tmp_var_name = grp + '_protein_ids_filtered'
            key_list.append(tmp_var_name)
            val_list.append(curr_exp_res_grp['Percent Coverage'].count())

            tmp_var_name = grp + '_peptide_ids_filtered'
            key_list.append(tmp_var_name)
            val_list.append(curr_exp_res_grp['Stripped Peptides'].sum())

    # convert paired lists to dict
    grp_dict = {key_list[idx3]: val_list[idx3] for idx3 in range(len(key_list))}

    # put grp_dict to exp_dict using fl as key and
    exp_dict[fl] = grp_dict

# calculate all scores
combination_ls = set(exp_contains_ls)
# extract trypsin dicts
exp_dict_try = dict()
for exp_fl in exp_ls:
    if exp_dict[exp_fl]['protease_combination'] == '00001':
        exp_dict_try[exp_fl] = exp_dict[exp_fl]

# calculate and store results:
for exp_fl in exp_ls:
    if not args.alternative_grouping:
        for idx2, grp in enumerate(args.grouping):
            if idx2 < len(args.grouping) - 1:
                for exp_try in exp_dict_try.keys():
                    curr_enz_code = exp_dict[exp_fl]['protease_combination']
                    key_str = 'group_' + str(eval(str(idx2))) + '_lower_bound'
                    grp_lower = exp_dict[exp_fl][key_str]
                    key_str = 'group_' + str(eval(str(idx2))) + '_upper_bound'
                    grp_upper = exp_dict[exp_fl][key_str]

                    key_str = 'group_' + str(eval(str(idx2))) + '_coverage_unfiltered'
                    cov = exp_dict[exp_fl][key_str]
                    cov_try = exp_dict_try[exp_try][key_str]

                    key_str = 'group_' + str(eval(str(idx2))) + '_protein_ids_unfiltered'
                    prot_id = exp_dict[exp_fl][key_str]
                    prot_id_try = exp_dict_try[exp_try][key_str]

                    key_str = 'group_' + str(eval(str(idx2))) + '_peptide_ids_unfiltered'
                    pept_id = exp_dict[exp_fl][key_str]
                    pept_id_try = exp_dict_try[exp_try][key_str]

                    score_unfiltered = (((prot_id/prot_id_try) ** args.prot_weight) * ((pept_id/pept_id_try) ** args.pep_weight) * ((cov/cov_try) ** args.cov_weight)) ** (1 / (args.prot_weight + args.pep_weight + args.cov_weight))

                    key_str = 'group_' + str(eval(str(idx2))) + '_coverage_filtered'
                    cov = exp_dict[exp_fl][key_str]
                    cov_try = exp_dict_try[exp_try][key_str]

                    key_str = 'group_' + str(eval(str(idx2))) + '_protein_ids_filtered'
                    prot_id = exp_dict[exp_fl][key_str]
                    prot_id_try = exp_dict_try[exp_try][key_str]

                    key_str = 'group_' + str(eval(str(idx2))) + '_peptide_ids_filtered'
                    pept_id = exp_dict[exp_fl][key_str]
                    pept_id_try = exp_dict_try[exp_try][key_str]

                    score_filtered = (((prot_id / prot_id_try) ** args.prot_weight) * (
                                (pept_id / pept_id_try) ** args.pep_weight) * (
                                                    (cov / cov_try) ** args.cov_weight)) ** (
                                                   1 / (args.prot_weight + args.pep_weight + args.cov_weight))

                    # grp_list idx must be corrected by +1 to fit indices from in-silico results
                    grp_lst.append('group_' + str(eval(str(idx2+1))))
                    grp_lower_lst.append(grp_lower)
                    grp_upper_lst.append(grp_upper)
                    combin_lst.append(curr_enz_code)
                    score_unfiltered_lst.append(score_unfiltered)
                    score_filtered_lst.append(score_filtered)
    else:
        for idx2, grp in enumerate(group_lst):
            for exp_try in exp_dict_try.keys():
                curr_enz_code = exp_dict[exp_fl]['protease_combination']
                key_str = grp + '_coverage_unfiltered'
                cov = exp_dict[exp_fl][key_str]
                cov_try = exp_dict_try[exp_try][key_str]

                key_str = grp + '_protein_ids_unfiltered'
                prot_id = exp_dict[exp_fl][key_str]
                prot_id_try = exp_dict_try[exp_try][key_str]

                key_str = grp + '_peptide_ids_unfiltered'
                pept_id = exp_dict[exp_fl][key_str]
                pept_id_try = exp_dict_try[exp_try][key_str]

                score_unfiltered = (((prot_id / prot_id_try) ** args.prot_weight) * (
                            (pept_id / pept_id_try) ** args.pep_weight) * ((cov / cov_try) ** args.cov_weight)) ** (
                                               1 / (args.prot_weight + args.pep_weight + args.cov_weight))

                key_str = grp + '_coverage_filtered'
                cov = exp_dict[exp_fl][key_str]
                cov_try = exp_dict_try[exp_try][key_str]

                key_str = grp + '_protein_ids_filtered'
                prot_id = exp_dict[exp_fl][key_str]
                prot_id_try = exp_dict_try[exp_try][key_str]

                key_str = grp + '_peptide_ids_filtered'
                pept_id = exp_dict[exp_fl][key_str]
                pept_id_try = exp_dict_try[exp_try][key_str]

                score_filtered = (((prot_id / prot_id_try) ** args.prot_weight) * (
                        (pept_id / pept_id_try) ** args.pep_weight) * (
                                          (cov / cov_try) ** args.cov_weight)) ** (
                                         1 / (args.prot_weight + args.pep_weight + args.cov_weight))

                # grp_list must be group string
                grp_lst.append(grp)
                combin_lst.append(curr_enz_code)
                score_unfiltered_lst.append(score_unfiltered)
                score_filtered_lst.append(score_filtered)
if not args.alternative_grouping:
    exp_result_df = pandas.DataFrame({'Protein group': grp_lst,
                              'group_lower_bond': grp_lower_lst,
                              'group_upper_bond': grp_upper_lst,
                              'combination_code': combin_lst,
                              'Protease score (unfiltered)': score_unfiltered_lst, 'Protease score (filtered)': score_unfiltered_lst})
else:
    exp_result_df = pandas.DataFrame({'Protein group': grp_lst,
                                      'combination_code': combin_lst,
                                      'Protease score (unfiltered)': score_unfiltered_lst,
                                      'Protease score (filtered)': score_unfiltered_lst})

exp_result_df.sort_values(by=['Protein group', 'combination_code', 'Protease score (filtered)', 'Protease score (unfiltered)'], ascending=[True,True,False,False], inplace=True)

# calculate stats,
# mad in pandas is deprecated as it was implemented as 'mean absolute deviation of the mean' see https://github.com/pandas-dev/pandas/issues/11787;
# thus use median_abs_deviation from scipy.stats instead
from scipy.stats import median_abs_deviation
exp_result_df_summary = exp_result_df.groupby(['Protein group', 'combination_code'])[['Protease score (filtered)', 'Protease score (unfiltered)']].aggregate(['mean', 'std', 'count', 'median', median_abs_deviation]).reset_index()

# remove multi-index from columns (level with e.g. 'Protease score' and level with e.g. 'mean')
exp_result_df_summary.columns = exp_result_df_summary.columns.get_level_values(0)

# assign new column header
exp_col_name_lst = ['Protein group',
                    'combination_code',
                    'mean experimental Protease score (filtered)',
                    'stdev experimental Protease score (filtered)',
                    'count experimental Protease score (filtered)',
                    'median experimental Protease score (filtered)',
                    'mad experimental Protease score (filtered)',
                    'mean experimental Protease score (unfiltered)',
                    'stdev experimental Protease score (unfiltered)',
                    'count experimental Protease score (unfiltered)',
                    'median experimental Protease score (unfiltered)',
                    'mad experimental Protease score (unfiltered)']
exp_result_df_summary.columns = exp_col_name_lst

#exp_columns_1 = [''.join(col) for col in exp_result_df_summary.columns.values[0:2]]
#exp_columns_2 = ['_'.join(col) for col in exp_result_df_summary.columns.values[2:]]
#exp_columns = exp_columns_1 + exp_columns_2
#exp_result_df_summary.columns = exp_columns

# proceed with the in-silico results
for idx4, fl in enumerate(in_silico_ls):
    print(f"\t current in-silico file: '{fl}'")
    # read in-silico file
    curr_insilico_res = pandas.read_csv(os.path.join(fl), sep = '\t')

    # handle 'trypsin-contained-in-chymotrypsin' ambiguity
    curr_insilico_res['protease_edit'] = curr_insilico_res['Protease combination'].str.replace('chymotrypsin', 'chymTrps')
    insilico_valid_enzymes_edit = ('chymTrps', 'glu-c', 'lys-c', 'lysarginase', 'trypsin')

    # convert protease column to list to use for loop
    insilico_contains_ls = list()
    tmp_protease_lst = curr_insilico_res['protease_edit'].to_list()

    # first loop through each result row
    for res in tmp_protease_lst:
        insilico_contains = list()
        # than loop through enzyme list
        for idx, enz in enumerate(insilico_valid_enzymes_edit):
            if enz in res:
                insilico_contains.append("1")
            else:
                insilico_contains.append("0")

        # append 5-char string as enzyme code to list
        insilico_contains_ls.append(''.join(insilico_contains))

    # return list to df, this df is similar to exp_result_df
    curr_insilico_res['combination_code'] = insilico_contains_ls
    curr_insilico_res.sort_values(by=['Protein group', 'combination_code', 'Protease score (filtered)', 'Protease score (unfiltered)'],
                              ascending=[True, True, False, False], inplace=True)
    # calculate stats; define mad to ignore nan's
    def mad(x):
        return median_abs_deviation(x.dropna(), nan_policy='omit')

    agg_funcs = ['mean', 'std', 'count', 'median', mad]
    # calculate stats
    insilico_result_df_summary = curr_insilico_res.groupby(['Protein group', 'combination_code']).agg({
        'Protease score (filtered)': agg_funcs,
        'Protease score (unfiltered)': agg_funcs
    }).reset_index()

    # Flatten the MultiIndex columns
    insilico_result_df_summary.columns = ['_'.join(col).strip() if isinstance(col, tuple) else col for col in
                                          insilico_result_df_summary.columns]
    '''
    # calculate stats
    insilico_result_df_summary = curr_insilico_res.groupby(['Protein group', 'combination_code'])[[
        'Protease score (filtered)', 'Protease score (unfiltered)']].aggregate(
        ['mean', 'std', 'count', 'median', median_abs_deviation]).reset_index()
    '''
    # remove multi-index from columns (level with e.g. 'Protease score' and level with e.g. 'mean')
    # first two columns do not have multi-index
    #insilico_columns_1 = [''.join(col) for col in insilico_result_df_summary.columns.values[0:2]]
    #insilico_columns_2 = ['_'.join(col) for col in insilico_result_df_summary.columns.values[2:]]
    #insilico_columns = insilico_columns_1 + insilico_columns_2
    #insilico_result_df_summary.columns = insilico_columns

    # use nested way to get parent folder name only;
    # requires usage of pathlibs relative_to combined with os.path.dirname
    # takes the path of fl relative to dirname(dirname(fl)) and returns the dirname w/o file name
    fl_path = Path(fl)
    curr_insilico_name = os.path.dirname(fl_path.relative_to(os.path.dirname(os.path.dirname(fl))))

    # remove multi-index from columns (level with e.g. 'Protease score' and level with e.g. 'mean')
    insilico_result_df_summary.columns = insilico_result_df_summary.columns.get_level_values(0)

    # assign new column header
    insilico_col_name_lst = ['Protein group',
                        'combination_code',
                        f'mean {curr_insilico_name} Protease score (filtered)',
                        f'stdev {curr_insilico_name} Protease score (filtered)',
                        f'count {curr_insilico_name} Protease score (filtered)',
                        f'median {curr_insilico_name} Protease score (filtered)',
                        f'mad {curr_insilico_name} Protease score (filtered)',
                        f'mean {curr_insilico_name} Protease score (unfiltered)',
                        f'stdev {curr_insilico_name} Protease score (unfiltered)',
                        f'count {curr_insilico_name} Protease score (unfiltered)',
                        f'median {curr_insilico_name} Protease score (unfiltered)',
                        f'mad {curr_insilico_name} Protease score (unfiltered)']
    insilico_result_df_summary.columns = insilico_col_name_lst

    # merge with other results or init results
    if idx4 == 0:
        insilico_result_df_summary_final = insilico_result_df_summary
    elif idx4 >=1:
        insilico_result_df_summary_final = pandas.merge(left=insilico_result_df_summary_final, right=insilico_result_df_summary,
                                             on=['combination_code', 'Protein group'], how='left')

# merge in-silico results into experimental ones
exp_result_df_summary = pandas.merge(left=exp_result_df_summary, right=insilico_result_df_summary_final, on=['combination_code', 'Protein group'], how='left')

# annotate final table
# read first in-silico file and repeat protease annotation
curr_insilico_res = pandas.read_csv(os.path.join(in_silico_ls[0]), sep = '\t')
# handle 'trypsin-contained-in-chymotrypsin' ambiguity
curr_insilico_res['protease_edit'] = curr_insilico_res['Protease combination'].str.replace('chymotrypsin', 'chymTrps')
insilico_valid_enzymes_edit = ('chymTrps', 'glu-c', 'lys-c', 'lysarginase', 'trypsin')
# convert protease column to list to use for loop
insilico_contains_ls = list()
tmp_protease_lst = curr_insilico_res['protease_edit'].to_list()
# first loop through each result row
for res in tmp_protease_lst:
    insilico_contains = list()
    # than loop through enzyme list
    for idx, enz in enumerate(insilico_valid_enzymes_edit):
        if enz in res:
            insilico_contains.append("1")
        else:
            insilico_contains.append("0")

    # append 5-char string as enzyme code to list
    insilico_contains_ls.append(''.join(insilico_contains))
# return list to df, this df is similar to exp_result_df
curr_insilico_res['combination_code'] = insilico_contains_ls

# generate annotation df
annot_res = curr_insilico_res[['Protease combination', 'combination_code']].drop_duplicates().reset_index(drop=True)
annot_res['Protease number'] = annot_res['Protease combination'].str.split(' - ').str.len()
# merge into results
exp_result_df_summary = pandas.merge(left=exp_result_df_summary, right=annot_res, on=['combination_code'], how='left')

# export
os.makedirs(os.path.dirname(args.out_file), exist_ok=True)
# use NA for missing values for compatibility during plotting in R
exp_result_df_summary.to_csv(args.out_file, sep='\t', index=False, na_rep='NA')

# ToDo: add formatted output with R2 and probably matrix-wise for different conditions (template table with file names?)
if args.format:
    def formatting_funct(res_df=exp_result_df_summary, formatting_file=args.format_file, what='mean', how='filtered'):
        # check options
        if what not in ('mean', 'count', 'median'):
            print(f"{colorama.Fore.RED}Error in formatting_funct: Parameter 'what' must be any of ('mean', 'count', 'median') but is currently {what}.{colorama.Style.RESET_ALL}")
            raise ValueError()
        if how not in ('filtered', 'unfiltered'):
            print(f"{colorama.Fore.RED}Error in formatting_funct: Parameter 'how' must be any of ('filtered', 'unfiltered') but is currently {how}.{colorama.Style.RESET_ALL}")
            raise ValueError()

        # construct column name
        exp_result_col = what +  " experimental Protease score (" + how + ")"

        # get insilico column names
        insilico_col_lst = list()
        # handle filtered is contained in unfiltered issue:
        how_modified = "(" + how + ")"
        for res_col in res_df.columns:
            if (what in res_col) and (how_modified in res_col) and (not 'experimental' in res_col):
                insilico_col_lst.append(res_col)

        grps = set(res_df['Protein group'])

        correl_lst = list()
        grp_lst = list()
        res_lst = list()
        for grp in grps:
            res_df_grp = res_df.loc[res_df['Protein group'] == grp]
            for res_col in insilico_col_lst:
                correl_lst.append(res_df_grp[exp_result_col].corr(res_df_grp[res_col], method='pearson'))
                grp_lst.append(grp)
                res_lst.append(res_col)

        correl_df = pandas.DataFrame({'Condition': res_lst,
                                      'Group': grp_lst,
                                      'Pearson Correlation': correl_lst})
        if formatting_file is None:
            return [correl_df], None

        else:
            # get formatting info and align data accordingly
            if os.path.isfile(os.path.join(formatting_file)):
                # header=0 and index_col=0 use first row/col as herder/index and keep folder names only as table content
                formatting_df = pandas.read_csv(formatting_file, sep="\t", header=0, index_col=0)
                import numpy

                formatted_res_df_lst = list()
                for grp in grps:
                    correl_df_grp = correl_df.loc[correl_df['Group'] == grp].reset_index(drop=True)
                    formatted_res_df = formatting_df.copy(deep=True)

                    res_lst_grp = correl_df_grp['Condition'].to_list()
                    correl_lst_grp = correl_df_grp['Pearson Correlation'].to_list()

                    for p_c, rs in zip(correl_lst_grp, res_lst_grp):
                        rs_lookup = rs.replace(what+" ", "")
                        rs_lookup = rs_lookup.replace(" Protease score (" + how + ")", "")

                        rw, cl = numpy.where(formatted_res_df == rs_lookup)
                        if (not rw is None) and (not cl is None):
                            formatted_res_df.iloc[rw, cl] = p_c
                    formatted_res_df_lst.append(formatted_res_df)

                return formatted_res_df_lst, grps

    res_lst_mean_filtered, grps_lst_mean_filtered = formatting_funct(res_df=exp_result_df_summary,
                                                                      formatting_file=args.format_file,
                                                                      what='mean', how='filtered')
    res_lst_mean_unfiltered, grps_lst_mean_unfiltered = formatting_funct(res_df=exp_result_df_summary,
                                                                          formatting_file=args.format_file,
                                                                          what='mean', how='unfiltered')

    res_lst_median_filtered, grps_lst_median_filtered = formatting_funct(res_df=exp_result_df_summary,
                                                                          formatting_file=args.format_file,
                                                                          what='median', how='filtered')
    res_lst_median_unfiltered, grps_lst_median_unfiltered = formatting_funct(res_df=exp_result_df_summary,
                                                                      formatting_file=args.format_file, what='median',
                                                                      how='unfiltered')

    def export_funct(res_lst, grp_lst, what='mean', how='filtered', res_path = os.path.dirname(args.out_file)):
        if grp_lst is None:
            f_name = os.path.join(res_path, f"{what}_{how}.tsv")
            res_lst[0].to_csv(f_name, sep='\t', index=False, header=True, na_rep='NA')
            print(res_lst[0])
            return None
        else:
            print("in else export")
            print(grp_lst)
            print(len(grp_lst))

            if not len(res_lst) == len(grp_lst):
                print(f"{colorama.Fore.RED}Error in export_funct: Length of result groups ({len(grp_lst)}) and number of result tables ({len(res_lst)}) for {what},{how} does not fit.{colorama.Style.RESET_ALL}")
                raise ValueError()

            for res, grp in zip(res_lst, grp_lst):
                print(grp)
                # result file name and path
                f_name = os.path.join(res_path, f"{what}_{how}_{grp}.tsv")
                res.to_csv(f_name, sep='\t', index=True, header= True, na_rep='NA')
            return None

    export_funct(res_lst_mean_filtered, grps_lst_mean_filtered, 'mean', 'filtered')
    export_funct(res_lst_mean_unfiltered, grps_lst_mean_unfiltered, 'mean', 'unfiltered')
    export_funct(res_lst_median_filtered, grps_lst_median_filtered, 'median', 'filtered')
    export_funct(res_lst_median_unfiltered, grps_lst_median_unfiltered, 'median', 'unfiltered')

print("finished")
