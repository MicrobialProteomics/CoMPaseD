import argparse
import colorama
from multiprocessing import cpu_count, Pool
from itertools import combinations
from os import environ, rename
from datetime import datetime
from time import perf_counter
from Bio import SeqIO
from pandas import read_csv, merge

try:
    from lib.CoMPaseD_gui_param_functions import *
except ModuleNotFoundError:
    from CoMPaseD_gui_param_functions import *

try:
    from lib.CoMPaseD_protein_class import *
except ModuleNotFoundError:
    from CoMPaseD_protein_class import *

# disable tensorflow warnings / info during import and reset to default
environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
try:
    from lib.CoMPaseD_DMSP import *
except ModuleNotFoundError:
    from CoMPaseD_DMSP import *
environ['TF_CPP_MIN_LOG_LEVEL'] = '0'


def main():
    parser = argparse.ArgumentParser(description="run CoMPaseD analysis functions")
    parser.add_argument('--param_file', required=True, help="path to parameter file")
    # cmd line arguments to perfom analysis only:
    parser.add_argument('--use_existing_sampling_output', required=False,
                        help="set to use pre-computed sampling output table, requires to set "
                             "sampling_output_path as well", action='store_true')
    parser.add_argument('--sampling_output_path', required=False, help="path to pre-computed sampling_output file",
                        default="This is not a path")
    parser.add_argument('--digestion_result', required=False, help="set unique peptides table from in-silico digestion", default="")

    # get start time
    time_0 = perf_counter()

    # output to progress tab:
    print(f"CoMPaseD - Comparison of Multiple-Protease Digestions", flush=True)
    print("", flush=True)
    print("", flush=True)
    print("---------------------------------------------------------------------------", flush=True)
    print("Analysis of In-silico digestion started", flush=True)

    # parse arguments
    args = parser.parse_args()

    # ensure absolute path in case relative path is provided on cmd line
    if not os.path.isabs(args.param_file):
        param_file = os.path.abspath(args.param_file)

    # check existing param file
    if not path.isfile(args.param_file):
        # raised errors are not displayed on progress tab when running from gui,
        # thus print and raise error afterward to output some useful error info
        print(f"Parameter file not found, please save file under {args.param_file}")
        raise FileNotFoundError(f"Parameter file not found, please save file under {args.param_file}")

    # load params
    params = CoMPaseD_Parameter()
    params.load_params(path.join(args.param_file))

    # check output dir
    if not path.isdir(path.join(params.Output_directory)):
        print("Output folder does not exist")
        raise NotADirectoryError("Output folder does not exist")

    # check protein weight file
    if not path.isfile(params.Protein_weight_file):
        print(f"Protein weight file not found, please save file.")
        raise FileNotFoundError(f"Protein weight file not found, please save file.")

    # check previous digestion
    if not path.join(args.digestion_result) == "":
        if path.isfile(path.join(args.digestion_result)):
            digest_file = path.join(args.digestion_result)

        elif path.isfile(path.abspath(args.digestion_result)):
            digest_file = path.abspath(args.digestion_result)
        else:
            print(f"Digestion result file ({args.digestion_result}) not found. \n"
                                    f"Did you forgot to digest?")
            raise FileNotFoundError(f"Digestion result file ({args.digestion_result}) not found. \n"
                                    f"Did you forgot to digest?")
    else:
        digest_file = path.join(params.Output_directory, "unique_peptides_table_filtered.tsv")
        if not path.isfile(digest_file):
            print(f"Digestion result file ({digest_file}) not found. \n"
                                    f"Did you forgot to digest?")
            raise FileNotFoundError(f"Digestion result file ({digest_file}) not found. \n"
                                    f"Did you forgot to digest?")

    # skip this part if --use_existing_sampling_output was set on cmd
    if args.use_existing_sampling_output:
        if path.isfile(path.join(args.sampling_output_path)):
            try:
                pep_df = read_csv(path.join(args.sampling_output_path), sep='\t')
            except Exception as e:
                print(f"{colorama.Fore.RED}ERROR: Could not open {path.join(args.sampling_output_path)} due to {e}. Stopping.{colorama.Style.RESET_ALL}")
                raise RuntimeError

    # execute random sampling when --use_existing_sampling_output was not set on cmd
    else:
        # import files
        pwf_df = read_csv(path.join(params.Protein_weight_file), sep='\t')
        pep_df = read_csv(digest_file, sep='\t')

        # merge raw sampling weights into pep_df
        pep_df = merge(left=pep_df, right=pwf_df,
                       left_on="protein", right_on="Identifier",
                       how="left").reset_index()

        # get protease / mc combinations
        protease_mc_df = protease_mc_expansion(params.Proteases, params.Max_MCs)
        # add mc freq
        freq_mc = get_numeric_list(params.Freq_MCs)

        if not len(freq_mc) == len(protease_mc_df):
            print(f"{colorama.Fore.RED}ERROR: Invalid protease table. Please ensure that in every row the number of values provided for 'MC frequency' equals 'Max MCs + 1' (there must be a frequency for peptides without any MC site). Stopping.{colorama.Style.RESET_ALL}")
            raise RuntimeError

        protease_mc_df['MC_Freq'] = freq_mc
        # correct mc freq sum to 1
        protease_mc_df = normalise_mc(protease_mc_df)
        # count peptides per mc
        protease_mc_df = get_pep_counts(pep_df, protease_mc_df)
        # calculate sampling size
        protease_mc_df = get_peps_required(protease_mc_df, params)  # Enzyme, MC, sampling_size

        # add DeepMSPeptide Prediction values
        # use_DeepMSPeptide_Predictions is string, not bool
        if params.Use_DeepMSPeptide_Predictions == "True":
            if path.isfile(path.join(params.Path_DeepMSPeptide_Model)):
                pep_df = predict_detectability(pep_df,
                                               path.join(params.Path_DeepMSPeptide_Model),
                                               float(params.Weights_DeepMSPeptide_Predictions))
            else:
                print(f"{colorama.Fore.RED}ERROR: No valid Deep-MS-Peptide model file provided. Please check: {path.join(params.Path_DeepMSPeptide_Model)}. Stopping.{colorama.Style.RESET_ALL}")
                raise FileNotFoundError
        else:
            pep_df['DeepMSPep_prediction'] = 1
        # modify random sampling columns by DeepMSPeptide prediction
        pep_df = multiply_dmsp(pep_df)

        # generate col with subset keys to allow subsetting by only one col
        pep_df["subset"] = pep_df["Enzyme"].astype(str) + "__" + pep_df["MC"].astype(str)

        # find random sampling columns
        rand_sampling_cols = list()
        for pep_df_col in pep_df.columns:
            if pep_df_col.lower().startswith("random_sampling_"):
                rand_sampling_cols.append(pep_df_col)

        # get pep_df index to assign sampled peptides later
        pep_df = pep_df.reset_index(drop=True)
        pep_df['ID'] = pep_df.index

        # loop through rand_sampling_cols and sample peptides
        # returning list and assigning samples peptides by ID requires less memory than merging within each loop
        smp_col_list = list()
        for smp_col in rand_sampling_cols:
            print(f"Started {smp_col}", flush=True)
            tmp_id_list = rand_smp(pep_df, protease_mc_df, column_to_sample=smp_col)
            smp_col_name = smp_col.lower().replace("random_sampling_", "sampling_")

            # init sampling_N col with 0 and replace with 1 where ID was obtained from rand_smp
            pep_df[smp_col_name] = 0
            pep_df.loc[pep_df["ID"].isin(tmp_id_list), smp_col_name] = 1
            smp_col_list.append(smp_col_name)

        # remove unused peptides to reduce file / df size
        pep_df['pep_used'] = pep_df[smp_col_list].sum(axis=1)
        pep_df = pep_df.loc[pep_df['pep_used'] > 0]
        pep_df.drop(['pep_used'], axis=1, inplace=True)
        pep_df.reset_index(drop=True, inplace=True)

        # convert 0/1 columns to int8 to save memory
        pep_df[smp_col_list] = pep_df[smp_col_list].astype('int8')

        # write output if required
        if params.Sampling_output:
            sampling_out_file = path.join(params.Output_directory, "RandomSampling.tsv")

            # if file exists, try to rename existing file with last modification date and time
            if path.isfile(sampling_out_file):
                mti = datetime.fromtimestamp(path.getmtime(sampling_out_file))
                rename_f_name = path.join(params.Output_directory, mti.strftime("%Y-%m-%d_%Hh%Mmin%Ssec_RandomSampling.tsv"))
                try:
                    rename(sampling_out_file, rename_f_name)
                except Exception as e:
                    print(f"{colorama.Fore.CYAN}WARNING: Could not rename existing file {sampling_out_file} due to {e}. \n File will be overwritten.{colorama.Style.RESET_ALL}")

            print("Started writing sampling output table", flush=True)
            pep_df.to_csv(sampling_out_file, sep="\t", index=False)
            print("Finished writing sampling output table", flush=True)

    # find possible protease combinations to analyse
    combin_list = get_protease_combinations(protease_list=params.Proteases,
                                            max_proteases=int(params.Number_of_Proteases))

    # get groups
    tmp_grouping_df = pep_df[["Group", "protein"]]
    tmp_grouping_df = tmp_grouping_df.drop_duplicates()
    groups_list = list(set(tmp_grouping_df["Group"].to_list()))
    groups_list.sort()

    # analyse results and store as nested list of CoMPaseD_results objs
    result_list = list()
    for curr_group in groups_list:

        print(f"Started sampling result analysis for group {curr_group}", flush=True)

        # subset pep_df by protein group
        group_pep_df = pep_df.loc[pep_df["Group"] == curr_group]

        # counter for protease combinations
        n = 1
        tot_n = len(combin_list)

        '''
        # variant w/o multiprocessing, used for development:
        result_list = list()
        for curr_group in groups_list:
            # subset pep_df by protein group
            group_pep_df = pep_df.loc[pep_df["Group"] == curr_group]
            n = 1
            tot_n = len(combin_list)
            for combin in combin_list:
                tmp_result_list = analyse_sampling(group_pep_df, combin, curr_group, params, n, tot_n)
                result_list.append(tmp_result_list)
                n+=1
            print(f"Finished sampling result analysis for group {curr_group}", flush=True)
        print(f"Finished sampling result analysis for all groups", flush=True)
        result_list = [result_list[0:63], result_list[63:126], result_list[126:189]]
        '''
        # generate analyse_sampling argument tuple for all protease combinations
        analysis_args_list = list()
        for combin in combin_list:
            analysis_args = (group_pep_df, combin, curr_group, params, n, tot_n)
            analysis_args_list.append(analysis_args)
            n += 1

        # use starmap from multiprocessing.Pool to run analysis in parallel
        # pool_n = multiprocessing.cpu_count() - 1
        pool_n = cpu_count() - 1
        # in case only one core is available, use this, otherwise keep one core free for other tasks
        if pool_n == 0:
            pool_n = 1
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
    final_res_df['Protein ID weight'] = float(params.Protein_IDs_weight)

    final_res_df['Peptide ID ratio (unfiltered)'] = final_res_df['Total number of peptides identified (unfiltered)'] / \
                                                    final_res_df[
                                                        'Total number of peptides identified (unfiltered) trypsin']
    final_res_df['Peptide ID ratio (filtered)'] = final_res_df['Total number of peptides identified (filtered)'] / \
                                                  final_res_df['Total number of peptides identified (filtered) trypsin']
    final_res_df['Peptide ID weight'] = float(params.Peptide_IDs_weight)

    final_res_df['Protein coverage ratio (unfiltered)'] = final_res_df['Mean protein coverage (unfiltered)'] / \
                                                          final_res_df['Mean protein coverage (unfiltered) trypsin']
    final_res_df['Protein coverage ratio (filtered)'] = final_res_df['Mean protein coverage (filtered)'] / final_res_df[
        'Mean protein coverage (filtered) trypsin']
    final_res_df['Protein coverage weight'] = float(params.Peptide_IDs_weight)

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
    final_res_df_file_name = path.join(params.Output_directory, "CoMPaseD_results.tsv")

    # if file exists, try to rename existing file with last modification date and time
    if path.isfile(final_res_df_file_name):
        mti = datetime.fromtimestamp(path.getmtime(final_res_df_file_name))
        rename_f_name = path.join(params.Output_directory, mti.strftime("%Y-%m-%d_%Hh%Mmin%Ssec_CoMPaseD_results.tsv"))
        try:
            rename(final_res_df_file_name, rename_f_name)
        except Exception as e:
            print(f"{colorama.Fore.CYAN}WARNING: Could not rename existing file {final_res_df_file_name} due to {e}. \n File will be overwritten.{colorama.Style.RESET_ALL}")

    print(f"Saved results to {final_res_df_file_name}")
    final_res_df.to_csv(final_res_df_file_name, sep='\t', index=False)

    # generate summary output table and save
    agg_res_df = final_res_df.groupby(['Protease combination', 'Protein group'], as_index=False).agg(
        Mean_score_unfiltered=("Protease score (unfiltered)", "mean"),
        SD_score_unfiltered=("Protease score (unfiltered)", "std"),
        Mean_score_filtered=("Protease score (filtered)", "mean"),
        SD_score_filtered=("Protease score (filtered)", "std"))

    # sort by increasing group names and decreasing unfiltered mean protease score
    agg_res_df = agg_res_df.sort_values(by=['Protein group', 'Mean_score_unfiltered'], ascending=[True, False])

    # output file name and save
    agg_res_df_file_name = path.join(params.Output_directory, "CoMPaseD_results_summary.tsv")

    # if file exists, try to rename existing file with last modification date and time
    if path.isfile(agg_res_df_file_name):
        mti = datetime.fromtimestamp(path.getmtime(agg_res_df_file_name))
        rename_f_name = path.join(params.Output_directory, mti.strftime("%Y-%m-%d_%Hh%Mmin%Ssec_CoMPaseD_results_summary.tsv"))
        try:
            rename(agg_res_df_file_name, rename_f_name)
        except Exception as e:
            print(f"{colorama.Fore.CYAN}WARNING: Could not rename existing file {agg_res_df_file_name} due to {e}. \n File will be overwritten.{colorama.Style.RESET_ALL}")

    print(f"Saved results summary to {agg_res_df_file_name}")
    agg_res_df.to_csv(agg_res_df_file_name, sep='\t', index=False)

    # output to progress tab:
    print("", flush=True)
    print("Analysis of In-silico digestion finished", flush=True)
    time_1 = perf_counter()
    print(f"Analysis took {time_1 - time_0:0.1f} seconds", flush=True)
    print("---------------------------------------------------------------------------", flush=True)


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


def analyse_sampling(pep_df, protease_combin, curr_group, params, n, tot_n) -> list:
    """Analyse results for one (combination of) protease(s)"""

    # convert potential protease list to string prior print
    combination_str = " - ".join(protease_combin)
    print(f"Started analysing protease combination {n} of {tot_n} ({combination_str})", flush=True)

    # result_list
    combin_result_list = list()
    # subset pep_df
    tmp_df = pep_df.loc[pep_df["Enzyme"].isin(protease_combin)]

    # number of random sampling is taken from params instead of df to allow fewer samplings than available
    sampling_col_list = list()
    for i in range(1, int(params.Sampling_Number) + 1):
        sampling_col_list.append("sampling_" + str(i))

    # generate list of "identified" proteins and fill protein objs with peptides for each random_sampling
    for sampling_col in sampling_col_list:
        protein_list = makeProteinList(SeqIO.parse(params.Fasta, "fasta"))
        tmp_df_smp = tmp_df.loc[tmp_df[sampling_col] == 1, ["peptide", "protein", "location"]].reset_index(drop=True)
        protein_list = fillProteinList(protein_list, tmp_df_smp)

        # after each sampling generate new result obj and put to list
        tmp_result = CoMPaseD_results(protease_combin, sampling_col, curr_group, min_peps_per_prot=2)
        tmp_result.get_results(protein_list, update_coverage=True)
        combin_result_list.append(tmp_result)

    print(f"\t Finished analysing protease combination {n} of {tot_n}", flush=True)

    return combin_result_list


def protease_mc_expansion(protease_list, mc_list):
    """Generate df with MC/protease pairs, can be accessed later by itertuples"""

    # convert mc_list to int
    mc_list = [int(i) for i in mc_list]
    # expand protease_list and MC_list by each other
    protease_list_expanded = list()
    MC_list_expanded = list()
    for actual_protease, actual_MC in zip(protease_list, mc_list):
        counter = 0
        while counter < (actual_MC + 1):
            protease_list_expanded.append(actual_protease)
            MC_list_expanded.append(counter)
            counter = counter + 1
    return_df = DataFrame()
    return_df['Enzyme'] = protease_list_expanded
    return_df['MC'] = MC_list_expanded
    return return_df


def predict_detectability(pep_df, dmsp_model, dmsp_weight):
    """Predict peptide detectability by DeepMSPeptide and annotate pep_df"""

    # get unique peptide sequences
    peptide_seqs = pep_df['peptide'].unique().tolist()

    print("Started peptide detectability prediction", flush=True)

    # predict
    peptide_detectability = run_DeepMSPep(peptide_seqs, path.join(dmsp_model))

    # modify by scaling to increase relative importance of this prediction, use power for realistic values,
    # factor 4 was tested empirically
    peptide_detectability['DeepMSPep_prediction'] = peptide_detectability['DeepMSPep_prediction'] ** (4*dmsp_weight)
    peptide_detectability['DeepMSPep_prediction'] = peptide_detectability['DeepMSPep_prediction'] / max(
        peptide_detectability['DeepMSPep_prediction'])

    # merge results into peptide_df
    if 'DeepMSPep_prediction' in pep_df.columns:
        pep_df = pep_df.drop('DeepMSPep_prediction', axis=1)
    pep_df = merge(left=pep_df, right=peptide_detectability, on='peptide', how='left')

    return pep_df


def multiply_dmsp(pep_df, rand_sampling_col="random_sampling_", dmsp_col="DeepMSPep_prediction"):
    """Multiply random sampling weights with DeepMSPep prediction"""

    # find random_sampling_ cols and multiply with dmsp prediction
    for df_col_name in pep_df.columns:
        if df_col_name.lower().startswith(rand_sampling_col):
            pep_df[df_col_name] = pep_df[df_col_name] * pep_df[dmsp_col]

    return pep_df


def get_numeric_list(mc_freq_string):
    """Convert missed-cleavage-frequency string to list of float"""

    # force mc_freq_string to string
    if isinstance(mc_freq_string, list):
        tmp = ','.join(mc_freq_string)
        mc_freq_string = tmp
        # and convert string to list of strings, separated by comma
        mc_freq_string = mc_freq_string.split(",")

    # find and replace brackets
    start_pos_list = list()
    end_pos_list = list()
    value_list = list()
    mc_freq_list = list()
    final_mc_freq_list = list()
    for count, mc_freq in enumerate(mc_freq_string):
        if mc_freq.startswith("["):
            start_pos_list.append(count)
            mc_freq = mc_freq.replace("[", "")
        if mc_freq.endswith("]"):
            end_pos_list.append(count)
            mc_freq = mc_freq.replace("]", "")
        value_list.append(float(mc_freq))
    for start, stop in zip(start_pos_list, end_pos_list):
        current_mc_list = value_list[start:stop + 1]
        mc_freq_list.append(current_mc_list)

    # remove white-spaces
    for sublist in mc_freq_list:
        SublistStr = str(sublist).replace(" ", "")
        final_mc_freq_list.append(SublistStr)

    # unlist nested lists
    mc_freq = [float(fl) for sublist in mc_freq_list for fl in sublist]

    return mc_freq


def normalise_mc(df):
    """Calculate and correct sums of missed cleavage site frequencies"""

    # group and sum
    grouped_df = df.groupby('Enzyme')
    summed_df = grouped_df['MC_Freq'].sum()

    # inform if deviation is more than one percent
    for idx, mc_sum in zip(summed_df.index, summed_df):
        if mc_sum < 0.99:
            print(f"WARNING: Frequency of all potential missed cleavage sites for {idx} sums to {mc_sum}. "
                  f"Sum will be normalised to 1.0 during sampling. Please check", flush=True)
        elif mc_sum > 1.01:
            print(f"WARNING: Frequency of all potential missed cleavage sites for {idx} sums to {mc_sum}. "
                  f"Sum will be normalised to 1.0 during sampling. Please check", flush=True)
    # convert to df for merging/naming
    summed_df = summed_df.to_frame('Norm_fact')

    # remove Norm_fact(or) column in case of multiple runs
    if "Norm_fact" in df.columns:
        df = df.drop("Norm_fact", axis=1)

    # merge and correct
    df = merge(left=df, right=summed_df, left_on='Enzyme', right_index=True, how='left')
    df['MC_Freq'] = df['MC_Freq'] / df['Norm_fact']

    return df


def get_pep_counts(df, protease_mc_df):
    """Count peptides per protease and mc in pep_df"""

    grouped_df = df.groupby(['Enzyme', 'MC']).agg({'peptide': ['count']})
    grouped_df.columns = ['pep_count']
    grouped_df = grouped_df.reset_index()
    # merge back
    protease_mc_df = merge(left=protease_mc_df, right=grouped_df, on=['Enzyme', 'MC'], how='left')

    return protease_mc_df


def get_peps_required(protease_mc_df, params: CoMPaseD_Parameter):
    """Calculate number of peptides required for each protease/mc combination"""

    sampling_base = params.Sampling_Size_Based_On

    # for numerical-based sampling sizes multiply with mc fraq
    if sampling_base == "number":
        # get peptide sampling size as numeric list
        pep_number = [int(i) for i in params.Peptides_Sampling_Size]
        proteases = params.Proteases

        # generate df from pep_number and proteases to merge into protease_mc_df
        merge_df = DataFrame()
        merge_df['Enzyme'] = proteases
        merge_df['total_peptides'] = pep_number

        protease_mc_df = merge(left=protease_mc_df, right=merge_df, on='Enzyme', how='left')

        # multiply and convert back to integer
        protease_mc_df['sampling_size'] = protease_mc_df['total_peptides'] * protease_mc_df['MC_Freq']
        protease_mc_df['sampling_size'] = protease_mc_df['sampling_size'].astype('int')
        # remove non-required column
        protease_mc_df = protease_mc_df.drop('total_peptides', axis=1)

        # use at most the available number of peptides
        missing_pep_len = ((protease_mc_df['sampling_size'] - protease_mc_df['pep_count']) > 0).sum()
        if missing_pep_len > 0:
            print("WARNING: Insufficient number of peptides for random sampling available. Please "
                  "check sampling size and fraction of missed cleavage sites.", flush=True)
        protease_mc_df['sampling_size'].where((protease_mc_df['sampling_size'] - protease_mc_df['pep_count']) < 0,
                                              protease_mc_df['pep_count'].astype('int'), inplace=True)

    elif sampling_base == "coverage":
        # get peptide sampling size as numeric list
        pep_frac = [float(f) for f in params.Pep_Level_Proteome_Cov]
        proteases = params.Proteases

        # generate df from pep_number and proteases to merge into protease_mc_df
        merge_df = DataFrame()
        merge_df['Enzyme'] = proteases
        merge_df['frac_peptides'] = pep_frac

        # calculate peptide count by protease
        grouped_df = protease_mc_df.groupby(['Enzyme']).agg({'pep_count': ['sum']})
        grouped_df.columns = ['pep_count']
        grouped_df = grouped_df.reset_index()

        # calculate actual total sampling size by protease
        merge_df = merge(left=merge_df, right=grouped_df, on='Enzyme', how='left')
        merge_df['total_peptides'] = merge_df['frac_peptides'] * merge_df['pep_count']
        merge_df['total_peptides'] = merge_df['total_peptides'].astype('int')
        merge_df = merge_df.drop(['pep_count', 'frac_peptides'], axis=1)

        protease_mc_df = merge(left=protease_mc_df, right=merge_df, on='Enzyme', how='left')

        # multiply and convert back to integer
        protease_mc_df['sampling_size'] = protease_mc_df['total_peptides'] * protease_mc_df['MC_Freq']
        protease_mc_df['sampling_size'] = protease_mc_df['sampling_size'].astype('int')
        # remove non-required column
        protease_mc_df = protease_mc_df.drop('total_peptides', axis=1)

        # use at most the available number of peptides
        missing_pep_len = ((protease_mc_df['sampling_size'] - protease_mc_df['pep_count']) > 0).sum()
        if missing_pep_len > 0:
            print("WARNING: Insufficient number of peptides for random sampling available. Please "
                  "check sampling size and fraction of missed cleavage sites.", flush=True)
        protease_mc_df['sampling_size'].where((protease_mc_df['sampling_size'] - protease_mc_df['pep_count']) < 0,
                                              protease_mc_df['pep_count'].astype('int'), inplace=True)

    # generate col with subset keys to allow subsetting by only one col
    protease_mc_df["subset"] = protease_mc_df["Enzyme"].astype(str) + "__" + protease_mc_df["MC"].astype(str)

    return protease_mc_df


def rand_smp(pep_df, protease_mc_df, column_to_sample="Random_sampling_1"):
    """Randomly sample peptides from pep_df"""

    # shuffle pep_df to ensure randomness
    pep_df = pep_df.sample(frac=1).reset_index(drop=True)

    # generate list for pep_df['ID']
    pep_id_list = list()

    # iterate through all protease / mc combinations and sample peptides
    for row in protease_mc_df.itertuples():
        # subset pep_df for current protease / mc combination
        tmp_df = pep_df.loc[(pep_df["subset"] == getattr(row, "subset")), :].reset_index(drop=True)

        # check number of available peptides
        tmp_sample_size = int(getattr(row, 'sampling_size'))
        tmp_df_size = len(tmp_df.index)
        if not tmp_df_size > tmp_sample_size:
            print(f"{colorama.Fore.RED}ERROR: Could not find enough peptides to sample for {getattr(row, 'subset')}. Please check parameters. Stopping.{colorama.Style.RESET_ALL}")
            print(f"{colorama.Fore.RED}\t Current settings require {tmp_sample_size} peptides to sample but there are less peptides for this category (you should not sample all peptides).{colorama.Style.RESET_ALL}")
            raise RuntimeError

        # sample according to sample size and remove column_to_sample
        tmp_df = tmp_df.sample(n=int(getattr(row, 'sampling_size')), replace=False,
                               weights=column_to_sample).reset_index(drop=True)

        # copy ID to pep_id_list
        pep_id_list.append(tmp_df['ID'])

    # flatten pep_id_list
    # use upper-case "ID" variable as lower-case is part of build-in namespace
    pep_id_list = [ID for sublist in pep_id_list for ID in sublist]

    return pep_id_list


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
        for protease_combination in tmp_combinations:
            combination_list.append(protease_combination)
    # convert tuple elements in combination_list to lists and return
    return [list(combin) for combin in combination_list]


if __name__ == "__main__":
    main()
