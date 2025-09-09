import argparse
import colorama
from os import path, makedirs, listdir, remove
from time import strftime
import subprocess
from sys import executable
import shutil

from lib import CoMPaseD_gui_param_functions, CoMPaseD_gui_export_functions

colorama.init()


# class from: https://stackoverflow.com/questions/29484443/python-argparse-insert-blank-line-between-help-entries
# which helps to format help text with additional empty lines between the arguments
# "formatter_class=BlankLinesHelpFormatter" in parser = argparse.ArgumentParser(... belongs to this class.
class BlankLinesHelpFormatter (argparse.HelpFormatter):
    def _split_lines(self, text, width):
        return super()._split_lines(text, width) + ['']


def main():
    parser = argparse.ArgumentParser(description="CoMPaseD - Comparison of Multiple-Protease Digestions",
                                     formatter_class=BlankLinesHelpFormatter)

    # main arguments
    main_args_param = parser.add_argument_group("Parameter File")
    main_args_param.add_argument('-p', '--param_file', required=True, help='path to CoMPaseD parameter file')

    main_args = parser.add_argument_group("CoMPaseD Mode. If none of these is provided, assuming full analysis, i.e. '-e -d -a'")
    main_args.add_argument('-e', '--export', help='export simulated protein abundance values', action='store_true')
    main_args.add_argument('-d', '--digest', help='perform in-silico digest using crux toolkit', action='store_true')
    main_args.add_argument('-a', '--analysis', help='perform analysis from simulated protein abundance and in-silico digestion', action='store_true')

    # digestion arguments
    digestion_args = parser.add_argument_group("Digestion Mode Arguments (required when no mode or -d is provided)")
    digestion_args.add_argument('--use_original_proteomapper', help='use original perl scripts for mapping in-silico digested peptides, this might be slower but requires less memory (try to use when large databases permit usage of python implementation)', action='store_true')
    digestion_args.add_argument('--differentiate_I_L', help='distinguish between peptide variants containing leucine or iso-leucine (default treat as identical)', action='store_false')
    digestion_args.add_argument('--indexing_key_len', help='length in amino acids of the indexing keys for mapping, shorter length result in longer mapping times while longer increase memory load (default = 5, min = 2, max = 6)', default=5, type=int)

    # analysis arguments
    analysis_args = parser.add_argument_group("Analysis Mode Arguments (required when no mode or -a is provided)")
    analysis_args.add_argument('--export_result', help="path to CoMPaseD export result file with simulated protein abundance values and protein group assignment", default="", type=str)
    analysis_args.add_argument('--digestion_result', help="path to CoMPaseD digestion result file ('unique_peptides_table_filtered')", default="", type=str)

    # parameter adjust arguments
    param_args = parser.add_argument_group("Overwrite Parameter File Content (all optional)")
    param_args.add_argument('--out_folder', help='change output directory', default="", type=str)
    param_args.add_argument('--fasta', help='change fasta file', default="", type=str)
    param_args.add_argument('--score_peptide', help="change weight of peptide IDs for protease score calculation, e.g. '--score_peptide 1.5", default=-1.0, type=float)
    param_args.add_argument('--score_protein', help="change weight of protein IDs for protease score calculation, e.g. '--score_protein 1.5", default=-1.0, type=float)
    param_args.add_argument('--score_coverage', help="change weight of protein coverage for protease score calculation, e.g. '--score_coverage 1.5", default=-1.0, type=float)
    param_args.add_argument('--enzymes', help="change enzymes, provide as comma-separated list, e.g. '--enzymes trypsin,lys-c,asp-n,lysarginase'", default="", type=str)
    param_args.add_argument('--mc', help="change maximal missed cleavage sites, provide as comma-separated list in same order as enzymes, e.g. '--mc 2,2,3,4'", default="", type=str)
    param_args.add_argument('--mc_freq', help="change missed cleavage sites frequency, e.g. '--mc_freq [0.7415,0.2090,0.0484],[0.9102,0.0836,0.0058],[0.5620,0.2753,0.1110,0.0419],...'", default="", type=str)
    param_args.add_argument('--num_peps', help="change number of peptides to sample, provide as comma-separated list in same order as enzymes, e.g. '--num_peps 10000,10000,10000,10000'; DO NOT use together with --frac_peps", default="", type=str)
    param_args.add_argument('--frac_peps', help="change fraction of generated peptides to sample, provide as comma-separated list in same order as enzymes, e.g. '--frac_peps 0.685,0.780,0.230'; DO NOT use together with --num_peps", default="", type=str)
    param_args.add_argument('--min_pep_mw', help="change minimal peptide mass (in Da), e.g. '--min_pep_mw 400', default = 400", default=-1, type=int)
    param_args.add_argument('--max_pep_mw', help="change maximal peptide mass (in Da), e.g. '--min_pep_mw 6000', default = 6000", default= -1, type=int)
    param_args.add_argument('--min_pep_len', help="change minimal peptide length (in amino acids), e.g. '--min_pep_len 6', default = 6", default=-1, type=int)
    param_args.add_argument('--max_pep_len', help="change maximal peptide length (in amino acids), e.g. '--min_pep_len 55', default = 55", default=-1, type=int)
    param_args.add_argument('--bins', help="change protein binning, provide as comma-separated list of protein length in amino acids, e.g. '--bins 0,200,999' will group all peptides between 0 and 200 aa into one group and all between 201 and 999 aa into another", default="", type=str)
    param_args.add_argument('--undetectable', help="change undetectable protein fraction in protein bins, provide as comma-separated list, e.g. '--undetectable 40,20' will assume 40 percent undetectable protein in first group and 20 percent in second", default="", type=str)
    param_args.add_argument('--DMSP_weight', help="change weighting factor of Deep-MS-Peptide prediction, set to zero to disable DMSP prediction, e.g. '--DMSP_weight 2.5'", default=-1, type=float)
    param_args.add_argument('--DMSP_model', help="change Deep-MS-Peptide prediction model path, e.g. '--DMSP_model path/to/prediction_model.h5'", default="", type=str)
    param_args.add_argument('--samplings', help="change number of random samplings, e.g. '--samplings 5", default=-1, type=int)
    param_args.add_argument('--dynamic_range', help="change dynamic range of protein abundance, e.g. '--dynamic_range 5.0", default=-1, type=float)
    param_args.add_argument('--use_unique_peptides_only', help="change whether to use only unique peptides or assemble protein groups and consider shared peptides, e.g. '--use_unique_peptides_only False' will assemble protein groups", default="", type=str)

    args = parser.parse_args()

    # check parameter file
    if not path.isfile(path.join(args.param_file)):
        print(f"{colorama.Fore.RED}ERROR: Parameter file not existing. "
              f"Please check: {args.param_file}{colorama.Style.RESET_ALL}")
        raise FileNotFoundError
    else:
        param_obj = CoMPaseD_gui_param_functions.CoMPaseD_Parameter()
        try:
            param_obj.load_params(path.join(args.param_file))
        except PermissionError:
            print(f"{colorama.Fore.RED}ERROR: Could not read parameter file. Check permission for: {args.param_file}{colorama.Style.RESET_ALL}")
            raise PermissionError
        except Exception:
            print(f"{colorama.Fore.RED}ERROR: Could not read parameter file. This is likely due to wrong formatting. Please check. You can generate valid parameter file by using CoMPaseD_gui.py{colorama.Style.RESET_ALL}")
            raise RuntimeError

    # change params according to cmd line options
    # # output folder
    if not args.out_folder == "":
        if not path.isdir(path.join(args.out_folder)):
            makedirs(path.join(args.out_folder))
        param_obj.Output_directory = path.join(args.out_folder)
    # # check that out_dir is empty and warn in case not
    if len(listdir(path.join(param_obj.Output_directory))) > 0:
        print(f"{colorama.Fore.CYAN}WARNING: Output directory {param_obj.Output_directory} is not empty. Files may be overwritten.{colorama.Style.RESET_ALL}")

    # # fasta file
    if not args.fasta == "":
        if path.isfile(path.join(args.fasta)):
            param_obj.Fasta = path.join(args.fasta)
        else:
            print(f"{colorama.Fore.CYAN}WARNING: Fasta file {args.fasta} does not exist. Keeping fasta file provided in parameter file ({param_obj.Fasta}).{colorama.Style.RESET_ALL}")
            if not path.isfile(path.join(param_obj.Fasta)):
                print(f"{colorama.Fore.RED}ERROR: Could not find any fasta file. Tried \n{param_obj.Fasta} and \n{args.fasta}\n Please check.{colorama.Style.RESET_ALL}")
                raise FileNotFoundError

    # # scores
    if not args.score_peptide < 0:
        param_obj.Peptide_IDs_weight = args.score_peptide

    if not args.score_protein < 0:
        param_obj.Protein_IDs_weight = args.score_protein

    if not args.score_coverage < 0:
        param_obj.Coverage_weight = args.score_coverage

    # # proteases
    if not args.enzymes == "":
        param_obj.Proteases = parse_enzyme_list(args.enzymes)

    # # missed cleavages
    if not args.mc == "":
        param_obj.Max_MCs = parse_mc_list(args.mc)

    # # missed cleavage frequencies
    if not args.mc_freq == "":
        param_obj.Freq_MCs = parse_mc_freq_list(args.mc_freq, param_obj.Max_MCs)

    # # number or fraction of peptides to sample
    if (not args.num_peps == "") & (not args.frac_peps == ""):
        print(f"{colorama.Fore.RED}ERROR: Provided sampling based on number of peptides (--num_peps) AND fraction of generated peptides (--frac_peps). Can only use one option at a time.{colorama.Style.RESET_ALL}")
        raise RuntimeError
    elif not args.num_peps == "":
        param_obj.Peptides_Sampling_Size = parse_num_peps(args.num_peps)
        param_obj.Sampling_Size_Based_On = "number"
    elif not args.frac_peps == "":
        param_obj.Pep_Level_Proteome_Cov = parse_num_peps(args.frac_peps)
        param_obj.Sampling_Size_Based_On = "coverage"

    # # peptide mw and length filter
    if not args.min_pep_mw < 0:
        param_obj.Min_Pep_MW = args.min_pep_mw
    if not args.max_pep_mw < 0:
        param_obj.Max_Pep_MW = args.max_pep_mw
    if not args.min_pep_len < 0:
        param_obj.Min_Pep_Len = args.min_pep_len
    if not args.max_pep_len < 0:
        param_obj.Max_Pep_Len = args.max_pep_len

    # # binning
    if not args.bins == "":
        param_obj.Bins = str(args.bins)
        if (not args.export) & (not((not args.export) & (not args.digest) & (not args.analysis))):
            print(f"{colorama.Fore.CYAN}WARNING: Defining new protein bins requires re-running export. This warning can be ignored, when export with 'bins = {param_obj.Bins}' was run before.{colorama.Style.RESET_ALL}")

    # # undetectable proteome fraction
    if not args.undetectable == "":
        param_obj.Not_expressed_fraction = str(args.undetectable)
        if (not args.export) & (not ((not args.export) & (not args.digest) & (not args.analysis))):
            print(f"{colorama.Fore.CYAN}WARNING: Defining new fraction of undetectable proteome requires re-running export. This warning can be ignored, when export with 'undetectable = {param_obj.Not_expressed_fraction}' was run before.{colorama.Style.RESET_ALL}")

    # # DMSP prediction weight
    if not args.DMSP_weight < 0:
        param_obj.Weights_DeepMSPeptide_Predictions = str(args.DMSP_weight)
        param_obj.Use_DeepMSPeptide_Predictions = "True"

    # # DMSP model path
    if not args.DMSP_model == "":
        if path.isfile(path.join(args.DMSP_model)):
            param_obj.Path_DeepMSPeptide_Model = path.join(args.DMSP_model)
            param_obj.Use_DeepMSPeptide_Predictions = "True"
        else:
            print(f"{colorama.Fore.CYAN}WARNING: File with Deep-MS-Peptide prediction model ({args.DMSP_model}) does not exist. Will use model file from parameters ({param_obj.Path_DeepMSPeptide_Model}) to predict peptide detectability.{colorama.Style.RESET_ALL}")

    # # check DMSP settings
    if param_obj.Use_DeepMSPeptide_Predictions == "True":
        if not path.isfile(path.join(param_obj.Path_DeepMSPeptide_Model)):
            print(f"{colorama.Fore.CYAN}WARNING: File with Deep-MS-Peptide prediction model ({param_obj.Path_DeepMSPeptide_Model}) does not exist. Will NOT predict peptide detectability.{colorama.Style.RESET_ALL}")
            param_obj.Use_DeepMSPeptide_Predictions = "False"
        if float(param_obj.Weights_DeepMSPeptide_Predictions) < 0:
            print(
                f"{colorama.Fore.CYAN}WARNING: Set negative weights for Deep-MS-Peptide prediction ({param_obj.Weights_DeepMSPeptide_Predictions}). Will NOT predict peptide detectability.{colorama.Style.RESET_ALL}")
            param_obj.Use_DeepMSPeptide_Predictions = "False"

    # # random samplings
    if not args.samplings < 0:
        param_obj.Sampling_Number = str(args.samplings)

    if not args.use_unique_peptides_only == "":
        if args.use_unique_peptides_only == "True":
            param_obj.Use_Unique_Peptides_Only = "True"
        else:
            param_obj.Use_Unique_Peptides_Only = "False"

    # # dynamic range
    if not args.dynamic_range < 0:
        param_obj.Protein_dynamic_range = str(args.dynamic_range)
        if (not args.export) & (not((not args.export) & (not args.digest) & (not args.analysis))):
            print(f"{colorama.Fore.CYAN}WARNING: Defining new dynamic range of protein expression requires re-running export. This warning can be ignored, when export with 'dynamic_range = {param_obj.Protein_dynamic_range}' was run before.{colorama.Style.RESET_ALL}")

    # # protein weight file
    if not args.export_result == "":
        if path.isfile(path.join(args.export_result)):
            param_obj.Protein_weight_file = path.join(args.export_result)

    # # set digestion result file to params
    if not args.digestion_result == "":
        param_obj.Digestion_result_file = path.join(args.digestion_result)
        try:
            if not path.isdir(path.dirname(path.join(args.digestion_result))):
                makedirs(path.dirname(path.join(args.digestion_result)))
        except Exception:
            param_obj.Digestion_result_file = path.join(param_obj.Output_directory,
                                                        "unique_peptides_table_filtered.tsv")
    else:
        param_obj.Digestion_result_file = path.join(param_obj.Output_directory, "unique_peptides_table_filtered.tsv")

    # save updated parameter file
    file_name_time = strftime("%Y-%m-%d_%Hh%Mmin%Ssec") # will be used to prefix all output
    param_file = path.join(param_obj.Output_directory, (file_name_time + "_Autosaved_CoMPaseD_Parameters.params"))

    param_obj.save_params_to_file_from_cli(param_file)
    print(f"Saved updated parameter to {param_file}")

    # select mode and run actual tasks
    if args.export:
        run_export(param_obj, file_name_time, param_file, args)

    if args.digest:
        run_digest(param_obj, args)

    if args.analysis:
        run_analysis(param_obj, param_file, args)

    # # can only be true when all other if statements before were false
    if (not args.export) & (not args.digest) & (not args.analysis):
        run_export(param_obj, file_name_time, param_file, args)
        run_digest(param_obj, args)
        run_analysis(param_obj, param_file, args)


def parse_enzyme_list(protease_string: str) -> list:
    protease_list = protease_string.split(",")
    return protease_list


def parse_mc_list(mc_string: str) -> list:
    mc_list = mc_string.split(",")
    # convert int to str to ensure compatibility with param_obj
    mc_list = [str(mc) for mc in mc_list]
    return mc_list


def parse_mc_freq_list(mc_freq_string: str, max_mc_list: list) -> list:
    mc_freq_string_original = mc_freq_string

    max_mc_list = [int(val) for val in max_mc_list]
    mc_len = 0
    start_pos = list()
    stop_pos = list()
    for max_mc in max_mc_list:
        start_pos.append(mc_len)
        for _ in range(max_mc + 1):
            mc_len+=1
        stop_pos.append(mc_len)

    # when provided normally, mc_freq_string should not be a list, conversion just to ensure this
    if isinstance(mc_freq_string, list):
        new_mc_freq_list = list()
        for el in mc_freq_string:
            el = el.replace("[", "")
            el = el.replace("]", "")
            new_mc_freq_list.append(el)
        mc_freq_string = ",".join(new_mc_freq_list)
        mc_freq_string = new_mc_freq_list

    # replace special chars to generate a list of values convertable to float
    mc_freq_string = mc_freq_string.replace("[", "")
    mc_freq_string = mc_freq_string.replace("]", "")

    mc_freq_list = mc_freq_string.split(",")
    mc_freq_list = [float(val) for val in mc_freq_list]

    # check list length to check for conversion-errors
    if not len(mc_freq_list) == mc_len:
        print(f"{colorama.Fore.RED}ERROR: MC_freq parameter has wrong format. Converted \n'{mc_freq_string_original}' to \n'{mc_freq_list}' \nwhich should have a length of {mc_len} values but has {len(mc_freq_list)} values.{colorama.Style.RESET_ALL}")
        raise RuntimeError

    # convert back to list of strings in proper format
    tmp_mc_list = list()
    for idx, start in enumerate(start_pos):
        tmp_mc_sublist = mc_freq_list[start_pos[idx]:stop_pos[idx]]
        tmp_mc_sublist = [str(val) for val in tmp_mc_sublist]
        tmp_mc_substr = ",".join(tmp_mc_sublist)
        tmp_mc_substr = "[" + tmp_mc_substr + "]"
        tmp_mc_list.append(tmp_mc_substr)

    return tmp_mc_list


def parse_num_peps(num_peps_string: str) -> list:
    num_peps_list = num_peps_string.split(",")
    # convert int to str to ensure compatibility with param_obj
    num_peps_list = [str(peps) for peps in num_peps_list]
    return num_peps_list


def run_export(param_obj: CoMPaseD_gui_param_functions.CoMPaseD_Parameter, file_name_time: str, param_file_name: str, args):
    # get list of fasta proteins and group by length-bins
    protein_df, group_list, result_check = CoMPaseD_gui_export_functions.load_proteins_cli(param_obj)

    # assign random abundance values
    protein_df = CoMPaseD_gui_export_functions.simulate_abundance_cli(param_obj, protein_df, group_list)

    # save to file
    # if user-provided export file name, export there
    if not args.export_result == "":
        pwf_file_name = path.join(args.export_result)
    # else export to file provided in parameter file
    else:
        pwf_file_name = param_obj.Protein_weight_file

    # if path does not exist, create it and if this does not work, save a pwf in output_directory with default name
    if not path.isdir(path.dirname(pwf_file_name)):
        try:
            makedirs(path.dirname(pwf_file_name))
        except Exception:
            # fallback to default file name
            pwf_file_name = path.join(param_obj.Output_directory, (file_name_time + "_Autosaved_CoMPaseD_ProteinAbundanceExport.tsv"))
    try:
        remove(pwf_file_name)
    except OSError:
        pass
    # set pwf to this file name
    param_obj.Protein_weight_file = pwf_file_name
    param_obj.save_params_to_file_from_cli(param_file_name)
    protein_df.to_csv(pwf_file_name, index=False, sep="\t")


def run_digest(param_obj: CoMPaseD_gui_param_functions.CoMPaseD_Parameter, args):

    # try to find current python executable by sys.executable, fallback to PATH values python3 or python else
    python_exec = ""
    if executable is not None:
        python_exec = executable
    elif shutil.which("python3") is not None:
        python_exec = shutil.which("python3")
    elif shutil.which("python") is not None:
        python_exec = shutil.which("python")

    if args.use_original_proteomapper:
        file_location = path.dirname(path.realpath(__file__))
        CoMPaseD_CruxScript = path.join(file_location, 'lib','CoMPaseD_crux_script.py')



        args_list = [python_exec,
                     CoMPaseD_CruxScript,
                     "--out_folder",
                     param_obj.Output_directory,
                     "--fasta",
                     param_obj.Fasta,
                     "--crux_path",
                     param_obj.Crux_path,
                     "--clips_path",
                     param_obj.Clips_path,
                     "--promast_path",
                     param_obj.Promast_path,
                     "--enzyme_list",
                     param_obj.Proteases,
                     "--max_mc_list",
                     param_obj.Max_MCs,
                     "--min_mass",
                     str(param_obj.Min_Pep_MW),
                     "--max_mass",
                     str(param_obj.Max_Pep_MW),
                     "--min_len",
                     str(param_obj.Min_Pep_Len),
                     "--max_len",
                     str(param_obj.Max_Pep_Len),
                     "--unique_peps_file",
                     param_obj.Digestion_result_file]
        if not python_exec == "":
            completed_process = subprocess.run(args_list)
    else:
        file_location = path.dirname(path.realpath(__file__))
        CoMPaseD_PeptideMapper = path.join(file_location, 'lib','CoMPaseD_PeptideMapper.py')

        args_list= [python_exec,
                    CoMPaseD_PeptideMapper,
                    "--out_folder",
                    param_obj.Output_directory,
                    "--fasta",
                    param_obj.Fasta,
                    "--crux_path",
                    param_obj.Crux_path,
                    "--enzyme_list",
                    param_obj.Proteases,
                    "--max_mc_list",
                    param_obj.Max_MCs,
                    "--min_mass",
                    str(param_obj.Min_Pep_MW),
                    "--max_mass",
                    str(param_obj.Max_Pep_MW),
                    "--min_len",
                    str(param_obj.Min_Pep_Len),
                    "--max_len",
                    str(param_obj.Max_Pep_Len),
                    "--unique_peps_file",
                    param_obj.Digestion_result_file,
                    "--indexing_key_len",
                    str(param_obj.Indexing_key_len)
                    ]
        if not param_obj.Differentiate_I_L:
            args_list.append("--differentiate_I_L")

        if not python_exec == "":
            completed_process = subprocess.run(args_list)


def run_analysis(param_obj: CoMPaseD_gui_param_functions.CoMPaseD_Parameter, param_file_name: str, args):
    # try to find current python executable by sys.executable, fallback to PATH values python3 or python else
    python_exec = ""
    if executable is not None:
        python_exec = executable
    elif shutil.which("python3") is not None:
        python_exec = shutil.which("python3")
    elif shutil.which("python") is not None:
        python_exec = shutil.which("python")

    file_location = path.dirname(path.realpath(__file__))
    CoMPaseD_analysis_script = path.join(file_location, 'lib', 'CoMPaseD_analysis_script.py')

    args_list = [python_exec,
                 CoMPaseD_analysis_script,
                 "--param_file",
                 path.join(param_file_name),
                 "--digestion_result",
                 param_obj.Digestion_result_file]

    if not python_exec == "":
        completed_process = subprocess.run(args_list)


if __name__ == "__main__":
    main()
