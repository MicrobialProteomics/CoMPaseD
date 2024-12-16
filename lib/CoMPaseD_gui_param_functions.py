import os.path
from os import path, getcwd, makedirs, remove
from time import strftime
from fileinput import input
from sys import stdout
import PyQt6.QtWidgets as QtW

# import CoMPaseD_gui_tabs


class ValidationClass:
    """Class to count errors and hold error massages during parameter validation"""

    def __init__(self):
        self.ErrCounter = int()
        self.ErrCounter = 0
        self.Errors = list()
        self.Errors = []

    def add_error(self, message):
        self.ErrCounter += 1
        self.Errors.append(message)

    def get_result(self):
        if self.ErrCounter == 0:
            return True
        else:
            return False


class CoMPaseD_Parameter:
    """
    Class holding CoMPaseD parameter values and related functions.

    This is to be used as an interface between set values on gui and values in param files.
    """

    def __init__(self, *args, **kwargs):
        """set default parameter values"""
        # find current CoMPaseD file location to obtain relative paths
        file_location = path.dirname(path.realpath(__file__))
        clips_location = path.join(file_location, '../bin', 'perl', 'clips2.pl')
        promast_location = path.join(file_location, '../bin', 'perl', 'promast_write2.pl')
        dmsp_model_location = path.join(file_location, '../bin', 'DeepMSPep_Confetti_Model.h5')
        # variable names as required in param file

        self.Crux_path = str(getcwd()) + str("\\crux.exe")
        self.Clips_path = str(clips_location)
        self.Promast_path = str(promast_location)
        self.Indexing_key_len = "5"
        self.Differentiate_I_L = "False"
        self.Use_perl_mapping = "False"
        self.Multi_Threads = "False"
        self.Sampling_output = "True"
        self.Fasta = str(getcwd()) + str("\\placeholder.fasta")
        self.Output_directory = str(getcwd())
        self.Proteases = ["trypsin", "lysarginase", "glu-c", "chymotrypsin", "lys-c"]
        self.Max_MCs = ["2", "2", "5", "5", "2"]
        self.Freq_MCs = ["[0.7415,0.2090,0.0484]",
                         "[0.5757,0.2899,0.1336]",
                         "[0.5620,0.2753,0.1110,0.0419,0.0086,0.0012]",
                         "[0.2002,0.3369,0.2648,0.1471,0.0498,0.0012]",
                         "[0.9102,0.0836,0.0058]"]
        self.Peptides_Sampling_Size = ["10000", "10000", "10000", "10000", "10000"]
        self.Pep_Level_Proteome_Cov = ["0.033051", "0.031973356", "0.014570424", "0.009367681", "0.053870932"]
        # can be coverage, number, constant_number
        self.Sampling_Size_Based_On = "number"
        self.Bins = "0,50,100,99999"
        self.Number_of_Proteases = "5"
        self.Sampling_Number = "10"
        self.Protein_dynamic_range = "6"
        self.Not_expressed_fraction = "40,30,20"
        self.Protein_IDs_weight = "1.0"
        self.Peptide_IDs_weight = "1.0"
        self.Coverage_weight = "1.0"
        self.Use_DeepMSPeptide_Predictions = "True"
        self.Weights_DeepMSPeptide_Predictions = "1.0000"
        self.Path_DeepMSPeptide_Model = str(dmsp_model_location)
        self.Protein_weight_file = str(getcwd()) + str("\\ProteinIdentifierList.tsv")
        self.Digestion_result_file = str(getcwd()) + str("\\unique_peptides_table_filtered.tsv")

    def get_params_from_file(self, param_file_path, tab_wdg):
        """get parameter values from param file and set param object accordingly"""
        param_import = list()
        with open(param_file_path, mode="r") as param_file:
            for param_row in param_file:
                param_row = param_row.rstrip()
                if param_row != "[MPD-config]":
                    if "=" in param_row:
                        param_import.append(param_row)
        param_import_list = list()
        for param_val in param_import:
            param_val_list = param_val.split("=")
            for index, element in enumerate(param_val_list):
                element = element.strip()
                param_val_list[index] = element
            param_import_list.append(tuple(param_val_list))
        param_import_dict = dict(param_import_list)
        param_key_list = ['Crux_path', 'Use_perl_mapping', 'Indexing_key_len', 'Differentiate_I_L', 'Use_perl_mapping',
                          'Sampling_output', 'Fasta', 'Output_directory',
                          'Proteases', 'Max_MCs', 'Freq_MCs', 'Peptides_Sampling_Size', 'Pep_Level_Proteome_Cov',
                          'Sampling_Size_Based_On', 'Bins', 'Number_of_Proteases', 'Sampling_Number',
                          'Protein_dynamic_range',
                          'Not_expressed_fraction', 'Protein_IDs_weight', 'Peptide_IDs_weight', 'Coverage_weight',
                          'Use_DeepMSPeptide_Predictions', 'Protein_weight_file']
        # check for presence of these keys by list comprehension and all()
        if not all(elem in param_import_dict.keys() for elem in param_key_list):
            # whenever not all keys are present do not update values in ParamClass object
            # valid_param_file determines if update_gui is called upon loading
            valid_param_file = False
            # show warning message when attempting to load invalid file
            invalid_param_file_msg = QtW.QMessageBox(parent=tab_wdg.ExportTab)
            invalid_param_file_msg.setWindowTitle("Warning")
            invalid_param_file_msg.setText("Selected file is not a valid CoMPaseD parameter file. \n "
                                           "Cannot load settings.")
            invalid_param_file_msg.setIcon(QtW.QMessageBox.Icon.Warning)
            invalid_param_file_msg.exec()
        else:
            # set switch to update gui with valid parameter file structure,
            # this does NOT check for validity of the actual values
            valid_param_file = True
            # set ParamClass obj to new values only if param file contains corresponding keys
            self.Crux_path = param_import_dict["Crux_path"]
            if self.Use_perl_mapping == "True":
                self.Clips_path = param_import_dict["Clips_path"]
                self.Promast_path = param_import_dict["Promast_path"]
            self.Indexing_key_len = param_import_dict['Indexing_key_len']
            self.Differentiate_I_L = param_import_dict['Differentiate_I_L']
            self.Use_perl_mapping = param_import_dict['Use_perl_mapping']
            self.Multi_Threads = "False"
            self.Sampling_output = param_import_dict["Sampling_output"]
            self.Fasta = param_import_dict["Fasta"]
            self.Output_directory = param_import_dict["Output_directory"]
            self.Proteases = param_import_dict["Proteases"].split(",")
            self.Max_MCs = param_import_dict["Max_MCs"].split(",")
            self.Freq_MCs = param_import_dict["Freq_MCs"].split(",")
            self.Peptides_Sampling_Size = param_import_dict["Peptides_Sampling_Size"].split(",")
            self.Pep_Level_Proteome_Cov = param_import_dict["Pep_Level_Proteome_Cov"].split(",")
            self.Sampling_Size_Based_On = param_import_dict["Sampling_Size_Based_On"]
            self.Bins = param_import_dict["Bins"]
            self.Number_of_Proteases = param_import_dict["Number_of_Proteases"]
            self.Sampling_Number = param_import_dict["Sampling_Number"]
            self.Protein_dynamic_range = param_import_dict["Protein_dynamic_range"]
            self.Not_expressed_fraction = param_import_dict["Not_expressed_fraction"]
            self.Protein_IDs_weight = param_import_dict["Protein_IDs_weight"]
            self.Peptide_IDs_weight = param_import_dict["Peptide_IDs_weight"]
            self.Coverage_weight = param_import_dict["Coverage_weight"]
            self.Use_DeepMSPeptide_Predictions = param_import_dict["Use_DeepMSPeptide_Predictions"]
            self.Protein_weight_file = param_import_dict["Protein_weight_file"]
            # handle optional settings
            if "Path_DeepMSPeptide_Model" in param_import_dict.keys():
                self.Path_DeepMSPeptide_Model = param_import_dict["Path_DeepMSPeptide_Model"]
            if "Weights_DeepMSPeptide_Predictions" in param_import_dict.keys():
                self.Weights_DeepMSPeptide_Predictions = param_import_dict["Weights_DeepMSPeptide_Predictions"]
            if "Digestion_result_file" in  param_import_dict.keys():
                self.Digestion_result_file = param_import_dict["Digestion_result_file"]
        # validate param values in ParamClass obj and correct typical formatting errors
        if valid_param_file:
            parameter_validation, parameter_error_count, parameter_errors = self.validate_params()
            # set tab_wdg param file path accordingly
            tab_wdg.param_file_name = param_file_path
        else:
            # don't do so when no valid param file was found
            # instead set param validation reporters to default values
            parameter_validation = False
            parameter_error_count = -1
            parameter_errors = ""

        return valid_param_file, parameter_validation, parameter_error_count, parameter_errors

    def set_params_to_gui(self, tab_wdg):
        """set current values from param obj to gui"""
        # self.get_params_from_file(self, param_file)
        tab_wdg.CruxPathField.setText(self.Crux_path)
        if self.Use_perl_mapping == "True":
            tab_wdg.ClipsPathField.setText(self.Clips_path)
            tab_wdg.PromastPathField.setText(self.Promast_path)

        if self.Sampling_output == "True":
            tab_wdg.SamplingOutputCheckbox.setChecked(True)
        if self.Sampling_output == "False":
            tab_wdg.SamplingOutputCheckbox.setChecked(False)

        tab_wdg.FastaPathField.setText(self.Fasta)
        tab_wdg.OutputPathField.setText(self.Output_directory)

        if self.Sampling_Size_Based_On == "number":
            tab_wdg.SamplingSizePeptides.setChecked(True)
        if self.Sampling_Size_Based_On == "coverage":
            tab_wdg.SamplingSizeFraction.setChecked(True)

        tab_wdg.ProteinBinsField.setText(self.Bins)
        tab_wdg.ProteinNotExprFracField.setText(self.Not_expressed_fraction)

        tab_wdg.MaxProteasesSpinbox.setValue(int(self.Number_of_Proteases))
        tab_wdg.RandomSamplingsSpinbox.setValue(int(self.Sampling_Number))
        tab_wdg.DynamicRangeSpinbox.setValue(float(self.Protein_dynamic_range))

        tab_wdg.ProtIDWeightField.setText(self.Protein_IDs_weight)
        tab_wdg.PepIDWeightField.setText(self.Peptide_IDs_weight)
        tab_wdg.CoverageWeightField.setText(self.Coverage_weight)

        if self.Protein_weight_file != "":
            self.export_run_counter = 1
            tab_wdg.Protein_weight_file = self.Protein_weight_file

        if self.Use_DeepMSPeptide_Predictions == "True":
            tab_wdg.DMSPBox.setChecked(True)
            tab_wdg.DMSPModelPathField.setText(self.Path_DeepMSPeptide_Model)
            tab_wdg.DMSPWeightSpinbox.setValue(float(self.Weights_DeepMSPeptide_Predictions))
        if self.Use_DeepMSPeptide_Predictions == "False":
            tab_wdg.DMSPBox.setChecked(False)

        # load protease table data
        # set row count by number of proteases
        tab_wdg.ProteaseTable.setRowCount(len(self.Proteases))
        # adjust column names
        if self.Sampling_Size_Based_On == "number":
            tab_wdg.ProteaseTable.setHorizontalHeaderLabels(
                ["Protease", "Max MCs", "MC Frequency", "Number of Peptides"])
        if self.Sampling_Size_Based_On == "coverage":
            tab_wdg.ProteaseTable.setHorizontalHeaderLabels(
                ["Protease", "Max MCs", "MC Frequency", "Fraction of Peptides"])

            # nested list of MC frequencies needs fixing as it was changed to a list of strings
            # e.g: input = ['[0.7415', '0.2090', '0.0484]', '[0.5757', '0.2899', '0.0058]']
            #     output = ['[0.7415,0.2090,0.0484]','[0.5757,0.2899,0.0058]']

        def mc_freq_string_to_list(mc_freq_string):
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

            for sublist in mc_freq_list:
                SublistStr = str(sublist).replace(" ", "")
                final_mc_freq_list.append(SublistStr)
            return mc_freq_list

        self.Freq_MCs = mc_freq_string_to_list(self.Freq_MCs)

        for table_row in range(len(self.Proteases)):
            tab_wdg.ProteaseTable.setItem(table_row, 0, QtW.QTableWidgetItem(self.Proteases[table_row]))
            tab_wdg.ProteaseTable.setItem(table_row, 1, QtW.QTableWidgetItem(self.Max_MCs[table_row]))
            tab_wdg.ProteaseTable.setItem(table_row, 2, QtW.QTableWidgetItem(str(self.Freq_MCs[table_row])))
            if self.Sampling_Size_Based_On == "number":
                tab_wdg.ProteaseTable.setItem(
                    table_row, 3, QtW.QTableWidgetItem(self.Peptides_Sampling_Size[table_row]))
            if self.Sampling_Size_Based_On == "coverage":
                tab_wdg.ProteaseTable.setItem(
                    table_row, 3, QtW.QTableWidgetItem(self.Pep_Level_Proteome_Cov[table_row]))
        return self

    def get_params_from_gui(self, tab_wdg):
        """set current values from gui to param obj"""
        self.Crux_path = tab_wdg.CruxPathField.text()
        if self.Use_perl_mapping == "True":
            self.Clips_path = tab_wdg.ClipsPathField.text()
            self.Promast_path = tab_wdg.PromastPathField.text()
        # set values also, if these are correct files but not in use
        elif self.Use_perl_mapping == "False":
            if path.isfile(tab_wdg.ClipsPathField.text()):
                self.Clips_path = tab_wdg.ClipsPathField.text()
            if path.isfile(tab_wdg.PromastPathField.text()):
                self.Promast_path = tab_wdg.PromastPathField.text()

        # keep Multi_Threads of (deprecated after replacement of perl scripts)
        self.Multi_Threads = "False"

        if tab_wdg.SamplingOutputCheckbox.isChecked():
            self.Sampling_output = "True"
        else:
            self.Sampling_output = "False"

        self.Fasta = tab_wdg.FastaPathField.text()
        self.Output_directory = tab_wdg.OutputPathField.text()

        # get table values into lists
        Protease_list = list()
        MC_list = list()
        Freq_list = list()
        Size_list = list()
        # loop through table and obtain cell values
        for row in range(tab_wdg.ProteaseTable.rowCount()):
            # use if condition to avoid adding empty row upon table clearing
            if tab_wdg.ProteaseTable.item(row, 0) is not None:
                Protease_list.append(tab_wdg.ProteaseTable.item(row, 0).text())
            if tab_wdg.ProteaseTable.item(row, 1) is not None:
                MC_list.append(tab_wdg.ProteaseTable.item(row, 1).text())
            if tab_wdg.ProteaseTable.item(row, 2) is not None:
                Freq_list.append(tab_wdg.ProteaseTable.item(row, 2).text())
            if tab_wdg.ProteaseTable.item(row, 3) is not None:
                Size_list.append(tab_wdg.ProteaseTable.item(row, 3).text())
        # set parameters to new values
        self.Proteases = Protease_list
        self.Max_MCs = MC_list
        self.Freq_MCs = Freq_list
        # for sampling size update correct parameter by selection
        if tab_wdg.SamplingSizePeptides.isChecked():
            self.Peptides_Sampling_Size = Size_list
            self.Sampling_Size_Based_On = "number"
        if tab_wdg.SamplingSizeFraction.isChecked():
            self.Sampling_Size_Based_On = "coverage"
            self.Pep_Level_Proteome_Cov = Size_list

        self.Bins = tab_wdg.ProteinBinsField.text()
        self.Not_expressed_fraction = tab_wdg.ProteinNotExprFracField.text()
        self.Number_of_Proteases = str(tab_wdg.MaxProteasesSpinbox.value())
        self.Sampling_Number = str(tab_wdg.RandomSamplingsSpinbox.value())
        self.Protein_dynamic_range = str(tab_wdg.DynamicRangeSpinbox.value())
        self.Protein_IDs_weight = str(tab_wdg.ProtIDWeightField.text())
        self.Peptide_IDs_weight = str(tab_wdg.PepIDWeightField.text())
        self.Coverage_weight = str(tab_wdg.CoverageWeightField.text())
        if tab_wdg.DMSPBox.isChecked():
            self.Use_DeepMSPeptide_Predictions = "True"
            self.Path_DeepMSPeptide_Model = tab_wdg.DMSPModelPathField.text()
            self.Weights_DeepMSPeptide_Predictions = str(tab_wdg.DMSPWeightSpinbox.value())
        if not tab_wdg.DMSPBox.isChecked():
            self.Use_DeepMSPeptide_Predictions = "False"

        self.Protein_weight_file = path.join(tab_wdg.Protein_weight_file)
        return self

    def save_params_to_file(self, param_file_path, tab_wdg):
        """save values from param obj to file"""
        # update param obj once
        self.get_params_from_gui(tab_wdg=tab_wdg)
        parameter_validation, parameter_error_count, parameter_errors = self.validate_params()
        # fix formatting issues
        self.fix_proteases()
        self.fix_max_missed_cleavages()
        self.fix_missed_cleavage_freq()
        self.fix_pep_level_cov()
        self.fix_pep_sampling_size()
        # write to param file in SeqIO.parse config format
        with open(param_file_path, mode="w") as param_file:
            param_file.write("[MPD-config]")
            param_file.write("\n")
            param_file.writelines(self.get_param_list())
        # update tab_wdg with new param_file_path
        tab_wdg.param_file_name = param_file_path

        return parameter_validation, parameter_error_count, parameter_errors


    def save_params_to_file_from_cli(self, param_file_path):
        """save values from param obj to file without gui"""
        parameter_validation, parameter_error_count, parameter_errors = self.validate_params()
        # fix formatting issues
        self.fix_proteases()
        self.fix_max_missed_cleavages()
        self.fix_missed_cleavage_freq()
        self.fix_pep_level_cov()
        self.fix_pep_sampling_size()
        # write to param file in SeqIO.parse config format
        with open(param_file_path, mode="w") as param_file:
            param_file.write("[MPD-config]")
            param_file.write("\n")
            param_file.writelines(self.get_param_list())

        return parameter_validation, parameter_error_count, parameter_errors

    # ####################
    # Helper functions:

    def get_param_list(self):
        """"generate parameter file-like string from CoMPaseD_Parameter obj"""
        parameter_file = list()
        for param_property in vars(self).keys():
            param_val = vars(self)[param_property]
            param_str = str(param_property) + " = " + str(param_val) + "\n"
            parameter_file.append(param_str)
        return parameter_file

    def get_pwf(self):
        '''return line with Protein_weight_file only'''
        param_str = str()
        for param_property in vars(self).keys():
            if param_property == "Protein_weight_file":
                param_val = vars(self)[param_property]
                param_str = str(param_property) + " = " + str(param_val) + "\n"
            return param_str

    # functions to convert sub-lists from get_param_list() to strings
    def fix_proteases(self):
        if isinstance(self.Proteases, list):
            tmp = ','.join(self.Proteases)
            self.Proteases = tmp

    def fix_max_missed_cleavages(self):
        if isinstance(self.Max_MCs, list):
            tmp = ','.join(self.Max_MCs)
            self.Max_MCs = tmp

    def fix_missed_cleavage_freq(self):
        if isinstance(self.Freq_MCs, list):
            tmp = ','.join(self.Freq_MCs)
            self.Freq_MCs = tmp

    def fix_pep_level_cov(self):
        if isinstance(self.Pep_Level_Proteome_Cov, list):
            tmp = ','.join(self.Pep_Level_Proteome_Cov)
            self.Pep_Level_Proteome_Cov = tmp

    def fix_pep_sampling_size(self):
        if isinstance(self.Peptides_Sampling_Size, list):
            tmp = ','.join(self.Peptides_Sampling_Size)
            self.Peptides_Sampling_Size = tmp

    # function to correct typical errors in parameter obj and check validity of the vals
    def validate_params(self, mode=1):
        # error counter
        Validation = ValidationClass()

        # check file path for existence
        if mode==1:
            if not path.isfile(self.Crux_path):
                Validation.add_error(message="Crux path does not exist.")
            # do NOT check clips and promast paths unless using perl scripts is explicitly on
            if self.Use_perl_mapping:
                if not path.isfile(self.Clips_path):
                    Validation.add_error(message="Clips path does not exist.")
                if not path.isfile(self.Promast_path):
                    Validation.add_error(message="Promast path does not exist.")
        if not path.isfile(self.Fasta):
            Validation.add_error(message="Fasta file does not exist.")
        # check output folder for existence and try to create if necessary
        if not path.isdir(self.Output_directory):
            try:
                makedirs(self.Output_directory)
            except Exception:
                Validation.add_error(message="Output folder does not exist and could not be created.")

        # function to remove invalid strings from parameter lists
        def clean_parameter_values(parameter_string, parameter_name, allowed_values='0123456789,.[]'):
            clean_parameter_string = ""
            for Chr in parameter_string:
                if Chr in allowed_values:
                    clean_parameter_string += Chr
                else:
                    if mode==1:
                        Validation.add_error(message=f"'{Chr}' was removed from {parameter_string} in {parameter_name}.")
                # convert underscore to hyphen
                if Chr == "_":
                    clean_parameter_string += "-"
                    if mode==1:
                        Validation.add_error(
                            message=f"'_' was replaced by '-' from {parameter_string} in {parameter_name}.")
                # convert semicolon to comma
                # should not effect input from protease table as it is handled as nested list
                if Chr == ";":
                    clean_parameter_string += ","
                    if mode==1:
                        Validation.add_error(
                            message=f"';' was replaced by ',' from {parameter_string} in {parameter_name}.")
            return clean_parameter_string

        # function to test if a string is a valid numeric value, accepts math expressions also
        def test_float(a):
            try:
                float(eval(a))
                return True
            except Exception:
                return False

        # correct protease table lists
        # limit correction by isinstance to cases when values have not been converted to one string
        if isinstance(self.Proteases, list):
            tmp_Protease = list()
            for Protease in self.Proteases:
                # convert protease automatically to lower as currently upper-case proteases are not defined in crux
                # ToDo: check at new crux releases
                Protease = Protease.lower()
                tmp_Protease.append(clean_parameter_values(parameter_string=Protease,
                                                           parameter_name="Proteases",
                                                           allowed_values='abcdefghijklmnopqrstuvwxyz-,/'))
            self.Proteases = tmp_Protease

        if isinstance(self.Max_MCs, list):
            tmp_MaxMCs = list()
            for MC in self.Max_MCs:
                tmp_MaxMCs.append(clean_parameter_values(parameter_string=MC,
                                                         parameter_name="Maximal missed cleavage sites",
                                                         allowed_values='0123456789'))
            self.Max_MCs = tmp_MaxMCs

        if isinstance(self.Freq_MCs, list):
            tmp_FreqMCs = list()
            for MC_Freq in self.Freq_MCs:
                tmp_FreqMCs.append(clean_parameter_values(parameter_string=MC_Freq,
                                                          parameter_name="Missed cleavage site frequencies"))
            self.Max_MCs = tmp_MaxMCs

        if isinstance(self.Peptides_Sampling_Size, list):
            tmp_SamplingSize = list()
            for SamplingSize in self.Peptides_Sampling_Size:
                tmp_SamplingSize.append(clean_parameter_values(parameter_string=SamplingSize,
                                                               parameter_name="Sampled peptides",
                                                               allowed_values="0123456789"))
            self.Peptides_Sampling_Size = tmp_SamplingSize

        if isinstance(self.Pep_Level_Proteome_Cov, list):
            tmp_CoverageSamplingSize = list()
            for SamplingSize in self.Pep_Level_Proteome_Cov:
                tmp_CoverageSamplingSize.append(clean_parameter_values(parameter_string=SamplingSize,
                                                                       parameter_name="Sampled peptide fraction",
                                                                       allowed_values="0123456789."))
            self.Pep_Level_Proteome_Cov = tmp_CoverageSamplingSize

        self.Bins = clean_parameter_values(parameter_string=self.Bins,
                                           parameter_name="Bin definition",
                                           allowed_values="0123456789,")
        for Bin in self.Bins.split(","):
            if not test_float(Bin):
                Validation.add_error(message="Bin definition contains non-numeric values.")

        self.Not_expressed_fraction = clean_parameter_values(parameter_string=self.Not_expressed_fraction,
                                                             parameter_name="Not expresseed fraction",
                                                             allowed_values="0123456789,.")

        if not test_float(self.Protein_IDs_weight):
            Validation.add_error(message="Protein ID weighting factor is not numeric.")

        if not test_float(self.Peptide_IDs_weight):
            Validation.add_error(message="Peptide ID weighting factor is not numeric.")

        if not test_float(self.Coverage_weight):
            Validation.add_error(message="Protein coverage weighting factor is not numeric.")

        # check DMSP parameters only if enabled
        if self.Use_DeepMSPeptide_Predictions == "True":
            if not test_float(self.Weights_DeepMSPeptide_Predictions):
                Validation.add_error(message="Deep MS-Peptide weighting factor is not numeric.")

            if not path.isfile(self.Path_DeepMSPeptide_Model):
                Validation.add_error(message="Deep MS-Peptide model path does not exist.")

        if Validation.ErrCounter > 0:
            pass

        if not Validation.get_result():
            return False, Validation.ErrCounter, Validation.Errors
        else:
            return True, "0", ""

    def param_load_file(self, tab_wdg):
        '''open parameter file'''

        # clear any existing protein_weight_file in tab_wdg
        tab_wdg.Protein_weight_file = ""

        # upon first param file save this will be true
        # but will be False again, when changing fasta file via the browse button
        if path.isdir(str(path.dirname(tab_wdg.param_file_name))):
            file_name_tmp, _ = QtW.QFileDialog.getOpenFileName(tab_wdg, "Load CoMPaseD Parameter File",
                                                           directory=(path.dirname(tab_wdg.param_file_name)),
                                                           filter="Parameter File (*.params)")
            if file_name_tmp != '':
                self.get_params_from_file(file_name_tmp, tab_wdg) # update param obj
                self.set_params_to_gui(tab_wdg) # set param obj vals to gui
                tab_wdg.param_file_name = file_name_tmp
                # set run_counter to 1 when parameter file does not contain a valid file name
                # and prevent new automatic abundance simulation
                if path.isfile(tab_wdg.Protein_weight_file):
                    tab_wdg.export_run_counter = 1
                else:
                    tab_wdg.export_run_counter = 0

        # handle other cases
        else:
            file_name_tmp, _ = QtW.QFileDialog.getOpenFileName(tab_wdg, "Load CoMPaseD Parameter File",
                                                               directory=getcwd(),
                                                               filter="Parameter File (*.params)")
            if file_name_tmp != '':
                self.get_params_from_file(file_name_tmp, tab_wdg)  # update param obj
                self.set_params_to_gui(tab_wdg)  # set param obj vals to gui
                tab_wdg.param_file_name = file_name_tmp
                tab_wdg.export_run_counter = 1
                # set run_counter to 1 when parameter file does not contain a valid file name
                # and prevent new automatic abundance simulation
                if path.isfile(tab_wdg.Protein_weight_file):
                    tab_wdg.export_run_counter = 1
                else:
                    tab_wdg.export_run_counter = 0

        return file_name_tmp


    def param_save_file(self, tab_wdg):
        '''save parameter file, ask user where to save'''
        self.get_params_from_gui(tab_wdg) # update param obj

        # open containing folder if param file was saved before
        if path.isdir(str(path.dirname(tab_wdg.param_file_name))):
            file_name_tmp, _ = QtW.QFileDialog.getSaveFileName(tab_wdg, "Save CoMPaseD Parameter File",
                                                               directory=path.dirname(tab_wdg.param_file_name),
                                                               filter="Parameter File (*.params)")
            # if user did not cancel, save
            if file_name_tmp != '':
                self.save_params_to_file(param_file_path=file_name_tmp, tab_wdg=tab_wdg)
                tab_wdg.param_file_name = file_name_tmp

        # open cwd for cases where no param file was saved before or this path is not valid
        else:
            file_name_tmp, _ = QtW.QFileDialog.getSaveFileName(tab_wdg, "Save CoMPaseD Parameter File",
                                                               directory=getcwd(),
                                                               filter="Parameter File (*.params)")
            if file_name_tmp != '':
                self.save_params_to_file(param_file_path=file_name_tmp, tab_wdg=tab_wdg)
                tab_wdg.param_file_name = file_name_tmp

        return file_name_tmp


    def param_silent_save_file(self, tab_wdg):
        '''save parameter file, without asking if possible'''

        # for existing param file, just prefix with Autosaved_
        if path.isfile(tab_wdg.param_file_name):
            out_path = path.dirname(path.join(tab_wdg.OutputPathField.text(), ""))
            param_file = path.basename(tab_wdg.param_file_name)
            new_param_file = "Autosaved_" + param_file
            file_name_tmp = path.join(out_path, new_param_file)
            # if "Autosaved_..." was not already created, rename accordingly
            if not path.isfile(file_name_tmp):
                self.save_params_to_file(param_file_path=file_name_tmp, tab_wdg=tab_wdg)
                # set file name to tab_wdg obj
                tab_wdg.param_file_name = file_name_tmp
                return file_name_tmp
            # if "Autosaved_..." was already created, try to replace
            else:
                try:
                    remove(file_name_tmp)
                    self.save_params_to_file(param_file_path=file_name_tmp, tab_wdg=tab_wdg)
                    # set file name to tab_wdg obj
                    tab_wdg.param_file_name = file_name_tmp
                    return file_name_tmp
                except PermissionError:
                    # show warning message when file can not be deleted
                    del_file_error_msg = QtW.QMessageBox(parent=tab_wdg.ExportTab)
                    del_file_error_msg.setWindowTitle("Error")
                    del_file_error_msg.setText("Could not auto-save parameter file due to missing "
                                               "permission to delete existing file: \n %s" % (str(file_name_tmp)))
                    del_file_error_msg.setIcon(QtW.QMessageBox.Icon.Critical)
                    del_file_error_msg.exec()
                except Exception:
                    # show warning message when file can not be deleted
                    del_file_error_msg = QtW.QMessageBox(parent=tab_wdg.ExportTab)
                    del_file_error_msg.setWindowTitle("Error")
                    del_file_error_msg.setText("Could not auto-save parameter file due to unknown reason.")
                    del_file_error_msg.setIcon(QtW.QMessageBox.Icon.Critical)
                    del_file_error_msg.exec()

        # for tab_wdg.param_file_name not pointing to  an existing param file, name it with timestamp_Autosaved_params
        else:
            param_path = path.dirname(path.join(tab_wdg.OutputPathField.text(), ""))
            param_file = strftime("%Y-%m-%d_%Hh%Mmin%Ssec_") + "Autosaved_params.params"
            file_name_tmp = path.join(param_path, param_file)

            # this should almost always be true
            if not path.isfile(file_name_tmp):
                self.save_params_to_file(param_file_path=file_name_tmp, tab_wdg=tab_wdg)
                tab_wdg.param_file_name = file_name_tmp
                return file_name_tmp

            # for rare cases of identical names, try to replace file
            else:
                try:
                    remove(file_name_tmp)
                    self.save_params_to_file(param_file_path=file_name_tmp, tab_wdg=tab_wdg)
                    tab_wdg.param_file_name = file_name_tmp
                    return file_name_tmp
                except PermissionError:
                    # show warning message when file can not be deleted
                    del_file_error_msg = QtW.QMessageBox(parent=tab_wdg.ExportTab)
                    del_file_error_msg.setWindowTitle("Error")
                    del_file_error_msg.setText("Could not auto-save parameter file due to missing "
                                               "permission to delete existing file: \n %s" % (str(file_name_tmp)))
                    del_file_error_msg.setIcon(QtW.QMessageBox.Icon.Critical)
                    del_file_error_msg.exec()
                except Exception:
                    # show warning message when file can not be deleted
                    del_file_error_msg = QtW.QMessageBox(parent=tab_wdg.ExportTab)
                    del_file_error_msg.setWindowTitle("Error")
                    del_file_error_msg.setText("Could not auto-save parameter file due to unknown reason.")
                    del_file_error_msg.setIcon(QtW.QMessageBox.Icon.Critical)
                    del_file_error_msg.exec()

    '''
    def update_pwf_param(self, tab_wdg):
        if path.isfile(tab_wdg.param_file_name):
            param_file_path = path.join(tab_wdg.param_file_name)
            try:
                for line in input(param_file_path, inplace=True):
                    if line.strip().startswith('Protein_weight_file = '):
                        line = self.get_pwf()
                    stdout.write(line)
            except PermissionError:
                perm_error_msg = QtW.QMessageBox(parent=tab_wdg.ExportTab)
                perm_error_msg.setWindowTitle("Error")
                perm_error_msg.setText("Could not save protein-weights file due to missing "
                                       "permission: \n %s" % (str(param_file_path)))
                perm_error_msg.setIcon(QtW.QMessageBox.Icon.Critical)
                perm_error_msg.exec()
        return None
    '''

    def load_params(self, param_file_path):
        """get parameter values from param file and set param object accordingly"""
        param_import = list()
        with open(param_file_path, mode="r") as param_file:
            for param_row in param_file:
                param_row = param_row.rstrip()
                if param_row != "[MPD-config]":
                    if "=" in param_row:
                        param_import.append(param_row)
        param_import_list = list()
        for param_val in param_import:
            param_val_list = param_val.split("=")
            for index, element in enumerate(param_val_list):
                element = element.strip()
                param_val_list[index] = element
            param_import_list.append(tuple(param_val_list))
        param_import_dict = dict(param_import_list)
        param_key_list = ['Crux_path', 'Use_perl_mapping', 'Indexing_key_len', 'Differentiate_I_L', 'Use_perl_mapping',
                          'Sampling_output', 'Fasta', 'Output_directory',
                          'Proteases', 'Max_MCs', 'Freq_MCs', 'Peptides_Sampling_Size', 'Pep_Level_Proteome_Cov',
                          'Sampling_Size_Based_On', 'Bins', 'Number_of_Proteases', 'Sampling_Number',
                          'Protein_dynamic_range',
                          'Not_expressed_fraction', 'Protein_IDs_weight', 'Peptide_IDs_weight', 'Coverage_weight',
                          'Use_DeepMSPeptide_Predictions', 'Protein_weight_file']
        # check for presence of these keys by list comprehension and all()
        if not all(elem in param_import_dict.keys() for elem in param_key_list):
            # whenever not all keys are present do not update values in ParamClass object but throw error
            for key in param_import_dict.keys():
                param_key_list.remove(key)
            raise ValueError(f"Parameter file does not contain all required information. Missing information: {param_key_list}")
        else:
            # set switch to update gui with valid parameter file structure,
            # this does NOT check for validity of the actual values
            valid_param_file = True
            # set ParamClass obj to new values only if param file contains corresponding keys
            self.Crux_path = param_import_dict["Crux_path"]
            if self.Use_perl_mapping == "True":
                self.Clips_path = param_import_dict["Clips_path"]
                self.Promast_path = param_import_dict["Promast_path"]
            self.Multi_Threads = "False"
            self.Sampling_output = param_import_dict["Sampling_output"]
            self.Fasta = param_import_dict["Fasta"]
            self.Output_directory = param_import_dict["Output_directory"]
            self.Proteases = param_import_dict["Proteases"].split(",")
            self.Max_MCs = param_import_dict["Max_MCs"].split(",")
            self.Freq_MCs = param_import_dict["Freq_MCs"].split(",")
            self.Peptides_Sampling_Size = param_import_dict["Peptides_Sampling_Size"].split(",")
            self.Pep_Level_Proteome_Cov = param_import_dict["Pep_Level_Proteome_Cov"].split(",")
            self.Sampling_Size_Based_On = param_import_dict["Sampling_Size_Based_On"]
            self.Bins = param_import_dict["Bins"]
            self.Number_of_Proteases = param_import_dict["Number_of_Proteases"]
            self.Sampling_Number = param_import_dict["Sampling_Number"]
            self.Protein_dynamic_range = param_import_dict["Protein_dynamic_range"]
            self.Not_expressed_fraction = param_import_dict["Not_expressed_fraction"]
            self.Protein_IDs_weight = param_import_dict["Protein_IDs_weight"]
            self.Peptide_IDs_weight = param_import_dict["Peptide_IDs_weight"]
            self.Coverage_weight = param_import_dict["Coverage_weight"]
            self.Use_DeepMSPeptide_Predictions = param_import_dict["Use_DeepMSPeptide_Predictions"]
            self.Protein_weight_file = param_import_dict["Protein_weight_file"]
            # handle optional settings
            if "Path_DeepMSPeptide_Model" in param_import_dict.keys():
                self.Path_DeepMSPeptide_Model = param_import_dict["Path_DeepMSPeptide_Model"]
            if "Weights_DeepMSPeptide_Predictions" in param_import_dict.keys():
                self.Weights_DeepMSPeptide_Predictions = param_import_dict["Weights_DeepMSPeptide_Predictions"]
            # validate param values in ParamClass obj and correct typical formatting errors

            parameter_validation, parameter_error_count, parameter_errors = self.validate_params()

        return valid_param_file, parameter_validation, parameter_error_count, parameter_errors


# functions, not directly related to parameter class but to parameter tab
def param_browse_fasta(tab_wdg):
    """browse for fasta to analyse"""
    if path.isdir(str(path.dirname(tab_wdg.FastaPathField.text()))):
        file_name, _ = QtW.QFileDialog.getOpenFileName(tab_wdg, "FASTA protein database",
                                                       directory=path.dirname(
                                                           tab_wdg.FastaPathField.text()),
                                                       filter="FASTA File (*.fasta *.fa *.fas *.faa);;All "
                                                              "files (*)")
    else:
        file_name, _ = QtW.QFileDialog.getOpenFileName(tab_wdg, "FASTA protein database",
                                                       directory=getcwd(),
                                                       filter="FASTA File (*.fasta *.fa *.fas *.faa);;All "
                                                              "files (*)")
    if file_name != "":
        tab_wdg.FastaPathField.setText(file_name)
        if not tab_wdg.param_file_name == "":
            tab_wdg.param_file_name = ""  # reset param file name when changing fasta
            tab_wdg.export_run_counter = 0
            tab_wdg.params_obj.Protein_weight_file = ""
            tab_wdg.Protein_weight_file = ""


def param_browse_out_folder(tab_wdg):
    """browse for out folder"""
    if path.isdir(str(path.dirname(tab_wdg.OutputPathField.text()))):
        folder_name = QtW.QFileDialog.getExistingDirectory(tab_wdg, "Output folder",
                                                           directory=path.dirname(
                                                               tab_wdg.OutputPathField.text()))
    else:
        folder_name = QtW.QFileDialog.getExistingDirectory(tab_wdg, "Output folder",
                                                           directory=getcwd())
    if folder_name != "":
        tab_wdg.OutputPathField.setText(folder_name)


def param_browse_dmsp_model(tab_wdg):
    """browse for deep-MS-peptide model file"""
    if path.isdir(str(path.dirname(tab_wdg.DMSPModelPathField.text()))):
        file_name, _ = QtW.QFileDialog.getOpenFileName(tab_wdg, "Deep MS-Peptide model",
                                                       directory=path.dirname(
                                                           tab_wdg.DMSPModelPathField.text()),
                                                       filter=".h5 Model (*.h5);;All files (*)")
    else:
        file_location = path.dirname(path.realpath(__file__))
        dmsp_model_loc = path.join(file_location, '../bin', 'DeepMSPep_Confetti_Model.h5')
        if path.isfile(dmsp_model_loc):
            file_name, _ = QtW.QFileDialog.getOpenFileName(tab_wdg, "Deep MS-Peptide model",
                                                           directory=dmsp_model_loc,
                                                           filter=".h5 Model (*.h5);;All files (*)")
        else:
            file_name, _ = QtW.QFileDialog.getOpenFileName(tab_wdg, "Deep MS-Peptide model",
                                                           directory=getcwd(),
                                                           filter=".h5 Model (*.h5);;All files (*)")
        if file_name != "":
            tab_wdg.DMSPModelPathField.setText(file_name)