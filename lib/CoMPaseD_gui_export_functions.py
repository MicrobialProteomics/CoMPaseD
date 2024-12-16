import colorama
import numpy as np
import PyQt6.QtWidgets as QtW
from os import path, getcwd, makedirs
from time import sleep
from pandas import DataFrame, read_csv, merge, to_numeric
from random import shuffle
from Bio import SeqIO
from lib import CoMPaseD_protein_class
from lib.CoMPaseD_tools import *
from lib.CoMPaseD_gui_param_functions import CoMPaseD_Parameter


def make_folders(tab_wdg, param_obj: CoMPaseD_Parameter) -> bool:
    '''check existence of out_folder and generate it if required
    does not check whether out directory is empty or if there are files to overwrite
    '''
    out_file = param_obj.Output_directory
    err_return = True
    if not path.dirname(out_file):
        try:
            makedirs(path.dirname(out_file))
        except ValueError:
            err_return = False
            # show error message when folders can not be created
            create_out_error_msg = QtW.QMessageBox(parent=tab_wdg.ExportTab)
            create_out_error_msg.setWindowTitle("Error")
            create_out_error_msg.setText("Could not create out folder. Please check permission.")
            create_out_error_msg.setIcon(QtW.QMessageBox.Icon.Critical)
            create_out_error_msg.exec()
    return err_return


def load_proteins(tab_wdg, param_obj: CoMPaseD_Parameter):
    err_return = True
    try:
        fasta_file = SeqIO.parse(param_obj.Fasta, "fasta")
        protein_list = CoMPaseD_protein_class.makeProteinList(fasta_file)
        # loop through fasta and store protein ID and len
        protein_names_list = list()
        protein_length_list = list()
        protein_group_list = list()
        grouping_column_list = list()

        bin_list = config_to_numeric_list(param_obj.Bins)
        # assume bin_list with len 3 or 4 indicates small, (medium), and large proteins; use 3 and 4 instead of 2 and 3 as there has to be one value more than the number of bins
        if len(bin_list) == 3:
            protein_group_list = ["small_proteins", "large_proteins"]
        elif len(bin_list) == 4:
            protein_group_list = ["small_proteins", "medium_proteins", "large_proteins"]
        else:
            for bin_pos in range(0, (len(bin_list)-1)):
                bin_list[bin_pos] = int(bin_list[bin_pos])
                protein_group_list.append("group_" + (str(bin_pos+1)))

        bin_list_array = np.array(bin_list)

        for record in protein_list:
            tmp_id = record.id
            tmp_length = record.length
            # test in which size group the protein belongs - the bins are ordered and calculated so that each return must be of length 1
            for lower_border in bin_list:
                lower_border_index = int(np.where(bin_list_array == lower_border)[0])
                # test only, if the protein is longer than the lower bound of the bin
                if tmp_length > lower_border:
                    # catch cases where the protein is longer than the provided max
                    if tmp_length > max(bin_list):
                        tmp_group = "unknown"
                        # correct protein_group_list in case, it was not corrected before, since it was was created in an earlier loop, the order should not be disturbed and 'unknown' will be the last element
                        if protein_group_list.count("unknown") == 0:
                            protein_group_list.append("unknown")
                    else:
                        # ensure, that the protein is shorter or as long as the upper bound;
                        if bin_list[lower_border_index + 1] >= tmp_length:
                            tmp_group = (protein_group_list[lower_border_index])

            protein_names_list.append(tmp_id)
            protein_length_list.append(tmp_length)
            grouping_column_list.append(tmp_group)

        # merge lists to pandas df
        out_df = DataFrame()
        out_df['Identifier'] = protein_names_list
        out_df['Sequence_Len[aa]'] = protein_length_list
        out_df['Group'] = grouping_column_list


        return out_df, protein_group_list, err_return
    except Exception:
        err_return = False
        return None, None, err_return

def load_proteins_cli(param_obj: CoMPaseD_Parameter):
    err_return = True
    try:
        fasta_file = SeqIO.parse(param_obj.Fasta, "fasta")
        protein_list = CoMPaseD_protein_class.makeProteinList(fasta_file)
        # loop through fasta and store protein ID and len
        protein_names_list = list()
        protein_length_list = list()
        protein_group_list = list()
        grouping_column_list = list()

        bin_list = config_to_numeric_list(param_obj.Bins)
        if len(bin_list) == 3:
            protein_group_list = ["small_proteins", "large_proteins"]
        elif len(bin_list) == 4:
            protein_group_list = ["small_proteins", "medium_proteins", "large_proteins"]
        for bin_pos in range(0, (len(bin_list)-1)):
            bin_list[bin_pos] = int(bin_list[bin_pos])
            protein_group_list.append("group_" + (str(bin_pos+1)))

        bin_list_array = np.array(bin_list)

        for record in protein_list:
            tmp_id = record.id
            tmp_length = record.length
            # test in which size group the protein belongs - the bins are ordered and calculated so that each return must be of length 1
            for lower_border in bin_list:
                lower_border_index = int(np.where(bin_list_array == lower_border)[0])
                # test only, if the protein is longer than the lower bound of the bin
                if tmp_length > lower_border:
                    # catch cases where the protein is longer than the provided max
                    if tmp_length > max(bin_list):
                        tmp_group = "unknown"
                        # correct protein_group_list in case, it was not corrected before, since it was was created in an earlier loop, the order should not be disturbed and 'unknown' will be the last element
                        if protein_group_list.count("unknown") == 0:
                            protein_group_list.append("unknown")
                    else:
                        # ensure, that the protein is shorter or as long as the upper bound;
                        if bin_list[lower_border_index + 1] >= tmp_length:
                            tmp_group = (protein_group_list[lower_border_index])

            protein_names_list.append(tmp_id)
            protein_length_list.append(tmp_length)
            grouping_column_list.append(tmp_group)

        # merge lists to pandas df
        out_df = DataFrame()
        out_df['Identifier'] = protein_names_list
        out_df['Sequence_Len[aa]'] = protein_length_list
        out_df['Group'] = grouping_column_list


        return out_df, protein_group_list, err_return
    except Exception:
        err_return = False
        return None, None, err_return


def simulate_abundance_cli(param_obj: CoMPaseD_Parameter, protein_df: DataFrame, protein_groups_list: list) -> DataFrame:
    '''set protein abundance for each round of testing to semi-random values'''
    dyn_range = float(param_obj.Protein_dynamic_range)
    leave_out_list = config_to_numeric_list(param_obj.Not_expressed_fraction)
    samplings = int(param_obj.Sampling_Number)

    if len(leave_out_list) != len(protein_groups_list):
        pass # err_handling_function

    # add column with sampling weights for each sampling round and init weight as 1
    for sampling_col in range(1, samplings+1):
        curr_col_name = "Random_sampling_" + str(sampling_col)
        protein_df[curr_col_name] = 1.0

        # define leave_out fraction for each round
        existing_proteins = list()
        for unique_group, unique_leave_out_fraction in zip(protein_groups_list, leave_out_list):
            # subset for current group
            subset_protein_df = protein_df[protein_df.Group == unique_group]
            # calculate the number of proteins to keep/leave_out and correct if all would be removed
            number_to_keep = round(len(subset_protein_df.Group)*((100-unique_leave_out_fraction)/100))
            if number_to_keep < 1:
                number_to_keep = 1
            number_to_discard = len(subset_protein_df.Group) - number_to_keep
            # use shuffled index values to select which proteins are left out
            current_group_index = protein_df[protein_df.Group == unique_group].index.to_list()
            shuffle(current_group_index)
            index_to_zero = list()
            # set number_to_discard indices from current_group_index to zero starting with the last
            for current_index in range(0,number_to_discard):
                index_to_zero.append(current_group_index[-1])
                # remove this index to avoid duplicated zero-setting
                current_group_index.pop(-1)
            protein_df.iloc[index_to_zero, protein_df.columns.get_loc(curr_col_name)] = 0.0
            # keep proteins that are not set to zero
            existing_proteins.extend(current_group_index)

        n_proteins_expressed = len(protein_df.loc[protein_df[curr_col_name] == 1])
        
        # load differentially sized pools by number of samples, small pool should be sufficient for vast majority of applications
        file_location = path.dirname(path.realpath(__file__))
        if n_proteins_expressed > 3000000:
            pool_file =  path.join(file_location, "abundance_pool_large.csv")
        elif n_proteins_expressed > 300000:
            pool_file =  path.join(file_location, "abundance_pool_medium.csv")
        else:
            pool_file =  path.join(file_location, "abundance_pool_small.csv")

        abundance_pool = load_abundance_pool(pool_file = pool_file)

        abundance_weights = get_abundance(n_proteins_expressed, abundance_pool, dyn_range)
        for index_to_modify, protein_weight in zip(existing_proteins, abundance_weights):
            protein_df.iloc[index_to_modify, protein_df.columns.get_loc(curr_col_name)] = protein_weight
    return protein_df


def load_abundance_pool(pool_file: str) -> list:
    '''load abundance sampling pool values from file'''
    pool = list()
    if path.isfile(pool_file):
        with open(pool_file, 'r') as p:
            for line in p:
                pool.append(int(line.rstrip()))
        return pool
    else:
        return list()

def get_abundance(N: int, pool: list, dyn_range: float):
    '''randomly samples a pool of abundance values'''

    # unlikely, but warn in case more than 3e7 proteins are in the fasta file
    if not (N < len(pool)):
        print(f"WARNING: Number of proteins for which abundance values should be generated ({N}) exceeds the number of available values {len(pool)}. Will assign identical abundance to some proteins.")

    # randomly select pool values
    tmp_list = np.random.choice(pool, size = N, replace = True, p = None)

    # define dynamic range cutoff; correct by +2 to keep realistic dynamic range values
    dyn_range_cutoff = 10**(12-dyn_range)
    # split tmp_list by dyn_range_cutoff
    keep_vals = [abund for abund in tmp_list if abund > dyn_range_cutoff]
    replace_vals = [abund for abund in tmp_list if not abund > dyn_range_cutoff]
    # randomise order
    shuffle(replace_vals)
    # select values equivalent to 2-sigma of a normal distribution (95 percent) for down-shifting
    dynamic_range_filter_size = int(round(len(replace_vals) * 0.95))
    # replace 95 percent of the replace_vals by their abundance divided by the number of expressed proteins,
    # using N should compensate for likely higher number of identified peptides for large proteomes
    replace_vals[0:dynamic_range_filter_size] = [int(round(abund/N)) for abund in replace_vals[0:dynamic_range_filter_size]]
    # join both lists again
    tmp_list = keep_vals + replace_vals
    # shuffle to ensure random order
    shuffle(tmp_list)
    # normalise to 1
    tmp_list = [x/sum(tmp_list) for x in tmp_list]

    return tmp_list


'''
Test:
all_test_vals = list()
for i in range(2, 8):
    curr_test_vals = get_abundance(3900, i)
    all_test_vals.append(curr_test_vals)
    sns.kdeplot(curr_test_vals)    

test_vals = get_abundance(3900, 6)
test_vals.sort()
CDF = list()
for i in range(0, len(test_vals)):
    CDF.append(sum(test_vals[0:i]))    
plt.plot(CDF)
'''

def df_to_table(tab_wdg, protein_df = DataFrame):
    # protein_df dims
    n_cols = protein_df.shape[1]
    n_rows = protein_df.shape[0]
    # reset and setup gui table
    tab_wdg.ProteinWeightTable.clear()
    tab_wdg.ProteinWeightTable.setColumnCount(n_cols)
    tab_wdg.ProteinWeightTable.setRowCount(n_rows)
    # set table headers
    headers = list()
    for ColName in protein_df.columns:
        headers.append(ColName)
    tab_wdg.ProteinWeightTable.setHorizontalHeaderLabels(headers)
    # fill table with values from df
    for DF_Col in range(n_cols):
        for DF_Row in range(n_rows):
            # numeric values are stored as np.float64 in pandas df, string conversion
            tab_wdg.ProteinWeightTable.setItem(DF_Row, DF_Col,
                                        QtW.QTableWidgetItem(str(protein_df.iat[DF_Row, DF_Col])))
    # adjust column widths
    table_header = tab_wdg.ProteinWeightTable.horizontalHeader()
    table_header.setSectionResizeMode(QtW.QHeaderView.ResizeMode.ResizeToContents)

    return None

def table_to_df(tab_wdg) -> DataFrame:
    n_cols = tab_wdg.ProteinWeightTable.columnCount()
    n_rows = tab_wdg.ProteinWeightTable.rowCount()
    # prepare empty df
    protein_df = DataFrame(np.empty((n_rows, n_cols), dtype=object))
    # list for column names
    protein_df_columns = list()

    # function to try to convert numeric values to proper floats
    def convert_to_float_type(value_str):
        try:
            return float(value_str)
        except ValueError:
            # if it fails return the value as string
            return value_str

    for DF_Col in range(n_cols):
        # horizontalHeaderItem returns an QTableWidgetItem which needs conversion to text
        protein_df_columns.append(str(tab_wdg.ProteinWeightTable.horizontalHeaderItem(DF_Col).text()))
        for DF_Row in range(n_rows):
            # here conversion to text is necessary as well
            cell_value_str = tab_wdg.ProteinWeightTable.item(DF_Row, DF_Col).text()
            protein_df.iat[DF_Row, DF_Col] = convert_to_float_type(cell_value_str)

    # check if all values in a column are integers and convert if possible
    for DF_Col in range(n_cols):
        all_integers = True
        for DF_Row in range(n_rows):
            value = protein_df.iat[DF_Row, DF_Col]
            if isinstance(value, float) and value.is_integer():
                continue
            elif isinstance(value, (int, float)):
                all_integers = False
                break
            else:
                all_integers = False
                break
        if all_integers:
            protein_df.iloc[:, DF_Col] = protein_df.iloc[:, DF_Col].astype(int)

    # set column names
    protein_df.columns = protein_df_columns
    # convert columns to numeric format where appropriate
    for col in protein_df.columns:
        try:
            protein_df[col] = to_numeric(protein_df[col])
        except ValueError:
            pass

    return protein_df

def run_export(tab_wdg, param_obj: CoMPaseD_Parameter, run_counter: int) -> int:
    # update current params with values from gui
    param_obj.get_params_from_gui(tab_wdg=tab_wdg)
    # make out folders if necessary
    result_check = make_folders(tab_wdg, param_obj)
    # upon error (e.g. from invalid params) exit function)
    if not result_check:
        return run_counter

    tab_wdg.ExportProgressWin = QtW.QProgressDialog("", "", 0, 10, tab_wdg)
    tab_wdg.ExportProgressWin.setCancelButton(None)
    tab_wdg.ExportProgressWin.setWindowTitle("Resampling")
    tab_wdg.ExportProgressWin.setLabelText("<b>Resample proteins from fasta</b><br><br>"
                                               "Please wait...")
    tab_wdg.ExportProgressWin.setStyleSheet("text-align: center")
    tab_wdg.ExportProgressWin.setGeometry(rel_pos(0, 0, 400, 150))
    center_point = QtGui.QGuiApplication.primaryScreen().availableGeometry().center()
    qt_rectangle = tab_wdg.ExportProgressWin.frameGeometry()
    qt_rectangle.moveCenter(center_point)
    tab_wdg.ExportProgressWin.move(qt_rectangle.topLeft())
    qt_rectangle = tab_wdg.ExportProgressWin.frameGeometry()
    tab_wdg.ExportProgressWin.setGeometry(qt_rectangle)
    tab_wdg.ExportProgressWin.setMinimumDuration(2)
    tab_wdg.ExportProgressWin.forceShow()
    # block main window during execution
    tab_wdg.ExportProgressWin.setWindowModality(QtCore.Qt.WindowModality.WindowModal)
    tab_wdg.ExportProgressWin.setValue(1)

    # get list of fasta proteins and group by length-bins
    protein_df, group_list, result_check = load_proteins(tab_wdg, param_obj)
    if not result_check:
        return run_counter
    tab_wdg.ExportProgressWin.setValue(5)

    # assign random abundance values
    #protein_df = simulate_abundance(tab_wdg, param_obj, protein_df, group_list)
    protein_df = simulate_abundance_cli(param_obj, protein_df, group_list)
    tab_wdg.ExportProgressWin.setValue(8)
    # update gui with table values
    df_to_table(tab_wdg, protein_df)
    tab_wdg.ExportProgressWin.setValue(9)
    # use run_counter to check for first or subsequent runs
    run_counter += 1
    sleep(0.3)
    tab_wdg.ExportProgressWin.setValue(10)
    return run_counter


def save_export_funct(tab_wdg, param_obj: CoMPaseD_Parameter):
    # update current params with values from gui
    param_obj.get_params_from_gui(tab_wdg=tab_wdg)
    # load table values to pd.DataFrame
    protein_df = table_to_df(tab_wdg=tab_wdg)

    if tab_wdg.param_file_name == "":
        # auto-save params when not done before
        param_file_name = param_obj.param_silent_save_file(tab_wdg)
        protein_df.to_csv(path.join(param_obj.Protein_weight_file), index=False, sep="\t")
        save_info_msg = "Saved protein weight file to:    " + str(path.join(param_obj.Protein_weight_file))
        tab_wdg.status.showMessage(save_info_msg, 3500)
    else:
        try:
            # try to update current protein weights file field only
            param_obj.Protein_weight_file = param_obj.get_pwf()
            protein_df.to_csv(path.join(param_obj.Protein_weight_file), index=False, sep="\t")
            save_info_msg = "Saved protein weight file to:    " + str(path.join(param_obj.Protein_weight_file))
            tab_wdg.status.showMessage(save_info_msg, 3500)
        except FileNotFoundError:
            # if params were deleted or moved during CoMPaseD exec, save again
            param_file_name = param_obj.param_silent_save_file(tab_wdg)
            protein_df.to_csv(path.join(param_obj.Protein_weight_file), index=False, sep="\t")
            save_info_msg = "Saved protein weight file to:    " + str(path.join(param_obj.Protein_weight_file))
            tab_wdg.status.showMessage(save_info_msg, 3500)


def open_export(tab_wdg, param_obj: CoMPaseD_Parameter, run_counter: int, use_gui_pwf: bool) -> int:
    if not use_gui_pwf:
        # open file dialogue
        if path.isdir(path.dirname(tab_wdg.Protein_weight_file)):
            pwf_name, _ = QtW.QFileDialog.getOpenFileName(tab_wdg, "Select Protein Weights File",
                                                          directory=path.dirname(tab_wdg.Protein_weight_file),
                                                          filter="Protein Weight File (*.tsv)")
        elif path.isdir(tab_wdg.params_obj.Output_directory):
            pwf_name, _ = QtW.QFileDialog.getOpenFileName(tab_wdg, "Select Protein Weights File",
                                                          directory=tab_wdg.params_obj.Output_directory,
                                                          filter="Protein Weight File (*.tsv)")
        else:
            pwf_name, _ = QtW.QFileDialog.getOpenFileName(tab_wdg, "Select Protein Weights File",
                                                          directory=getcwd(),
                                                          filter="Protein Weight File (*.tsv)")
        if pwf_name != "":

            tab_wdg.Protein_weight_file = pwf_name
            tab_wdg.params_obj.Protein_weight_file = tab_wdg.Protein_weight_file
            # update current params with values from gui
            param_obj.get_params_from_gui(tab_wdg=tab_wdg)
            # read abundance values from file
            protein_df = read_csv(param_obj.Protein_weight_file, sep = "\t")
            # update gui with table values
            df_to_table(tab_wdg, protein_df)
            # use run_counter to check for first or subsequent runs
            run_counter += 1

        return run_counter
    else:
        try:
            # update current params with values from gui
            param_obj.get_params_from_gui(tab_wdg=tab_wdg)
            # read abundance values from file
            protein_df = read_csv(param_obj.Protein_weight_file, sep="\t")
            # update gui with table values
            df_to_table(tab_wdg, protein_df)
            # use run_counter to check for first or subsequent runs
            run_counter += 1
            return run_counter
        except Exception as e:
            # should only be printed to console but GUI should remain unaffected
            print(f"ERROR: {e}")
            pass


def open_export_funct(tab_wdg):
    tab_wdg.export_run_counter = open_export(tab_wdg, tab_wdg.params_obj, tab_wdg.export_run_counter, False)
    if not tab_wdg.locked:
        tab_wdg.LoadAnnotButton.setEnabled(True)


def run_first_export(tab_wdg, export_run_counter):
    '''init export according to current params when clicked on export tab
        ensure that this will only be done when params are good by using try-except
    '''
    if path.isfile(tab_wdg.params_obj.Protein_weight_file):
        try:
            tab_wdg.export_run_counter = open_export(tab_wdg, tab_wdg.params_obj, tab_wdg.export_run_counter, True)
            if not tab_wdg.locked:
                tab_wdg.LoadAnnotButton.setEnabled(True)
        except Exception:
            try:
                if (tab_wdg.tabs.currentIndex() == 2) or (tab_wdg.tabs.currentIndex() == 3):
                    tab_idx_backup = tab_wdg.tabs.currentIndex()
                    if export_run_counter == 0:
                        try:
                            tab_wdg.export_run_counter = run_export(tab_wdg, tab_wdg.params_obj, export_run_counter)
                            if not tab_wdg.locked:
                                tab_wdg.LoadAnnotButton.setEnabled(True)
                            tab_wdg.tabs.setCurrentIndex(tab_idx_backup)
                        except FileNotFoundError:
                            pass
                else:
                    pass
            except Exception:
                pass
    else:
        try:
            if (tab_wdg.tabs.currentIndex() == 2) or (tab_wdg.tabs.currentIndex() == 3):
                tab_idx_backup = tab_wdg.tabs.currentIndex()
                if export_run_counter == 0:
                    try:
                        tab_wdg.export_run_counter = run_export(tab_wdg, tab_wdg.params_obj, export_run_counter)
                        if not tab_wdg.locked:
                            tab_wdg.LoadAnnotButton.setEnabled(True)
                        tab_wdg.tabs.setCurrentIndex(tab_idx_backup)
                    except FileNotFoundError:
                        pass
            else:
                pass
        except Exception as e:
            # should only be printed to console but GUI should remain unaffected
            print(f"ERROR: {e}")
            pass


def run_further_export(tab_wdg, export_run_counter):
    if export_run_counter > 0:
        try:
            tab_wdg.export_run_counter = run_export(tab_wdg, tab_wdg.params_obj, export_run_counter)
            if not tab_wdg.locked:
                tab_wdg.LoadAnnotButton.setEnabled(True)
        except FileNotFoundError:
            pass

    else:
        pass


def annotate_pwf_table(tab_wdg):
    if path.isdir(path.dirname(tab_wdg.params_obj.Protein_weight_file)):
        annot_file_name, _ = QtW.QFileDialog.getOpenFileName(tab_wdg, "Load protein annotation file",
                                                             directory=path.dirname(tab_wdg.params_obj.Protein_weight_file),
                                                             filter="Annotation File (*.tsv)")
    else:
        annot_file_name, _ = QtW.QFileDialog.getOpenFileName(tab_wdg, "Load protein annotation file",
                                                             directory=getcwd(),
                                                             filter="Annotation File (*.tsv)")
    if not path.isfile(annot_file_name):
        return
    else:
        annot_df = read_csv(annot_file_name, sep="\t")
        pwf_df = table_to_df(tab_wdg=tab_wdg)
        out_df = merge(left=pwf_df, right=annot_df,
                              on="Identifier",
                              how='left',
                              suffixes=('', '_AnnotDF'))

        out_df.drop(out_df.filter(regex='_AnnotDF$').columns.tolist(), axis=1, inplace=True)
        # get column names in df for grouping col selection, omit Sampling cols and ID
        new_group_column_list = list()
        for ColName in out_df.columns:
            if not str(ColName).startswith("Random_sampling_"):
                if str(ColName) != "Identifier":
                    new_group_column_list.append(ColName)
        # if 'Group' column is present, locate it in the list and set start_item to this index
        start_item = 0
        try:
            start_item = new_group_column_list.index("Group")
        except Exception:
            pass
        # input dialog to return new name and True for ok clicked, False for cancel clicked
        # construction: (parent, window title, message, item list, auto-selection, editability)
        tab_wdg.GroupingColSelection, ok = QtW.QInputDialog.getItem(tab_wdg.ExportTab,
                                                                 "Protein grouping",
                                                                 "Please select proteome sub-set definition "
                                                                 "column",
                                                                 new_group_column_list,
                                                                 start_item,
                                                                 False)
        # change selected column to new Group column if ok was clicked
        if ok and (not tab_wdg.GroupingColSelection == "Group"):
            if "Group" in out_df.columns:
                # change old grouping column name
                if "OldGroup" not in out_df.columns:
                    out_df.rename(columns={'Group': "OldGroup"}, inplace=True)
                # handle multiple annotate tries with changing columns
                else:
                    # generate increasing _number for OldGroup and check existence in out_df
                    counter = 1
                    new_col_name = "OldGroup_" + str(counter)
                    while new_col_name in out_df.columns:
                        counter += 1
                        new_col_name = "OldGroup_" + str(counter)
                    out_df.rename(columns={'Group': new_col_name}, inplace=True)

            out_df.rename(columns={tab_wdg.GroupingColSelection: 'Group'}, inplace=True)
            # sort Group column to third position
            col_names = out_df.columns.to_list()
            col_names.insert(1, col_names.pop(col_names.index('Group')))
            out_df = out_df.loc[:, col_names]

        # update gui table
        df_to_table(tab_wdg, out_df)