import PyQt6.QtWidgets as QtW
# import PyQt6.QtGui as QtGui
# import PyQt6.QtCore as QtCore


# initialise / reset table
def load_protease_table_data(tab_wdg):
    """reset protease table to it's defaults"""
    if tab_wdg.SamplingSizePeptides.isChecked():
        based_on = "Peptides"

    if tab_wdg.SamplingSizeFraction.isChecked():
        based_on = "Fraction"

    tab_wdg.ProteaseTable.clearContents()
    tab_wdg.ProteaseTable.setColumnCount(4)
    tab_wdg.ProteaseTable.setRowCount(5)
    row_number_warning_label = ""
    tab_wdg.RowNumberWarning.setText(row_number_warning_label)
    if based_on == "Peptides":
        tab_wdg.ProteaseTable.setHorizontalHeaderLabels(
            ["Protease", "Max MCs", "MC Frequency", "Number of Peptides"])
        default_bases_on_peptides = [
            {'Protease': "trypsin", 'Max_MCs': "2",
             'MC_Frequencies': "[0.7415,0.2090,0.0484]", '# Peptides': "10000"},
            {'Protease': "lysarginase", 'Max_MCs': "2",
             'MC_Frequencies': "[0.5757,0.2899,0.1336]", '# Peptides': "10000"},
            {'Protease': "glu-c", 'Max_MCs': "5",
             'MC_Frequencies': "[0.5620,0.2753,0.1110,0.0419,0.0086,0.0012]", '# Peptides': "10000"},
            {'Protease': "chymotrypsin", 'Max_MCs': "5",
             'MC_Frequencies': "[0.2002,0.3369,0.2648,0.1471,0.0498,0.0012]", '# Peptides': "10000"},
            {'Protease': "lys-c", 'Max_MCs': "2",
             'MC_Frequencies': "[0.9102,0.0836,0.0058]", '# Peptides': "10000"}]
        row_counter = 0
        for Protease in default_bases_on_peptides:
            tab_wdg.ProteaseTable.setItem(row_counter, 0, QtW.QTableWidgetItem(Protease["Protease"]))
            tab_wdg.ProteaseTable.setItem(row_counter, 1, QtW.QTableWidgetItem(Protease["Max_MCs"]))
            tab_wdg.ProteaseTable.setItem(row_counter, 2, QtW.QTableWidgetItem(Protease["MC_Frequencies"]))
            tab_wdg.ProteaseTable.setItem(row_counter, 3, QtW.QTableWidgetItem(Protease["# Peptides"]))
            row_counter += 1

    else:
        tab_wdg.ProteaseTable.setHorizontalHeaderLabels(
            ["Protease", "Max MCs", "MC Frequency", "Fraction of Peptides"])
        default_bases_on_fraction = [
            {'Protease': "trypsin", 'Max_MCs': "2",
             'MC_Frequencies': "[0.7415,0.2090,0.0484]", 'ID Fraction': "0.033051"},
            {'Protease': "lysarginase", 'Max_MCs': "2",
             'MC_Frequencies': "[0.5757,0.2899,0.1336]", 'ID Fraction': "0.031973356"},
            {'Protease': "glu-c", 'Max_MCs': "5",
             'MC_Frequencies': "[0.5620,0.2753,0.1110,0.0419,0.0086,0.0012]", 'ID Fraction': "0.014570424"},
            {'Protease': "chymotrypsin", 'Max_MCs': "5",
             'MC_Frequencies': "[0.2002,0.3369,0.2648,0.1471,0.0498,0.0012]", 'ID Fraction': "0.009367681"},
            {'Protease': "lys-c", 'Max_MCs': "2",
             'MC_Frequencies': "[0.9102,0.0836,0.0058]", 'ID Fraction': "0.053870932"}]
        row_counter = 0
        for Protease in default_bases_on_fraction:
            tab_wdg.ProteaseTable.setItem(row_counter, 0, QtW.QTableWidgetItem(Protease["Protease"]))
            tab_wdg.ProteaseTable.setItem(row_counter, 1, QtW.QTableWidgetItem(Protease["Max_MCs"]))
            tab_wdg.ProteaseTable.setItem(row_counter, 2, QtW.QTableWidgetItem(Protease["MC_Frequencies"]))
            tab_wdg.ProteaseTable.setItem(row_counter, 3, QtW.QTableWidgetItem(Protease["ID Fraction"]))
            row_counter += 1


# add rows to protease table
def add_row(tab_wdg):
    tab_wdg.ProteaseTable.setRowCount(tab_wdg.ProteaseTable.rowCount() + 1)
    row_number_warning_label = ""
    if tab_wdg.ProteaseTable.rowCount() > 10:
        row_number_warning_label = "WARNING: Comparing a large number of proteases results in long " \
                                   "processing times. Consider to reduce the proteases of interest."
    tab_wdg.RowNumberWarning.setText(row_number_warning_label)

# remove rows from protease table
def subtract_row(tab_wdg):
    tab_wdg.ProteaseTable.setRowCount(tab_wdg.ProteaseTable.rowCount() - 1)
    row_number_warning_label = ""
    if tab_wdg.ProteaseTable.rowCount() > 10:
        row_number_warning_label = "WARNING: Comparing a large number of proteases results in long " \
                                   "processing times. Consider to reduce the proteases of interest."
    tab_wdg.RowNumberWarning.setText(row_number_warning_label)


# set table with trypsin only and based on peptide number
def clear_table(tab_wdg):
    """empty protease table and set trypsin as the only protease, set sampling based on peptide number"""
    tab_wdg.ProteaseTable.clearContents()
    tab_wdg.ProteaseTable.setColumnCount(4)
    tab_wdg.ProteaseTable.setRowCount(2)
    tab_wdg.ProteaseTable.setHorizontalHeaderLabels(
        ["Protease", "Max MCs", "MC Frequency", "Number of Peptides"])
    trypsin_data = [
        {'Protease': "trypsin", 'Max_MCs': "2",
         'MC_Frequencies': "[0.7415,0.2090,0.0484]", '# Peptides': "10000"}]
    row_counter = 0
    for Protease in trypsin_data:
        tab_wdg.ProteaseTable.setItem(row_counter, 0, QtW.QTableWidgetItem(Protease["Protease"]))
        tab_wdg.ProteaseTable.setItem(row_counter, 1, QtW.QTableWidgetItem(Protease["Max_MCs"]))
        tab_wdg.ProteaseTable.setItem(row_counter, 2, QtW.QTableWidgetItem(Protease["MC_Frequencies"]))
        tab_wdg.ProteaseTable.setItem(row_counter, 3, QtW.QTableWidgetItem(Protease["# Peptides"]))
        row_counter += 1
    tab_wdg.SamplingSizePeptides.setChecked(True)
    row_number_warning_label = ""
    tab_wdg.RowNumberWarning.setText(row_number_warning_label)


def header_based_on_peptides(tab_wdg):
    tab_wdg.ProteaseTable.setHorizontalHeaderLabels(
        ["Protease", "Max MCs", "MC Frequency", "Number of Peptides"])


def header_based_on_fraction(tab_wdg):
    tab_wdg.ProteaseTable.setHorizontalHeaderLabels(
        ["Protease", "Max MCs", "MC Frequency", "Fraction of Peptides"])
