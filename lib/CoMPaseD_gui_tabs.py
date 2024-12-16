import shutil
import colorama
import re
from os import rename, stat, path
from datetime import datetime
from sys import executable
from PyQt6.QtSvgWidgets import QSvgWidget

from lib.CoMPaseD_gui_protease_table_functions import *
from lib.CoMPaseD_gui_config_functions import *
from lib.CoMPaseD_gui_param_functions import *
from lib.CoMPaseD_gui_export_functions import *
from lib.CoMPaseD_gui_result_plot import *


class CoMPaseD_Tabs(QtW.QWidget):
    """class to define tab layout"""

    def __init__(self, parent):
        # init as Qwidget and add basic layout
        super(QtW.QWidget, self).__init__(parent)
        super().__init__(parent)
        self.layout = QtW.QVBoxLayout(self)

        # prepare QTabWidget as tabs and add QWidgets as container for the actual tabs
        self.tabs = QtW.QTabWidget(self)
        self.tabs.setGeometry(rel_pos(0, 0, 1200, 900))
        self.tabs.setTabShape(QtW.QTabWidget.TabShape.Triangular)
        self.tabs.setIconSize(QtCore.QSize(100, 100))

        # prepare tabs
        self.ConfigTab = QtW.QWidget(self.tabs)
        self.ParamTab = QtW.QWidget(self.tabs)
        self.ExportTab = QtW.QWidget(self.tabs)
        self.ProgressTab = QtW.QWidget(self.tabs)
        self.ResultsTab = QtW.QWidget(self.tabs)

        # add status bar
        self.status = QtW.QStatusBar(parent)
        parent.setStatusBar(self.status)
        self.status.setSizeGripEnabled(False)
        # # setup status bar
        self.locked_msg = QtW.QLabel()
        self.locked_msg.setText("")
        self.status.addWidget(self.locked_msg)

        ##############################
        # setup configuration tab

        # # set external programs box
        self.FilePathBox = QtW.QGroupBox(parent=self.ConfigTab)
        self.FilePathBox.setTitle("")
        self.FilePathBox.setGeometry(rel_pos(30, 40, 800, 50))

        # # set DeepMSPeptide box
        self.PerlFilePathBox = QtW.QGroupBox(parent=self.ConfigTab)
        self.PerlFilePathBox.setTitle(" ")
        self.PerlFilePathBox.setCheckable(True)
        self.PerlFilePathBox.setChecked(False)
        self.PerlFilePathBox.setToolTip("Try this option if memory issues occur during the <em> in-silico </em><br>"
                                        "digestion of a very large fasta database.<br>"
                                        "The original perl-based implementation of ProteoMapper <br>"
                                        "is more memory efficient but slower. A working installation <br>"
                                        "of the Trans-Proteomics-Pipeline may be required for the <br>"
                                        "perl scripts to run properly. <br>"
                                        "<b>Default: Do NOT use this option!</b>")
        self.PerlFilePathBox.setGeometry(rel_pos(30, 140, 800, 120))
        self.PerlFilePathBoxTitle = QtW.QLabel(self.PerlFilePathBox)
        # self.PerlFilePathBoxTitle.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.LinksAccessibleByMouse)
        self.PerlFilePathBoxTitle.setText("Use original Perl-based <a style='color: #031059; font: bold 14px' href='http://www.tppms.org/tools/pm/'>ProteoMapper</a> for peptide mapping.")
        self.PerlFilePathBoxTitle.setOpenExternalLinks(True)
        self.PerlFilePathBoxTitle.setProperty("class", "large_lable_title")
        self.PerlFilePathBoxTitle.setGeometry(rel_pos(30, -3, 600, 25))

        # # # generate all labels
        self.CruxPathLabel = QtW.QLabel(self.FilePathBox)
        self.ClipsPathLabel = QtW.QLabel(self.PerlFilePathBox)
        self.PromastPathLabel = QtW.QLabel(self.PerlFilePathBox)
        # # # # setup labels
        self.CruxPathLabel.setText("<a  style='color: #031059; font: bold 14px' href='https://www.crux.ms'>Crux</a> path")
        self.CruxPathLabel.setOpenExternalLinks(True)
        self.ClipsPathLabel.setText("Clips path")
        self.PromastPathLabel.setText("Promast path")
        self.CruxPathLabel.setProperty("class", "large_lable")
        self.ClipsPathLabel.setProperty("class", "normal_lable")
        self.PromastPathLabel.setProperty("class", "normal_lable")
        # # # # add tooltips
        self.CruxPathLabel.setToolTip(
            "Path to Crux executable. Should be version >= 3.2.\n"
            "See 'https://crux.ms/' for further information.")
        self.ClipsPathLabel.setToolTip(
            "Path to Clips perl script from ProteoMapper.\n"
            "A suitable version can be found in the /bin folder\n"
            "of CoMPaseD.\nSee 'http://www.tppms.org/tools/pm/'\n"
            "for further information.")
        self.PromastPathLabel.setToolTip(
            "Path to Promast perl script from ProteoMapper.\n"
            "A suitable version can be found in the /bin folder\n"
            "of CoMPaseD.\nSee 'http://www.tppms.org/tools/pm/'\n"
            "for further information.")
        # # # # position label
        self.CruxPathLabel.setGeometry(rel_pos(8, 20, 90, 25))
        self.ClipsPathLabel.setGeometry(rel_pos(50, 30, 90, 25))
        self.PromastPathLabel.setGeometry(rel_pos(50, 70, 90, 25))

        # # # generate all entry boxes
        self.CruxPathField = QtW.QLineEdit(self.FilePathBox)
        self.ClipsPathField = QtW.QLineEdit(self.PerlFilePathBox)
        self.PromastPathField = QtW.QLineEdit(self.PerlFilePathBox)
        # # # # setup entry boxes with default text that will be displayed when no configuration is found
        self.CruxPathField.setText("Crux path")
        self.ClipsPathField.setText("Clips path")
        self.PromastPathField.setText("Promast path")
        # # # # position entry boxes
        self.CruxPathField.setGeometry(rel_pos(120, 20, 600, 25))
        self.ClipsPathField.setGeometry(rel_pos(150, 30, 570, 25))
        self.PromastPathField.setGeometry(rel_pos(150, 70, 570, 25))

        # # # generate browse buttons
        self.CruxBrowseButton = QtW.QPushButton(self.FilePathBox)
        self.ClipsBrowseButton = QtW.QPushButton(self.PerlFilePathBox)
        self.PromastBrowseButton = QtW.QPushButton(self.PerlFilePathBox)
        # # # # setup browse buttons
        self.CruxBrowseButton.setText("Browse")
        self.ClipsBrowseButton.setText("Browse")
        self.PromastBrowseButton.setText("Browse")
        self.CruxBrowseButton.setProperty("class", "small_btn")
        self.ClipsBrowseButton.setProperty("class", "small_btn")
        self.PromastBrowseButton.setProperty("class", "small_btn")
        # # # # position buttons, move 1 pt up to compensate for shadow
        self.CruxBrowseButton.setGeometry(rel_pos(730, 19, 60, 25))
        self.ClipsBrowseButton.setGeometry(rel_pos(730, 29, 60, 25))
        self.PromastBrowseButton.setGeometry(rel_pos(730, 69, 60, 25))

        # # set further settings box
        self.SettingsBox = QtW.QGroupBox(parent=self.ConfigTab)
        self.SettingsBox.setTitle("")
        self.SettingsBox.setGeometry(rel_pos(30, 255, 800, 150))

        # # # generate checkboxes
        self.SamplingOutputCheckbox = QtW.QCheckBox(self.SettingsBox)
        # # # # setup checkboxes
        self.SamplingOutputCheckbox.setText("Output peptide lists from random sampling?")
        self.SamplingOutputCheckbox.setChecked(True)
        self.SamplingOutputCheckbox.setProperty("class", "large_chkbx")
        # # # # add checkbox tooltips
        self.SamplingOutputCheckbox.setToolTip(
            "<b>Default: On</b><br><br>May be disabled for very large sampling numbers.")
        # # # # position checkboxes
        self.SamplingOutputCheckbox.setGeometry(rel_pos(8, 35, 800, 25))

        # # # add save and load buttons
        self.SaveConfigButton = QtW.QPushButton(parent=self.ConfigTab)
        self.LoadConfigButton = QtW.QPushButton(parent=self.ConfigTab)
        # # # # setup save and load buttons
        self.SaveConfigButton.setText("Save configuration")
        self.LoadConfigButton.setText("Load configuration")
        # # # # add tooltips to load and save button
        self.SaveConfigButton.setToolTip(
            "Save current settings as default.\n\n"
            "By default the configuration file is located \nin the current "
            "working directory and called 'CoMPaseD_Config.txt'.")
        self.LoadConfigButton.setToolTip(
            "Load last saved settings.\n\n"
            "By default the configuration file is located \nin the current "
            "working directory and called 'CoMPaseD_Config.txt'.")
        # # # # position save and load button
        self.SaveConfigButton.setGeometry(rel_pos(490, 400, 150, 35))
        self.LoadConfigButton.setGeometry(rel_pos(680, 400, 150, 35))

        # # # add placeholder for error message
        self.ErrorLabelCorruptedFile = QtW.QLabel(parent=self.ConfigTab)
        self.ErrorLabelCorruptedFile.setGeometry(rel_pos(500, 440, 400, 35))
        self.ErrorLabelCorruptedFile.setText("")
        self.ErrorLabelCorruptedFile.setProperty("class", "err_lable")

        # add CoMPaseD logo
        # ToDo: fix svg rendering shadows!
        file_location = path.dirname(path.realpath(__file__))
        bg_img = path.join(file_location, "../bin", "CoMPaseD_logo.svg")
        if path.isfile(bg_img):
            self.bg_img = QSvgWidget(bg_img, parent=self.ConfigTab)
            self.bg_img.setGeometry(rel_pos(30, 400, 450, 400))
            self.bg_img.show()


        ##############################
        # setup parameter tab

        # # # set fasta file and output folder labels
        self.FastaPathLabel = QtW.QLabel(parent=self.ParamTab)
        self.OutputPathLabel = QtW.QLabel(parent=self.ParamTab)
        # # # # setup fasta file and output folder labels
        self.FastaPathLabel.setText("Fasta file")
        self.OutputPathLabel.setText("Output folder")
        self.FastaPathLabel.setProperty("class", "normal_lable")
        self.OutputPathLabel.setProperty("class", "normal_lable")
        # # # # add tooltips to fasta file and output folder labels
        self.FastaPathLabel.setToolTip(
            "Select database file.<BR><BR>Do <b>not</b> use a database with decoy entries.")
        self.OutputPathLabel.setToolTip("Select result output folder.")
        # # # # position fasta file and output folder labels
        self.FastaPathLabel.setGeometry(rel_pos(30, 40, 90, 22))
        self.OutputPathLabel.setGeometry(rel_pos(30, 80, 90, 22))

        # # # generate entry boxes for fasta file and output folder
        self.FastaPathField = QtW.QLineEdit(self.ParamTab)
        self.OutputPathField = QtW.QLineEdit(self.ParamTab)
        # # # # setup entry boxes with default text
        self.FastaPathField.setText("Fasta path")
        self.OutputPathField.setText("Output folder")
        # # # # position entry boxes
        self.FastaPathField.setGeometry(rel_pos(120, 40, 600, 22))
        self.OutputPathField.setGeometry(rel_pos(120, 80, 600, 22))

        # # # generate browse buttons for fasta file and output folder
        self.FastaPathButton = QtW.QPushButton(self.ParamTab)
        self.OutputPathButton = QtW.QPushButton(self.ParamTab)
        # # # # setup browse buttons
        self.FastaPathButton.setText("Browse")
        self.OutputPathButton.setText("Browse")
        self.FastaPathButton.setProperty("class", "small_btn")
        self.OutputPathButton.setProperty("class", "small_btn")
        # # # # position buttons, move 2 pt up to compensate for shadow
        self.FastaPathButton.setGeometry(rel_pos(730, 38, 60, 25))
        self.OutputPathButton.setGeometry(rel_pos(730, 78, 60, 25))

        # # set protease table box
        self.ProteaseBox = QtW.QGroupBox(parent=self.ParamTab)
        self.ProteaseBox.setTitle("Proteases")
        self.ProteaseBox.setGeometry(rel_pos(30, 130, 760, 350))
        # # # # add label for warning message in case protease table has more than 10 rows
        self.RowNumberWarning = QtW.QLabel(parent=self.ProteaseBox)
        self.RowNumberWarning.setGeometry(rel_pos(650, 45, 85, 300))
        self.RowNumberWarning.setText("")
        self.RowNumberWarning.setWordWrap(True)
        self.RowNumberWarning.setProperty("class", "err_lable")

        # # + generate table
        self.ProteaseTable = QtW.QTableWidget(parent=self.ProteaseBox)

        # # set box for sampling size basis
        self.SamplingSizeBasisBox = QtW.QGroupBox(parent=self.ProteaseBox)
        self.SamplingSizeBasisBox.setTitle("Sampling Size Based On")
        self.SamplingSizeBasisBox.setToolTip(
            "The number of peptides randomly sampled for each protease <br>"
            "can be based on the fraction of all <b>unique peptides</b> (Peptidefraction) or <br>"
            "be an <b>absolute number</b> (Peptides). <br><br>Provide 'Peptidefraction' as factor, "
            "<em>i.e.</em> a value less than 1. "
        )
        self.SamplingSizeBasisBox.setGeometry(rel_pos(585, 25, 150, 100))

        # # # add radiobuttons for selection
        self.SamplingSizePeptides = QtW.QRadioButton(parent=self.SamplingSizeBasisBox)
        self.SamplingSizeFraction = QtW.QRadioButton(parent=self.SamplingSizeBasisBox)
        # # # # setup radiobuttons for selection
        self.SamplingSizePeptides.setText("Peptides")
        self.SamplingSizeFraction.setText("Peptidefraction")
        self.SamplingSizePeptides.setChecked(True)
        # # # # position radiobuttons
        self.SamplingSizePeptides.setGeometry(rel_pos(10, 30, 140, 22))
        self.SamplingSizeFraction.setGeometry(rel_pos(10, 60, 140, 22))

        # # # set buttons to add/subtract rows
        self.AddRowButton = QtW.QPushButton(parent=self.ProteaseBox)
        self.SubtractRowButton = QtW.QPushButton(parent=self.ProteaseBox)
        # # # # setup add/subtract row buttons
        self.AddRowButton.setText("+ Row")
        self.SubtractRowButton.setText("- Row")
        self.AddRowButton.setProperty("class", "small_btn")
        self.SubtractRowButton.setProperty("class", "small_btn")
        # # # # add tooltips to add/subtract row buttons
        self.AddRowButton.setToolTip(
            "Add row to protease table. <br><br>Possible protease names are:<br> "
            "trypsin', 'trypsin/p', 'chymotrypsin', 'elastase', 'clostripain', <br>"
            "'cyanogen-bromide', 'iodosobenzoate', 'proline-endopeptidase',<br>"
            "'staph-protease', 'asp-n', 'lys-c', <br>'lys-n', 'arg-c', 'glu-c',<br>"
            "'pepsin-a', 'elastase-trypsin-chymotrypsin', 'lysarginase'.<br><br>"
            "Please see crux documentation for further details.")
        self.SubtractRowButton.setToolTip(
            "Remove last row from protease table.")
        # # # # position add/subtract row buttons
        self.AddRowButton.setGeometry(rel_pos(585, 170, 60, 25))
        self.SubtractRowButton.setGeometry(rel_pos(585, 200, 60, 25))

        # # # add reset and load default buttons
        self.ProteaseClearButton = QtW.QPushButton(parent=self.ProteaseBox)
        self.ProteaseResetButton = QtW.QPushButton(parent=self.ProteaseBox)
        # # # # setup reset and load default buttons
        self.ProteaseClearButton.setText("Clear table")
        self.ProteaseResetButton.setText("Reset table")
        self.ProteaseClearButton.setProperty("class", "small_btn")
        self.ProteaseResetButton.setProperty("class", "small_btn")
        # # # # add tooltips to  reset and load default buttons
        self.ProteaseClearButton.setToolTip(
            "Empty table.<br><br>Trypsin must not be removed from the comparison and thus is kept.")
        self.ProteaseResetButton.setToolTip(
            "Reset table to default values.")
        # # # # position buttons
        self.ProteaseClearButton.setGeometry(rel_pos(585, 265, 120, 25))
        self.ProteaseResetButton.setGeometry(rel_pos(585, 305, 120, 25))

        # ## load initial data to protease table and setup column width's
        load_protease_table_data(self)
        self.ProteaseTable.setColumnWidth(0, 80)
        self.ProteaseTable.setColumnWidth(1, 70)
        self.ProteaseTable.setColumnWidth(2, 250)
        self.ProteaseTable.setColumnWidth(3, 120)
        self.ProteaseTable.setGeometry(rel_pos(20, 30, 537, 300))

        # # set Score weights box
        self.ScoreWeightsBox = QtW.QGroupBox(parent=self.ParamTab)
        self.ScoreWeightsBox.setTitle("Protease score weights")
        self.ScoreWeightsBox.setToolTip("Protease score is calculated as the weighted geometric mean of <br><br>"
                                        "a) the number of identified proteins;<br>"
                                        "b) the number of identified peptides;<br>"
                                        "c) the average sequence coverage<br><br>"
                                        "of all proteins within a bin of proteins. Each<br>"
                                        "value is normalised to the result observed for trypsin.")
        self.ScoreWeightsBox.setGeometry(rel_pos(30, 500, 760, 90))

        # # # set weight parameter labels
        self.ProtIDWeightLabel = QtW.QLabel(parent=self.ScoreWeightsBox)
        self.PepIDWeightLabel = QtW.QLabel(parent=self.ScoreWeightsBox)
        self.CoverageWeightLabel = QtW.QLabel(parent=self.ScoreWeightsBox)
        # # # # setup weight parameter labels
        self.ProtIDWeightLabel.setText("# Protein IDs")
        self.PepIDWeightLabel.setText("# Peptide IDs")
        self.CoverageWeightLabel.setText("Sequence coverage")
        self.ProtIDWeightLabel.setProperty("class", "normal_lable")
        self.PepIDWeightLabel.setProperty("class", "normal_lable")
        self.CoverageWeightLabel.setProperty("class", "normal_lable")
        # # # # position weight parameter labels
        self.ProtIDWeightLabel.setGeometry(rel_pos(30, 25, 150, 22))
        self.PepIDWeightLabel.setGeometry(rel_pos(320, 25, 150, 22))
        self.CoverageWeightLabel.setGeometry(rel_pos(580, 25, 150, 22))

        # # # set weight parameter boxes
        self.ProtIDWeightField = QtW.QLineEdit(parent=self.ScoreWeightsBox)
        self.PepIDWeightField = QtW.QLineEdit(parent=self.ScoreWeightsBox)
        self.CoverageWeightField = QtW.QLineEdit(parent=self.ScoreWeightsBox)
        # # # # setup weight parameter boxes
        self.ProtIDWeightField.setText("1.0")
        self.PepIDWeightField.setText("1.0")
        self.CoverageWeightField.setText("1.0")
        self.ProtIDWeightField.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        self.PepIDWeightField.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        self.CoverageWeightField.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        self.ProtIDWeightField.setProperty("class", "score_weight_field")
        self.PepIDWeightField.setProperty("class", "score_weight_field")
        self.CoverageWeightField.setProperty("class", "score_weight_field")
        # # # # position weight parameter labels
        self.ProtIDWeightField.setGeometry(rel_pos(30, 50, 100, 22))
        self.PepIDWeightField.setGeometry(rel_pos(320, 50, 100, 22))
        self.CoverageWeightField.setGeometry(rel_pos(580, 50, 100, 22))

        # # set protein binning box
        self.ProteinBinningBox = QtW.QGroupBox(parent=self.ParamTab)
        self.ProteinBinningBox.setTitle("Protein binning")
        self.ProteinBinningBox.setGeometry(rel_pos(825, 130, 300, 120))
        # # # set protein binning labels
        self.ProteinBinsLabel = QtW.QLabel(parent=self.ProteinBinningBox)
        self.ProteinNotExprFracLabel = QtW.QLabel(parent=self.ProteinBinningBox)
        # # # # setup protein binning labels
        self.ProteinBinsLabel.setText("Protein size binning (aa)")
        self.ProteinNotExprFracLabel.setText("Undetectable fraction (%)")
        self.ProteinBinsLabel.setProperty("class", "normal_lable")
        self.ProteinNotExprFracLabel.setProperty("class", "normal_lable")
        # # # # add tooltips to protein binning labels
        self.ProteinBinsLabel.setToolTip("Define protein binning by number<br>of amino acids.<br><br>"
                                         "Protease scores will be calculated<br>"
                                         "separately for each bin. Provide a comma-separated<br>"
                                         "list of bin stops, <em>e.g.</em><br>"
                                         "'0,50,100,99999'<br>"
                                         "will create three bins containing proteins:<br>"
                                         "- up to 50 aa length (bin_1)<br>"
                                         "- 51 to 100 aa length (bin_2)<br>"
                                         "- 101 to 99999 aa length (bin_3)")
        self.ProteinNotExprFracLabel.setToolTip("Set fraction of undetectable proteins<br>"
                                                "for each bin [%].<br><br>"
                                                "The number of values must fit the number<br>"
                                                "of bins defined above. Adjust this value<br>"
                                                "carefully in case a targeted enrichment of <em>e.g.</em><br>"
                                                "small proteins was applied or when the<br>"
                                                "database contains mostly protein entries<br>"
                                                "that are never translated (<em>in-silico</em> 6-frame translations).")
        # # # # position protein binning labels
        self.ProteinBinsLabel.setGeometry(rel_pos(145, 25, 215, 22))
        self.ProteinNotExprFracLabel.setGeometry(rel_pos(145, 75, 215, 22))

        # # # set protein binning fields
        self.ProteinBinsField = QtW.QLineEdit(parent=self.ProteinBinningBox)
        self.ProteinNotExprFracField = QtW.QLineEdit(parent=self.ProteinBinningBox)
        # # # # setup protein binning fields
        self.ProteinBinsField.setText("0,50,100,99999")
        self.ProteinNotExprFracField.setText("40,30,20")
        # # # # position protein binning field styles
        self.ProteinBinsField.setGeometry(rel_pos(20, 25, 120, 22))
        self.ProteinNotExprFracField.setGeometry(rel_pos(20, 75, 120, 22))

        # # set DeepMSPeptide box
        self.DMSPBox = QtW.QGroupBox(parent=self.ParamTab)
        self.DMSPBox.setTitle("Deep MS-Peptide prediction")
        self.DMSPBox.setCheckable(True)
        self.DMSPBox.setChecked(False)
        self.DMSPBox.setToolTip("Enable Deep-MS-Peptide detectability prediction<br>"
                                "for random peptide sampling.<br><br>"
                                "See DOI: '10.1093/bioinformatics/btz708' for details.")
        self.DMSPBox.setGeometry(rel_pos(825, 280, 300, 200))

        # # # set DeepMSPeptide model path label, field and button
        self.DMSPModelPathField = QtW.QLineEdit(parent=self.DMSPBox)
        self.DMSPModelPathLabel = QtW.QLabel(parent=self.DMSPBox)
        self.DMSPModelPathButton = QtW.QPushButton(parent=self.DMSPBox)
        self.DMSPModelPathLabel.setText("Deep MS-Peptide model")
        self.DMSPModelPathLabel.setProperty("class", "normal_lable")
        self.DMSPModelPathButton.setText("Browse")
        self.DMSPModelPathButton.setProperty("class", "small_btn")
        # # # # add tooltip to DeepMSPeptide model path label
        dmsp_path_tooltip = "Select Deep-MS-Peptide model.<br><br>" \
                            "An alternative prediction model,<br>" \
                            "trained on the <em>Confetti dataset</em><br>" \
                            "(Pride accession number: PXD000900,<br>" \
                            "see DOI: '10.1074/mcp.M113.035170' for details)<br>" \
                            "is available at CoMPaseD’s GitHub repository."
        self.DMSPModelPathField.setToolTip(dmsp_path_tooltip)
        self.DMSPModelPathLabel.setToolTip(dmsp_path_tooltip)
        self.DMSPModelPathButton.setToolTip(dmsp_path_tooltip)
        # # # # position DeepMSPeptide model path label, field and button
        self.DMSPModelPathLabel.setGeometry(rel_pos(20, 40, 200, 25))
        self.DMSPModelPathField.setGeometry(rel_pos(20, 70, 190, 25))
        self.DMSPModelPathButton.setGeometry(rel_pos(225, 70, 60, 25))

        # # # set DeepMSPeptide weight spinbox and label
        self.DMSPWeightSpinbox = QtW.QDoubleSpinBox(parent=self.DMSPBox)
        self.DMSPWeightLabel = QtW.QLabel(parent=self.DMSPBox)
        # # # # setup DeepMSPeptide weight spinbox and label
        self.DMSPWeightSpinbox.setDecimals(1)
        self.DMSPWeightSpinbox.setValue(1.0)
        self.DMSPWeightSpinbox.setRange(0.1, 8)
        self.DMSPWeightSpinbox.setSingleStep(0.1)
        self.DMSPWeightSpinbox.setAccelerated(True)
        self.DMSPWeightLabel.setText("Prediction weight")
        self.DMSPWeightLabel.setProperty("class", "normal_lable")
        # add tooltip to DeepMSPeptide model path label
        self.DMSPWeightLabel.setToolTip("Set Deep-MS-Peptide weight during random sampling.<br><br>"
                                        "Deep-MS-Peptide score of each peptide will be potentiated<br>"
                                        "by this value. (Default = 4.0, Minimum = 1.0, Maximum = 8.0)")
        self.DMSPWeightLabel.setGeometry(rel_pos(20, 100, 200, 25))
        self.DMSPWeightSpinbox.setGeometry(rel_pos(20, 130, 80, 25))

        # # set Further settings box
        self.FurtherSettingsBox = QtW.QGroupBox(parent=self.ParamTab)
        self.FurtherSettingsBox.setTitle("Further settings")
        self.FurtherSettingsBox.setGeometry(rel_pos(825, 500, 300, 200))
        # # # set labels for further settings
        self.MaxProteasesLabel = QtW.QLabel(parent=self.FurtherSettingsBox)
        self.RandomSamplingsLabel = QtW.QLabel(parent=self.FurtherSettingsBox)
        self.DynamicRangeLabel = QtW.QLabel(parent=self.FurtherSettingsBox)
        # # # # setup labels for further settings
        self.MaxProteasesLabel.setText("Maximal proteases")
        self.RandomSamplingsLabel.setText("Random samplings")
        self.DynamicRangeLabel.setText("Dynamic range")
        self.MaxProteasesLabel.setProperty("class", "normal_lable")
        self.RandomSamplingsLabel.setProperty("class", "normal_lable")
        self.DynamicRangeLabel.setProperty("class", "normal_lable")
        # # # # position labels for further settings
        self.MaxProteasesLabel.setGeometry(rel_pos(110, 40, 200, 25))
        self.RandomSamplingsLabel.setGeometry(rel_pos(110, 90, 200, 25))
        self.DynamicRangeLabel.setGeometry(rel_pos(110, 140, 200, 25))
        # # # set spin-boxes for further settings
        self.MaxProteasesSpinbox = QtW.QSpinBox(parent=self.FurtherSettingsBox)
        self.RandomSamplingsSpinbox = QtW.QSpinBox(parent=self.FurtherSettingsBox)
        self.DynamicRangeSpinbox = QtW.QDoubleSpinBox(parent=self.FurtherSettingsBox)
        # # # # setup spin-boxes for further settings
        self.MaxProteasesSpinbox.setValue(5)
        self.MaxProteasesSpinbox.setMinimum(1)
        self.MaxProteasesSpinbox.setMaximum(20)
        self.MaxProteasesSpinbox.setWrapping(True)
        self.RandomSamplingsSpinbox.setValue(10)
        self.RandomSamplingsSpinbox.setMinimum(1)
        self.RandomSamplingsSpinbox.setMaximum(100)
        self.RandomSamplingsSpinbox.setWrapping(True)
        self.DynamicRangeSpinbox.setValue(6.5)
        self.DynamicRangeSpinbox.setDecimals(1)
        self.DynamicRangeSpinbox.setSingleStep(0.1)
        self.DynamicRangeSpinbox.setMinimum(0.0)
        self.DynamicRangeSpinbox.setMaximum(10.0)
        self.DynamicRangeSpinbox.setWrapping(True)
        # # # # add tooltips to spin-boxes and labels of further settings
        max_proteases_tooltip = "Maximal number of proteases to use during the experiment.<br><br>" \
                                "Larger numbers will usually improve the protease score<br>" \
                                "but at the cost of increased measurement time and<br>" \
                                "sample consumption. Do not select unrealistically<br>" \
                                "large values to avoid long running time of CoMPaseD.<br>" \
                                "Must not be larger than the number of proteases<br>" \
                                "provided in the table.<br>(Default = 5, Maximum = 20)"
        random_sampling_tooltip = "Number of random samplings to perform by CoMPaseD.<br><br>" \
                                  "Protein abundances are randomly initialised for<br>" \
                                  "each run. Averaged scores are reported together with<br>" \
                                  "their standard deviation in the main output table<br>" \
                                  "for each protease combination.<br>" \
                                  "Higher values will increase running time linearly. (Default = 10, Maximum = 100)"
        dynamic_range_tooltip = "Set the dynamic range of the abundance of<br>" \
                                "detectable proteins in orders of magnitude.<br><br>" \
                                "Higher values result in clustering of the randomly sampled peptides<br>" \
                                "around fewer proteins. Select 0 to disable. (Default = 6, Maximum = 10)"
        self.MaxProteasesLabel.setToolTip(max_proteases_tooltip)
        self.MaxProteasesSpinbox.setToolTip(max_proteases_tooltip)
        self.RandomSamplingsLabel.setToolTip(random_sampling_tooltip)
        self.RandomSamplingsSpinbox.setToolTip(random_sampling_tooltip)
        self.DynamicRangeLabel.setToolTip(dynamic_range_tooltip)
        self.DynamicRangeSpinbox.setToolTip(dynamic_range_tooltip)
        # # # # position spin-boxes for further settings
        self.MaxProteasesSpinbox.setGeometry(rel_pos(20, 40, 80, 25))
        self.RandomSamplingsSpinbox.setGeometry(rel_pos(20, 90, 80, 25))
        self.DynamicRangeSpinbox.setGeometry(rel_pos(20, 140, 80, 25))

        # # set save and load parameter buttons
        self.SaveParamsButton = QtW.QPushButton(parent=self.ParamTab)
        self.LoadParamsButton = QtW.QPushButton(parent=self.ParamTab)
        # # # # setup save and load parameter buttons
        self.SaveParamsButton.setText("Save parameters")
        self.LoadParamsButton.setText("Load parameters")
        # # # # position save and load button
        self.SaveParamsButton.setGeometry(rel_pos(788, 725, 150, 35))
        self.LoadParamsButton.setGeometry(rel_pos(976, 725, 150, 35))

        # # set placeholder lable to show errors in param file
        self.ParamErrorLabel = QtW.QLabel(parent=self.ParamTab)
        self.ParamErrorLabel.setText("")
        self.ParamErrorLabel.setProperty("class", "err_lable")
        self.ParamErrorLabel.setWordWrap(True)
        self.ParamErrorLabel.setGeometry(rel_pos(788, 760, 300, 50))

        ##############################
        # setup export tab

        # # set box for protein weight file
        self.ProteinWeightBox = QtW.QGroupBox(parent=self.ExportTab)
        self.ProteinWeightBox.setGeometry(rel_pos(30, 95, 1140, 750))

        # + set table for protein weight file
        self.ProteinWeightTable = QtW.QTableWidget(parent=self.ProteinWeightBox)
        header = self.ProteinWeightTable.horizontalHeader()
        header.setSectionResizeMode(QtW.QHeaderView.ResizeMode.ResizeToContents)
        self.ProteinWeightTable.setGeometry(rel_pos(0, 00, 1120, 730))
        self.ProteinWeightTable.setProperty("class", "large_tbl")

        # # set protein weight file box
        self.ProteinWeightFileBox = QtW.QLabel(parent=self.ExportTab)
        self.ProteinWeightFileBox.setText("Protein weight file")
        self.ProteinWeightFileBox.setGeometry(rel_pos(1000, 830, 180, 30))
        self.ProteinWeightFileBox.setToolTip("Select or save a file containing protein<br>"
                                             "identifier, protein bins and abundance<br>"
                                             "values for each round of random sampling.<br>")

        # # # set buttons to re-sample proteins, save values, load annotation
        self.LoadAnnotButton = QtW.QPushButton(parent=self.ExportTab)
        self.ResampleButton = QtW.QPushButton(parent=self.ExportTab)
        self.OpenPWFileButton = QtW.QPushButton(parent=self.ExportTab)
        self.SavePWFileButton = QtW.QPushButton(parent=self.ExportTab)
        self.LoadAnnotButton.setEnabled(False)
        # # # # setup buttons to re-sample proteins, save values, load annotation
        self.LoadAnnotButton.setText("Load annotation")
        self.ResampleButton.setText("Sample proteins")
        self.OpenPWFileButton.setText("Open")
        self.SavePWFileButton.setText("Save")
        # # # # add tooltips for buttons to re-sample proteins, save values, load annotation
        self.LoadAnnotButton.setToolTip(
            "Replace current group assignment.<br><br>"
            "Load a tab-separated text file that contains<br>"
            "at least the column 'Identifier' and a new grouping column.")
        self.ResampleButton.setToolTip(
            "Repeat random protein abundance assignment.<br><br><em>Current values will be lost if not "
            "saved in advance.</em>")
        self.SavePWFileButton.setToolTip(
            "Open existing file with random protein abundance values.")
        self.SavePWFileButton.setToolTip(
            "Save current random protein abundance values.")
        # position buttons for buttons to re-sample proteins, save values, load annotation
        self.LoadAnnotButton.setGeometry(rel_pos(30, 30, 120, 35))
        self.ResampleButton.setGeometry(rel_pos(160, 30, 120, 35))
        self.SavePWFileButton.setGeometry(rel_pos(730, 30, 100, 35))
        self.OpenPWFileButton.setGeometry(rel_pos(600, 30, 100, 35))

        ##############################
        # setup progress tab

        # # set buttons for start, stop, pause
        self.ProgressStartPipelineButton = QtW.QPushButton(parent=self.ProgressTab)
        self.ProgressStartDigestButton = QtW.QPushButton(parent=self.ProgressTab)
        self.ProgressStartAnalysisButton = QtW.QPushButton(parent=self.ProgressTab)
        self.ProgressClearButton = QtW.QPushButton(parent=self.ProgressTab)
        self.ProgressSaveAsButton = QtW.QPushButton(parent=self.ProgressTab)
        # # # # setup save and load parameter buttons
        self.ProgressStartPipelineButton.setText("Start pipeline")
        self.ProgressStartDigestButton.setText("Start digest")
        self.ProgressStartAnalysisButton.setText("Start analysis")
        self.ProgressClearButton.setText("Clear")
        self.ProgressSaveAsButton.setText("Save as")
        # # # # set tooltips for all progress tab buttons
        self.ProgressStartPipelineButton.setToolTip("Start in-silico digestion and result analysis.")
        self.ProgressStartDigestButton.setToolTip("Start in-silico digestion and perform analysis later.")
        self.ProgressStartAnalysisButton.setToolTip("Start result analysis based on previous in-silico digestion.")
        self.ProgressSaveAsButton.setToolTip("Save log file.")
        self.ProgressClearButton.setToolTip("Clear log display.")
        # # # # position save and load button
        self.ProgressStartPipelineButton.setGeometry(rel_pos(30, 30, 100, 35))
        self.ProgressStartDigestButton.setGeometry(rel_pos(160, 30, 100, 35))
        self.ProgressStartAnalysisButton.setGeometry(rel_pos(290, 30, 100, 35))
        self.ProgressClearButton.setGeometry(rel_pos(600, 30, 100, 35))
        self.ProgressSaveAsButton.setGeometry(rel_pos(730, 30, 100, 35))

        # # add output window in cmd style
        self.ProgressReportFrame = QtW.QTextEdit(parent=self.ProgressTab)
        self.ProgressReportFrame.setReadOnly(True)
        self.ProgressReportFrame.setUndoRedoEnabled(False)  # reduce memory footprint a bit?
        self.ProgressReportFrame.setAcceptRichText(True)  # same purpose
        self.ProgressReportFrame.setGeometry(rel_pos(30, 95, 1140, 700))

        ##############################
        # setup Results tab

        # # set buttons for load and save
        self.ResultLoadFileButton = QtW.QPushButton(parent=self.ResultsTab)
        self.ResultSaveFigureButton = QtW.QPushButton(parent=self.ResultsTab)
        # # # # setup save and load buttons
        self.ResultLoadFileButton.setText("Load results")
        self.ResultSaveFigureButton.setText("Save figure")
        # # # # position save and load button
        self.ResultLoadFileButton.setGeometry(rel_pos(830, 800, 150, 35))
        self.ResultSaveFigureButton.setGeometry(rel_pos(1000, 800, 150, 35))

        # # set scroll area for plot-wdg
        self.ResultScrollArea = QtW.QScrollArea(parent=self.ResultsTab)
        # # # # setup scroll area for plot-wdg
        self.ResultScrollArea.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        self.ResultScrollArea.setWidgetResizable(False) # required!
        self.ResultScrollArea.setGeometry(rel_pos(30, 35, 890, 760))

        # # set buttons for load and save
        self.ResultIncreasePltSizeButton = QtW.QPushButton(parent=self.ResultsTab)
        self.ResultDecreasePltSizeButton = QtW.QPushButton(parent=self.ResultsTab)
        self.ResultResetPltSizeButton = QtW.QPushButton(parent=self.ResultsTab)
        # # # # setup save and load buttons
        self.ResultIncreasePltSizeButton.setText("+")
        self.ResultDecreasePltSizeButton.setText("-")
        self.ResultResetPltSizeButton.setText("=")
        # # # # set tooltips for size adjustment buttons
        self.ResultIncreasePltSizeButton.setToolTip("Increase plot size.")
        self.ResultDecreasePltSizeButton.setToolTip("Decrease plot size.")
        self.ResultResetPltSizeButton.setToolTip("Reset plot size.")
        self.ResultIncreasePltSizeButton.setProperty("class", "round_btn")
        self.ResultDecreasePltSizeButton.setProperty("class", "round_btn")
        self.ResultResetPltSizeButton.setProperty("class", "round_btn")
        # # # # position save and load button
        self.ResultIncreasePltSizeButton.setGeometry(rel_pos(30, 800, 35, 35))
        self.ResultDecreasePltSizeButton.setGeometry(rel_pos(70, 800, 35, 35))
        self.ResultResetPltSizeButton.setGeometry(rel_pos(110, 800, 35, 35))

        # # set radio buttons for filtered / unfiltered data
        self.ResultFilteredButton = QtW.QRadioButton(parent=self.ResultsTab)
        self.ResultUnfilteredButton = QtW.QRadioButton(parent=self.ResultsTab)
        # # # # setup radio buttons for filtered / unfiltered data
        self.ResultFilteredButton.setText("\u22652 Peptides")
        self.ResultUnfilteredButton.setText("Unfiltered")
        self.ResultFilteredButton.setChecked(True)
        # # # # set tooltips for size adjustment buttons
        self.ResultFilteredButton.setToolTip("Show results filtered for at least two unique peptides per protein.")
        self.ResultUnfilteredButton.setToolTip("Show unfiltered results")
        # # # # position save and load button
        self.ResultFilteredButton.setGeometry(rel_pos(950, 35, 150, 35))
        self.ResultUnfilteredButton.setGeometry(rel_pos(1050, 35, 150, 35))

        # # set empty plot
        self.new_plt = CoMPaseD_Boxplot(parent=self.ResultScrollArea, dpi = 96)
        self.new_plt.setGeometry(rel_pos(0, 0, 0, 0))
        self.ResultScrollArea.setWidget(self.new_plt)

        # Button to open information window with score formulae
        self.ResultFormulaInfoButton = QtW.QPushButton(parent=self.ResultsTab)
        self.ResultFormulaInfoButton.setText("Show Score Definition")
        self.ResultFormulaInfoButton.setToolTip("Show protease score definition.")
        self.ResultFormulaInfoButton.setGeometry(rel_pos(150, 800, 150, 35))

        '''
        score_file_location = path.dirname(path.realpath(__file__))
        score_img = path.join(score_file_location, "bin", "CoMPaseD_score.svg")
        if path.isfile(score_img):
            self.bg_img = QSvgWidget(score_img, parent=self.ResultsTab)
            self.bg_img.setGeometry(rel_pos(925, 600, 250, 200))
            self.bg_img.show()
        '''


        ##############################
        # setup objects and processes
        self.Protein_weight_file = ""
        self.param_file_name = ""  # init parameter file location as empty string within CoMPaseD_tabs object
        self.params_obj = CoMPaseD_Parameter()  # init parameter object instance within CoMPaseD_tabs object
        self.digest_proc = None  # init QProcess as None before testing 'is None' in function is important
        self.analysis_proc = None # init QProcess as None before testing 'is None'
        self.pipeline = None # init pipeline as None and set during pipeline or separate exec
        self.locked = False # init with non-locked gui, switch when running lock_params
        self.chkbx_lst = list() # init empty list of checkboxes for result filtering
        self.plt_width = 850 # plot width
        self.plt_height = 32 # plot height per bar
        self.plt_width_initial = self.plt_width # copy for reset
        self.plt_height_initial = self.plt_height # copy for reset
        self.finished_analysis = False
        self.digest_counter = 0 # counter during digest process to supress crux warning on existing output folder


        # # add tabs
        self.tabs.addTab(self.ConfigTab, "Configuration")
        self.tabs.addTab(self.ParamTab, "Parameter")
        self.tabs.addTab(self.ExportTab, "Export")
        self.tabs.addTab(self.ProgressTab, "Progress")
        self.tabs.addTab(self.ResultsTab, "Results")

        # add QTabWidget with tab containers to arbitrary QWidget and set this layout active
        self.layout.addWidget(self.tabs)
        self.setLayout(self.layout)

        # set functionalities for config tab
        config_initialise(self)
        self.LoadConfigButton.clicked.connect(lambda: config_load_defaults(self))
        self.SaveConfigButton.clicked.connect(lambda: config_save_defaults(self))
        self.CruxBrowseButton.clicked.connect(lambda: config_browse_crux(self))
        self.ClipsBrowseButton.clicked.connect(lambda: config_browse_clips(self))
        self.PromastBrowseButton.clicked.connect(lambda: config_browse_promast(self))

        # set functionalities for param tab
        self.FastaPathButton.clicked.connect(lambda: param_browse_fasta(self))
        self.OutputPathButton.clicked.connect(lambda: param_browse_out_folder(self))
        self.DMSPModelPathButton.clicked.connect(lambda: param_browse_dmsp_model(self))
        self.ProteaseResetButton.clicked.connect(lambda: load_protease_table_data(self))
        self.ProteaseClearButton.clicked.connect(lambda: clear_table(self))
        self.SubtractRowButton.clicked.connect(lambda: subtract_row(self))
        self.AddRowButton.clicked.connect(lambda: add_row(self))
        self.LoadParamsButton.clicked.connect(lambda: self.params_obj.param_load_file(tab_wdg=self))
        self.SaveParamsButton.clicked.connect(lambda: self.params_obj.param_save_file(tab_wdg=self))

        # set functionalities for export tab
        self.export_run_counter = 0
        self.tabs.currentChanged.connect(lambda: run_first_export(self, self.export_run_counter))
        self.ResampleButton.clicked.connect(lambda: run_further_export(self, self.export_run_counter))
        self.SavePWFileButton.clicked.connect(lambda: self.save_pwf_funct())
        self.OpenPWFileButton.clicked.connect(lambda: open_export_funct(self))
        self.LoadAnnotButton.clicked.connect(lambda: annotate_pwf_table(self))

        # set functionalities for progress tab
        self.ProgressStartPipelineButton.clicked.connect(self.start_pipeline)
        self.ProgressStartDigestButton.clicked.connect(self.start_digest)
        self.ProgressStartAnalysisButton.clicked.connect(self.start_analysis)
        self.ProgressSaveAsButton.clicked.connect(self.save_progress_as)
        self.ProgressClearButton.clicked.connect(self.clear_progress)

        # set functionalities for result tab
        self.ResultLoadFileButton.clicked.connect(lambda: self.init_new_plt())
        self.ResultFilteredButton.clicked.connect(lambda: self.set_filtered(self.new_plt))
        self.ResultUnfilteredButton.clicked.connect(lambda: self.set_filtered(self.new_plt))
        self.ResultIncreasePltSizeButton.clicked.connect(lambda: self.increase_plot_size())
        self.ResultDecreasePltSizeButton.clicked.connect(lambda: self.decrease_plot_size())
        self.ResultResetPltSizeButton.clicked.connect(lambda: self.reset_plot_size())
        self.ResultFormulaInfoButton.clicked.connect(lambda: self.show_score_formulae())
        self.ResultSaveFigureButton.clicked.connect(lambda: self.save_plt())

    def save_pwf_funct(self):
        """save protein abundance file"""

        # save as dialogue
        if path.isdir(path.join(self.params_obj.Output_directory)):
            file_name, _ = QtW.QFileDialog.getSaveFileName(self, "Save protein weight file",
                                                           directory=path.join(self.params_obj.Output_directory),
                                                           filter="ProteinAbundance (*.tsv);; ProteinAbundance * (*)")
        else:
            file_name, _ = QtW.QFileDialog.getSaveFileName(self, "Save protein weight file",
                                                       directory=getcwd(),
                                                       filter="ProteinAbundance (*.tsv);; ProteinAbundance * (*)")

        # do not change existing name when cancel was clicked
        if file_name != "":
            # if file exists, try to rename existing file with last modification date and time
            if path.isfile(file_name):
                mti = datetime.fromtimestamp(path.getmtime(file_name))
                rename_f_time = str(mti.strftime('%Y-%m-%d_%Hh%Mmin%Ssec_autosaved_'))
                rename_f_basename = str(path.basename(file_name))
                rename_f_name = rename_f_time  + rename_f_basename
                rename_f_name = path.join(path.dirname(file_name), rename_f_name)
                try:
                    rename(file_name, rename_f_name)
                except Exception as e:
                    print(f"{colorama.Fore.CYAN}WARNING: Could not rename existing protein abundance file {file_name} due to {e}. \n File will be overwritten.{colorama.Style.RESET_ALL}")

            # set protein weight file to chosen name
            self.Protein_weight_file = file_name

            # run export where all params are updated and
            save_export_funct(self, self.params_obj)


    def save_pwf_funct_silent(self):
        """save protein abundance file silently"""

        generic_file_name = path.join(self.params_obj.Output_directory, "ProteinAbundance.tsv")

        if path.isfile(generic_file_name):
            file_name_prefix = strftime("%Y-%m-%d_%Hh%Mmin%Ssec_") + "ProteinAbundance.tsv"
            generic_file_name = path.join(self.params_obj.Output_directory, file_name_prefix)

        # set protein weight file to chosen name
        self.Protein_weight_file = generic_file_name

        # run export where all params are updated and
        save_export_funct(self, self.params_obj)

    def save_plt(self):
        """save results boxplot"""
        if path.isdir(path.join(self.params_obj.Output_directory)):
            file_name, _ = QtW.QFileDialog.getSaveFileName(self, "Save plot as",
                                                           directory=path.join(self.params_obj.Output_directory),
                                                           filter="jpg (*.jpg);; tiff (*.tiff)")
        else:
            file_name, _ = QtW.QFileDialog.getSaveFileName(self, "Save plot as",
                                                           directory=getcwd(),
                                                           filter="jpg (*.jpg);; tiff (*.tiff)")
        if file_name != "":
            self.new_plt.print_figure(file_name, dpi=600)


    def init_new_plt(self):
        """load results and display as boxplot"""
        if path.isdir(path.join(self.params_obj.Output_directory)):
            file_name, _ = QtW.QFileDialog.getOpenFileName(self, "Result file",
                                                           directory=path.join(self.params_obj.Output_directory),
                                                           filter="tsv (*.tsv);; Result file (*)")
        else:
            file_name, _ = QtW.QFileDialog.getOpenFileName(self, "Result file",
                                                       directory=getcwd(),
                                                       filter="tsv (*.tsv);; Result file (*)")

        # load results temporary and count number of protein groups
        if path.isfile(file_name):
            tmp_res_df = read_csv(path.join(file_name), sep="\t")

            # check presence of required columns and remind user if columns are wrong and stop loading
            col_list = ['Protein group', 'Protease combination', "Protease score (filtered)", "Protease score (unfiltered)"]
            if not set(col_list).issubset(tmp_res_df.columns):

                # show error message when folders can not be created
                self.error_msg = QtW.QDialog(parent=self.ResultsTab)
                self.error_msg.setWindowTitle("Error - Selected file seems not to be a valid CoMPaseD results file.")
                # move to center
                self.error_msg.setGeometry(rel_pos(0, 0, 480, 165))
                center_point = QtGui.QGuiApplication.primaryScreen().availableGeometry().center()
                qt_rectangle = self.error_msg.frameGeometry()
                qt_rectangle.moveCenter(center_point)
                self.error_msg.move(qt_rectangle.topLeft())

                self.error_msg.setProperty("class", "wrong_file_warn")

                self.error_msg_txt = QtW.QLabel(parent = self.error_msg)
                self.error_msg_txt.setProperty("class", "wrong_file_warn")
                self.error_msg_txt.setText("Please select a valid file with at least the following columns: 'Protein group', 'Protease combination', 'Protease score (filtered)', 'Protease score (unfiltered)'.")
                self.error_msg_txt.setGeometry(rel_pos(20, 0, 440, 75))
                self.error_msg_txt.setWordWrap(True)

                self.error_msg_txt2 = QtW.QLabel(parent = self.error_msg)
                self.error_msg_txt2.setProperty("class", "wrong_file_warn2")
                self.error_msg_txt2.setText("Did you select CoMPaseD_results_<b>summary</b>.tsv instead of CoMPaseD_results.tsv?")
                self.error_msg_txt2.setGeometry(rel_pos(20, 75, 440, 35))
                self.error_msg_txt2.setWordWrap(True)

                self.error_msg_ok = QtW.QPushButton(parent = self.error_msg)
                self.error_msg_ok.setProperty("class", "wrong_file_warn_btn")
                self.error_msg_ok.setText("Close")

                self.error_msg_ok.setGeometry(rel_pos(350, 115, 100, 35))

                self.error_msg_ok.clicked.connect(lambda: self.error_msg.close())
                self.error_msg.exec()



                return


            n_plots = len(set(tmp_res_df['Protein group']))
            # init plot with n subplots equal to number of protein groups
            self.new_plt = CoMPaseD_Boxplot(parent=self.ResultScrollArea, n_groups=n_plots, dpi=96)
            self.new_plt.setGeometry(rel_pos(0, 0, 0, 0))  # hide while loading
            # assign plot to scroll area
            self.ResultScrollArea.setWidget(self.new_plt)

        # load result file again and fill plot
        if path.isfile(file_name):
            self.new_plt.read_results(file_name)
            # self.new_plt.combinations is list with n proteases and list is updated by (un)checking buttons in gui
            self.new_plt.get_bp_data(self.new_plt.combinations)
            # store combination checkboxes (changed, now QPushbutton subclass) to list for later access
            self.chkbx_lst = self.new_plt.make_combination_checkboxes(self)
            self.new_plt.set_bp_data(filtered=True)
            # calculate plot height and width by number of boxes
            n_plts = len(self.new_plt.grps)
            n_bars = self.new_plt.height_n_bars
            self.new_plt.setGeometry(rel_pos(0, 0, self.plt_width, ((n_plts * self.plt_height * n_bars)+(n_plts * 150))))

        # set functionality for protease number buttons
        for i, chkbx in enumerate(self.chkbx_lst):
            chkbx.clicked.connect(lambda _, i=i: self.adjust_combin_list(self.new_plt, i))

    def increase_plot_size(self):
        if self.plt_width < 5000:
            self.plt_width = self.plt_width * 1.1
            self.plt_height = self.plt_height * 1.1
            self.set_filtered(plt_wdg=self.new_plt)

    def decrease_plot_size(self):
        if self.plt_width > 10:
            self.plt_width = self.plt_width * 0.9
            self.plt_height = self.plt_height * 0.9
            self.set_filtered(plt_wdg=self.new_plt)

    def reset_plot_size(self):
        self.plt_width = self.plt_width_initial
        self.plt_height = self.plt_height_initial
        self.set_filtered(plt_wdg=self.new_plt)

    def show_score_formulae(self):
        # locate score formulae file
        score_file_location = path.dirname(path.realpath(__file__))
        score_img = path.join(score_file_location, "../bin", "CoMPaseD_score.svg")
        if path.isfile(score_img):
            # add dialog window
            self.ScoreWindow = QtW.QDialog(parent = self.ResultsTab)
            self.ScoreWindow.setWindowTitle('Help')
            # move to center
            self.ScoreWindow.setGeometry(rel_pos(0, 0, 900, 300))
            center_point = QtGui.QGuiApplication.primaryScreen().availableGeometry().center()
            qt_rectangle = self.ScoreWindow.frameGeometry()
            qt_rectangle.moveCenter(center_point)
            self.ScoreWindow.move(qt_rectangle.topLeft())

            # add formulae as 'background' image
            self.ScoreWindow_bg = QSvgWidget(score_img, parent=self.ScoreWindow)
            #self.ScoreWindow_bg.setGeometry(rel_pos(-60, 30, 1050, 90))
            self.ScoreWindow_bg.setGeometry(rel_pos(30, 30, 840, 90))
            self.ScoreWindow_bg.show()

            # add figure description
            self.ScoreWindow_txt = QtW.QLabel(parent = self.ScoreWindow)
            self.ScoreWindow_txt.setProperty("class", "info_lable")
            '''
            self.ScoreWindow_txt.setText("The <b>Protease Score</b> of a particular "
                                         "combination (<i>comb</i>) of proteases for a group of proteins "
                                         "is calculated as the "
                                         "weighted geometric mean of the number of identified "
                                         "proteins (<i>N<sub>prot</sub></i>), the number of identified "
                                         "peptides (<i>N<sub>pep</sub></i>) and the average sequence coverage "
                                         "of the identified proteins (<i>cov</i>). Each of these values is normalised "
                                         "to the same value in a tryptic digest (<i>try</i>). The weighting factors (<i>w</i>) "
                                         "can be adjusted for each of these values.")
            '''
            self.ScoreWindow_txt.setText("The <b>Protease Score</b> <i>S</i> for a group <i>G</i> of proteins and protease "
                                         "combination <i>A</i> is calculated as the weighted geometric mean of the "
                                         "number of identified proteins (<i>N<sub>Prot-G-A</sub></i>), the number of "
                                         "identified peptides (<i>N<sub>Pep-G-A</sub></i>) and the average sequence "
                                         "coverage of the identified proteins (<i>COV<sub>G-A</sub></i>). Each of these "
                                         "values is normalised to the same value in a tryptic digest (combination <i>C</i>)."
                                         "The weighting factors (<i>W<sub>Prot</sub></i>, <i>W<sub>Pep</sub></i>, "
                                         "<i>W<sub>COV</sub></i>) may be adjusted.")



            self.ScoreWindow_txt.setWordWrap(True)
            self.ScoreWindow_txt.setGeometry(rel_pos(20, 130, 860, 190))
            self.ScoreWindow.exec()

    def set_filtered(self, plt_wdg: CoMPaseD_Boxplot):
        """run when switching radiobutton selection"""

        # run only when data were loaded
        if not plt_wdg.res_df.empty:
            plt_wdg.get_bp_data(plt_wdg.combinations)
            plt_wdg.set_bp_data(filtered=self.ResultFilteredButton.isChecked())
            n_plts = len(plt_wdg.grps)
            n_bars = plt_wdg.height_n_bars
            plt_wdg.setGeometry(rel_pos(0, 0, self.plt_width, ((n_plts * self.plt_height * n_bars)+(n_plts * 150))))
            # set position of scroll bar to avoid showing uppermost plot upon change at any position
            tmp_pos = self.ResultScrollArea.verticalScrollBar().value()
            self.ResultScrollArea.verticalScrollBar().setValue(0)
            self.ResultScrollArea.verticalScrollBar().setValue(tmp_pos)

    def adjust_combin_list(self, plt_wdg: CoMPaseD_Boxplot, i: int):

        # disable unchecking last combination to avoid crashes
        if plt_wdg.combinations.count(None) == (len(plt_wdg.combinations) - 1):
            self.chkbx_lst[i].setChecked(True)

        # check whether checkbox is unchecked or checked and set corresponding value in
        # plt_wdg.combinations to None or original value
        if self.chkbx_lst[i].isChecked():
            plt_wdg.combinations[i] = self.chkbx_lst[i].combination_val
        if not self.chkbx_lst[i].isChecked():
            plt_wdg.combinations[i] = None

        # re-draw figure
        self.set_filtered(plt_wdg=self.new_plt)


    def clear_progress(self):
        """Clear text in progress cmd-style window"""

        self.ProgressReportFrame.clear()


    def save_progress(self, log_file_path):
        """Save progress report frame content as log_file_path"""

        log_text = self.ProgressReportFrame.toPlainText()

        with open(path.join(log_file_path), 'w') as f:
            for line in log_text:
                f.writelines(line)


    def save_progress_as(self):
        """Open save as dialog and save progress frame content"""

        if path.isdir(self.params_obj.Output_directory):
            file_name, _ = QtW.QFileDialog.getSaveFileName(self, "Save logfile as",
                                                           directory=self.params_obj.Output_directory,
                                                           filter="Text file (*.txt);;Log file (*.log)")
        else:
            file_name, _ = QtW.QFileDialog.getSaveFileName(self, "Save logfile as",
                                                           directory=getcwd(),
                                                           filter="Text file (*.txt);;Log file (*.log)")
        if file_name != "":
            self.save_progress(log_file_path=file_name)


    def lock_params(self):
        """lock editable fields to prevent parameter changes during analysis"""

        self.locked = True
        self.CruxPathField.setEnabled(False)
        self.PerlFilePathBox.setEnabled(False)
        self.ClipsPathField.setEnabled(False)
        self.PromastPathField.setEnabled(False)
        self.CruxBrowseButton.setEnabled(False)
        self.ClipsBrowseButton.setEnabled(False)
        self.PromastBrowseButton.setEnabled(False)
        self.SamplingOutputCheckbox.setEnabled(False)
        self.SaveConfigButton.setEnabled(False)
        self.LoadConfigButton.setEnabled(False)

        self.FastaPathField.setEnabled(False)
        self.OutputPathField.setEnabled(False)
        self.FastaPathButton.setEnabled(False)
        self.OutputPathButton.setEnabled(False)
        self.ProteaseTable.setEnabled(False)
        self.SamplingSizeBasisBox.setEnabled(False)
        self.AddRowButton.setEnabled(False)
        self.SubtractRowButton.setEnabled(False)
        self.ProteaseResetButton.setEnabled(False)
        self.ProteaseClearButton.setEnabled(False)
        self.ProtIDWeightField.setEnabled(False)
        self.PepIDWeightField.setEnabled(False)
        self.CoverageWeightField.setEnabled(False)
        self.ProteinBinsField.setEnabled(False)
        self.ProteinNotExprFracField.setEnabled(False)
        self.DMSPBox.setEnabled(False)
        self.DMSPModelPathField.setEnabled(False)
        self.DMSPWeightSpinbox.setEnabled(False)
        self.MaxProteasesSpinbox.setEnabled(False)
        self.RandomSamplingsSpinbox.setEnabled(False)
        self.DynamicRangeSpinbox.setEnabled(False)
        self.SaveParamsButton.setEnabled(False)
        self.LoadParamsButton.setEnabled(False)

        self.LoadAnnotButton.setEnabled(False)
        self.ResampleButton.setEnabled(False)
        self.OpenPWFileButton.setEnabled(False)
        self.SavePWFileButton.setEnabled(False)

        self.ProgressClearButton.setEnabled(False)
        # self.ProgressSaveAsButton.setEnabled(False) # save as should be possible at any time
        self.ProgressStartAnalysisButton.setEnabled(False)
        self.ProgressStartDigestButton.setEnabled(False)
        self.ProgressStartPipelineButton.setEnabled(False)

        # inform user on locking in status bar
        self.locked_msg.setText("Analysis running, parameter can not be changed.")


    def unlock_params(self):
        """unlock editable fields after analysis"""
        self.locked_msg.setText("")
        self.locked = False

        self.CruxPathField.setEnabled(True)
        self.PerlFilePathBox.setEnabled(True)
        self.ClipsPathField.setEnabled(True)
        self.PromastPathField.setEnabled(True)
        self.CruxBrowseButton.setEnabled(True)
        self.ClipsBrowseButton.setEnabled(True)
        self.PromastBrowseButton.setEnabled(True)
        self.SamplingOutputCheckbox.setEnabled(True)
        self.SaveConfigButton.setEnabled(True)
        self.LoadConfigButton.setEnabled(True)

        self.FastaPathField.setEnabled(True)
        self.OutputPathField.setEnabled(True)
        self.FastaPathButton.setEnabled(True)
        self.OutputPathButton.setEnabled(True)
        self.ProteaseTable.setEnabled(True)
        self.SamplingSizeBasisBox.setEnabled(True)
        self.AddRowButton.setEnabled(True)
        self.SubtractRowButton.setEnabled(True)
        self.ProteaseResetButton.setEnabled(True)
        self.ProteaseClearButton.setEnabled(True)
        self.ProtIDWeightField.setEnabled(True)
        self.PepIDWeightField.setEnabled(True)
        self.CoverageWeightField.setEnabled(True)
        self.ProteinBinsField.setEnabled(True)
        self.ProteinNotExprFracField.setEnabled(True)
        self.DMSPBox.setEnabled(True)
        self.DMSPModelPathField.setEnabled(True)
        self.DMSPWeightSpinbox.setEnabled(True)
        self.MaxProteasesSpinbox.setEnabled(True)
        self.RandomSamplingsSpinbox.setEnabled(True)
        self.DynamicRangeSpinbox.setEnabled(True)
        self.SaveParamsButton.setEnabled(True)
        self.LoadParamsButton.setEnabled(True)

        self.LoadAnnotButton.setEnabled(True)
        self.ResampleButton.setEnabled(True)
        self.OpenPWFileButton.setEnabled(True)
        self.SavePWFileButton.setEnabled(True)

        self.ProgressClearButton.setEnabled(True)
        # self.ProgressSaveAsButton.setEnabled(True) # save as should be possible at any time
        self.ProgressStartAnalysisButton.setEnabled(True)
        self.ProgressStartDigestButton.setEnabled(True)
        self.ProgressStartPipelineButton.setEnabled(True)

    def report_frame_print_progress(self, curr_print):
        self.ProgressReportFrame.append(curr_print)

    def progress_digest_finished(self):
        self.digest_proc = None
        self.unlock_params()
        self.digest_counter = 0

        # set result to param_obj
        self.params_obj.Digestion_result_file = path.join(self.params_obj.Output_directory, "unique_peptides_table_filtered.tsv")
        # if not run in pipeline mode, save logfile
        if self.pipeline == False:
            file_name = strftime("%Y-%m-%d_%Hh%Mmin%Ssec_CoMPaseD_digestion_log.txt")
            log_file_path = path.join(self.params_obj.Output_directory, file_name)
            self.save_progress(log_file_path=log_file_path)
            # reset pipeline if not in pipeline mode
            self.pipeline = None

        # if run in pipeline mode, start analysis
        elif self.pipeline == True:
            # ensure that digestion and peptide pooling are finished
            # file size should change within 3 secs
            def test_file_writing_finished(filename):
                size_0 = stat(filename).st_size
                sleep(3)
                size_1 = stat(filename).st_size
                # if file size does not change / is equal, return true as writing appears to be finished
                if size_0 == size_1:
                    return True
                # if there was a change, return false and repeat testing in while below
                else:
                    return False
            # wait until filtered table has been created
            if not path.isfile(path.join(self.params_obj.Output_directory, 'unique_peptides_table_filtered.tsv')):
                #
                while not test_file_writing_finished(path.join(self.params_obj.Output_directory, 'unique_peptides_table_unfiltered.tsv')):
                    sleep(1)
                # wait another three sec to ensure creation of filtered file
                sleep(3)
            else:
                # if filtered file already exists, directly wait until it does not change anymore
                while not test_file_writing_finished(path.join(self.params_obj.Output_directory, 'unique_peptides_table_filtered.tsv')):
                    sleep(1)
                # start analysis only when filtered file was finished
                self.start_analysis(pipeline=True)

            # test again if filtered file now exists to avoid errors when this file is not created by some reason
            if path.isfile(path.join(self.params_obj.Output_directory, 'unique_peptides_table_filtered.tsv')):
                while not test_file_writing_finished(path.join(self.params_obj.Output_directory, 'unique_peptides_table_filtered.tsv')):
                    sleep(1)
                # start analysis only when filtered file was finished
                self.start_analysis(pipeline=True)
            else:
                self.ProgressReportFrame.append(f"\t Error: File 'unique_peptides_table_filtered' does not exist in output directory. Please check")

    # formatting for print statements
    def format_byte_str(self, byte_str):
        '''
        use html formatting for output to report window
        '''
        # split into lines
        lines = byte_str.split("\r\n")
        formatted_lines = []
        previous_line_blank = False
        current_indent = 0

        for line in lines:
            # remove leading/trailing spaces for checks
            stripped_line = line.strip()

            # find empty lines
            # if line does not exist, i.e. if it is empty
            if not stripped_line:
                # add a new line only, if the previous line was not blank
                if not previous_line_blank:
                    formatted_lines.append("<br>")
                    # set checkpoint
                    previous_line_blank = True
                continue  # skip further processing for blank lines

            # reset if the last line was not an empty one
            previous_line_blank = False

            # count the indentation level
            indent = len(line) - len(line.lstrip())
            current_indent = max(current_indent, indent)

            # wrap all lines in <span...> to avoid mixed html and normal character parsing,
            # this could otherwise result in &nbsp to be printed literally
            html_line = line.replace('\t', '&nbsp;' * (current_indent + 4))

            # remove ANSI formatting information in case this is present; remove this line if the output gets affected
            html_line = re.sub(r'\x1b\[([0-?]*[ -/]*[@-~])', '', html_line)

            # count digestions
            if stripped_line.startswith("Started digestion with"):
                self.digest_counter += 1

            # list of line starts that indicate info lines
            info_keywords = ['INFO:',
                             'Protease:',
                             'Missed cleavages:',
                             'Indexed fasta file:',
                             'Index creation date:',
                             'Protein entries:',
                             'Index len:',
                             'Indexed amino acids:',
                             'Existing files will be overwritten.']

            info_keywords_indent = ['INFO:','Existing files will be overwritten.' ]

            if not any(char in line for char in "<>"):
                html_line = f"<span style='white-space: pre;'>{html_line}</span>"

            # format header line larger and bold
            if stripped_line.startswith("CoMPaseD - Comparison of Multiple-Protease Digestions"):
                formatted_lines.append(f"<b><span style='font-size: 20px;'>{stripped_line}</span></b>")

            # format module headers (like "In-silico digestion started")
            # this captures both; 'In-silico digestion started' and 'Analysis of In-silico digestion started'
            elif "digestion started" in stripped_line.lower():
                formatted_lines.append(f"<b><u>{stripped_line}</u></b>")

            # format warnings in orange
            elif stripped_line.startswith("WARNING:"):
                # check if this warning is related to the existing crux output order
                if stripped_line.startswith("WARNING: The output directory 'crux-output' already exists."):
                    if self.digest_counter > 1:
                        pass
                    else:
                        formatted_lines.append(f"<span style='color: #FFA500'>{html_line}</span>")
                else:
                    formatted_lines.append(f"<span style='color: #FFA500'>{html_line}</span>")

            # format errors in red and italic
            elif stripped_line.startswith("ERROR:"):
                formatted_lines.append(
                    f"<span style='color: #FF1059'><i>{html_line}</i></span>"
                )

            # keep error messages from raise errorXY separate and remove the automatically added formatting information
            elif 'ERROR:' in stripped_line:
                # remove the ANSI escape sequences from error messages
                html_line = re.sub(r'\x1b\[([0-?]*[ -/]*[@-~])', '', html_line)
                formatted_lines.append(
                    f"<span style='color: #FF1059'><i>{html_line}</i></span>"
                )

            # format information lines from crux to default blue
            elif any(stripped_line.startswith(keyword) for keyword in info_keywords):
                if any(stripped_line.startswith(kw) for kw in info_keywords_indent):
                    formatted_lines.append(
                        f"<span style='color: #031059'>&nbsp;&nbsp;&nbsp;&nbsp;{html_line}</span>"
                    )
                else:

                    formatted_lines.append(
                        f"<span style='color: #031059'>{html_line}</span>"
                    )

            elif stripped_line.startswith('Write results to'):
                split_html_line = re.sub('Write results to', '', html_line)
                formatted_lines.append(
                    f"Write results to: <b><u> {path.normpath(split_html_line)}</b></u><br>"
                )

            elif stripped_line.startswith('Write results summary to'):
                split_html_line = re.sub('Write results summary to', '', html_line)
                formatted_lines.append(
                    f"Write results summary to: <b><u> {path.normpath(split_html_line)}</b></u>"
                )

            # Default case for unformatted lines
            else:
                formatted_lines.append(html_line)

        # Join formatted lines without trailing `<br>` to prevent extra gaps
        formatted_lines = [re.sub(r'<br>+', '', f) for f in formatted_lines]
        formatted_lines = [f for f in formatted_lines if f]

        result = "<br>".join(formatted_lines)
        result = re.sub(r'<br>+', '<br>', result)
        return result

    def handle_digest_stdout(self):

        data = self.digest_proc.readAllStandardOutput()
        err = self.digest_proc.readAllStandardError()
        stdout = bytes(data).decode("iso-8859-1")
        stderr = bytes(err).decode("iso-8859-1")
        # process only if stdout is non-empty
        if stdout.strip():
            # format non-empty outputs
            stdout = self.format_byte_str(stdout)
            # print(f"stdout: {repr(stdout)}")
            self.report_frame_print_progress(stdout)

        if stderr.strip():
            stderr = self.format_byte_str(stderr)
            # print(f"stderr: {repr(stderr)}")
            self.report_frame_print_progress(stderr)
        # print("###################")



    def start_digest(self, pipeline=False):

        # check validity of parameters
        self.params_obj.get_params_from_gui(tab_wdg=self)
        valid_state, valid_err_count, valid_err = self.params_obj.validate_params(mode=2)

        err_number_backup = 1
        if not valid_state:
            self.ProgressReportFrame.append("Parameter error(s) can not start analysis:")
            for err, err_number in zip(valid_err, range(valid_err_count+1)):
                self.ProgressReportFrame.append(f"\t Error #{err_number}: {err}")
                err_number_backup = err_number
            return

        # adjust pipeline flag
        self.pipeline = pipeline

        if self.digest_proc is None:
            self.digest_proc = QtCore.QProcess()
            self.digest_proc.readyReadStandardOutput.connect(self.handle_digest_stdout)
            self.digest_proc.finished.connect(self.progress_digest_finished)
            # save parameter file under the output_directory
            param_file = strftime("%Y-%m-%d_%Hh%Mmin%Ssec_CoMPaseD_parameter.params")
            param_file_path = path.join(self.params_obj.Output_directory, param_file)
            try:
                # for safety, due to the time string it should be impossible to have a param file with identical name
                remove(param_file_path)
            except OSError:
                pass

            try:
                # ensure that output directory exists
                if not path.isdir(self.params_obj.Output_directory):
                    makedirs(self.params_obj.Output_directory)
                # auto-save parameters
                self.params_obj.save_params_to_file(param_file_path=param_file_path,tab_wdg=self)
                save_info_msg = "Auto-saved parameter file to " + param_file
                self.status.showMessage(save_info_msg, 3500)
            except Exception as e:
                # should only be printed to console but GUI should remain unaffected
                print(f"ERROR: {e}")
                pass

            # lock gui to prevent changes during execution
            self.lock_params()

            # parse params to cmd line args
            # try to find current python executable by sys.executable, fallback to PATH values python3 or python else
            python_exec = ""
            if executable is not None:
                python_exec = executable
            elif shutil.which("python3") is not None:
                python_exec = shutil.which("python3")
            elif shutil.which("python") is not None:
                python_exec = shutil.which("python")

            if self.PerlFilePathBox.isChecked():

                out_dir = self.params_obj.Output_directory
                fasta = self.params_obj.Fasta
                crux_path = self.params_obj.Crux_path
                clips_path = self.params_obj.Clips_path
                promast_path = self.params_obj.Promast_path
                enzyme_string = self.params_obj.Proteases
                max_mc_str = self.params_obj.Max_MCs

                file_location = path.dirname(path.realpath(__file__))
                CoMPaseD_crux_script = path.join(file_location, 'CoMPaseD_crux_script.py')

                args_list = [CoMPaseD_crux_script,
                             "--out_folder",
                             out_dir,
                             "--fasta",
                             fasta,
                             "--crux_path",
                             crux_path,
                             "--clips_path",
                             clips_path,
                             "--promast_path",
                             promast_path,
                             "--enzyme_list",
                             enzyme_string,
                             "--max_mc_list",
                             max_mc_str]
            else:

                out_dir = self.params_obj.Output_directory
                fasta = self.params_obj.Fasta
                crux_path = self.params_obj.Crux_path
                enzyme_string = self.params_obj.Proteases
                max_mc_str = self.params_obj.Max_MCs
                digestion_result_file = self.params_obj.Digestion_result_file
                indexing_key_len = self.params_obj.Indexing_key_len
                differentiate_I_L = self.params_obj.Differentiate_I_L


                file_location = path.dirname(path.realpath(__file__))
                CoMPaseD_PeptideMapper = path.join(file_location, 'CoMPaseD_PeptideMapper.py')

                args_list = [CoMPaseD_PeptideMapper,
                             "--out_folder",
                             out_dir,
                             "--fasta",
                             fasta,
                             "--crux_path",
                             crux_path,
                             "--enzyme_list",
                             enzyme_string,
                             "--max_mc_list",
                             max_mc_str,
                             "--unique_peps_file",
                             digestion_result_file,
                             "--indexing_key_len",
                             str(indexing_key_len)
                             ]
                if not differentiate_I_L:
                    args_list.append("--differentiate_I_L")

            if not python_exec == "":
                self.digest_proc.start(python_exec, args_list)


    def progress_analysis_finished(self):
        self.analysis_proc = None
        self.unlock_params()
        self.finished_analysis = True

        # if not run in pipeline mode, save logfile
        if self.pipeline == False:
            file_name = strftime("%Y-%m-%d_%Hh%Mmin%Ssec_CoMPaseD_analysis_log.txt")
            log_file_path = path.join(self.params_obj.Output_directory, file_name)
            self.save_progress(log_file_path=log_file_path)
            # reset pipeline
            self.pipeline = None
        # if running in pipeline mode export only after analysis and use different log file name
        elif self.pipeline == True:
            file_name = strftime("%Y-%m-%d_%Hh%Mmin%Ssec_CoMPaseD_pipeline_log.txt")
            log_file_path = path.join(self.params_obj.Output_directory, file_name)
            self.save_progress(log_file_path=log_file_path)
            # reset pipeline
            self.pipeline = None


    def handle_analysis_stdout(self):
        data = self.analysis_proc.readAllStandardOutput()
        err = self.analysis_proc.readAllStandardError()
        stdout = bytes(data).decode("iso-8859-1")
        stderr = bytes(err).decode("iso-8859-1")
        #print(f"stdout-pre: {repr(stdout)}")
        #print(f"stderr-pre: {repr(stderr)}")
        # process only if stdout is non-empty
        if stdout.strip():
            # format non-empty outputs
            stdout = self.format_byte_str(stdout)
            #print(f"stdout: {repr(stdout)}")
            self.report_frame_print_progress(stdout)

        if stderr.strip():
            stderr = self.format_byte_str(stderr)
            #print(f"stderr: {repr(stderr)}")
            self.report_frame_print_progress(stderr)
        #print("###################")


    def start_analysis(self, pipeline=False):

        # adjust pipeline flag
        self.pipeline = pipeline

        # check validity of parameters
        self.params_obj.get_params_from_gui(tab_wdg=self)
        valid_state, valid_err_count, valid_err = self.params_obj.validate_params(mode=2)

        err_number_backup = 1
        if not valid_state:
            self.ProgressReportFrame.append("Parameter error(s) can not start analysis:")
            for err, err_number in zip(valid_err, range(valid_err_count+1)):
                self.ProgressReportFrame.append(f"\t Error #{err_number}: {err}")
                err_number_backup = err_number
            if not path.isfile(path.join(self.params_obj.Digestion_result_file)):
                self.ProgressReportFrame.append(f"\t Error #{err_number_backup}: Digestion result file does not exist ({self.params_obj.Digestion_result_file}).")
            return
        if not path.isfile(path.join(self.params_obj.Digestion_result_file)):
            self.ProgressReportFrame.append(
                f"\t Error #{err_number_backup}: Digestion result file does not exist ({self.params_obj.Digestion_result_file}).")
            return


        if not path.isfile(self.Protein_weight_file):
            if self.export_run_counter == 0:
                run_first_export(self, self.export_run_counter)
            self.save_pwf_funct_silent()

        if self.analysis_proc is None:
            self.analysis_proc = QtCore.QProcess()
            self.analysis_proc.readyReadStandardOutput.connect(self.handle_analysis_stdout)
            self.analysis_proc.readyReadStandardError.connect(self.handle_analysis_stdout)
            self.analysis_proc.finished.connect(self.progress_analysis_finished)

            # do not repeatedly auto-save params in pipeline mode,
            # usage of pipeline instead of self.pipeline as the latter could change from finished signal
            if not pipeline == True:
                # save parameter file under the output_directory
                param_file = strftime("%Y-%m-%d_%Hh%Mmin%Ssec_CoMPaseD_parameter.params")
                param_file_path = path.join(self.params_obj.Output_directory, param_file)
                try:
                    # for safety, try to delete file,
                    # due to the time string it should be impossible to have a param file with identical name
                    remove(param_file_path)
                except OSError:
                    pass

                try:
                    # ensure that output directory exists
                    if not path.isdir(self.params_obj.Output_directory):
                        makedirs(self.params_obj.Output_directory)
                    # auto-save parameters
                    self.params_obj.save_params_to_file(param_file_path=param_file_path,tab_wdg=self)
                    save_info_msg = "Auto-saved parameter file to " + param_file
                    self.status.showMessage(save_info_msg, 3500)
                except Exception as e:
                    # should only be printed to console but GUI should remain unaffected
                    print(f"ERROR: {e}")
                    pass
            # lock gui to prevent changes during execution
            self.lock_params()

            # parse params to cmd line args

            file_location = path.dirname(path.realpath(__file__))
            CoMPaseD_analysis_script = path.join(file_location, 'CoMPaseD_analysis_script.py')

            args_list = [CoMPaseD_analysis_script,
                         "--param_file",
                         self.param_file_name]

            # try to find current python executable by sys.executable, fallback to PATH values python3 or python else
            python_exec = ""
            if executable is not None:
                python_exec = executable
            elif shutil.which("python3") is not None:
                python_exec = shutil.which("python3")
            elif shutil.which("python") is not None:
                python_exec = shutil.which("python")

            if not python_exec == "":
                self.analysis_proc.start(python_exec, args_list)

        else:
            pass


    def start_pipeline(self):

        # start only digest from here and use pipeline=True to modify
        # QProcess.finished signal behavior, analysis is started from
        # progress_digest_finished()
        self.start_digest(pipeline=True)

