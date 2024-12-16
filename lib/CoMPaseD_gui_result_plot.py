import matplotlib as mpl

mpl.use('QtAgg')

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
# from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import PyQt6.QtWidgets as QtW

from os import path
import numpy
from pandas import read_csv, DataFrame
from lib.CoMPaseD_tools import *


def sort_extract(tmp_res_df, score_col="Protease score (filtered)"):
    """sort results by mean score and return score arrays, combination list and protease count list"""

    # .transform adds column directly
    tmp_res_df["mean_" + score_col] = tmp_res_df.groupby(["Protease combination"])[score_col].transform('mean')
    # sort by score mean, then by Number of proteases and combination to break ties and finally by Random sampling
    mean_score_col = "mean_" + score_col
    tmp_res_df = tmp_res_df.sort_values([mean_score_col,
                                         "Number of proteases",
                                         "Protease combination",
                                         "Random sampling"], ascending=True)
    tmp_res_df = tmp_res_df.reset_index(drop=True)

    # do NOT use set() as the order is not preserved by set()
    combinations = tmp_res_df["Protease combination"].drop_duplicates().to_list()

    # lists to hold vals for each combination
    score = list()
    protease_count = list()  # in protease_combinations

    for comb in combinations:
        tmp_df = tmp_res_df.loc[tmp_res_df["Protease combination"] == comb, ]

        tmp_score = numpy.asarray(tmp_df[score_col].to_list())
        tmp_score = tmp_score[~numpy.isnan(tmp_score)]
        tmp_protease_count = tmp_df["Number of proteases"].to_list()

        score.append(tmp_score)
        # append only first element of protease count to obtain list with len == len(combinations)
        protease_count.append(tmp_protease_count[0])

    return combinations, score, protease_count


def set_bp_layout(ax, grp=""):
    """ define CoMPaseD boxplot layout"""

    if not grp == "":
        ax.set_xlabel(f'Protease score for {grp}')

    ax.tick_params(axis='both', which='major', labelsize=8, labelcolor="#031059")
    plt.subplots_adjust(left=0.4, top=0.97, bottom=0.05, right=0.97)
    ax.xaxis.set_tick_params(which='both', labelbottom=True)

    for spine in ax.spines.values():
        spine.set_edgecolor('#031059')

    return ax


class QCheckBox_Data(QtW.QPushButton):
    """QCheckbox class with associated data container for number of proteases"""

    def __init__(self, data, bg_col,*args, **kwargs):
        super().__init__(*args, **kwargs)
        self.combination_val = int(data)
        self.setCheckable(True)

        self.setStyleSheet(
            "QPushButton{"
            f"background-color: {bg_col};"
            "border-color: #0D0D0D;"
            "border-radius: 15px;"
            "color: 0D0D0D;"
            "font: bold 12px;"
            "}"
            "QPushButton:!checked{"
            "background-color: #A6A6A6;"
            "border-color: #0D0D0D;"
            "border-radius: 15px;"
            "color: #0D0D0D;"
            "font: bold 12px;"
            "}")


class CoMPaseD_Boxplot(FigureCanvasQTAgg):
    """boxplot"""

    def __init__(self, parent=None, n_groups=3, width=300, height=300, dpi=96):
        """plots are initialised only after loading data"""

        self.fig, self.axes = plt.subplots(n_groups, 1, figsize=rel_pos_fig([width, height]),sharex=False)
        # convert axes to 1d array in case of single-group experiments
        self.axes = numpy.array(self.axes)
        self.axes = numpy.reshape(self.axes, n_groups)

        self.fig.set_dpi(dpi)

        super(CoMPaseD_Boxplot, self).__init__(self.fig)
        self.setParent(parent)

        # set subplot for each group
        for n in range(n_groups):
            self.axes[n].boxplot([[], []], vert=False, showmeans=True, patch_artist=True)
            self.axes[n] = set_bp_layout(self.axes[n])
            self.axes[n].axvline(x=1.0, color='#A6A6A6')

        # init complete results
        self.res_df = DataFrame()
        self.combinations = list()

        # init data as lists of lists
        self.grps = list()
        self.protease_combinations_filtered = list()
        self.score_filtered = list()
        self.protease_count_filtered = list()  # in protease_combinations
        # as sorting might be different, repeat for all three fields
        self.protease_combinations_unfiltered = list()
        self.score_unfiltered = list()
        self.protease_count_unfiltered = list()  # in protease_combinations

        # init size properties
        self.rel_height = 1
        self.height_n_bars = 1
        self.rel_width = width


    def read_results(self, res_f_name):
        """read CoMPaseD result file
        default name is 'CoMPaseD_results.tsv' but might be prefixed with a time stamp or renamed
        """
        res_f_name = path.join(res_f_name)
        self.res_df = read_csv(path.join(res_f_name), sep="\t")
        self.res_df['Number of proteases'] = self.res_df['Protease combination'].str.split(" - ").str.len()
        self.combinations = sorted(set(self.res_df['Number of proteases']))


    def get_bp_data(self, n_proteases: list):
        """filter res_df, sort by score and fill data fields"""

        # clear old data
        self.grps = list()
        self.protease_combinations_filtered = list()
        self.score_filtered = list()
        self.protease_count_filtered = list()
        self.protease_combinations_unfiltered = list()
        self.score_unfiltered = list()
        self.protease_count_unfiltered = list()

        # filter for number of proteases and protein group
        subset_res_df = self.res_df.loc[self.res_df["Number of proteases"].isin(n_proteases), ]
        subset_res_df = subset_res_df.reset_index(drop=True)

        # loop through all protein_groups alphabetically sorted
        grps = sorted(set(subset_res_df["Protein group"]))

        # make sure that small proteins are the first plot in case default names are used
        default_2_grps = ["small_proteins", "large_proteins"]
        default_3_grps = ["small_proteins", "medium_proteins", "large_proteins"]
        if grps == sorted(default_2_grps):
            grps = default_2_grps
        elif grps == sorted(default_3_grps):
            grps = default_3_grps

        for i, grp in enumerate(grps):
            tmp_subset_res_df = subset_res_df.loc[subset_res_df["Protein group"] == grp, ]
            # reset index before putting df to list
            tmp_subset_res_df = tmp_subset_res_df.reset_index(drop=True)

            combinations_filtered, score_filtered, protease_count_filtered = sort_extract(tmp_subset_res_df,
                                                                                          "Protease score (filtered)")
            combinations_unfiltered, score_unfiltered, protease_count_unfiltered = sort_extract(tmp_subset_res_df,
                                                                                                "Protease score ("
                                                                                                "unfiltered)")

            # set new data
            self.grps.insert(i, str(grp))
            self.protease_combinations_filtered.insert(i, combinations_filtered)
            self.score_filtered.insert(i, score_filtered)
            self.protease_count_filtered.insert(i, protease_count_filtered)
            self.protease_combinations_unfiltered.insert(i, combinations_unfiltered)
            self.score_unfiltered.insert(i, score_unfiltered)
            self.protease_count_unfiltered.insert(i, protease_count_unfiltered)


    def make_combination_checkboxes(self, tab_wdg):
        """ generate n checkboxes for number of protease filtering"""

        # get list of combinations from df, sort set increasing
        n_combin = sorted(set(self.res_df["Number of proteases"]))

        # list to save all checkbox widgets
        combin_checkbox_lst = list()

        # list with button colors
        color_lst = ['#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB']


        # loop through combinations, create checkbox, set name and properties
        for i, combin in enumerate(n_combin):

            # change active tab to enable widget adding
            tab_wdg.tabs.setCurrentIndex(0)
            new_chkbx = QCheckBox_Data(parent=tab_wdg.ResultsTab, data=int(combin), bg_col=color_lst[(i+1) % len(color_lst)])


            if str(combin) == "1":
                new_chkbx.setText("1 Protease")
            else:
                new_chkbx.setText(f"{combin} Proteases")
            new_chkbx.setChecked(True)
            new_chkbx.setToolTip(f"Uncheck to hide combinations with {combin} protease(s).")
            # calculate position and place on gui
            new_chkbx.setGeometry(rel_pos(950, (60 + (40*(i+1))), 200, 35))
            # append checkbox to list for later access
            combin_checkbox_lst.append(new_chkbx)

            tab_wdg.tabs.setCurrentIndex(4)
        # list must be returned to calling tab_wdg
        return combin_checkbox_lst


    def set_bp_data(self, filtered=True):
        """plot actual data"""

        # check if data is avail, this will only be true after running get_bp_data
        if len(self.grps) > 0:

            # get filtered or unfiltered data
            if filtered:
                scores = self.score_filtered
                lables = self.protease_combinations_filtered
                lable_grp = self.protease_count_filtered
            else:
                scores = self.score_unfiltered
                lables = self.protease_combinations_unfiltered
                lable_grp = self.protease_count_unfiltered

            # init dict with boxplots and list with colors
            plt_dict_lst = list()
            color_lst = ['#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB']

            # set number of bars
            bar_idx_lst = list()

            for n, grp in enumerate(self.grps):
                # clear init boxplot or boxplot created by earlier set_bp_data
                self.axes[n].clear()

                # plot with current data and append plot to plt_dict_list for later access
                plt_dict_lst.append(self.axes[n].boxplot(scores[n], vert=False, labels=lables[n] , showmeans=True,
                                                         meanprops=dict(marker='.', markeredgecolor='#031059',
                                                                        markerfacecolor='#031059'),
                                                         medianprops=dict(color="#36454F", linestyle=':', linewidth=1.5),
                                                         boxprops=dict(linestyle=':', linewidth=1, color="black"),
                                                         whiskerprops =dict(linestyle='-', linewidth=1, color="#36454F"),
                                                         capprops=dict(linestyle='-', linewidth=1, color="#36454F"),
                                                         showfliers=False,
                                                         patch_artist=True))

                # set layout
                self.axes[n] = set_bp_layout(self.axes[n], grp)

                # add vertical line at score = 1.0
                self.axes[n].axvline(x=1.0, color='#A6A6A6')

                # update subplot
                self.axes[n].figure.canvas.draw()

                # extract lable_grp sublist
                curr_lable = lable_grp[n]

                # set box color by number of proteases
                for i, bx in enumerate(plt_dict_lst[n]['boxes']):
                    # set color as color_lst index of the i'th curr_lable item; use modulo len(color_lst) to 'recycle' color_lst)
                    bx.set_facecolor(color_lst[curr_lable[i] % len(color_lst)])
                    bar_idx_lst.append(i)


            # set n_bars to maximum number of bars per group
            self.height_n_bars = max(bar_idx_lst)
            if self.height_n_bars < 1:
                self.height_n_bars = 1
