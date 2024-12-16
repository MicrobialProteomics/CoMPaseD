from os import path, getcwd
from sys import platform
import PyQt6.QtWidgets as QtW


def config_read_values(tab_wdg):
    """read config values from gui"""
    crux_path = tab_wdg.CruxPathField.text()
    if tab_wdg.PerlFilePathBox.isChecked():
        use_perl = "True"
        clips_path = tab_wdg.ClipsPathField.text()
        promast_path = tab_wdg.PromastPathField.text()
    else:
        use_perl = "False"
        if path.exists(tab_wdg.ClipsPathField.text()):
            clips_path = tab_wdg.ClipsPathField.text()
        else:
            clips_path = " "

        if path.exists(tab_wdg.PromastPathField.text()):
            promast_path = tab_wdg.PromastPathField.text()
        else:
            promast_path = " "

    settings_multi = 'False' # keep always off
    settings_output = str(tab_wdg.SamplingOutputCheckbox.isChecked())
    return crux_path, use_perl, clips_path, promast_path, settings_multi, settings_output


def error_corrupted_config(value_list):
    """check config values submitted as value_list"""
    if path.exists(value_list[0]) and (value_list[1] == "True"):
        if path.exists(value_list[2]) and path.exists(value_list[3]) and (value_list[5] == "True" or value_list[5] == "False"):
            ErrorLabelCorruptedFile = ""
        else:
            ErrorLabelCorruptedFile = "Error. Configuration file is corrupted. Please save again."
    elif path.exists(value_list[0]) and (value_list[1] == "False"):
        if value_list[5] == "True" or value_list[5] == "False":
            ErrorLabelCorruptedFile = ""
        else:
            ErrorLabelCorruptedFile = "Error. Configuration file is corrupted. Please save again."
    else:
        ErrorLabelCorruptedFile = "Error. Configuration file is corrupted. Please save again."
    return ErrorLabelCorruptedFile


def config_write_values(tab_wdg, value_list):
    """write config values from file to gui"""
    current_error_state = error_corrupted_config(value_list=value_list)
    tab_wdg.ErrorLabelCorruptedFile.setText(current_error_state)

    if current_error_state == "":
        tab_wdg.CruxPathField.setText(value_list[0])
        if value_list[1] == "True":
            tab_wdg.PerlFilePathBox.setChecked(True)
            tab_wdg.ClipsPathField.setText(value_list[2])
            tab_wdg.PromastPathField.setText(value_list[3])
        else:
            tab_wdg.PerlFilePathBox.setChecked(False)
            if path.exists(value_list[2]):
                tab_wdg.ClipsPathField.setText(value_list[2])
            if path.exists(value_list[3]):
                tab_wdg.PromastPathField.setText(value_list[3])

        if value_list[5] == "True":
            tab_wdg.SamplingOutputCheckbox.setChecked(True)
        else:
            tab_wdg.SamplingOutputCheckbox.setChecked(False)
    return None


def config_save_defaults(tab_wdg):
    """save current config values from gui to file"""
    # get values for config and put into list
    crux_path, use_perl, clips_path, promast_path, settings_multi, settings_output = config_read_values(tab_wdg)

    config_list = [crux_path, use_perl, clips_path, promast_path, settings_multi, settings_output]
    # check if config_list is a valid configuration and adjust error message accordingly
    current_error_state = error_corrupted_config(config_list)
    tab_wdg.ErrorLabelCorruptedFile.setText(current_error_state)
    # try open an existing config file in the current wd and overwrite if possible or create
    try:
        file_location = path.dirname(path.realpath(__file__))
        config_location = path.join(file_location, 'CoMPaseD_Config.txt')

        with open(config_location, "w+") as f:
            for items in config_list:
                f.writelines(items)
                f.writelines("\n")
    # if overwrite / create is not possible due to missing permissions, ask where to save
    except PermissionError:
        config_file, _ = QtW.QFileDialog.getSaveFileName(tab_wdg, 'Save File')
        f = open(config_file, mode='w+')
        for items in config_list:
            f.writelines(items)
            f.writelines("\n")
    return None


def config_load_defaults(tab_wdg):
    """output config values from file to gui"""
    # try to open config file by default name and location
    file_location = path.dirname(path.realpath(__file__))
    config_location = path.join(file_location, 'CoMPaseD_Config.txt')

    if not path.isfile(config_location):
        tab_wdg.ErrorLabelCorruptedFile.setText("Error. File 'CoMPaseD_Config.txt' was not found. Please save again.")
    else:
        try:
            with open(config_location, "r") as f:
                config_list = list()
                for line in f:
                    config_list.append(line.rstrip())
        except FileNotFoundError:
            config_file, _ = QtW.QFileDialog.getOpenFileName(tab_wdg, 'Open File')
            try:
                f = open(config_file, mode='r')
                if f is not None:
                    config_list = list()
                    for line in f:
                        config_list.append(line.rstrip())
            except Exception:
                config_list = list()
        if len(config_list) > 0:
            config_write_values(tab_wdg=tab_wdg, value_list=config_list)
    return None


def config_browse_crux(tab_wdg):
    """browse for location of crux, clips and promast"""
    # check current entry and use this path as starting point when it looks like a path
    if path.isdir(str(path.dirname(tab_wdg.CruxPathField.text()))):
        if 'linux' in platform:
            file_name, _ = QtW.QFileDialog.getOpenFileName(tab_wdg, "Crux Executable",
                                                           directory=path.dirname(tab_wdg.CruxPathField.text()),
                                                           filter="Crux (*);;Application (*.exe)")
        else:
            file_name, _ = QtW.QFileDialog.getOpenFileName(tab_wdg, "Crux Executable",
                                                           directory=path.dirname(tab_wdg.CruxPathField.text()),
                                                           filter="Application (*.exe);;Crux (*)")
    # if there is not a path in the text field, use working directory as starting point
    else:
        if 'linux' in platform:
            file_name, _ = QtW.QFileDialog.getOpenFileName(tab_wdg, "Crux Executable",
                                                           directory=getcwd(),
                                                           filter="Crux (*);;Application (*.exe)")
        else:
            file_name, _ = QtW.QFileDialog.getOpenFileName(tab_wdg, "Crux Executable",
                                                           directory=getcwd(),
                                                           filter="Application (*.exe);;Crux (*)")
    # if the user cancels and does not select a file, do not change text field content
    if file_name != "":
        tab_wdg.CruxPathField.setText(file_name)
    return None


def config_browse_clips(tab_wdg):
    """browse for location of clips"""
    if path.isdir(str(path.dirname(tab_wdg.ClipsPathField.text()))):
        file_name, _ = QtW.QFileDialog.getOpenFileName(tab_wdg, "Clips perl script",
                                                       directory=path.dirname(tab_wdg.ClipsPathField.text()),
                                                       filter="Perl script (*.pl);;All files (*)")
    else:
        file_location = path.join(path.dirname(path.realpath(__file__)), '../bin', 'perl')
        if path.exists(file_location):
            file_name, _ = QtW.QFileDialog.getOpenFileName(tab_wdg, "Clips perl script",
                                                           directory=file_location,
                                                           filter="Perl script (*.pl);;All files (*)")
        else:
            file_name, _ = QtW.QFileDialog.getOpenFileName(tab_wdg, "Clips perl script",
                                                           directory=getcwd(),
                                                           filter="Perl script (*.pl);;All files (*)")
    if file_name != "":
        tab_wdg.ClipsPathField.setText(file_name)
    return None


def config_browse_promast(tab_wdg):
    """browse for location of promast"""
    if path.isdir(str(path.dirname(tab_wdg.PromastPathField.text()))):
        file_name, _ = QtW.QFileDialog.getOpenFileName(tab_wdg, "Promast perl script",
                                                       directory=path.dirname(tab_wdg.PromastPathField.text()),
                                                       filter="Perl script (*.pl);;All files (*)")
    else:
        file_location = path.join(path.dirname(path.realpath(__file__)), '../bin', 'perl')
        if path.exists(file_location):
            file_name, _ = QtW.QFileDialog.getOpenFileName(tab_wdg, "Promast perl script",
                                                           directory=file_location,
                                                           filter="Perl script (*.pl);;All files (*)")
        else:
            file_name, _ = QtW.QFileDialog.getOpenFileName(tab_wdg, "Promast perl script",
                                                       directory=getcwd(),
                                                       filter="Perl script (*.pl);;All files (*)")
    if file_name != "":
        tab_wdg.PromastPathField.setText(file_name)
    return None


def config_initialise(tab_wdg):
    """try loading values from default config file at startup"""
    file_location = path.dirname(path.realpath(__file__))
    config_location = path.join(file_location, 'CoMPaseD_Config.txt')

    if not path.isfile(config_location):
        tab_wdg.ErrorLabelCorruptedFile.setText("Error. File 'CoMPaseD_Config.txt' was not found. Please save again.")
    else:
        try:
            with open(config_location, "r") as f:
                config_list = list()
                for line in f:
                    config_list.append(line.rstrip())
            current_error_state = error_corrupted_config(value_list=config_list)
            tab_wdg.ErrorLabelCorruptedFile.setText(current_error_state)
            config_write_values(tab_wdg=tab_wdg, value_list=config_list)
            if current_error_state == "":
                tab_wdg.tabs.setCurrentIndex(1)
        except FileNotFoundError:
            # do not force the user to select a config file when there is none available
            pass
    return None
