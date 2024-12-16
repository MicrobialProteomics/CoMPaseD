import colorama
from sys import stdout
import PyQt6.QtGui as QtGui
import PyQt6.QtCore as QtCore


def rel_pos(old_x, old_y, old_width, old_height, return_type='QRect'):
    """convert QRect positions by screen resolution, returns QRect obj"""

    # layout was optimised for 1200x900 pt on a 1920x1040 pt monitor
    avail_height = QtGui.QGuiApplication.primaryScreen().availableGeometry().height()
    avail_width = QtGui.QGuiApplication.primaryScreen().availableGeometry().width()
    new_width = int((old_width / 1920) * avail_width)
    new_height = int((old_height / 1040) * avail_height)

    new_x = int((old_x / 1920) * avail_width)
    new_y = int((old_y / 1040) * avail_height)

    if return_type == 'QRect':
        return QtCore.QRect(new_x, new_y, new_width, new_height)
    else:
        return new_x, new_y, new_width, new_height


# may not be required for param_obj but argparse parameter file needs conversion
def config_to_numeric_list(config_list, invalid_characters=" []"):
    """Convert config values to numeric"""

    for invalid_character in invalid_characters:
        config_list = config_list.replace(invalid_character, "")
    config_list = list(config_list.split(","))
    numeric_config_list = list()
    for config_element in config_list:
        try:
            numeric_config_list.append(float(config_element))
        except ValueError:
            print(f"{colorama.Fore.RED}ERROR: Configuration corrupted. Please check. Maybe confused 0 and \'o\' resulting in non-numeric? {colorama.Style.RESET_ALL}", flush=True)
            print("", flush=True)
            stdout.flush()

    config_list = numeric_config_list
    return config_list

def rel_pos_fig(x_y_pts: list):
    x_dpi = QtGui.QGuiApplication.primaryScreen().logicalDotsPerInchX()
    y_dpi = QtGui.QGuiApplication.primaryScreen().logicalDotsPerInchY()

    x_y_inch = list()
    x_y_inch.append(x_y_pts[0] / x_dpi)
    x_y_inch.append(x_y_pts[1] / y_dpi)

    return x_y_inch
