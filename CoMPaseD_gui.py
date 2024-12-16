import sys
import os
import colorama
from lib.CoMPaseD_gui_tabs import *
from sys import argv
from os import path
from sys import platform

try:
    from PyQt6 import QtWidgets as QtW
except Error as e:
    # only for debugging purposes
    print(f"{colorama.Fore.RED}ERROR: PyQt6 import failed due to error: {e}. Please check.{colorama.Style.RESET_ALL}")
    raise ImportError


# MainWindow class
class CoMPaseD(QtW.QMainWindow):
    """CoMPaseD main window class"""
    def __init__(self):
        super().__init__()
        self.setWindowTitle('CoMPaseD - Comparison of Multiple-Protease Digestions')
        # set to static 1200x900 pt resolution and place at center of the main-screen
        x, y, w, h = rel_pos(0, 0, 1200, 930, return_type='list')
        self.setGeometry(x, y, w, h)
        center_point = QtGui.QGuiApplication.primaryScreen().availableGeometry().center()
        qt_rectangle = self.frameGeometry()
        qt_rectangle.moveCenter(center_point)
        self.move(qt_rectangle.topLeft())
        qt_rectangle = self.frameGeometry()
        self.setGeometry(qt_rectangle)

        # tabs need to be inserted as a widget and
        # collection of them is created as an own class
        self.table_widget = CoMPaseD_Tabs(self)
        self.setCentralWidget(self.table_widget)

        # change window item
        file_location = path.dirname(path.realpath(__file__))
        item_img = path.join(file_location, "bin", "CoMPaseD_logo.png")
        if path.isfile(item_img):
            self.setWindowIcon(QtGui.QIcon(item_img))

# create a QApplication and a MainWindow object to execute
CoMPasedApplication = QtW.QApplication(argv)
MainWindow = CoMPaseD()

# find relative location of css styles
file_location = path.dirname(path.realpath(__file__))
css_location = path.join(file_location, 'lib/CoMPaseD_gui_style.css')
with open(css_location, 'r') as f:
    MainWindow.setStyleSheet(f.read())

# for win os explicitly state UID to indicate that CoMPaseD is its own python instance
# cygwin is untested
if ("win32" in platform) or ("cygwin" in platform):
    import ctypes
    CoMPasedApplicationID = 'CoMPaseD'
    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(CoMPasedApplicationID)

MainWindow.show()
CoMPasedApplication.exec()
