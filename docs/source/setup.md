(setup-installation)=
# Installation

CoMPaseD is written in Python 3 and has been tested on Windows (≥ Windows 10) and Linux machines. Installation and program start should be done via the command prompt while configuring and running CoMPaseD analyses is possible via the command prompt or the graphical interface.

## Dependencies

### Crux Toolkit

Binaries for the crux mass spectrometry analysis toolkit are available via the official [crux webpage](https://crux.ms/download.html).
Download the binaries suitable for your operating system (Windows/Linux) (version ≥ 3.2).

### Python

To run CoMPaseD python Version 3.5–3.12 (tested) is required.

Install python using:

````{tabs}

```{group-tab} Windows

   1. Download python for your [Windows distribution](https://www.python.org/downloads/windows/).

   2. Install python by executing the downloaded file and following the instructions.

```

```{group-tab} Linux

   - Ubuntu/Debian:

       sudo apt update
       sudo apt install -y python3 python3-pip

   - Fedora

       sudo dnf install -y python3 python3-pip

   - Arch Linux

       sudo pacman -Syu python python-pip

```
````

### Python packages

The Python packages with suitable version numbers are listed in the provided `requirements.txt`.

```{dropdown} requirements.txt

| Package name                | Version    |
|-----------------------------|------------|
| bio                         | &ge;1.5.9  |
| biopython                   | &ge;1.81   |
| colorama                    | &ge;0.4    |
| keras                       | &ge;2.13.1 |
| matplotlib                  | &#61;3.7.2 |
| numpy                       | &ge;1.22   |
| pandas                      | &ge;2.0    |
| PyQt6-Qt6                   | &ge;6.5    |
| PyQt6-sip                   | &ge;13.5   |
| PyQt6                       | &ge;6.5    |
| python-dateutil             |            |
| tensorflow-estimator        |            |
| tensorflow-intel            |            |
| tensorflow-io-gcs-filesystem|            |
| tensorflow                  | &ge;2.12.0 |

```

```{tip} We recommend using a virtual environment and pip to install these dependencies.
```

Install the python requirements using `pip` and `requirements.txt`:


````{tabs}

```{group-tab} Windows

1. **Set Installation Directory**

   Open a command-line window and navigate to the directory where you want to install CoMPaseD:

       mkdir C:\Programs\CoMPaseD
       cd C:\Programs\CoMPaseD


2. **Clone the Repository**

   Copy the CoMPaseD repository into the chosen directory:

       git clone https://github.com/MicrobialProteomics/CoMPaseD.git

   ```{tip}
   You can use [Git Bash](https://git-scm.com/downloads) for a better command-line experience on Windows.

3. **Create and Activate a Virtual Environment**

   Set up a virtual environment to manage dependencies:

       C:\Programs\Python\python.exe -m venv venv
       C:\Programs\CoMPaseD\venv\Scripts\activate.bat

4. **Install Dependencies**

    Install the required Python packages:

       python -m pip install -r requirements.txt

5. **Run CoMPaseD**

    Launch the graphical interface:

       CoMPaseD_gui.py

6. **Configure Crux Path**

    Set the path to the Crux toolkit within the GUI and save the configuration. CoMPaseD is now ready to use.

    You can start it at any time by activating the virtual environment and running:


    - `CoMPaseD_gui.py` (Graphical Interface)
    - `CoMPaseD_cli.py` (Command-Line Interface)

```

```{group-tab} Linux

1. **Set Installation Directory**

   Open a command-line window and navigate to the directory where you want to install CoMPaseD:

       mkdir -p  ~/Programs/CoMPaseD
       cd ~/Programs/CoMPaseD


2. **Clone the Repository**

   Copy the CoMPaseD repository into the chosen directory:

       git clone https://github.com/MicrobialProteomics/CoMPaseD.git

   ```{tip}
   If you do not have `git` installed. Install it using ```sudo apt install git``` or your prefered package manager.

3. **Create and Activate a Virtual Environment**

   Set up a virtual environment to manage dependencies:

       python -m venv venv
       source venv/bin/activate

4. **Install Dependencies**

    Install the required Python packages:

       python -m pip install -r requirements.txt

5. **Run CoMPaseD**

    Launch the graphical interface:

       CoMPaseD_gui.py

6. **Configure Crux Path**

    Set the path to the Crux toolkit within the GUI and save the configuration. CoMPaseD is now ready to use.

    You can start it at any time by activating the virtual environment and running:


    - `CoMPaseD_gui.py` (Graphical Interface)
    - `CoMPaseD_cli.py` (Command-Line Interface)

```

````




# Quickstart

CoMPaseD provides both a command-line interface (CLI) and a graphical user interface (GUI).

```{note} The following command-line calls might slightly differ depending on your python setup.
```

## Commandline Interface (CLI)

To use the CLI run the `CoMPaseD_CLI.py` command.
All values for the parameters can be adjusted in the `para.meters` file.

Running a full analysis can be done by calling the CLI file together wih the parameter file.

```bash
python CoMPaseD_CLI.py -p path/to/parameter.param
```

For running individual steps of the tool, the export (-e), digestion (-d) or result analysis (-a) flags.

```bash
python CoMPaseD_CLI.py -p path/to/parameter.param [-e] [-d] [-a]
```

```{note} The [-a] flag requires you to run the rest of the analysis beforehand [-e] and [-d].
```

For a complete listing of command-line options use the help flag or refer to the [main documentation](documentation-cli).

```bash
python CoMPaseD_CLI.py -h
```

## Graphical User Interface (GUI)

The GUI can be accessed using `CoMPaseD_gui.py` command.

```bash
python CoMPaseD_gui.py
```

```{important} Upon the first program start an error message in the configuration tab indicates errors in the configuration file. Adjust the path to the crux executable and save the altered configuration. For subsequent program starts, the configuration tab should be skipped and CoMPaseD shows the analysis parameter tab directly.
When the checkbox to use original Perl-based ProteoMapper scripts is checked, paths to these scripts must be selected as well and the configuration saved again.
```

![configuration](../assets/images/configuration-tab.png)

For a detailed description of the GUI, please refer to the [main documentation](documentation-gui).


