(setup-installation)=
# Installation

CoMPaseD is written in Python 3 and has been tested on Windows (≥ Windows 10) and Linux machines. Installation and program start should be done via the command prompt while configuring and running CoMPaseD analyses is possible via the command prompt or the graphical interface.

## Dependencies

### Crux Toolkit

Binaries for the crux mass spectrometry analysis toolkit are available via the official [crux webpage](https://crux.ms/download.html).
Download the binaries suitable for your operating system (Windows/Linux) (version ≥ 3.2).

### Python

To run CoMPaseD python (versions 3.5–3.12 tested) is required.
```{note} While newer python versions may be compatible, some dependencies may not be immediately available for the latest releases. It is recommended to use a tested version to ensure full functionality.
```

Install python using:

````{tabs}

```{group-tab} Windows

   1. Download Python for your [Windows distribution](https://www.python.org/downloads/windows/).

   2. Install Python by executing the downloaded file and following the instructions.

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

Install the Python requirements using `pip` and `requirements.txt`:


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
   The Git command may not be available by default in many Windows installations. You may need to install [Git](https://git-scm.com/downloads) in advance.

3. **Create and Activate a Virtual Environment**

   Set up a virtual environment to manage dependencies:

       C:\Programs\python\python.exe -m venv venv
       C:\Programs\CoMPaseD\venv\Scripts\activate.bat

4. **Install Dependencies**

    Install the required Python packages:

       python -m pip install -r requirements.txt

5. **Run CoMPaseD**

    Launch the graphical interface:

       CoMPaseD_gui.py

6. **Configure Crux Path**

    Set the path to the Crux toolkit within the GUI and save the configuration. CoMPaseD is now ready to use.

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

```

````

```{figure} ../assets/images/configuration-tab.png
:width: 75 %
:name: configuration_initial
Configuration tab upon first start of CoMPaseD. The red error message will disapear when the path to the crux executable is properly set and the configuration is saved.
```  

You can now start CoMPaseD at any time by activating the virtual environment and running:


- `CoMPaseD_gui.py` (Graphical Interface)
- `CoMPaseD_cli.py` (Command-Line Interface)