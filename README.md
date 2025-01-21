[![build-sphinx-documentation](https://github.com/MicrobialProteomics/CoMPaseD/actions/workflows/build-documentation.yaml/badge.svg?branch=main)](https://github.com/MicrobialProteomics/CoMPaseD/actions/workflows/build-documentation.yaml)

# CoMPaseD
**Comparison of Multiple Protease Digestions**

CoMPaseD is a bioinformatics tool to select suitable protease combinations in proteomics experiments involving multiple proteases.

---

## Requirements
To run CoMPaseD, you will need:
- **Python**: Version 3.5–3.12 (tested). [Download Python](https://www.python.org/downloads/)
- **Crux Toolkit**: Download and install from the [Crux website](https://www.crux.ms).

---

## Installation
Follow these steps to install CoMPaseD:

1. **Set Installation Directory**  
   Open a command-line window and navigate to the directory where you want to install CoMPaseD:
   ```bash
   cd C:\Programs\CoMPaseD
   ```

2. **Clone the Repository**  
   Copy the CoMPaseD repository into the chosen directory:
   ```bash
   git clone https://github.com/MicrobialProteomics/CoMPaseD.git
   ```

3. **Create and Activate a Virtual Environment**  
   Set up a virtual environment to manage dependencies:
   ```bash
   C:\Programs\Python\python.exe -m venv venv
   C:\Programs\CoMPaseD\venv\Scripts\activate.bat
   ```

4. **Install Dependencies**  
   Install the required Python packages:
   ```bash
   python -m pip install -r requirements.txt
   ```

5. **Run CoMPaseD**  
   Launch the graphical interface:
   ```bash
   CoMPaseD_gui.py
   ```

6. **Configure Crux Path**  
   Set the path to the Crux toolkit within the GUI and save the configuration. CoMPaseD is now ready to use.  

   You can start it at any time by activating the virtual environment and running:
   - `CoMPaseD_gui.py` (Graphical Interface)
   - `CoMPaseD_cli.py` (Command-Line Interface)
 
---
    
## License
CoMPaseD is distributed under the [MIT License](LICENSE).

For questions or contributions, please open an issue or contact us via the repository.

