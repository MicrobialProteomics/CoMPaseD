# Configuration File  

The configuration file defines default paths and settings for CoMPaseD analyses. If no graphical interface is available, the file can be manually created. It must be saved in the `lib` folder of CoMPaseD as **CoMPaseD_Config.txt**.  

The file consists of six lines, each specifying a particular setting:  

- **Line 1:** Absolute path to the Crux executable (*e.g.*, `C:\Programs\crux\bin\crux.exe`).  
- **Line 2:** Indicator of whether to use the original Perl-based peptide mapping scripts. For large *FASTA* files, the Python-based implementation may encounter memory errors. In such cases, activating the Perl-based version is recommended. To ensure proper configuration of Perl, it is advisable to install the **Trans-Proteomics Pipeline** (*Default: False; Options: True or False*).  
- **Line 3:** Absolute path to the `clips.pl` Perl script (*default location: `bin/Perl` folder of CoMPaseD*).  
- **Line 4:** Absolute path to the `promast.pl` Perl script (*default location: `bin/Perl` folder of CoMPaseD*).  
- **Line 5:** *Deprecated*.  
- **Line 6:** Indicator of whether to save peptide lists generated during random sampling as output (*Default: True; Options: True or False*).  
  
---
  
## Default Configuration File
```
C:/Programs/crux/bin/crux.exe
False
C:/Programs/CoMPaseD/bin/Perl/clips.pl
C:/Programs/CoMPaseD/bin/Perl/promast.pl
False
True

```