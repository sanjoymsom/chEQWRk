#######################################################################
# READ.me last update: March 17, 2025 - Sanjoy M. Som                 |
# (for study-specific READ.me files, see study_files/)                |
#                                                                     |
#######################################################################
# INTRODUCTION                                                        |
#######################################################################
chEQWRk is a series of wrappers code around EQ3/6 designed to perform |
Water:Rock reactions. "Study files" are built that call these wrappers|
to solve specific problems. chEQWRk is designed to be as robust and   |
flexible as possible, but you may encounter errors/bugs. Sorry. Thanks|
for understanding that anticipating every contingency is next to      |
impossible. EQ3/6 needs to be installed in the root directory.        |
See installation guidelines.                                          |
#######################################################################
# VERSION MAPPING TO PUBLICATIONS                                     |
#######################################################################
V0.1.1 --> Som et al. 2024                                            |
v0.2.0 --> Boden et al. 2025 (back-compatible                         |
#######################################################################
# KEY FILES                                                           |
#######################################################################
executor.py          : Program launcher.                              |
evaporation.py       : Creates brines.                                | 
serpentinization.py  : Runs the W:R rock reaction.                    |
saturation.py        : Saturates water file with salt.                |
plotcsv.py	         : Generates all the plots.                       | 
eq36python.py        : Function library that act on EQ3/6 input and   |
                       output files.                                  |
Important            : Study-specific .py files are in study_files/   |
                       including study-specific READ.me files         |
                                                                      |
#######################################################################
# DEPENDENCIES                                                        |
#######################################################################
The following python libraries are needed to successfull run chEQWRk: |
  os, sys, time, pandas 1.0.3, numpy 1.19.5, scipy 1.5.4, itertools,  |
  shutil                                                              |
                                                                      |
#######################################################################
# INSTALLATION GUIDELINES                                             |
#######################################################################
To install WSL (Windows System Linux) on Windows 10/11, open the      |
Microsoft Store app in Windows and search for 'Ubuntu'.               |
  Click 'Get' and follow the prompts.                                 |
                                                                      |
Follow this Gist to install EQ3/6 on linux:                           |
  https://gist.github.com/sanjoymsom/529d9100b75d5803a4be54903a759167 |
                                                                      |
On August 31, 2023, EQ3/6 could be found at                            \-|
  https://seaborg.llnl.gov/resources/geochemical-databases-modeling-codes|
                                                                       /-|
An X-window is required to run the code in WSL. XMing is good:        | 
  https://sourceforge.net/projects/xming/                             |
                                                                      |
The code is launched by typing in the terminal: python3 executor.py   |
                                                                      |
#######################################################################
# INITIAL CONSIDERATIONS                                              |
#######################################################################
1. If you're simulating water:rock reactions, do you allow cabonates? |
   If so, put a star in front of SIDERITE, HYDROMAGNESITE, ARTINITE,  |
   NESQUEHONITE, CALCITE, Calcite-SS, MAGNESITE, ARAGONITE, HUNTITE,  |
   DOLOMITE,ORDERED, DOLOMITE,DISORDERED, DOLOMITE in the 3i file     |
   located in the input_files/ folder (this example assumes use of the|
   mbn database).                                                     |
2. It's a good habit to clear everything (option 0) prior to running a|
   set of computations. You can save your results at any time with    |
   option 9, and load those saved files with option 1. Saved file are |
   stored in the /saved_files folder.                                 |
3. Feel free to mess with the source code to improve things. If you   | 
   screw up and need to revert, the original .py files released at    |
   manuscript publication are also saved in /python_backup of the     |
   release corresponding to the manuscript                .           |
                                                                      |
#######################################################################
# START A NEW STUDY                                                   |
#######################################################################
One of the tenets of chEQWRk is that it can store "studies" that can  | 
be easily loaded again. For example, a study can be loaded to recreate| 
data or reproduce figures from a publication that used chEQWRk.       |
                                                                      |
All studies are stored in the study_files/ folder and each have their |
own READ.me. To start a new study, create a new folder in the         |
study_files/ folder and the python file of the study needs to have the|
same name as the folder. It is recommended to use a starting file from|
a previous study as a template.                                       |  
                                                                      |
#######################################################################
# FAQ                                                                 |
#######################################################################
                                                                      |
  - Why is the code called chEQWRk? Well, it has a double-meaning.    |
    First, it is read as "checkwork", implying that it can be used to | 
    check the science that is presented by using it. A nod to open    |
    science. Second, it also can mean 'C'onsiderable 'h'elp when      |
    running 'EQ'3/6 for 'W'ater:'R'ock 'k'alculations.                |      
                                                                      |
#######################################################################
# Troubleshooting                                                     |
#######################################################################
                                                                      |
Setting up Xwindow within Linux in Windows 11 is a bit tricky.        | 
Here is one solution using XMing                                      |
                                                                      |
WSL1: in ~/.bashrc, ensure that DISPLAY is exported to localhost:     |
      export DISPLAY=$localhost:0.0                                   |
                                                                      |
      in /etc/ssh/sshd_config                                         |
      ensure that X11Fowarding is set to yes                          |
                                                                      |
WSL2: in ~/.bashrc, export the DISPLAY as follows:                     \-------------
      export DISPLAY=$(cat /etc/resolv.conf | grep nameserver | awk '{print $2}'):0 |
                                                                       /-------------
      in /etc/ssh/sshd_config                                         |
      ensure that X11Forarding is set to yes, and                     |
      ensure that X11UseLocalHost is set to yes                       |
                                                                      |
#######################################################################
# END of READ.me                                                      |
#######################################################################

