#######################################################################
# chEQWRk launcher 
#
# This code looks for study files and invites user to select a study
# 
# Sanjoy Som - started January 2024
#######################################################################
import os 
import sys
import pandas as pd
#######################################################################
# I/O init
#######################################################################
if sys.version_info[0] < 3:
    print("Python3 is required to run this code.");sys.exit()
cwdpath      = os.getcwd()+'/'
studyfolder  = cwdpath + 'study_files/'
inputfolder  = cwdpath + 'input_files/'
#######################################################################
# I/O functions
#######################################################################
def load_study_files(studyfolder):
  print(51*'-')
  print('    chEQWRk - EQ3/6 Water:Rock Reactions wrapper ')
  print('       (disable this launcher in executor.py) ')
  print(51*'-')
  available_study_files = os.listdir(studyfolder)
  df = pd.DataFrame(available_study_files,columns=['Available studies'])
  print(df) #part of code, don't delete
  ans = input('Enter index of archive to load: ')
  study_to_load = df.at[int(ans),'Available studies']
  return study_to_load

def launch_study(study):
  os.system('cp '+studyfolder+study+'/'+study+'.py '+cwdpath)
  os.system('python3 '+study+'.py')

def clean(study):
  os.system('rm *study*.py') #removes study executors after runs (declutter)
#######################################################################
# Support functions
#######################################################################
def perform_checks(inputfolder,studyfolder,study):
  #check that study folder exists in input_files
  itExists = os.path.exists(inputfolder+study)
  if itExists:
    pass
  else:
    print(51*'-')
    print('Error: input files for the study do not exist.')
    print('       Create study folder in input_files/')
    sys.exit()
  #check that study .py file has same name as study
  itExists = os.path.exists(studyfolder+study+'/'+study+'.py')
  if itExists:
    pass
  else:
    print(51*'-')
    print('Error: study .py file does not exist or does not have')
    print('       the same name as study folder.')
    print('       Create or rename .py file in study_files/study')
    sys.exit()
#######################################################################
# Main function
#######################################################################
def main(studyfolder):
  #study_to_load = 'som_2024_brines_study'       #uncomment to set default
  #study_to_load = 'boden_2025_phosphorus_study' #uncomment to set default
  study_to_load = load_study_files(studyfolder) #comment to set default
  perform_checks(inputfolder,studyfolder,study_to_load)
  launch_study(study_to_load)
  clean(study_to_load)
#######################################################################
# MAIN PROGRAM
#######################################################################
if __name__ == '__main__':
  main(studyfolder)
