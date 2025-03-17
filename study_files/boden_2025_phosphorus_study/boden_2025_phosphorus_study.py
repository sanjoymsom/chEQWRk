#####################################################################
# executor.py
#
# Program Launcher. Performs serpentinization with brines
#
# Sanjoy Som, started December 2022
#######################################################################
import os
import sys
import time
import pandas as pd
import numpy as np
from os.path import expanduser
import matplotlib.pyplot as plt
from sys import platform
#custom libraries
import serpentinization as serp
import eq36python       as eq
import plotcsv          as pc
start=time.time()
#######################################################################
# I/O INIT
#######################################################################
studyfolder  = 'boden_2025_phosphorus_study/'
home         = expanduser("~")
cwdpath      = os.getcwd()+'/'
eqpath       = home+'/EQ3_6v8.0a/'
brine6opath   = '' #not used
raw6opath    = 'plotting_files/raw6o/'
csv_destination = 'plotting_files/csv/'
saved_folder = 'saved_files/'+studyfolder 
platform     = 'linux'
serp_db      = 'tde'   # <---- database selection
#serp_db      = 'mbn'  # <---- database selection
evap_db      = '' #not used
if serp_db == 'tde':
  filename3i   = 'swmajp_tde.3i'
  tophalf_rock = 'tophalf_mbn_oliss_wP.rock'
elif serp_db == 'mbn':
  filename3i   = 'swmajp_mbn.3i'
  tophalf_rock = 'tophalf_mbn_oliss.rock'
filename6i   = 'serp.6i'
tophalf_evap = '' #not used
tophalf_NaCl = '' #not used
path_files   = [cwdpath, eqpath, platform, brine6opath]
wr_files     = [filename3i, tophalf_rock, filename6i, serp_db, 
                tophalf_evap, evap_db, tophalf_NaCl]

#load input files
os.system(' '.join(['cp',cwdpath+'input_files/'+studyfolder+'*.*',cwdpath]))

def ask_questions():
  print(51*'-')
  print('      chEQWRk v0.2 - Boden et al. 2025 ')
  print(51*'-')
  print('Choose from the setup options:')
  print('0:  clear files')
  print('1:  load data')
  print(51*'-')
  print('Or choose from these other data creation options:')
  print('2:  seawater            + serp')
  print('3:  ternary of seawater + serp')
  print(51*'-')
  print('Or from the following data consolidation options:')
  print('4:  (re)create  csv(s)')
  print('5:  compile csvs to one file (use only for ternary)')
  print('6:  zip & save')
  print(51*'-')
  print('Or from the following plotting options:')
  print('7:  plot 1 variable vs another')
  print('8:  plot 2 variables vs another')
  print('9:  plot the ratio of 2 vars vs another')
  print('10: plot minerals')
  print('11: plot ternary diagram')
  print(51*'-')
  print('Water file used: '+ filename3i)
  print(51*'-')
  option=input('Enter choice: ')
  option=CheckOption(option)
  return option

def CheckOption(option):
  #make sure option is an integer
  try: option=int(option)
  except:  raise TypeError('Option must be an integer')
  #catch if options are irrelevant
  if (option > 11) or (option < 0):
    raise Exception('Not a valid option.')
  return option

def read_colorfile(colorfile):
  dfm = pd.read_csv(colorfile)
  dfm = dfm.set_index('Variables')
  dfm = dfm.drop(dfm.columns[0],axis=1)
  return dfm

def update_wr_files(tophalf_evap):
  wr_files2 = [filename3i, tophalf_rock, filename6i, serp_db,
              tophalf_evap, evap_db, tophalf_NaCl]
  return wr_files2

def cleanup():
  os.system('mv *.temp *.txt *.6* *.3* output_files/ 2>/dev/null')
  os.system('mv *.evap *.rock *.halite *.heat output_files/ 2>/dev/null')
  os.system('mv *.csv output_files/ 2>/dev/null')

def duration():
  print('That took {:.2f} seconds'.format(time.time()-start))
#######################################################################
# EXEC FUNCTIONS
#######################################################################
def zip_and_save():
  import shutil
  _checkfolder('saved_files')
  zip_format = 'zip'
  output_filename=input('Enter filename (w/o extension): ')
  shutil.make_archive(output_filename, zip_format, 'plotting_files')
  os.system('mv '+output_filename+'.'+zip_format+' saved_files/'+studyfolder)

def unzip_and_load(saved_data):
  import shutil
  available_zips = os.listdir(saved_data)
  df = pd.DataFrame(available_zips,columns=['Available archives'])
  print(df) #part of code, don't delete
  ans = input('Enter index of archive to load: ')
  archive_to_load = df.at[int(ans),'Available archives']
  #remove existing data
  os.system('rm -r '+cwdpath+'plotting_files/raw6o/')
  os.system('rm -r '+cwdpath+'plotting_files/brine6o/')
  os.system('rm -r '+cwdpath+'plotting_files/csv/')
  #unpack
  zip_format = 'zip'
  shutil.unpack_archive(saved_data+archive_to_load, 'plotting_files', zip_format)
  print('Data from '+archive_to_load+' is loaded.') 
  print('Run code again to start plotting.')

def _checkfolder(foldername):
  check = os.path.isdir(foldername)
  if not check:
    os.makedirs(foldername)
    print('Created folder: '+foldername)

def check_output_folders():
  _checkfolder('output_files')
  _checkfolder('saved_files')
  _checkfolder('saved_files/'+studyfolder)
  _checkfolder('plotting_files/raw6o')
  _checkfolder('plotting_files/brine6o')
  _checkfolder('plotting_files/csv')

def check_database():
  itexists = os.path.isfile(eqpath+'db/data1.'+serp_db)
  if itexists: pass
  else: 
    print(eqpath+'db/data1.'+serp_db+' is missing.')
    print('Error: EQ3/6 database is not installed properly.');sys.exit()

def add_WR_to_compil(cwdpath,csv_destination,compilfile='Compil.csv'):
  data = pd.read_csv(cwdpath+csv_destination+compilfile)
  df   = pd.DataFrame(data)
  df['WR']=df['Variables'].str.split('-').str[4].astype(float)
  df.to_csv(cwdpath+csv_destination+compilfile, index=False)

def execute_query(option):
  if option == 0:
    os.system('rm -r '+cwdpath+'plotting_files/raw6o/')
    os.system('rm -r '+cwdpath+'plotting_files/brine6o/')
    os.system('rm -r '+cwdpath+'plotting_files/csv/')
    os.system('mv '+cwdpath+'minerals.colors'+' 2>/dev/null')
    os.system('mkdir plotting_files/raw6o/')
    os.system('mkdir plotting_files/brine6o/')
    os.system('mkdir plotting_files/csv/')
    print('output files have been deleted')
  elif option == 1:
    unzip_and_load(saved_folder)
  elif option == 2: #seawater + serp
    T_start  = 10
    T_array  = [400]
    WR_array = [0.2,1]
    OLIwt    = [95] #95
    serp.react_heating(T_start,T_array,WR_array,OLIwt,path_files,wr_files)
  elif option == 3: #ternary of seawater + serp
    T_start  = 300
    T_array  = [T_start]
    WR_array = [0.2]
    density = [10000] #10000 gives good coverage --> 5134 actual density
    serp.react_heating(T_start,T_array,WR_array,density,path_files,wr_files,ternary=True)
  elif option == 4:
    eq.generate_tablefrom6o('waterrock',cwdpath,eqpath,raw6opath,\
      csv_destination)
  elif option == 5:
    eq.compile_csvs_into_one_file(cwdpath, csv_destination)
    add_WR_to_compil(cwdpath,csv_destination,'Compil.csv')
  elif option == 6:
    zip_and_save()    
  elif option == 7:
    plot_lit   = True
    plot_field = True
    pc.setup_and_plot_vars(0,cwdpath,csv_destination,plot_lit,plot_field,default_label='Boden2025')
  elif option == 8:
    plot_lit   = False
    plot_field = False
    pc.setup_and_plot_vars(1,cwdpath,csv_destination,plot_lit,plot_field,default_label='Boden2025')
  elif option == 9:
    plot_lit   = False
    plot_field = False
    pc.setup_and_plot_vars(2,cwdpath,csv_destination,plot_lit,plot_field,default_label='Boden2025')
  elif option == 10: #plot minerals
    visualize = False # visualized the color palette used to display minerals
    colors_exist = os.path.isfile('plotting_files/mineral.colors')
    if colors_exist:
      print('----------------------------------------------------------------------')
      print('Warning: mineral.colors file exists in plotting_files/. Unless you are')
      print('         comparing figures from a previous run, it should be removed.')
      print('         mineral.colors exists to preserve the colors chosen to')
      print('         represent minerals between runs.')
      print('----------------------------------------------------------------------')
      dfm = read_colorfile('plotting_files/mineral.colors')
    else:
      dfm = pc.mineral_superset(cwdpath)
    pc.plot_minerals(dfm,cwdpath,visualize)
  elif option == 11:
    mcb = True #map compositional baseline
    pc.csv2ternary(cwdpath,csv_destination,'Compil.csv',mcb)
#######################################################################
# MAIN FUNCTION
#######################################################################
def main():
  #check database
  check_database()
  #create storage folders
  check_output_folders()
  #query the user
  option = ask_questions()
  #execute user's wishes
  execute_query(option)
#######################################################################
# MAIN PROGRAM
#######################################################################
if __name__ == '__main__':
  main()
  cleanup()
  duration()
