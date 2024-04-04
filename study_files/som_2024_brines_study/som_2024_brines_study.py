#####################################################################
# som_2024_brines_study.py
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
import evaporation      as evap
import saturation       as sat
import eq36python       as eq
import plotcsv          as pc
start=time.time()
#######################################################################
# I/O INIT
#######################################################################
studyfolder  = 'som_2024_brines_study/'
home         = expanduser("~")
cwdpath      = os.getcwd()+'/'
eqpath       = home+'/EQ3_6v8.0a/'
brine6opath   = 'plotting_files/brine6o/'
raw6opath    = 'plotting_files/raw6o/'
csv_destination = 'plotting_files/csv/'
saved_folder = 'saved_files/'+studyfolder
platform     = 'linux'
serp_db      = 'mbn'
evap_db      = 'ypf'
NaCl_db      = 'ypf'
filename3i   = 'swmajp_mbn.3i'
filename6i   = 'serp.6i'
tophalf_rock = 'tophalf_mbn_oliss.rock'
tophalf_evap = 'tophalf_mbn_pCO2.evap'#'tophalf_mbn.evap'
tophalf_NaCl = 'tophalf_mbn.halite'
path_files   = [cwdpath, eqpath, platform, brine6opath]
wr_files     = [filename3i, tophalf_rock, filename6i, serp_db, 
                tophalf_evap, evap_db, tophalf_NaCl]
#load input files
os.system(' '.join(['cp',cwdpath+'input_files/'+studyfolder+'*.*',cwdpath]))

def ask_questions():
  print(51*'-')
  print('READ.me is in study_files/'+studyfolder)
  print(51*'-')
  print('Choose from the setup options:')
  print('0: clear files')
  print('1: load data')
  print(51*'-')
  print('Or choose from the following data creation options:')
  print('2: seawater evaporation')
  print(51*'-')
  print('Or choose from these other data creation options:')
  print('3: pure water         + NaCl to saturation')
  print('4: seawater           + serp')
  print('5: NaCl(sat) seawater + serp')
  print('6: seawater brine     + serp')
  print(51*'-')
  print('Or from the following data verification options:')
  print('7: test brines against evaporation curve')
  print(49*'-')
  print('Or from the following data consolidation options:')
  print('8: (re)create  csv')
  print('9: zip & save')
  print(51*'-')
  print('Or from the following plotting options:')
  print('10: plot 1 variable vs another')
  print('11: plot minerals')
  print('12: plot H2 generated by minerals')
  print('13: plot Fe-content of minerals')
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
  if (option > 13) or (option < 0):
    raise Exception('Not a valid option.')
  return option

def obtain_aw_array_from_user():
  Aw_array = input('Target Aw of brine or [a]ll or [s]et: ')
  if Aw_array == 'a' or Aw_array == 'all':
    Aw_array = np.arange(98,30,-1)
    Aw_array = [str(x/100.) for x in Aw_array]
  elif Aw_array == 's' or Aw_array == 'set':
    Aw_array = [0.9,0.8,0.73,0.6,0.5,0.4,0.31]
  else:
    Aw_array = [Aw_array]
  return Aw_array

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
  os.system('mv *.evap *.rock *.halite output_files/ 2>/dev/null')
  os.system('mv *.csv output_files/ 2>/dev/null')

def duration():
  print('That took {:.2f} seconds'.format(time.time()-start))

def check_EQ36():
  #check that EQ3/6 exists
  does_eq36_exist = os.path.isdir(eqpath)
  if does_eq36_exist:
    print(51*'-')
    print('EQ3/6 folder detected.')
  else:
    print('   chEQWRk needs EQ3/6 installed at root.  It was not found.  See READ.me.')
    print('If it is installed, edit PATH of eqpath variable in '+studfolder+'.py (l.29).')
    sys.exit()
#######################################################################
# EXEC FUNCTIONS
#######################################################################
def clear_files():
  os.system('rm -r '+cwdpath+'plotting_files/raw6o/')
  os.system('rm -r '+cwdpath+'plotting_files/brine6o/')
  os.system('rm -r '+cwdpath+'plotting_files/csv/')
  os.system('rm -r '+cwdpath+'plotting_files/rock6i/')
  os.system('rm '+cwdpath+'output_files/*.*')
  os.system('rm '+cwdpath+'plotting_files/mineral.colors  2>/dev/null')
  os.system('mkdir plotting_files/raw6o/')
  os.system('mkdir plotting_files/brine6o/')
  os.system('mkdir plotting_files/csv/')
  os.system('mkdir plotting_files/rock6i/')
  print('output files have been deleted')

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
  os.system('rm -r '+cwdpath+'plotting_files/rock6i/')
  #unpack
  zip_format = 'zip'
  shutil.unpack_archive(saved_data+archive_to_load, 'plotting_files', zip_format)
  print('Data from '+archive_to_load+' is loaded.')

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
  _checkfolder('plotting_files/rock6i')

def execute_query(option):
  if option == 0:
    clear_files()
  elif option == 1:
    unzip_and_load(saved_folder)
  elif option == 2: #evaporate surface seawater
    tophalf_evap = 'tophalf_mbn_pCO2.evap'
    wr_files2 = update_wr_files(tophalf_evap)
    pH = 8.2 #seawater surface pH
    T_start  = 33.
    WR_array = [1]
    evap.evaporation_and_nothing_else(pH,T_start,WR_array[0],path_files,wr_files2)
    eq.generate_tablefrom6o('evaporation',cwdpath,eqpath,raw6opath,\
      csv_destination)
  elif option == 3: #purewater + NaCl to saturation
    drop_minerals = False # drop minerals from the brine system after formation
    move_files    = True 
    filename3i = 'purew_ypf.3i'
    T_start  = 25.
    WR_array = [1]
    sat.saturation(T_start,WR_array[0],filename3i,tophalf_NaCl,platform,NaCl_db,cwdpath,\
      eqpath,drop_minerals,move_files)
  elif option == 4: #seawater + serp
    T_start  = 10
    T_array  = [400]
    WR_array = [1]
    OLIwt    = [80]
    serp.react_heating(T_start,T_array,WR_array,OLIwt,path_files,wr_files)
  elif option == 5: #NaCl(sat) seawater + serp
    drop_minerals = True # drop minerals from the brine system after formation
    move_files    = False
    T_start  = 10
    T_array  = [400]
    WR_array = [1]
    OLIwt    = [80]
    sat.saturation(T_start,WR_array[0],filename3i,tophalf_NaCl,platform,NaCl_db,cwdpath,\
      eqpath,drop_minerals,move_files)
    serp.react_brine_heating('NaCl(sat)', T_start, T_array, WR_array, OLIwt,\
      path_files, wr_files)
  elif option == 6: #seawater-derived brines + serp
    drop_minerals = True # drop minerals from the brine system after formation
    T_start  = 10
    T_array  = [400]
    WR_array = [1]
    OLIwt    = [80]
    Aw_array = obtain_aw_array_from_user()
    for Aw in Aw_array:
      Aw = evap.make_brine(T_start,WR_array[0],Aw,path_files,wr_files,drop_minerals)
      serp.react_brine_heating(Aw, T_start, T_array, WR_array, OLIwt,\
        path_files, wr_files)
  elif option == 7:
    T = 10
    WR = 1
    evap.test_brines(T,WR,cwdpath,csv_destination,path_files,wr_files,brine6opath)
    pc.setup_and_plot_brines(0,cwdpath,csv_destination,plot_lit=True)
  elif option == 8:
    eq.generate_tablefrom6o('waterrock',cwdpath,eqpath,raw6opath,\
      csv_destination)
  elif option == 9:
    zip_and_save()    
  elif option == 10:
    plot_lit   = True
    plot_field = True
    pc.setup_and_plot_vars(0,cwdpath,csv_destination,plot_lit=True,plot_field=True)
  elif option == 11: #plot minerals
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
  elif option == 12:
    pc.plot_h2_by_mineral(cwdpath,csv_destination)
  elif option == 13:
    pc.plot_Fe_content_by_mineral(cwdpath,csv_destination)
#######################################################################
# MAIN FUNCTION
#######################################################################
def main():
  #check that EQ3/6 exists on system
  check_EQ36()
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