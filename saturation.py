
#####################################################################
# saturation.py
#
# Saturation initializer. Saturates waters with salts
#
# Sanjoy Som, December 2022
#######################################################################
# Initialization
#######################################################################
import os
import sys
import itertools
import eq36python as eq
from itertools import islice
from scipy import interpolate
#######################################################################
# Common Functions
#######################################################################
def extract_path_files(file):
  cwdpath, eqpath, platform = [x for x in file]
  return cwdpath, eqpath, platform

def extract_wr_files(file):
  f3i, _, f6i, _, _, db, th_NaCl = [x for x in file]
  return f3i, th_evap, f6i, db
#######################################################################
# Saturation
#######################################################################
def saturation(T,WR,filename3i,tophalf,platform,db,cwdpath,\
      eqpath,drop_minerals=False,move_files=False):
  '''
  All the steps needed to saturate water using EQ3/6
  '''
  #first, set the T and WR in the 3i file
  eq.edit_in3i('|Temperature',T,cwdpath,'Temp',filename3i)
  eq.edit_in3i('|Aq.',WR,cwdpath,'W:R',filename3i)
  #Now we are ready to run EQ3
  eq.run_eq(3,eqpath,cwdpath,filename3i,db,platform)
  #Now update the temperature in the rock file in LINEAR TRACKING
  eq.edit_in6i('|'+13*' '+'Value (C',T,cwdpath,'Temp',tophalf)
  #now add the output of EQ3 under the rock file
  filename6i = 'brine.6i'
  eq.stitch(filename3i,tophalf,filename6i)
  if drop_minerals:
    #drop minerals at the end of the run
    eq.edit_iopt_in_6i(filename6i,cwdpath,7,'Do it')
  #And run EQ6
  conv=eq.run_eq(6,eqpath,cwdpath,filename6i,db,platform)
  if move_files:
    #Move files
    #rename and move to output folder
    newbrinefile = 'r6-purew-sat'
    filename6o = 'brine.6o'
    os.system('cp '+filename6o+' '+newbrinefile+'.6o')
    os.system('mv '+newbrinefile+'.6o'+' plotting_files/raw6o')
  return conv
