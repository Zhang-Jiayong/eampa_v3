################################################################
#    Main Program
#
#
#
#
################################################################


#!/bin/python3
########################################################################
import os
import time
import datetime
import re
import sys
import shutil
import numpy
import matplotlib.pyplot as plt
import copy
from f2py_lib.f_fnc import fnc
from f2py_lib.f_interp import interp
from f2py_lib.f_efs import efs
from f2py_lib.f_bp import bp
from f2py_lib.f_spline import spline
from f2py_lib.f_bp import polyfit
from f2py_lib.f_es import es
from f2py_lib.f_sorting import sort
import random

###########################################
#  CLASS
###########################################
class g: 

  run_type = 'bp'
  
  dirs = {
         'wd': 'wd',
         'log': 'wd/log',  
         'output': 'wd/output',   
         'results': 'wd/results',   
         'plots': 'wd/plots',  
         'eos': 'wd/plots/eos', 
         'ec': 'wd/plots/ec', 
         'pots': 'wd/plots/pots',  
         'fitting': 'wd/fitting',  
         }
  
  times = {
          'start' : 0.0,
          'end' : 0.0,
          'duration' : 0.0,
          }
          
  pot_functions = {
                  'pot_name': '',
                  'pot_dir': '',
                  'zbl_file': '',
                  'functions': [],
                  'zbl': [],
                  'functions_original': [],
                  }
                  
  configs = {
            'config_files': [],
            'configs': [],
            }  
            
  fit = {
  'pop_size': 16,
  'fresh_size': 8,
  'generations': 10, 
  'spline_cycles': 1, 
  'spline_generations': 5, 
  'extinction_percentage': 50,  
  'extinction_frequency': 3,
  }
            
  dft_energy_adjustments = {}
  
  bulk_properties = {}
  
  labels = {}
  
  tab_size = 1001
  tab_width = 4

  outputs = True
  results_fh = None
  log_fh = None
         
  file_counter = 0 
         
  def file_name():
    globals.file_counter = globals.file_counter + 1
    name = "file_"
    file_counter_str = str(globals.file_counter)
    while(len(file_counter_str) < 6):
      file_counter_str = '0' + file_counter_str
    name = name + file_counter_str    
    return name
         
###########################################
#  CLASS st
###########################################
class std:

  @staticmethod
  def file_to_list(file_name, clean=False):
# Init variable
    file_data = []
# Read it in line by line
    fh = open(file_name, "r")
    for line in fh:
      if(clean):
        line = line.strip()
        if(line != ""):
          file_data.append(line)          
      else:
        file_data.append(line[0:-1])
# Return
    return file_data
    
  @staticmethod
  def split_fields(line, sep=" "):
    out = line.split(sep)
    key = out[0]
    value = out[1]
    value_out = ''    
    indata = False
    for char in value:
      if(indata and char != '"'):
        value_out = value_out + char
      elif(indata and char == '"'):
        indata = False
      elif(not indata and char == '"'):
        indata = True
    return key, value_out
    
  @staticmethod
  def one_space(line, sep=" "):
    out = ''   
    indata = 0
    last_char = None
    for char in line:
      if(indata == 1 and char != "'" and last_char != "\\"):
        out = out + char
      elif(indata == 1 and char == "'" and last_char != "\\"):
        out = out + char
        indata = 0
      elif(indata == 2 and char != '"' and last_char != "\\"):
        out = out + char
      elif(indata == 2 and char == '"' and last_char != "\\"):
        out = out + char
        indata = 0
      elif(indata == 0 and not (char == " " and last_char == " ")):
        out = out + char
      last_char = char
    return out   
    
  @staticmethod
  def to_fields(line, sep=" "):
    out = []
    temp = ''
    indata = 0
    last_char = None
    for char in line:
      if(indata == 1 and char != "'" and last_char != "\\"):
        temp = temp + char
      elif(indata == 1 and char == "'" and last_char != "\\"):
        temp = temp + char
        indata = 0
      elif(indata == 2 and char != '"' and last_char != "\\"):
        temp = temp + char
      elif(indata == 2 and char == '"' and last_char != "\\"):
        temp = temp + char
        indata = 0
      elif(indata == 0 and not (char == sep and last_char == sep)):
        if(char == sep):
          temp = temp.strip()
          if(temp != ""):
            out.append(temp)
            temp = ''
        else:
          temp = temp + char
    
    temp = temp.strip()
    if(temp != ""):
      out.append(temp)      
    return out    
    
  @staticmethod
  def make_dir(dir):
    try:
      if(not os.path.exists(dir) and dir.strip() != ''):
        os.mkdir(dir) 
        return True
      return False
    except:
      return False
    
  @staticmethod
  def remove_comments(content):
    data = ''
    i = 0
    for line in content:
      if(i > 0):
        data += '\n'
      data += line
      i = i + 1
    out = ''
    indata = 0
    incomment = 0
    for i in range(len(data)):
# Get char and next char
      char = data[i]
      next = None
      prev = None
      if(i < len(data)-1):
        next = data[i + 1]
      if(i > 0):
        prev = data[i - 1]
# If in '  '
      if(indata == 1 and char != "'" and last_char != "\\"):
        out = out + char
      elif(indata == 1 and char == "'" and last_char != "\\"):
        out = out + char
        indata = 0
# If in "  "
      elif(indata == 2 and char != '"' and last_char != "\\"):
        out = out + char
      elif(indata == 2 and char == '"' and last_char != "\\"):
        out = out + char
        indata = 0
      elif(indata == 0):
        if(incomment == 0 and char == "/" and next == "/"):
          incomment = 1
        elif(incomment == 1 and char == "\n"):
          incomment = 0
        if(incomment == 0 and char == "!"):
          incomment = 2
        elif(incomment == 2 and char == "\n"):
          incomment = 0
        if(incomment == 0 and char == "/" and next == "*"):
          incomment = 3
        elif(incomment == 3 and prev == "*" and char == "/"):
          incomment = 0
        elif(incomment == 0):
          out = out + char  
    return out.split("\n")    
    
# Remove comments from a block of data/text
  @staticmethod
  def remove_comments_data(data):
    out = ""
    n = 0
    inquotes = 0
    incomment = 0
    while n < len(data):
# Get char and next char
      char = data[n]
      next = None
      prev = None
      if(n < len(data)-1):
        next = data[n + 1]
      if(n > 0):
        prev = data[n - 1]
        
# If in '  '
      if(inquotes == 1 and char != "'" and last_char != "\\"):
        out = out + char
      elif(inquotes == 1 and char == "'" and last_char != "\\"):
        out = out + char
        inquotes = 0
# If in "  "
      elif(inquotes == 2 and char != '"' and last_char != "\\"):
        out = out + char
      elif(inquotes == 2 and char == '"' and last_char != "\\"):
        out = out + char
        inquotes = 0
# If not inside quotes
      elif(inquotes == 0):
# Comment on a line
        if(incomment == 0 and char == "/" and next == "/"):
          incomment = 1
        elif(incomment == 0 and char == "!"):
          incomment = 1
        elif(incomment == 0 and char == "#"):
          incomment = 1    
# Comment on line close
        elif(incomment == 1 and char == "\n"):
          incomment = 0
# Comment block
        elif(incomment == 0 and char == "/" and next == "*"):
          incomment = 3
        elif(incomment == 3 and prev == "*" and char == "/"):
          incomment = 0
        elif(incomment == 0):
          out = out + char  
# Increment counter
      n = n + 1
    return out        

# Single spaces, tabs to spaces
  @staticmethod
  def prep_data(content):
    out = []
    for line in content:
      line_new = std.prep_data_line(line)
      if(line_new != ''):
        out.append(line_new)
    return out  
      
  @staticmethod
  def prep_data_line(line): 
    temp = ''
    indata = 0
    last_char = None
    for char in line:
      if(char == '\t'):
        char = ' '
      if(indata == 1 and char != "'" and last_char != "\\"):
        temp = temp + char
      elif(indata == 1 and char == "'" and last_char != "\\"):
        temp = temp + char
        indata = 0
      elif(indata == 2 and char != '"' and last_char != "\\"):
        temp = temp + char
      elif(indata == 2 and char == '"' and last_char != "\\"):
        temp = temp + char
        indata = 0
      elif(indata == 0 and not (char == ' ' and last_char == ' ')):
        temp = temp + char       
      last_char = char  
    return temp.strip()    
    
  @staticmethod
  def remove_quotes(inp): 
    if(isinstance(inp, list)):    
      for i in range(len(inp)):
        inp[i] = std.remove_quotes(inp[i])        
      return inp
    else:
      inp = inp.strip()
      if(inp[0] == '"' and inp[-1] == '"'):
        return inp[1:-1]
      if(inp[0] == "'" and inp[-1] == "'"):
        return inp[1:-1]
      return inp
      
  @staticmethod
  def config_file_to_list(file_name):
# Init variable
    file_data = []
# Read it in line by line
    fh = open(file_name, "r")
    for line in fh:
      if(line.strip() != ""):
        line = line.strip()
        line = std.remove_comments(line)
        line = std.prep_data_line(line)
        fields = std.to_fields(line)
        file_data.append(fields)         
# Return
    file_data = std.remove_quotes(file_data)
    return file_data
    
  @staticmethod
  def get_dir(file_path):
    directory = ''
    read = False
    for i in range(len(file_path)):
      if(read):
        directory = file_path[-1-i] + directory
      if(file_path[-1-i] == "/"):
        read = True
    return directory
  
  @staticmethod
  def read_csv(filename, sep=","):
    data = []
    if(os.path.isfile(filename)):
# Read from file into memory
      fh = open(filename, 'r')
      file_data = ""
      for line in fh:
        file_data = file_data + line
      fh.close()
# Remove comments
      file_data = std.remove_comments_data(file_data)
# Read Data
      lines = file_data.split("\n")
      for line in lines:
        line = line.strip()
        if(line != ""):
          data.append(line.split(sep))  
    return data
     
  @staticmethod
  def read_csv_array(filename, sep=","):
    data = []
    if(os.path.isfile(filename)):
# Read from file into memory
      fh = open(filename, 'r')
      file_data = ""
      for line in fh:
        file_data = file_data + line.strip() + '\n'
      fh.close()
# Remove double spaces
      file_data = std.one_space(file_data)
# Remove comments
      file_data = std.remove_comments_data(file_data)
# Read Data
      lines = file_data.split("\n")
      for line in lines:
        line = line.strip()
        if(line != ""):
          data.append(line.split(sep))  
          
      lst = []
      for i in range(len(data)):
        row = []
        for j in range(len(data[0])):
          try:
            row.append(float(data[i][j]))
          except:
            pass
        lst.append(row)
        
      arr = numpy.zeros((len(lst),len(lst[0]),),)
      for i in range(len(lst)):
        for j in range(len(lst[0])):
          arr[i,j] = float(lst[i][j])
      return arr
    return None
    
  @staticmethod
  def write_csv(filename, arr):  
    fh = open(filename, 'w')
    for i in range(len(arr)):
      for j in range(len(arr[i])):
        fh.write(std.float_padded(arr[i,j],8) + " ")
      fh.write("\n")
    fh.close()
  
  @staticmethod
  def option(input):
    input = input.strip().upper()
    if(input[0:1] == "Y"):
      return True
    elif(input[0:2] == "ON"):
      return True
    elif(input[0:1] == "T"):
      return True
    else:
      return False
    
  @staticmethod
  def float_padded(inp, pad=7):
    out = float(inp)
    out = round(out, pad-3)
    out = str(out)  
    while(len(out)<pad):
      out = out + " "      
    return out[0:pad]
    
  @staticmethod
  def write_file_line(fh, title, title_pad, fields, field_pad):
    if(type(fields) == numpy.ndarray):
      t = fields
      fields = []
      for ti in t:
        fields.append(ti)    
    elif(type(fields) != list):
      fields = [fields]
    
    line = str(title)
    while(len(line)<title_pad):
      line = line + ' '
    for f in fields:
      f_str = str(f)
      while(len(f_str)<field_pad):
        f_str = f_str + ' '
      line = line + f_str + ' '
    line = line + '\n'
    fh.write(line)  
  
  @staticmethod
  def print_file_line(title, title_pad, fields, field_pad):
    if(type(fields) == numpy.ndarray):
      t = fields
      fields = []
      for ti in t:
        fields.append(ti)    
    elif(type(fields) != list):
      fields = [fields]
    
    line = str(title)
    while(len(line)<title_pad):
      line = line + ' '
    for f in fields:
      f_str = str(f)
      while(len(f_str)<field_pad):
        f_str = f_str + ' '
      line = line + f_str + ' '
    line = line + '\n'
    print(line,end='')  

###########################################
#  CLASS read_confi
###########################################
class read_config:
  
  @staticmethod
  def read_file(file_path):
  
# Input dictionary
    input = {}
  
# READ DATA
    d = []
    fh = open(file_path, 'r')
    for line in fh:
      line = line.strip()
      if(len(line) > 0 and line[0] != "#"):
        d.append(line)      
    fh.close()
    
# Count commands
    commands = {}
    for line in d:
      fields = read_config.split_by(line, ' ')
      c = fields[0].lower()
      if(c in commands.keys()):
        commands[c] = commands[c] + 1
      else:
        commands[c] = 1
        
# Prepare input dictionary
    for k in commands.keys():
      if(commands[k] == 1):
        input[k] = None
      else:
        input[k] = []
    
# Read Data into input
    for line in d:
      fields = read_config.split_by(line, ' ')
      fkey = fields[0].lower()
      
      fd_size = {}
      for i in range(1, len(fields)):
        f = fields[i]
        fs = f.split("=")
        fc = fs[0].lower()
        if(fc in fd_size.keys()):
          fd_size[fc] = fd_size[fc] + 1
        else:
          fd_size[fc] = 1
          
# Prepare dictionary
      fd = {} 
      for k in fd_size.keys():
        if(fd_size[k] == 1):
          fd[k] = None
        else:
          fd[k] = []        
        
      for i in range(1, len(fields)):
        f = fields[i]
        fs = f.split("=")     
        fc = fs[0].lower()        
        fs = read_config.split_by(fs[1], ',')         
        fs = read_config.store(fs)
        
        if(fd_size[fc] == 1):
          if(len(fs) == 1):
            fd[fc] = read_config.store(fs[0])
          else:
            fd[fc] = read_config.store(fs)
        else:
          if(len(fs) == 1):
            fd[fc].append(read_config.store(fs[0]))
          else:
            fd[fc].append(read_config.store(fs))
            
      if(commands[fkey] == 1):
        input[fkey] = fd
      else:
        input[fkey].append(fd)  

    return input
        
  @staticmethod  
  def split_by(line, sep=' ', ignore_double_sep=True):
    last_char = None
    in_quotes = 0
    fields = []
    temp_line = ""
    
    for char in line:
      if(char == "'" and in_quotes == 0 and last_char != "\\"):
        in_quotes = 1
      elif(char == "'" and in_quotes == 1 and last_char != "\\"):
        in_quotes = 0
      elif(char == '"' and in_quotes == 0 and last_char != "\\"):
        in_quotes = 2
      elif(char == '"' and in_quotes == 2 and last_char != "\\"):
        in_quotes = 0
      elif(in_quotes > 0):
        temp_line = temp_line + char
      elif(in_quotes == 0 and char != sep):
        temp_line = temp_line + char
      elif(char == sep and last_char == sep and ignore_double_sep):
        pass
      elif(char == sep):
        fields.append(temp_line)
        temp_line = "" 
    if(temp_line != ""):
      fields.append(temp_line)
    
    return fields
    
  @staticmethod
  def store(inp):  
    if(isinstance(inp, list)):
      for i in range(len(inp)):
        try:
          if('.' in inp[i]  or 'e' in inp[i]):
            inp[i] = float(inp[i])
          else:
            inp[i] = int(inp[i])
        except:
          pass
    else:
      try:
        if('.' in inp or 'e' in inp):
          inp = float(inp)
        else:
          inp = int(inp)
      except:
        pass
    return inp
      
###########################################
#  CLASS eamp
###########################################
class eampa:
 
  def run():
    print("RUNNING")
    
# Read Type
    eampa.run_type()
    
# Load potentials
    potential.load()
    
# Load configs
    configs.load()
    
# Energy Adjustments (dft)
    e_adjust.load()
    
# Bulk Properties
    b_props.load()
    
# Convert to ev/ang etc and adjust energies
    configs.complete()

    labels.output()    
    configs.output()
    print(g.run_type)
    if(g.run_type == 'e'):
      efs_calc.run_energy()
    elif(g.run_type == 'ef'):
      efs_calc.run_energy_force()
    elif(g.run_type == 'efs'):
      efs_calc.run_energy_force_stress()
    elif(g.run_type == 'bp'):
      bp_calc.run()
    elif(g.run_type == 'es'):
      es_calc.run()
    elif(g.run_type == 'rss'):
      rss_calc.run()
    elif(g.run_type == 'fit'): 
      pf.run()
    elif(g.run_type == 'plot'): 
      potential.run()
      
  def run_type():
# Types
# config       just calculate energy/forces/stress of configs
# bp           calculate bulk properties (bulk modulus, elastic constants etc)
# rss
# fit
  
# Default
    g.run_type = 'bp'
    
    try:
      if(g.inp['run']['type'].lower() == 'e'):
        g.run_type = 'e'
      elif(g.inp['run']['type'].lower() == 'ef'):
        g.run_type = 'ef'
      elif(g.inp['run']['type'].lower() == 'efs'):
        g.run_type = 'efs'
      elif(g.inp['run']['type'].lower() == 'bp'):
        g.run_type = 'bp'
      elif(g.inp['run']['type'].lower() == 'es'):
        g.run_type = 'es'
      elif(g.inp['run']['type'].lower() == 'rss'):
        g.run_type = 'rss'
      elif(g.inp['run']['type'].lower() == 'fit'):
        g.run_type = 'fit'
      elif(g.inp['run']['type'].lower() == 'plot'):
        g.run_type = 'plot'
    except:
      pass
  
###########################################
#  CLASS label
###########################################
class labels:

  @staticmethod
  def add(label):
    label = label.upper()
    if(label in g.labels.keys()):
      return label, g.labels[label]
    else: 
      g.labels[label] = len(g.labels) + 1
      return label, g.labels[label]
     
  @staticmethod
  def output(): 
    if(g.outputs):     
      fh = open(g.dirs['output'] + '/' + 'labels.dat', 'w')   
      for k in g.labels.keys():
        fh.write(str(k) + '  ' + str(g.labels[k]) + '\n')
      fh.close()  
        
  @staticmethod
  def get(l_id):     
    for k in g.labels.keys():
      if(l_id == g.labels[k]):
        return k
    return None
        
###########################################
#  CLASS potentia
###########################################
class potential:

  def run():

    print("Potential") 
  
    efs.init()                           # Initialise (allocate arrays)
    potential.efs_add_potentials()       # Load potentials
    potential.plot_fortran_potentials()
    potential.plot_python_potentials()
    
  def load():
    pot_file = g.inp['potential']['file']
    if(not os.path.isfile(pot_file)):
      return False
# Read potential index
    potential.read_potential(pot_file)
    potential.load_tabulated()
    potential.make_tabulated_points()
    potential.load_analytic()
    potential.make_analytic_points()
    potential.load_fit_data()
    potential.pf_output()
    potential.make_copies()
       
    return True
  
  @staticmethod
  def pot_function():
    return {      
    'f_on': 1,  
    'a_text': '',
    'b_text': '',
    'a': 0,
    'b': 0,
    'f_type': '',             # PAIR, EMBE, DENS
    'f_type_id': 0,           # 1=PAIR, 2=EMBE, 3=DENS
    'f_group': 1,
    'r_cut': 6.5,
    'file': None,
    'function_type': 0,       # 1 tab, 2 analytic
    'f_points': None,         # READ IN TO PYTHON
    'a_type': '',
    'f': None,
    'a_params': None,
    'a_l': 0.0,
    'a_u': 10.0,
    'zoor': 1,
    'points': numpy.zeros((g.tab_size,g.tab_width,),),         # THESE ARE USED BY FORTRAN
    'points_original': numpy.zeros((g.tab_size,g.tab_width,),),         # THESE ARE USED BY FORTRAN
    'fit_file': None,
    'fit_type': None,         # 1 spline, 2 analytic
    'fit_parameters': None,
    'fit_size': None,
    'fit_mult': None,
    }

# LOAD FROM FILE
  @staticmethod
  def read_potential(file_name):
    g.pot_functions['pot_dir'] = ''
    if('/' in file_name):
      lst = file_name.split('/')
      for i in range(len(lst) - 1):
        if(i > 0):
          g.pot_functions['pot_dir'] += '/'
        g.pot_functions['pot_dir'] += lst[i]
        
    index = std.config_file_to_list(file_name)  
    pot = potential.pot_function()
    for row in index:    
      if(len(row) > 1 and row[0].upper() == "POTNAME"):
        g.pot_functions['pot_name'] = row[1]
      elif(len(row) > 1 and row[0].upper() == "ZBLFILE"):
        g.pot_functions['zbl_file'] = row[1]
      elif(len(row) > 1 and row[0].upper() == "F_ON"):
        value = row[1].upper()
        if(value[0].upper() == "N"):
          pot['f_on'] = 0
        if(value[0].upper() == "F"):
          pot['f_on'] = 0
      elif(len(row) > 1 and row[0].upper() == "LABEL"):
        label_str, label_id = labels.add(row[1])
        pot['a_text'] = label_str
        pot['a'] = label_id
        if(len(row) > 2):
          label_str, label_id = labels.add(row[2])
          pot['b_text'] = label_str
          pot['b'] = label_id
      elif(len(row) > 1 and row[0].upper() == "FILE"):
        pot['file'] = row[1]
      elif(len(row) > 1 and row[0].upper() == "FIT"):
        pot['fit_file'] = row[1]
      elif(len(row) > 1 and row[0].upper() == "F_TYPE"):
        value = row[1].upper()
        if(value[0] == "P"):
          pot['f_type'] = 'PAIR'
          pot['f_type_id'] = 1
        elif(value[0] == "D"):
          pot['f_type'] = 'DENS'
          pot['f_type_id'] = 2
        elif(value[0] == "E"):
          pot['f_type'] = 'EMBE'
          pot['f_type_id'] = 3
        else:
          pot['f_type'] = 'NONE'
      elif(len(row) > 1 and row[0].upper() == "F_GROUP"):
        pot['f_group'] = int(row[1])  
      elif(len(row) > 1 and row[0].upper() == "ZOOR"):
        val = row[1].upper() 
        pot['zoor'] = 0  
        if(val[0] == "T" or val[0] == "1" or val[0] == "Y"):
          pot['zoor'] = 1  
      elif(len(row) > 0 and row[0].upper() == "END"):
        g.pot_functions['functions'].append(pot)
        pot = potential.pot_function()
    
# READ ZBL DATA
    
    read_zbl = False
    for row in index:  
      if(read_zbl == False and row[0].upper() == "ZBLSTART"):  
        read_zbl = True
      elif(read_zbl == True and row[0].upper() == "ZBLEND"):  
        read_zbl = False
      elif(read_zbl):
        label_str, id_1 = labels.add(row[0])
        label_str, id_2 = labels.add(row[1])
        on = True
        if(row[2].upper()[0:1] == "N" or row[2].upper()[0:1] == "F" or row[2].upper()[0:3] == "OFF"):
          on = False
        z1 = float(row[3])
        z2 = float(row[4])
        ra = float(row[5])
        rb = float(row[6])
        
        spline_type = 1
        if(row[7].upper()[0:5] == 'POLY3'):
          spline_type = 1  
        elif(row[7].upper()[0:5] == 'POLY5'):
          spline_type = 2      
        elif(row[7].upper()[0:4] == 'EXP3'):
          spline_type = 3       
        elif(row[7].upper()[0:4] == 'EXP5'):
          spline_type = 4        
        z = {
            'id_1': id_1,
            'id_2': id_2,
            'on': on ,
            'z1': z1 ,
            'z2': z2 ,
            'ra': ra ,
            'rb': rb ,
            'spline_type': spline_type ,
            }
        g.pot_functions['zbl'].append(z)
        
  @staticmethod
  def load_tabulated():      
    for i in range(len(g.pot_functions['functions'])): 
      pf_file = g.pot_functions['pot_dir'] + '/' + g.pot_functions['functions'][i]['file']      
      if(pf_file is not None and os.path.isfile(pf_file) and potential.pf_file_type(pf_file) == 'T'):
        g.pot_functions['functions'][i]['f_points'] = std.read_csv_array(pf_file, ' ')
        g.pot_functions['functions'][i]['function_type'] = 1
        
# Fill in tabulated points using interp.fill
  @staticmethod
  def make_tabulated_points():  
    for i in range(len(g.pot_functions['functions'])):    
      if(g.pot_functions['functions'][i]['function_type'] == 1):
        g.pot_functions['functions'][i]['points'] = interp.fill(g.pot_functions['functions'][i]['f_points'][:,0], g.pot_functions['functions'][i]['f_points'][:,1], g.tab_size, g.tab_width)   

# Spline and vary accordingly
  @staticmethod
  def vary_tabulated_points(fn, yvar=None):  
    if(type(yvar) != numpy.ndarray):
      yvar = numpy.zeros((10,),)
    g.pot_functions['functions'][fn]['points'] = spline.vary(g.pot_functions['functions'][fn]['points'][:,0], 
                                                             g.pot_functions['functions'][fn]['points'][:,1], yvar)

  @staticmethod
  def load_analytic():  
    for i in range(len(g.pot_functions['functions'])): 
      pf_file = g.pot_functions['pot_dir'] + '/' + g.pot_functions['functions'][i]['file']      
      if(pf_file is not None and os.path.isfile(pf_file) and potential.pf_file_type(pf_file) == 'A'):
        g.pot_functions['functions'][i]['function_type'] = 2
        fd = std.config_file_to_list(pf_file)  
        param = []
        for l in fd:
          if(l[0].upper() == '#TYPE'):
            g.pot_functions['functions'][i]['a_type'] = l[1].lower()
          elif(l[0].upper()[0:2] == '#P'):
            param.append([int(l[0][2:]), float(l[1])])
          elif(l[0].upper()[0:2] == '#L'):
            g.pot_functions['functions'][i]['a_l'] = float(l[1])
          elif(l[0].upper()[0:2] == '#U'):
            g.pot_functions['functions'][i]['a_u'] = float(l[1])
        g.pot_functions['functions'][i]['a_params'] = numpy.zeros((len(param),),)

        for p in range(len(param)):
          g.pot_functions['functions'][i]['a_params'][p] = float(param[p][1])
          
# Save function
        g.pot_functions['functions'][i]['f'] = getattr(potential_functions, g.pot_functions['functions'][i]['a_type'])  
        
  @staticmethod
  def make_analytic_points():  
    for fn in range(len(g.pot_functions['functions'])):  
      if(g.pot_functions['functions'][fn]['function_type'] == 2):  
        potential.make_analytic_points_inner(fn)
        
  @staticmethod
  def make_analytic_points_inner(fn):  
# Temp x,y array
    temp = numpy.zeros((g.tab_size,2,),)
    temp[:,0] = numpy.linspace(g.pot_functions['functions'][fn]['a_l'], g.pot_functions['functions'][fn]['a_u'], g.tab_size)
    temp[:,1] = g.pot_functions['functions'][fn]['f'](temp[:,0],  g.pot_functions['functions'][fn]['a_params'])
       
# Interpfill
    g.pot_functions['functions'][fn]['points'] = interp.fill(temp[:,0],temp[:,1], g.tab_size, g.tab_width)
    
  @staticmethod
  def load_fit_data(): 
# READ FIT DATA
    for fn in range(len(g.pot_functions['functions'])): 
      if(g.pot_functions['functions'][fn]['fit_file'] != None):
        params = None
        fit_file = g.pot_functions['pot_dir'] + '/' + g.pot_functions['functions'][fn]['fit_file']
        if(os.path.isfile(fit_file)):
          fh = open(fit_file, 'r')
          for line in fh:
            line = std.one_space(line.strip().upper())
            f = line.split(" ")
            if(line[0:4] == "#FIT"):
              if(f[-1] == 'S'):
                g.pot_functions['functions'][fn]['fit_type'] = 1
              elif(f[-1] == 'A'):
                g.pot_functions['functions'][fn]['fit_type'] = 2                
            if(line[0:3] == "#PL"):
              params_lower = f[1:]                     
            if(line[0:3] == "#PU"):
              params_upper = f[1:]     
            if(line[0:2] == "#M"):
              g.pot_functions['functions'][fn]['fit_mult'] = numpy.zeros((2,),)
              """
              try:
                g.pot_functions['functions'][fn]['fit_mult'][0] = float(f[1])
                g.pot_functions['functions'][fn]['fit_mult'][1] = float(f[2])
              except:  
                g.pot_functions['functions'][fn]['fit_mult'][0] = 0.1
                g.pot_functions['functions'][fn]['fit_mult'][1] = 10.0
              """
                
          fh.close()
        if(g.pot_functions['functions'][fn]['fit_type'] != None and type(params_lower) == list):
        
# Analytic fitting
          if(g.pot_functions['functions'][fn]['fit_type'] == 2 and
             g.pot_functions['functions'][fn]['function_type'] == 2):
            g.pot_functions['functions'][fn]['fit_size'] = len(g.pot_functions['functions'][fn]['a_params'])
            g.pot_functions['functions'][fn]['fit_parameters'] = numpy.zeros((2,len(g.pot_functions['functions'][fn]['a_params']),),)     
            
            for i in range(g.pot_functions['functions'][fn]['fit_size']):
              g.pot_functions['functions'][fn]['fit_parameters'][0,i] = float(params_lower[i])
              g.pot_functions['functions'][fn]['fit_parameters'][1,i] = float(params_upper[i])
              
# Spline fitting
          else:
            s = len(params_lower)
            g.pot_functions['functions'][fn]['fit_type'] = 1
            g.pot_functions['functions'][fn]['fit_size'] = s            
            g.pot_functions['functions'][fn]['fit_parameters'] = numpy.zeros((2,s,),) 
            
            for i in range(len(params_lower)):
              g.pot_functions['functions'][fn]['fit_parameters'][0,i] = float(params_lower[i])
              g.pot_functions['functions'][fn]['fit_parameters'][1,i] = float(params_upper[i])

  @staticmethod
  def make_copies(): 
    for fn in range(len(g.pot_functions['functions'])): 
      g.pot_functions['functions'][fn]['points_original'] = numpy.copy(g.pot_functions['functions'][fn]['points'])
  
    g.pot_functions['functions_original'] = copy.deepcopy(g.pot_functions['functions'])

  @staticmethod
  def spline_prep(): 
    for fn in range(len(g.pot_functions['functions'])): 
      if(g.pot_functions['functions'][fn]['fit_type'] == 2):   # If spline fit
        yvar = numpy.zeros((len(g.pot_functions['functions'][fn]['fit_parameters']),),)
        potential.vary_tabulated_points(fn, yvar)
    
  @staticmethod
  def pf_file_type(file): 
    fh = open(file, 'r')
    for row in fh:
      if(row.strip().upper() == '#A'):
        return 'A'
      return 'T'
      
  @staticmethod
  def pf_output(): 
    if(g.outputs):     
      fh = open(g.dirs['output'] + '/' + 'pot.dat', 'w')
      
      n = 0
      for pf in g.pot_functions['functions']:
        n = n + 1
        fh.write('#############################################################\n')
        fh.write('#############################################################\n')
        fh.write('POTENTIAL ' + str(n) + '\n')
        fh.write('#############################################################\n')
        fh.write('#############################################################\n')
        fh.write('f_on: ' + str(pf['f_on']) + '\n')
        fh.write('a_text: ' + str(pf['a_text']) + '\n')
        fh.write('b_text: ' + str(pf['b_text']) + '\n')
        fh.write('a: ' + str(pf['a']) + '\n')
        fh.write('b: ' + str(pf['b']) + '\n')
        fh.write('f_type: ' + str(pf['f_type']) + '\n')
        fh.write('f_group: ' + str(pf['f_group']) + '\n')
        fh.write('r_cut: ' + str(pf['r_cut']) + '\n')
        fh.write('file: ' + str(pf['file']) + '\n')
        fh.write('function_type: ' + str(pf['function_type']) + '\n')
        fh.write('a_type: ' + str(pf['a_type']) + '\n')
        fh.write('f: ' + str(pf['f']) + '\n')
        fh.write('a_l: ' + str(pf['a_l']) + '\n')
        fh.write('a_u: ' + str(pf['a_u']) + '\n')
        fh.write('zoor: ' + str(pf['zoor']) + '\n')
        fh.write('a_params: \n')
        try:
          for i in range(len(pf['a_params'])):
            fh.write(str(pf['a_params'][i]) + '\n')    
        except:
          fh.write(str('None\n')) 
          
        fh.write('f_points: \n')
        try:
          for i in range(len(pf['f_points'])):
            for j in range(len(pf['f_points'][i])):
              fh.write(str(pf['f_points'][i, j]) + ' ')  
            fh.write('\n')  
        except:
          fh.write(str('None\n')) 
          
        fh.write('points: \n')
        try:
          for i in range(len(pf['points'])):
            for j in range(len(pf['points'][i])):
              fh.write(str(pf['points'][i, j]) + ' ')  
            fh.write('\n') 
        except:
          fh.write(str('None\n'))         
      fh.close()
      
  """
  def pot_function():
    return {      
    'f_on': 1,  
    'a_text': '',
    'b_text': '',
    'a': 0,
    'b': 0,
    'f_type': '',             # PAIR, EMBE, DENS
    'f_type_id': 0,           # 1=PAIR, 2=EMBE, 3=DENS
    'f_group': 1,
    'r_cut': 6.5,
    'file': None,
    'function_type': 0,       # 1 tab, 2 analytic
    'f_points': None,         # READ IN TO PYTHON
    'a_type': '',
    'f': None,
    'a_params': None,
    'a_l': 0.0,
    'a_u': 10.0,
    'zoor': 1,
    'points': numpy.zeros((g.tab_size,g.tab_width,),),         # THESE ARE USED BY FORTRAN
    }
  """    
      
  def plot_python_potentials(dir = None):  
    if(dir == None):
      dir = g.dirs['pots']
      
    for i in range(len(g.pot_functions['functions'])): 
      pot_name = 'py_pot'
      pot_count = i + 1
      
      pot_type = g.pot_functions['functions'][i]['f_type_id']
      label_a = g.pot_functions['functions'][i]['a']
      
      if(pot_type == 1):
        label_b = g.pot_functions['functions'][i]['b']
      else:
        label_b = g.pot_functions['functions'][i]['f_group']
        
      potential.plot_potential(dir, pot_name, pot_count, pot_type, label_a, label_b, 
                               g.pot_functions['functions'][i]['points'][:,0], 
                               g.pot_functions['functions'][i]['points'][:,1], 
                               g.pot_functions['functions'][i]['points'][:,2], 
                               g.pot_functions['functions'][i]['points'][:,3])
      
  def plot_fortran_potentials(dir = None):  
    
    if(dir == None):
      dir = g.dirs['pots']
    
    if(efs.pc > 0):
      for i in range(efs.pc):
    
        a = efs.pkey[i,0] - 1
        b = efs.pkey[i,1]
      
        pot_name = 'fort_pot'
        pot_count = i + 1
       
        pot_type = efs.pkey[i,2]
        label_a = efs.pkey[i,3]
        label_b = efs.pkey[i,4]
      
        potential.plot_potential(dir, pot_name, pot_count, pot_type, label_a, label_b, 
                                 efs.pot[a:b,0], efs.pot[a:b,1], 
                                 efs.pot[a:b,2], efs.pot[a:b,3])
    
    elif(bp.pc > 0):
      for i in range(bp.pc):
    
        a = bp.pkey[i,0] - 1
        b = bp.pkey[i,1]
      
        pot_name = 'fort_pot'
        pot_count = i + 1
       
        pot_type = bp.pkey[i,2]
        label_a = bp.pkey[i,3]
        label_b = bp.pkey[i,4]
      
        potential.plot_potential(dir, pot_name, pot_count, pot_type, label_a, label_b, 
                                 bp.pot[a:b,0], bp.pot[a:b,1], 
                                 bp.pot[a:b,2], bp.pot[a:b,3])
                                 
  def plot_potential(dir, pot_name, pot_count, pot_type, label_a, label_b, x, y, yp, ypp):  
      
      if(pot_type == 1):
        plot_title = "Pair - " + labels.get(label_a) + ' ' + labels.get(label_b)    
        x_axis = "Seperation (Ang)"
        y_axis = "Potential Energy (eV)"
      elif(pot_type == 2):   
        plot_title = "Density - " + labels.get(label_a) + ' Group ' + str(label_b)    
        x_axis = "Seperation (Ang)"
        y_axis = "Electron Density"   
      elif(pot_type == 3):
        plot_title = "Embedding - " + labels.get(label_a) + ' Group ' + str(label_b)     
        x_axis = "Electron Density"
        y_axis = "Potential Energy (eV)"
      
      file_name = str(pot_count)
      while(len(file_name)<3):
        file_name = '0' + file_name
      file_name = pot_name + '_' + file_name + '.eps'
      
      plt.clf()    
      plt.rc('font', family='serif')
      plt.rc('xtick', labelsize='x-small')
      plt.rc('ytick', labelsize='x-small')
      fig, axs = plt.subplots(1, 1, figsize=(12,9))
      fig.tight_layout(pad=5.0)   
      
      fig.suptitle(plot_title)
      axs.plot(x, y, color='k', ls='solid')
      axs.plot(x, yp, color='k', ls='dashed')
      axs.plot(x, ypp, color='k', ls='dotted')
      axs.set_title('')
      axs.set_xlabel(x_axis)
      axs.set_ylabel(y_axis)
      
      min_y = min(y)
      max_y = max(y)
      r = max_y - min_y
      
      min_y = min_y - 0.05 * r
      max_y = max_y + 0.05 * r
      
      if(min_y < -1.0E03):
        min_y = -1.0E03
      if(max_y > 1.0E04):
        max_y = 1.0E04
         
      axs.set_ylim(min_y, max_y)
      axs.set_yscale('symlog', linthreshy=10)
      plt.savefig(dir + '/' + file_name, format='eps')
      
###########################################################
# F2PY functions
###########################################################

  def efs_add_potentials():
  
# Clear
    efs.clear_potentials()
    
# ADD POTENTIALS
    for pf in g.pot_functions['functions']:
      if(pf['f_on']):
        if(pf['f_type_id'] == 1):
          efs.add_potential(
                            pf['f_type_id'],
                            pf['a'], 
                            pf['b'],  
                            pf['r_cut'], 
                            pf['points'] 
                           )
        elif(pf['f_type_id'] > 1):
          efs.add_potential(
                            pf['f_type_id'],
                            pf['a'], 
                            pf['f_group'],  
                            pf['r_cut'], 
                            pf['points'] 
                           )    
# ADD ZBL
    for zbl in g.pot_functions['zbl']:
      efs.add_zbl(
                  zbl['id_1'],
                  zbl['id_2'], 
                  zbl['on'],  
                  zbl['z1'], 
                  zbl['z2'] , 
                  zbl['ra'] , 
                  zbl['rb'] , 
                  zbl['spline_type'] 
                 )    

# SET POTENTIALS
    efs.set_potentials()   
    
  def es_add_potentials():
  
# Clear
    es.clear_potentials()
    
# ADD POTENTIALS
    for pf in g.pot_functions['functions']:
      if(pf['f_on']):
        if(pf['f_type_id'] == 1):
          es.add_potential(
                            pf['f_type_id'],
                            pf['a'], 
                            pf['b'],  
                            pf['r_cut'], 
                            pf['points'] 
                           )
        elif(pf['f_type_id'] > 1):
          es.add_potential(
                            pf['f_type_id'],
                            pf['a'], 
                            pf['f_group'],  
                            pf['r_cut'], 
                            pf['points'] 
                           )    
# ADD ZBL
    for zbl in g.pot_functions['zbl']:
      es.add_zbl(
                  zbl['id_1'],
                  zbl['id_2'], 
                  zbl['on'],  
                  zbl['z1'], 
                  zbl['z2'] , 
                  zbl['ra'] , 
                  zbl['rb'] , 
                  zbl['spline_type'] 
                 )    

# SET POTENTIALS
    es.set_potentials()     

  def bp_add_potentials():
  
# Clear
    bp.clear_potentials()
    
# ADD POTENTIALS
    for pf in g.pot_functions['functions']:
      if(pf['f_on']):
        if(pf['f_type_id'] == 1):
          bp.add_potential(
                            pf['f_type_id'],
                            pf['a'], 
                            pf['b'],  
                            pf['r_cut'], 
                            pf['points'] 
                           )
        elif(pf['f_type_id'] > 1):
          bp.add_potential(
                            pf['f_type_id'],
                            pf['a'], 
                            pf['f_group'],  
                            pf['r_cut'], 
                            pf['points'] 
                           )    
      
# ADD ZBL
    for zbl in g.pot_functions['zbl']:
      bp.add_zbl(
                  zbl['id_1'],
                  zbl['id_2'], 
                  zbl['on'],  
                  zbl['z1'], 
                  zbl['z2'] , 
                  zbl['ra'] , 
                  zbl['rb'] , 
                  zbl['spline_type'] 
                 )    
                                               
# SET POTENTIALS
    bp.set_potentials()
      
###########################################
#  CLASS potential_function
###########################################
class potential_functions:

##############################################
# PAIR FUNCTIONS
##############################################

# Lennard Jones Potential
# p[0] = e
# p[1] = rm
# f(x) = A * ((B / r)**12 - 2 * (B/r)**6)
  @staticmethod
  def lennard_jones(r, p):
    return fnc.lennard_jones_v(r, p)
    
# Morse Potential
# p[0] = d
# p[1] = a
# p[2] = re
# f(x) = A * (exp(-2.0D0 * B * (r - C)) - 2.0D0 * exp(-B*(r - C)))
  @staticmethod
  def morse(r, p):
    return fnc.morse_v(r, p)
#return

# Buckingham Potential
# p[0] = A
# p[1] = B
# p[2] = C
# f(x) = A * exp(-1 * B * r) - C / r**6
  @staticmethod
  def buckingham(r, p):
    return fnc.buckingham_v(r, p)

##############################################
# DENSITY FUNCTIONS
##############################################

# Embedding Finnis-Sinclair
  @staticmethod
  def quadratic_density(r, p):
    return (r - p[0])**2
    
# Density Finnis-Sinclair
  @staticmethod
  def density_fs(r, p):
    return (r - p[0])**2

# Embedding Finnis-Sinclair
# f(x) = -A sqrt(rho)
  @staticmethod
  def fs_embedding(r, p):
    return fnc.fs_embedding(r, p)

# Embedding Mendelev
# f(x) = -sqrt(rho) + A*rho**2
  @staticmethod
  def mendelev_embedding_v(r, p):
    return fnc.mendelev_embedding_v(r, p)

# Triple Embedding
# f(x) = A * sqrt(r) + B * r + C * r**2
  @staticmethod
  def triple_embedding(r, p):
    return fnc.triple_embedding_v(r, p)

# Embedding Mendelev, Han, Srolovitz, Ackland, Sun, Asta
  @staticmethod
  def embedding_1(r, p):
    return -1 * numpy.sqrt(r) + p[0] * r**2
    
##############################################
# SIMPLE SPLINES
##############################################

  @staticmethod
  def summed_spline(r, p):
    return fnc.summed_spline_v(r, p)
    
  @staticmethod
  def simple_spline(r, p):
    return fnc.simple_spline_v(r, p)
    
##############################################
# NODE SPLINES
##############################################
   
  @staticmethod
  def spline_n_node(r, p):
    return fnc.spline_n_node_v(r, p)
    
###########################################
#  CLASS rescale_densit
###########################################
class rescale_density:
  
  def run():
    
    for fn in range(len(g.pot_functions['functions'])): 
      if(g.pot_functions['functions'][fn]['f_type_id'] == 2):
        min_rho = min(g.pot_functions['functions'][fn]['points'][:,1])
        g.pot_functions['functions'][fn]['points'][:,1] = g.pot_functions['functions'][fn]['points'][:,1] - min_rho 
        rho = rescale_density.estimate_density(fn)        
        g.pot_functions['functions'][fn]['points'][:,1:3] = (0.5 / rho) * g.pot_functions['functions'][fn]['points'][:,1:3]
        rho = rescale_density.estimate_density(fn)

  def estimate_density(fn):
    r = numpy.zeros((7,),)
    rn = numpy.zeros((7,),)
    r[0] = 7.48332e0
    r[1] = 6.32456e0
    r[2] = 6.92821e0
    r[3] = 4.89898e0
    r[4] = 5.65686e0
    r[5] = 2.82843e0
    r[6] = 4.0e0
    rn[0] = 48
    rn[1] = 24
    rn[2] = 8
    rn[3] = 24
    rn[4] = 12
    rn[5] = 12
    rn[6] = 6    
    rho = 0.0    
    for i in range(7):
      y = interp.search_x(r[i], g.pot_functions['functions'][fn]['points'][:,0], g.pot_functions['functions'][fn]['points'][:,1])
      rho = rho + rn[i] * y
    return rho
    
  def max_densities():    
    density_list = None
    for fn in range(len(g.pot_functions['functions'])): 
      if(g.pot_functions['functions'][fn]['f_type_id'] == 2):
        rho = rescale_density.estimate_density(fn) 
        if(density_list is None):
          density_list = str(rho)
        else:
          density_list = density_list + "/" + str(rho)
    return density_list
        
###########################################
#  CLASS potential_var
###########################################
class potential_vary:
 
  def vary_all():
    for pn in range(len(g.pot_functions['functions'])):
      if(g.pot_functions['functions'][pn]['function_type'] == 1):
        potential_vary.vary_tabulated(pn)
      elif(g.pot_functions['functions'][pn]['function_type'] == 2):
        potential_vary.vary_analytic(pn)
  
  def vary_tabulated(pn):
    pass
  
  def vary_analytic(pn):
    print(g.pot_functions['functions'][pn])
    
###########################################
#  CLASS config
###########################################
class configs:

  @staticmethod
  def load():
    try:
      configs_dir = g.inp['configs']['dir']
    except:  
      return False
    if(not os.path.isdir(configs_dir)):
      return False      
    
# Make file list
    g.configs['config_files'] = []
    g.configs['config_files'] = configs.load_files(configs_dir, g.configs['config_files'])
    
# Read files in
    configs.read()
    
#print(g.configs)
    
    return True

  @staticmethod
  def load_files(path_in, files):
    for path in os.listdir(path_in):
      path_new = path_in + "/" + path
      if(os.path.isdir(path_new)):
        files = configs.load_files(path_new, files)
      else:
        type = configs.file_type(path_new)
        files.append([path_new,type])
    return files
    
  @staticmethod
  def read():
# Loop through files and read in
    for file in g.configs['config_files']:      
      if(file[1] == 'std'):        
        configs.std(file[0])
      elif(file[1] == 'qe'):
        configs.qe(file[0])
        
  @staticmethod
  def std(file_path):   
  
# Read content from file
    content = std.file_to_list(file_path)
    content = std.prep_data(content)
    content = std.remove_comments(content)
  
# Split up into individual configs
    config_list = []    
    temp = []
    last_line = None
    for line in content:
      if(last_line != None and line != '' and last_line != '' and last_line[0] != "#" and line[0] == "#"):
        config_list.append(temp)
        temp = []
      if(line != ''):  
        temp.append(line)      
        last_line = line
    if(len(temp)>0):
      config_list.append(temp)
      
# Read each config
    for i in range(len(config_list)):
      configs.add_config(config_list[i], file_path, i)
  
  @staticmethod
  def make_config():
    return {
    'file_type': '',  
    'file_path': '',  
    'file_part': 0,  
    'alat': 0.0,  
    'uv_prim': numpy.zeros((3,3,),),
    'uv': numpy.zeros((3,3,),),
    'stress': numpy.zeros((3,3,),),
    'c': numpy.zeros((3,),dtype=numpy.int32,),
    'rcut': 0.0,
    'rverlet': 0.0,
    'mtemp': {},
    'm': {},
    'coord_count_prim': 0,
    'coords_label_prim': None,
    'coords_label_id_prim': None,
    'coords_prim': None,
    'forces_prim': None,
    'coord_count': 0,
    'coords_label': None,
    'coords_label_id': None,
    'coords': None,
    'forces': None,
    'energy_per_atom': None,
    'energy': None,
    'e': 0,
    'f': 0,
    's': 0,
    'n_atoms_prim': 0,
    'n_atoms': 0,
    'l_units': 'ANG',
    'e_units': 'EV',
    'f_units': 'EV/ANG',
    's_units': 'EV/ANG3',
    }    
  
  @staticmethod
  def add_config(content, file_path, i, config_type='STANDARD'):
    fd = configs.make_config()
    fd['file_type'] = config_type
    fd['file_path'] = file_path
    fd['file_part'] = i
    
# FIRST READ
    for line in content:
      line = line.strip() 
      f = std.to_fields(line, ' ')
      epa_set = False

      if(len(f) > 1 and f[0].upper() == "#ALAT"):
        fd['alat'] = float(f[1])
      if(len(f) > 3 and f[0].upper() == "#X"):
        fd['uv_prim'][0,0] = float(f[1])
        fd['uv_prim'][0,1] = float(f[2])
        fd['uv_prim'][0,2] = float(f[3])
      if(len(f) > 3 and f[0].upper() == "#Y"):
        fd['uv_prim'][1,0] = float(f[1])
        fd['uv_prim'][1,1] = float(f[2])
        fd['uv_prim'][1,2] = float(f[3])
      if(len(f) > 3 and f[0].upper() == "#Z"):
        fd['uv_prim'][2,0] = float(f[1])
        fd['uv_prim'][2,1] = float(f[2])
        fd['uv_prim'][2,2] = float(f[3])
      if(len(f) > 3 and f[0].upper() == "#SX"):
        fd['s'] = 1
        fd['stress'][0,0] = float(f[1])
        fd['stress'][0,1] = float(f[2])
        fd['stress'][0,2] = float(f[3])
      if(len(f) > 3 and f[0].upper() == "#SY"):
        fd['s'] = 1
        fd['stress'][1,0] = float(f[1])
        fd['stress'][1,1] = float(f[2])
        fd['stress'][1,2] = float(f[3])
      if(len(f) > 3 and f[0].upper() == "#SZ"):
        fd['s'] = 1
        fd['stress'][2,0] = float(f[1])
        fd['stress'][2,1] = float(f[2])
        fd['stress'][2,2] = float(f[3])
      if(len(f) > 2 and f[0].upper() == "#M"):
        fd['mtemp'][f[1]] = f[2]      
      if(len(f) > 3 and f[0].upper() == "#C"):
        fd['c'][0] = int(f[1])
        fd['c'][1] = int(f[2])
        fd['c'][2] = int(f[3])
      if(len(f) > 1 and f[0].upper() == "#E"):
        fd['energy'] = float(f[1])
        fd['e'] = 1    
      if(len(f) > 1 and f[0].upper() == "#EPA"):
        fd['energy_per_atom'] = float(f[1])
        fd['e'] = 1   
      if(len(f) > 1 and f[0].upper() == "#RCUT"):
        fd['rcut'] = f[1]
        if(fd['rverlet'] == 0):
          fd['rverlet'] = f[1]        
      if(len(f) > 1 and f[0].upper() == "#RVERLET"):
        fd['rverlet'] = f[1]      
      if(len(f) > 1 and f[0].upper() == "#L_UNITS"):
        fd['l_units'] = f[1]  
      if(len(f) > 1 and f[0].upper() == "#E_UNITS"):
        fd['e_units'] = f[1]  
      if(len(f) > 1 and f[0].upper() == "#F_UNITS"):
        fd['f_units'] = f[1]
      if(len(f) > 1 and f[0].upper() == "#S_UNITS"):
        fd['s_units'] = f[1]
      if(len(f) >= 4 and f[0][0] != "#"):    
        count_coord = False
        if(len(f) >= 4):
          count_coord = True
        if(len(f) >= 7):
          fd['f'] = 1
        if(count_coord):
          fd['coord_count_prim'] = fd['coord_count_prim'] + 1
          
# COORD SIZE
    c = numpy.identity(3)
    
# EXPAND
    c[0,0] = float(fd['c'][0]) 
    c[1,1] = float(fd['c'][1]) 
    c[2,2] = float(fd['c'][2]) 
    fd['coord_count'] = (fd['c'][0] * fd['c'][1] * fd['c'][2]) * fd['coord_count_prim']
    fd['uv'] = numpy.matmul(c, fd['uv_prim'])
    
# MAKE ARRAYS
    fd['coords_label_prim'] = [None] * fd['coord_count_prim']
    fd['coords_label_id_prim'] = numpy.zeros((fd['coord_count_prim'],), dtype=numpy.int32,)
    fd['coords_prim'] = numpy.zeros((fd['coord_count_prim'],3,),)
    if(fd['f'] == 1):
      fd['forces_prim'] = numpy.zeros((fd['coord_count_prim'],3,),)
    
    fd['coords_label'] = [None] * fd['coord_count']
    fd['coords_label_id'] = numpy.zeros((fd['coord_count'],), dtype=numpy.int32,)
    fd['coords'] = numpy.zeros((fd['coord_count'],3,),)
    if(fd['f'] == 1):
      fd['forces'] = numpy.zeros((fd['coord_count'],3,),)
        
    fd['n_atoms_prim'] = fd['coord_count_prim']
    fd['n_atoms'] = fd['coord_count']   
    
# READ COORDS
    n = 0
    for line in content:
      line = line.strip() 
      f = std.to_fields(line, ' ')
      if(len(f) >= 4 and f[0][0] != "#"):    
        if(len(f) >= 4):
          label_str, label_id = labels.add(f[0])
          fd['coords_label_prim'][n] = label_str
          fd['coords_label_id_prim'][n] = label_id
          fd['coords_prim'][n,0] = float(f[1])
          fd['coords_prim'][n,1] = float(f[2])
          fd['coords_prim'][n,2] = float(f[3])
          
        if(len(f) >= 7):
          fd['forces_prim'][n,0] = float(f[4])
          fd['forces_prim'][n,1] = float(f[5])
          fd['forces_prim'][n,2] = float(f[6])
        n = n + 1
        
# EXPAND COORDS
    m = 0
    for i in range(fd['c'][0]):
      for j in range(fd['c'][1]):
        for k in range(fd['c'][2]):
          for n in range(fd['coord_count_prim']):
            fd['coords_label'][m] = fd['coords_label_prim'][n]
            fd['coords_label_id'][m] = fd['coords_label_id_prim'][n]
            fd['coords'][m,0] = (i + fd['coords_prim'][n,0]) / fd['c'][0]
            fd['coords'][m,1] = (j + fd['coords_prim'][n,1]) / fd['c'][1]
            fd['coords'][m,2] = (k + fd['coords_prim'][n,2]) / fd['c'][2]
            if(fd['f'] == 1):
              fd['forces'][m,0] = fd['forces_prim'][n,0]
              fd['forces'][m,1] = fd['forces_prim'][n,1]
              fd['forces'][m,2] = fd['forces_prim'][n,2]
            m = m + 1
            
# SORT ENERGY
    if(fd['e'] == 1 and fd['energy_per_atom'] != None and fd['energy'] == None):
      try:
        fd['energy'] = fd['energy_per_atom'] * fd['n_atoms']
      except:
        pass            
    if(fd['e'] == 1 and fd['energy'] != None and fd['energy_per_atom'] == None):
      try:
        fd['energy_per_atom'] = fd['energy'] / fd['n_atoms_prim']
        fd['energy'] = fd['energy_per_atom'] * fd['n_atoms']
      except:
        pass

# APPEND
    g.configs['configs'].append(fd)

#
#  FILE TYPE (standard, quantum espresso)
#
        
  @staticmethod
  def file_type(file_path):
    content = std.file_to_list(file_path)
    
# Check if standard file
    count = 0
    for line in content:
      fields = line.split(" ")
      if(fields[0].upper() == "#ALAT"):
        count = count + 1
      if(fields[0].upper() == "#X"):
        count = count + 1
      if(fields[0].upper() == "#Y"):
        count = count + 1
      if(fields[0].upper() == "#Z"):
        count = count + 1
    if(count >= 4):
      return 'std'
  
# Check if pwscf/qe file
    count = 0
    for line in content:
      if(line.strip()[0:13] == "Program PWSCF"):
        count = count + 1
      if(line.strip()[0:27] == "bravais-lattice index     ="):
        count = count + 1
      if(line.strip()[0:27] == "kinetic-energy cutoff     ="):
        count = count + 1
      if(line.strip()[0:27] == "mixing beta               ="):
        count = count + 1
      if(line.strip()[0:27] == "Exchange-correlation      ="):
        count = count + 1
      if(line.strip()[0:9] == "JOB DONE."):
        count = count + 1
    if(count >= 4):
      return 'qe'   
      
# Espresso files

  @staticmethod
  def qe(file_path):   
    
    qe = pwscf_output(file_path)
    xyz = qe.make_xyz()    
    n = 0
    for xyz_inner in xyz:
      n = n + 1
      configs.add_config(xyz_inner, file_path, n, config_type='QE')

  @staticmethod
  def output(): 
    if(g.outputs):     
      fh = open(g.dirs['output'] + '/' + 'configs.dat', 'w')
      
      n = 0
      for c in g.configs['configs']:
              
        n = n + 1
        fh.write('#############################################################\n')
        fh.write('#############################################################\n')
        fh.write('CONFIG ' + str(n) + '\n')
        fh.write('#############################################################\n')
        fh.write('#############################################################\n')
        fh.write('file_type  ' + str(c['file_type']) + ' \n')
        fh.write('file_path  ' + str(c['file_path']) + ' \n')
        fh.write('file_part  ' + str(c['file_part']) + ' \n')
        fh.write('alat  ' + str(c['alat']) + ' \n')
        fh.write('uv_prim   ' + str(c['uv_prim'][0,0]) + ' ' + str(c['uv_prim'][0,1]) + ' ' + str(c['uv_prim'][0,2]) + ' ' + ' \n')
        fh.write('          ' + str(c['uv_prim'][1,0]) + ' ' + str(c['uv_prim'][1,1]) + ' ' + str(c['uv_prim'][1,2]) + ' ' + ' \n')
        fh.write('          ' + str(c['uv_prim'][2,0]) + ' ' + str(c['uv_prim'][2,1]) + ' ' + str(c['uv_prim'][2,2]) + ' ' + ' \n')
        fh.write('uv        ' + str(c['uv'][0,0]) + ' ' + str(c['uv'][0,1]) + ' ' + str(c['uv'][0,2]) + ' ' + ' \n')
        fh.write('          ' + str(c['uv'][1,0]) + ' ' + str(c['uv'][1,1]) + ' ' + str(c['uv'][1,2]) + ' ' + ' \n')
        fh.write('          ' + str(c['uv'][2,0]) + ' ' + str(c['uv'][2,1]) + ' ' + str(c['uv'][2,2]) + ' ' + ' \n')
        fh.write('stress    ' + str(c['stress'][0,0]) + ' ' + str(c['stress'][0,1]) + ' ' + str(c['stress'][0,2]) + ' ' + ' \n')
        fh.write('          ' + str(c['stress'][1,0]) + ' ' + str(c['stress'][1,1]) + ' ' + str(c['stress'][1,2]) + ' ' + ' \n')
        fh.write('          ' + str(c['stress'][2,0]) + ' ' + str(c['stress'][2,1]) + ' ' + str(c['stress'][2,2]) + ' ' + ' \n')
        fh.write('c         ' + str(c['c'][0])  + ' ' + str(c['c'][1]) + ' ' + str(c['c'][2]) + ' ' + ' \n')
        fh.write('rcut      ' + str(c['rcut']) + ' \n')
        fh.write('rverlet   ' + str(c['rverlet']) + ' \n')
        fh.write('e         ' + str(c['e']) + ' \n')
        fh.write('f         ' + str(c['f']) + ' \n')
        fh.write('s         ' + str(c['s']) + ' \n')
        
        fh.write('\n')        
        fh.write('COORDS PRIM (' + str(c['coord_count_prim']) + ')\n')
        for k in range(c['coord_count_prim']):        
          fh.write(str(c['coords_label_prim'][k]) + '  ')  
          fh.write(str(c['coords_prim'][k,0]) + '  ')  
          fh.write(str(c['coords_prim'][k,1]) + '  ')  
          fh.write(str(c['coords_prim'][k,2]) + '  ')             
          if(c['f'] == 1):     
            fh.write(str(c['forces_prim'][k,0]) + '  ')  
            fh.write(str(c['forces_prim'][k,1]) + '  ')  
            fh.write(str(c['forces_prim'][k,2]) + '  ')           
          fh.write('\n')  
        fh.write('\n')    
        
        fh.write('\n')        
        fh.write('COORDS (' + str(c['coord_count']) + ')\n')
        for k in range(c['coord_count']):        
          fh.write(str(c['coords_label'][k]) + '  ')  
          fh.write(str(c['coords'][k,0]) + '  ')  
          fh.write(str(c['coords'][k,1]) + '  ')  
          fh.write(str(c['coords'][k,2]) + '  ')             
          if(c['f'] == 1):     
            fh.write(str(c['forces'][k,0]) + '  ')  
            fh.write(str(c['forces'][k,1]) + '  ')  
            fh.write(str(c['forces'][k,2]) + '  ')           
          fh.write('\n')  
        fh.write('\n')    
             
        fh.write('\n')               
        fh.write('\n')  
         
    fh.close()

# Run after configs read in
# After dft_energy_adjustments
    
  @staticmethod
  def complete(): 
#print("Complete Configs")
#print(g.configs['configs'])
    
# CONVERT
    for i in range(len(g.configs['configs'])):
    
      g.configs['configs'][i]['alat'] = units.convert(g.configs['configs'][i]['l_units'], 'ang', g.configs['configs'][i]['alat'])
      g.configs['configs'][i]['l_units'] = 'ANG'
      
#print(g.configs['configs'][i]['energy'], g.configs['configs'][i]['energy_per_atom'])
      
      g.configs['configs'][i]['energy_per_atom'] = units.convert(g.configs['configs'][i]['e_units'], 'ev', g.configs['configs'][i]['energy_per_atom'])
      g.configs['configs'][i]['energy'] = units.convert(g.configs['configs'][i]['e_units'], 'ev', g.configs['configs'][i]['energy'])
      g.configs['configs'][i]['e_units'] = 'EV'
#print(g.configs['configs'][i]['energy'], g.configs['configs'][i]['energy_per_atom'])
      
# print(g.configs['configs'][i]['f_units'])
      if(g.configs['configs'][i]['f'] == 1):
        for na in range(len(g.configs['configs'][i]['forces_prim'])):
          for nb in range(len(g.configs['configs'][i]['forces_prim'][na])):
            g.configs['configs'][i]['forces_prim'][na,nb] = units.convert(g.configs['configs'][i]['f_units'], 'EV/ANG', g.configs['configs'][i]['forces_prim'][na,nb])
      
      if(g.configs['configs'][i]['f'] == 1):
        for na in range(len(g.configs['configs'][i]['forces'])):
          for nb in range(len(g.configs['configs'][i]['forces'][na])):
            g.configs['configs'][i]['forces'][na,nb] = units.convert(g.configs['configs'][i]['f_units'], 'EV/ANG', g.configs['configs'][i]['forces'][na,nb])
     
      if(g.configs['configs'][i]['f'] == 1):
        for na in range(len(g.configs['configs'][i]['stress'])):
          for nb in range(len(g.configs['configs'][i]['stress'][na])):
            g.configs['configs'][i]['stress'][na,nb] = units.convert(g.configs['configs'][i]['s_units'], 'EV/ANG3', g.configs['configs'][i]['stress'][na,nb])
    
# ADJUST ENERGY QE
    for i in range(len(g.configs['configs'])):
      if(g.configs['configs'][i]['file_type'] == 'QE'):
#print(g.configs['configs'][i]['coord_count'], g.configs['configs'][i]['coord_count_prim'])
        for k in range(g.configs['configs'][i]['coord_count']):      
          e_adj = g.dft_energy_adjustments[g.configs['configs'][i]['coords_label_id'][k]]          
          g.configs['configs'][i]['energy'] = g.configs['configs'][i]['energy'] + e_adj['calc_apaev']
        g.configs['configs'][i]['energy_per_atom'] = g.configs['configs'][i]['energy'] / g.configs['configs'][i]['coord_count'] 
#print(g.configs['configs'][i]['energy'], g.configs['configs'][i]['energy_per_atom'])
          
#print(e_adj['calc_apaev'])

#c['coords_label'][k]
#print('QE')
        
#  g.dft_energy_adjustments[label_id]
  
###########################################################
# F2PY functions
###########################################################
  
  @staticmethod
  def efs_add_config():  
    for c in g.configs['configs']:
      efs.add_config( 
                     6.5,
                     c['alat'], 
                     c['uv'], 
                     c['coords_label_id'], 
                     c['coords'], 
                     c['energy'], 
                     c['forces'], 
                     c['stress'], 
                     c['f'], 
                     c['s']
                    )                 
    efs.make_nl()
  
###########################################
#  CLASS pwscf_outpu
###########################################
class pwscf_output:

  def __init__(self, file_in=None):
    self.reset()
    if(file_in != None):
      self.load(file_in)
      self.calcs()

  def reset(self):
  
    self.z = numpy.zeros((3,3))
    
# Important so store in it's own variable
    self.atom_count = 1   

# Control
    self.data = {
      "ok": False,
      "job_done": False,
      "error": False,
      "error_code": None,
      "converged": False,
      "converged_in": None,
      
      "type": None,
      "summary": None,
      "mpi_processes": None,
      "threads_per_mpi_process": None,
      
      "species": {},
      "scf_settings": None,      
      "crystals": [],
      "results": [],
      
      "initial_positions": None,   
      "total_energy": None,
      "density_full": None,
      "stress": numpy.zeros((3,3)),
      "stress_sum": None,      
      "cpu_time": None,
      "wall_time": None,   
      "xyz": [],
      
      "mass_per_crystal": 0.0,
      "density": [],
    }
    
# Defaults
    self.xyz_units = 'evang'
    self.stress_units = 'gpa'
    
# Constants
    self.avogadro = 6.02214086E23

  def scf_settings(self):
    return {
    "bravais_lattice_index": None,
    "alat": None,
    "volume": None,
    "electrons": None, 
    "electrons_up": None, 
    "electrons_down": None, 
    "ecut_wfc": None, 
    "ecut_rho": None, 
    "convergence_threshold": None, 
    "mixing_beta": None, 
    "atomic_species": {},
    }
    
  def scf_crystal(self):
    return {
    "alat": 0.0,
    "cell_parameters": numpy.zeros((3,3)),
    "position_units": None,
    "atomic_labels": [],
    "atomic_positions": numpy.zeros((self.atom_count,3)),    
    "alat_adj": 0.0,
    "cell_parameters_adj": numpy.zeros((3,3)),
    "crystal_positions": numpy.zeros((self.atom_count,3)),
    "cell_volume": 0.0,
    "cell_density": 0.0,
    }

  def scf_results(self):
    return {
    "energy": 0.0,
    "total_force": 0.0,
    "stress": numpy.zeros((3,3)),
    "forces": numpy.zeros((self.atom_count,3)),
    "f_on": False,
    "s_on": False,
    }

#  Load, and use in another program
  def load(self, file_name): 
  
# Load data from file
    data_file = self.load_from_file(file_name)
    self.d = data_file.split("\n") 
  
#print(self.d)
  
# Reset data store
    self.reset()
       
# Load
    self.load_status()
    self.load_type()
    self.load_times()
    self.load_count()
    self.load_cpuinfo()
    self.load_crystal()
    self.load_scf_settings()
    self.load_results()
    self.load_species()
    
#print(len(self.data['crystals']))
#print(len(self.data['results']))
#self.load_scf('final_scf')
#self.load_results('initial_scf')
#self.load_results('final_scf')
    
  def load_status(self):    
# OK
###################################
    self.data['ok'] = False
    counter = 0

    for line in self.d:
      line = line.strip()
      if(pwscf_output.compare(line, "JOB DONE.")):
        self.data['job_done'] = True
      if(pwscf_output.compare(line, "Exit code:")):
        self.data['error'] = True              
        try:
          self.data['error_code'] = int(line[13:-1])
        except:
          pass
      if(pwscf_output.compare(line, "convergence has been achieved")):
        self.data['converged'] = True
        try:
          self.data['converged_in'] = int(line[39:42])
        except:
          pass
    if(self.data['job_done']):
      self.data['ok'] = True
  
  def load_type(self):
# Calc Type
###################################
    self.data['type'] = "SCF"
    for line in self.d:
      line = line.strip()
      if(line[0:23] == "A final scf calculation"):
        self.data['type'] = "VC-RELAX"
    
  def load_times(self):
# Calc Type
###################################
    for line in self.d:
      line = line.strip()
      if(line[0:14] == "PWSCF        :"):
        fa = line.split(':')
        fb = fa[1].split('CPU')
        cpu = fb[0].strip()
        fc = fb[1].split('WALL')
        wall = fc[0].strip()   
        self.data['cpu_time'] = pwscf_output.time_in_seconds(cpu)  
        self.data['wall_time'] = pwscf_output.time_in_seconds(wall)  
#print(self.data['cpu_time'])
#print(self.data['wall_time'])
      
  @staticmethod
  def time_in_seconds(inp): 
    h = 0.0
    m = 0.0
    s = 0.0
    if('h' in inp):
      f = inp.split('h')
      try:  
        h = float(f[0])
      except:
        pass
      inp = f[1]
    if('m' in inp):
      f = inp.split('m')
      try:  
        m = float(f[0])
      except:
        pass
      inp = f[1]
    if('s' in inp):
      f = inp.split('s')
      try:  
        s = float(f[0])
      except:
        pass
    return 3600 * h + 60 * m + s
        
    """    
PWSCF        : 15m15.00s CPU    16m33.12s WALL 
PWSCF        :    28.10s CPU        29.31s WALL 
PWSCF        :     1h29m CPU        1h32m WALL         
    """
   
  def load_count(self):
    n = 0
    while(n < len(self.d)):
      n, line, line_uc = self.next_line(n, self.d)
      if(pwscf_output.compare(line, "number of atoms/cell      =")):
        count = pwscf_output.extract(line, "=", "", "i")  
        try:
          self.atom_count = int(count)
          return self.atom_count 
        except:
          return 0
          
  def load_cpuinfo(self):
    counter = 0
    n = 0
    while(n < len(self.d)):
      n, line, line_uc = self.next_line(n, self.d)
      if(line != ""):
        counter += 1
        if(counter == 1):
          self.data['summary'] = line
        else:
          if(pwscf_output.compare(line, "Number of MPI processes:")):
            self.data['mpi_processes'] = pwscf_output.extract(line, ":", "", "i") 
          elif(pwscf_output.compare(line, "Threads/MPI process:")):
            self.data['threads_per_mpi_process'] = pwscf_output.extract(line, ":", "", "i") 
   
  def load_species(self):
    n = 0
    while(n < len(self.d)):
      n, line, line_uc = self.next_line(n, self.d)
      if(line[0:14] == "atomic species"):
        n, line, line_uc = self.next_line(n, self.d)
        while(line.strip() != ""):      
          line = pwscf_output.single_spaces(line).strip()
          f = line.split(" ")
          self.data['scf_settings'][f[0]] = [float(f[1]), float(f[2])]   # Valence electrons, atomic mass
          n, line, line_uc = self.next_line(n, self.d)
        break
      
#self.data['scf_settings']
#
#  atomic species   valence    mass     pseudopotential
#  Al             3.00    26.98200     Al( 1.00)
 
  def load_scf_settings(self):
    self.data['scf_settings'] = self.scf_settings()  
    
    n = 0
    while(n < len(self.d)):
      n, line, line_uc = self.next_line(n, self.d)
      if(pwscf_output.compare(line, "bravais-lattice index     =")):
        self.data['scf_settings']['bravais_lattice_index'] = pwscf_output.extract(line, "=", "", "i") 
      elif(pwscf_output.compare(line, "lattice parameter (alat)  =")):
        self.data['scf_settings']['alat'] = pwscf_output.extract(line, "=", "a.u.", "f")  
      elif(pwscf_output.compare(line, "unit-cell volume          =")):
        self.data['scf_settings']['volume'] = pwscf_output.extract(line, "=", "(a.u.)^3", "f")     
      elif(pwscf_output.compare(line, "number of atoms/cell      =")):
        self.data['scf_settings']['nat'] = pwscf_output.extract(line, "=", "", "i")  
      elif(pwscf_output.compare(line, "number of atomic types    =")):
        self.data['scf_settings']['types'] = pwscf_output.extract(line, "=", "", "i")  
      elif(pwscf_output.compare(line, "number of electrons       =")):
        str_e = pwscf_output.extract(line, "=", "", "s")
        e, eu, ed = pwscf_output.electron_string(str_e)
        self.data['scf_settings']['electrons'] = e
        self.data['scf_settings']['electrons_up'] = eu
        self.data['scf_settings']['electrons_down'] = ed
      elif(pwscf_output.compare(line, "kinetic-energy cutoff     =")):
        self.data['scf_settings']['ecut_wfc'] = pwscf_output.extract(line, "=", "Ry", "f")  
      elif(pwscf_output.compare(line, "charge density cutoff     =")):
        self.data['scf_settings']['ecut_rho'] = pwscf_output.extract(line, "=", "Ry", "f")  
      elif(pwscf_output.compare(line, "convergence threshold     =")):
        self.data['scf_settings']['convergence_threshold'] = pwscf_output.extract(line, "=", "", "f")  
      elif(pwscf_output.compare(line, "mixing beta               =")):
        self.data['scf_settings']['mixing_beta'] = pwscf_output.extract(line, "=", "", "f")  
      elif(("atomic species" in line) and ("valence" in line) and ("mass" in line) and ("pseudopotential" in line)):
        loop = True
        while(loop):
          n, line, line_uc = self.next_line(n, self.d)
          if(line.strip() == ""):
            loop = False
          else:  
            line = pwscf_output.single_spaces(line).strip()
            line_arr = line.split(" ")
            self.data['scf_settings']['atomic_species'][line_arr[0]] = {}
            self.data['scf_settings']['atomic_species'][line_arr[0]]['valence'] = line_arr[1]
            self.data['scf_settings']['atomic_species'][line_arr[0]]['mass'] = line_arr[2]

# End of file/loop
        n = len(self.d)

###################################
# LOAD CRYSTALS FROM OUTPUT FILE
###################################

  def load_crystal(self):
# Make new list for crystals
    self.data['crystals'] = []
    
# FIRST
    n = 0
    crystal = self.scf_crystal()
    while(n < len(self.d)):
      n, line, line_uc = self.next_line(n, self.d)  
      if(line[0:10] == "celldm(1)="):
        crystal['alat'] = float(line[10:21].strip())
      elif(pwscf_output.compare(line.strip(), "crystal axes:")):
        for j in range(3):              
          n, line, line_uc = self.next_line(n, self.d)
          fields = pwscf_output.extract(line, "= (", ")", "s", " ")
          for i in range(len(fields)):
            crystal['cell_parameters'][j, i] = float(fields[i])
      elif(pwscf_output.compare(line.strip(), "Cartesian axes")):
        n, line, line_uc = self.next_line(n, self.d)
        n, line, line_uc = self.next_line(n, self.d)
        
# Unit
        crystal['position_units'] = "alat"
        
        loop = True 
        k = 0
        while(loop):
          n, line, line_uc = self.next_line(n, self.d)   
          if(line.strip() == ""):
            loop = False
          else:
            line_arr = line.split("tau(")
            label = line_arr[0][-15:].strip()
            crystal['atomic_labels'].append(label)
            
            coords = line_arr[1]
            x = float(coords[9:21])
            y = float(coords[22:33])
            z = float(coords[34:44])
            crystal['atomic_positions'][k, 0] = x
            crystal['atomic_positions'][k, 1] = y
            crystal['atomic_positions'][k, 2] = z 
            
# Increment
            k = k + 1
        
# Cell Volume
        crystal['cell_volume'] = pwscf_output.cell_volume(crystal['alat'], crystal['cell_parameters'][:,:])    
        
# Add/Save
        self.data['crystals'].append(crystal)
        n = len(self.d)
    
# MIDDLE
    n = 0
    while(n < len(self.d)):
      n, line, line_uc = self.next_line(n, self.d)   
      if(line[0:15] == "CELL_PARAMETERS"):
# Create
        crystal = self.scf_crystal()
        
# Get alat
        line_arr = line.split("=")
        line_arr = line_arr[1].split(")")
        crystal['alat'] = float(line_arr[0].strip())
        
#Cell Parameters
        for j in range(3):              
          n, line, line_uc = self.next_line(n, self.d)
          line = pwscf_output.single_spaces(line)
          fields = line.split(" ")
          for i in range(len(fields)):
            crystal['cell_parameters'][j, i] = float(fields[i])
      elif(line[0:9] == "density ="):  
        crystal['density'] = float(line[10:23])      
            
      elif(line[0:16] == "ATOMIC_POSITIONS"):     
        
# Unit
        crystal['position_units'] = "crystal"
        
# Read Coords
        loop = True 
        k = 0
        while(loop):
          n, line, line_uc = self.next_line(n, self.d)   
          if(line.strip() == ""):
            loop = False
          elif(line.strip() == "End final coordinates"):
            loop = False
          else:
            line = pwscf_output.single_spaces(line)
            line_arr = line.split(" ")
            crystal['atomic_labels'].append(line_arr[0])
            
            crystal['atomic_positions'][k, 0] = float(line_arr[1])
            crystal['atomic_positions'][k, 1] = float(line_arr[2])
            crystal['atomic_positions'][k, 2] = float(line_arr[3]) 
            
# Increment
            k = k + 1
            
# Cell Volume
        crystal['cell_volume'] = pwscf_output.cell_volume(crystal['alat'], crystal['cell_parameters'][:,:])    
            
# Add/Save
        self.data['crystals'].append(crystal)
    
# END
    n = 0
    crystal = self.scf_crystal()
    d = 0
    while(n < len(self.d)):
      n, line, line_uc = self.next_line(n, self.d)   
      if(line[0:10] == "celldm(1)="):
        d = d + 1
        if(d == 2):
          crystal['alat'] = float(line[10:21].strip())
      elif(d == 2 and pwscf_output.compare(line.strip(), "crystal axes:")):
        for j in range(3):              
          n, line, line_uc = self.next_line(n, self.d)
          fields = pwscf_output.extract(line, "= (", ")", "s", " ")
          for i in range(len(fields)):
            crystal['cell_parameters'][j, i] = float(fields[i])
     
      elif(d == 2 and pwscf_output.compare(line.strip(), "Cartesian axes")):      
        
# Unit
        crystal['position_units'] = "alat"
        
# Read coords
        n, line, line_uc = self.next_line(n, self.d)
        n, line, line_uc = self.next_line(n, self.d)
        
        loop = True 
        k = 0
        while(loop):
          n, line, line_uc = self.next_line(n, self.d)   
          if(line.strip() == ""):
            loop = False
          else:
            line_arr = line.split("tau(")
            label = line_arr[0][-15:].strip()
            crystal['atomic_labels'].append(label)
            
            coords = line_arr[1]
            x = float(coords[9:21])
            y = float(coords[22:33])
            z = float(coords[34:44])
            crystal['atomic_positions'][k, 0] = x
            crystal['atomic_positions'][k, 1] = y
            crystal['atomic_positions'][k, 2] = z 
            
# Increment
            k = k + 1
            
# Cell Volume
        crystal['cell_volume'] = pwscf_output.cell_volume(crystal['alat'], crystal['cell_parameters'][:,:])    
            
# Add/Save
        self.data['crystals'].append(crystal)
        n = len(self.d)
    
# Loop through crystals
    for i in range(len(self.data['crystals'])):
# Adjust alat and cell_parameters so celldm(1)=1.0
      factor = 1.0 / self.data['crystals'][i]['cell_parameters'][0, 0]
      
      self.data['crystals'][i]['alat_adj'] = self.data['crystals'][i]['alat'] * self.data['crystals'][i]['cell_parameters'][0, 0]
      self.data['crystals'][i]['cell_parameters_adj'][:, :] = factor * self.data['crystals'][i]['cell_parameters'][:, :]
    
# Make crystal_positions
      if(self.data['crystals'][i]['position_units'] == 'crystal'):
        self.data['crystals'][i]['crystal_positions'][:,:] = self.data['crystals'][i]['atomic_positions'][:,:] 
      elif(self.data['crystals'][i]['position_units'] == 'alat'):  
        minv = numpy.linalg.inv(self.data['crystals'][i]['cell_parameters'][:, :])
        for j in range(len(self.data['crystals'][i]['atomic_positions'])):
          self.data['crystals'][i]['crystal_positions'][j, :] = numpy.matmul(minv[:,:], self.data['crystals'][i]['atomic_positions'][j, :])

  def load_results(self):
  
# Make new list for results
    self.data['results'] = []
    
    n = 0
    while(n < len(self.d)):
      n, line, line_uc = self.next_line(n, self.d)
      
# READ ENERGY
      if(pwscf_output.compare(line, "!    total energy")):
# Create dictionary
        results = self.scf_results()  
        results['energy'] = pwscf_output.extract(line, "=", "Ry", "f") 
        
# READ FORCES
      elif(pwscf_output.compare(line, "Forces acting on atoms")):
        n, line, line_uc = self.next_line(n, self.d)  
        loop = True
        f = 0
        while(loop):
          n, line, line_uc = self.next_line(n, self.d)  
          if(line.strip() == ""):
            loop = False
          else:
            line_arr = line.split("force =")
            fields = pwscf_output.single_spaces(line_arr[1].strip()).split(" ")
            results['forces'][f,0] = float(fields[0])
            results['forces'][f,1] = float(fields[1])
            results['forces'][f,2] = float(fields[2])
            f = f + 1
        if(f>0):
          results['f_on'] = True
        
# READ TOTAL FORCE
      elif(pwscf_output.compare(line, "Total force =")):
        results['total_force'] = pwscf_output.extract(line, "=", "T", "f")
        
# READ STRESS
      elif(pwscf_output.compare(line, "total   stress  (Ry/bohr**3)")):  
        results['s_on'] = True      
        for j in range(3):              
          n, line, line_uc = self.next_line(n, self.d)  
          fields = pwscf_output.extract(line, "", "", "f", " ", True)  
          results['stress'][j,0] = fields[3] 
          results['stress'][j,1] = fields[4] 
          results['stress'][j,2] = fields[5]

#SAVE
        self.data['results'].append(results)

  def calcs(self):
  
# MASS PER CELL
    self.data['mass_per_crystal'] = 0.0    
    for l in self.data['crystals'][0]['atomic_labels']:
      self.data['mass_per_crystal'] = self.data['mass_per_crystal'] + self.data['scf_settings'][l][1]    
      
# CELL VOLUMES
    for i in range(len(self.data['crystals'])):
      self.data['crystals'][i]['cell_volume'] = pwscf_output.cell_volume(self.data['crystals'][i]['alat'], self.data['crystals'][i]['cell_parameters'][:,:])  
      self.data['crystals'][i]['cell_density'] = pwscf_output.cell_density(self.data['crystals'][i]['cell_volume'], self.data['mass_per_crystal'])
      
  def next_line(self, n, data):
    if(n < len(data)):
      line = data[n].strip()
      line_uc = line.upper()
      n = n + 1
      return n, line, line_uc
    else:
      n = n + 1
      return n, None, None
    
  def store(self, store, line, field, n=0):
    l, f = pwscf_output.read_line(line, field)  
    if(l != False):
      self.data[store] = f[n]

#  Run as it's own program
  def run(self):
    self.reset()

    option = ""
    file_name = ""

    if(len(sys.argv) > 1 and sys.argv[1] is not None):
      option = sys.argv[1]

    if(len(sys.argv) > 2 and sys.argv[2] is not None):
      file_name = sys.argv[2]

    if(option.lower().strip() == "" or option.lower().strip() == "interactive"):
      self.menu()
      exit()
    elif(option.lower().strip() == "quiet"):
      print("Quiet")
    else:
      return 0

#################################
# READ/LOAD input file
#################################

  def load_from_file(self, file_name):
# Init variable
    file_data = ""

# Read it in line by line
    fh = open(file_name, "r")
    for file_row in fh:
      file_data = file_data + file_row.strip() + '\n'

    return file_data

#################################
# Get
#################################

  def get_nat(self): 
    return self.atom_count

  def get_alat(self):
    return self.data['alat']
    
  def get_volume(self):
    return self.data['scf_settings']['volume'] 
    
  def get_volume_per_atom(self):
    return self.data['scf_settings']['volume'] / self.atom_count
  
  def get_total_energy(self, n=None):
    if(n == None):
      return self.data['results'][-1]['energy']  
    else:
      return self.data['results'][n]['energy']  
    
  def get_energy_per_atom(self, n=None):
    if(n == None):
      return self.data['results'][-1]['energy'] / self.atom_count  
    else:
      return self.data['results'][n]['energy'] / self.atom_count 
    
  def get_total_force(self):
    if(n == None):
      return self.data['results'][-1]['total_force']  
    else:
      return self.data['results'][n]['total_force']
    
  def get_force_per_atom(self, n=None):
    if(n == None):
      return self.data['results'][-1]['total_force'] / self.atom_count  
    else:
      return self.data['results'][n]['total_force'] / self.atom_count  
  
  def get_density(self):
    return self.data['density']  
    
  def get_cell_parameters(self):
    cp = ['alat', 
          [str(self.data['crystal_calc'][0,0]), str(self.data['crystal_calc'][0,1]), str(self.data['crystal_calc'][0,2])], 
          [str(self.data['crystal_calc'][1,0]), str(self.data['crystal_calc'][1,1]), str(self.data['crystal_calc'][1,2])], 
          [str(self.data['crystal_calc'][2,0]), str(self.data['crystal_calc'][2,1]), str(self.data['crystal_calc'][2,2])]]
    return cp

# Return relaxed unit vector
  def get_cell_array(self):
    return self.data['crystal_calc']

# return alat and normalised unit vector
  def get_norm_relaxed(self):
    alat = self.data['crystals'][-2]['alat_adj']
    cp = self.data['crystals'][-2]['cell_parameters_adj']
    return alat, cp

# Get stress
  def get_stress(self):
    return self.data['stress']
    
  def get_stress_sum(self):
    return self.data['stress_sum']
  
  def get_job_done(self):
    return self.data['job_done']
    
  def get_job_error(self):
    return self.data['error']
    
  def get_job_converged(self):
    return self.data['converged']

  def get_ok(self):
    return self.data['ok']

# GET - VC-RELAXED

  def get_vc_relax_alat(self):
    if(self.data['type']=="VC-RELAX"):
      return self.data['crystals'][-1]['alat_adj']
    else:
      return None

  def get_vc_relax_cp(self):
    if(self.data['type']=="VC-RELAX"):
      return self.data['crystals'][-1]['cell_parameters_adj']
    else:
      return None

  def get_vc_relax_volume(self):
    if(self.data['type']=="VC-RELAX"):
      return self.data['crystals'][-1]['cell_volume']
    else:
      return None

  def get_vc_relax_density(self):
    if(self.data['type']=="VC-RELAX"):
      return self.data['crystals'][-1]['cell_density']
    else:
      return None
      
  def get_mass_per_crystal(self):
    return self.data['mass_per_crystal']
      
  def get_times(self):
    return self.data['cpu_time'], self.data['wall_time']
      
#################################
# Interactive
#################################

  def menu(self):
    while(True):
      choice = self.print_menu().upper()
      print(choice)
      if(choice == "X"):
        exit()
      elif(choice == "1"):
        self.i_load()
      elif(choice == "2"):
        self.i_display()

  def print_menu(self):
    pwscf_output.header("Menu")
    print("1. Load File")
    print("2. Display File")
    print("X. Exit")
    return input("Choice: ")

  def i_load(self):
    pwscf_output.header("Load Output File")
    file_name = input("Enter file name: ")
    self.load(file_name)
    print("File loaded.")
    input()

  def i_display(self):
    pwscf_output.header("Display File")
    self.output_details()
    input()

  def output_details(self):
    print("Output")
    print("=======================================================================")
    for key in sorted(self.data.keys()):
      value = self.data[key]
      print(key, ":  ", value)
    print("=======================================================================")
    print()
    
  def xyz_evang(self):
    self.xyz_units = 'ev/ang'
    self.stress_units = 'gpa'

  def xyz_stress_gpa(self):
    self.stress_units = 'gpa'

  def make_xyz(self, option=None):
    self.xyz = []
    if(option == None):
      for rn in range(len(self.data['results'])):
        option = rn + 1
    elif(option == -1):
      option = len(self.data['results'])
    else:
      option = (option - 1) % len(self.data['results']) + 1
    self.make_xyz_inner(option)
    return self.xyz
    
  def make_xyz_inner(self, option):
    if(len(self.data['results'])==0):
      return False
    if(len(self.data['crystals'])==0):
      return False  
   
    rn = (option - 1) % len(self.data['results'])
    cn = rn
    if(rn == len(self.data['results']) - 1):
      cn = len(self.data['crystals']) - 1
      
    crystal = self.data['crystals'][cn]
    result = self.data['results'][rn]
    settings = self.data['scf_settings']
    species = settings['atomic_species']
    
# Add new list and set counter n
    self.xyz.append([])
    n = len(self.xyz) - 1
    
# Add data
    self.xyz[n].append("#ALAT " + str(crystal['alat_adj']))
    self.xyz[n].append("#X " + str(crystal['cell_parameters_adj'][0][0]) + " " + str(crystal['cell_parameters_adj'][0][1]) + " " + str(crystal['cell_parameters_adj'][0][2]))
    self.xyz[n].append("#Y " + str(crystal['cell_parameters_adj'][1][0]) + " " + str(crystal['cell_parameters_adj'][1][1]) + " " + str(crystal['cell_parameters_adj'][1][2]))
    self.xyz[n].append("#Z " + str(crystal['cell_parameters_adj'][2][0]) + " " + str(crystal['cell_parameters_adj'][2][1]) + " " + str(crystal['cell_parameters_adj'][2][2]))
    
# Just use 1 1 1
    self.xyz[n].append("#C 2 2 2")
    self.xyz[n].append("#RCUT 6.5")
    self.xyz[n].append("#L_UNITS bohr")
    self.xyz[n].append("#E_UNITS ry")
    self.xyz[n].append("#F_UNITS ry/bohr")
    self.xyz[n].append("#S_UNITS kbar")
    
    self.xyz[n].append("#E " + str(result['energy']))
     
    if(result['s_on']):
      self.xyz[n].append("#SX " + str(result['stress'][0,0]) + " " +  str(result['stress'][0,1]) + " " + str(result['stress'][0,2]))
      self.xyz[n].append("#SY " + str(result['stress'][1,0]) + " " +  str(result['stress'][1,1]) + " " + str(result['stress'][1,2]))
      self.xyz[n].append("#SZ " + str(result['stress'][2,0]) + " " +  str(result['stress'][2,1]) + " " + str(result['stress'][2,2]))
      
    for label in species.keys():
      self.xyz[n].append("#M " + label + " " + str(species[label]['mass']))
    
    for i in range(self.atom_count):
      line = crystal['atomic_labels'][i]
      line = line + " " + str(crystal['crystal_positions'][i,0]) + " " + str(crystal['crystal_positions'][i,1]) + " " + str(crystal['crystal_positions'][i,2])      
      if(result['f_on']):   
        line = line + " " + str(result['forces'][i,0]) + " " + str(result['forces'][i,1]) + " " + str(result['forces'][i,2])
#"s_on": False,
      self.xyz[n].append(line)
      
# Return
    return self.xyz[n]
    
  def get_data(self, file=None): 
  
    """
# Important so store in it's own variable
    self.atom_count = 1    

# Control
    self.data = {
      "ok": False,
      "job_done": False,
      "error": False,
      "type": None,
      "summary": None,
      "mpi_processes": None,
      "threads_per_mpi_process": None,
      
      "scf_settings": None,      
      "crystals": [],
      "results": [],
      
      "initial_positions": None,   
      "total_energy": None,
      "density_full": None,
      "density": None,
      "stress": numpy.zeros((3,3)),
      "stress_sum": None,      
      "cpu_time": None,
      "wall_time": None,   
      "xyz": [],
    }
    """
    
    out = "##############################################################################################################\n"
    out = out + "atom count:                  " + str(self.atom_count) + "\n"
    out = out + "ok:                          " + str(self.data['ok']) + "\n"
    out = out + "job_done:                    " + str(self.data['job_done']) + "\n"
    out = out + "error:                       " + str(self.data['error']) + "\n"
    out = out + "type:                        " + str(self.data['type']) + "\n"
    out = out + "summary:                     " + str(self.data['summary']) + "\n"
    out = out + "mpi_processes:               " + str(self.data['mpi_processes']) + "\n"
    out = out + "threads_per_mpi_process:     " + str(self.data['threads_per_mpi_process']) + "\n"
    out = out + "scf_settings:                " + str(len(self.data['scf_settings'])) + "\n"
    
    for k in self.data['scf_settings'].keys():    
      out = out + "                             " + k + '  '+ str(self.data['scf_settings'][k]) + "\n"  
    
    out = out + "crystals (count):            " + str(len(self.data['crystals'])) + "\n"
    out = out + "results (count):             " + str(len(self.data['results'])) + "\n"
    i = 0
    for i in range(len(self.data['crystals'])-1):
      out = out + "##############################################################################################################\n"
      c = self.data['crystals'][i]
      n = i + 1
      if(n == 1):
        out = out + "crystal " + str(n) + " (input):\n"         
      elif(n == len(self.data['crystals']) - 1):
        out = out + "crystal " + str(n) + " (relaxed):\n" 
      else:
        out = out + "crystal " + str(n) + ":\n"  
      out = out + "##############################################################################################################\n"
      out = out + "alat:                        " + str(c['alat']) + "\n"  
      out = out + "cp:                          " + str(c['cell_parameters'][0,0]) + " " + str(c['cell_parameters'][0,1]) + " " + str(c['cell_parameters'][0,2]) + "\n" 
      out = out + "                             " + str(c['cell_parameters'][1,0]) + " " + str(c['cell_parameters'][1,1]) + " " + str(c['cell_parameters'][1,2]) + "\n"
      out = out + "                             " + str(c['cell_parameters'][2,0]) + " " + str(c['cell_parameters'][2,1]) + " " + str(c['cell_parameters'][2,2]) + "\n" 
      out = out + "alat adjusted:               " + str(c['alat_adj']) + "\n"  
      out = out + "cp adjusted:                 " + str(c['cell_parameters_adj'][0,0]) + " " + str(c['cell_parameters_adj'][0,1]) + " " + str(c['cell_parameters_adj'][0,2]) + "\n" 
      out = out + "                             " + str(c['cell_parameters_adj'][1,0]) + " " + str(c['cell_parameters_adj'][1,1]) + " " + str(c['cell_parameters_adj'][1,2]) + "\n"
      out = out + "                             " + str(c['cell_parameters_adj'][2,0]) + " " + str(c['cell_parameters_adj'][2,1]) + " " + str(c['cell_parameters_adj'][2,2]) + "\n" 
      if(i<len(self.data['crystals'])-2):
        s = self.data['results'][i]
        out = out + "energy:                      " + str(s['energy']) + "\n"  
        out = out + "total_force:                 " + str(s['total_force']) + "\n"  
        out = out + "stress:                      " + str(s['stress'][0,0]) + " " + str(s['stress'][0,1]) + " " + str(s['stress'][0,2]) + "\n"  
        out = out + "stress:                      " + str(s['stress'][1,0]) + " " + str(s['stress'][1,1]) + " " + str(s['stress'][1,2]) + "\n"  
        out = out + "stress:                      " + str(s['stress'][2,0]) + " " + str(s['stress'][2,1]) + " " + str(s['stress'][2,2]) + "\n"   
        out = out + "alat positions,  crystal positions,  forces:   " + "\n" 
        for m in range(self.atom_count):
          out = out + "                             " + str(c['atomic_positions'][m,0]) + " " + str(c['atomic_positions'][m,1]) + " " + str(c['atomic_positions'][m,2]) + "    " + str(c['crystal_positions'][m,0]) + " " + str(c['crystal_positions'][m,1]) + " " + str(c['crystal_positions'][m,2]) + "    " + str(s['forces'][m,0]) + " " + str(s['forces'][m,1]) + " " + str(s['forces'][m,2]) + "\n"
      
      out = out +     "cell volume:                 " + str(c['cell_volume']) + "\n"  
      out = out +     "density:                     " + str(c['cell_density']) + "\n" 
       
    c = self.data['crystals'][-1]
    s = self.data['results'][-1]
    n = len(self.data['crystals'])
    
    out = out + "##############################################################################################################\n"
    out = out + "crystal " + str(n) + " (final): \n"     
    out = out + "##############################################################################################################\n"
    out = out + "alat:                        " + str(c['alat']) + "\n" 
    out = out + "cp:                          " + str(c['cell_parameters'][0,0]) + " " + str(c['cell_parameters'][0,1]) + " " + str(c['cell_parameters'][0,2]) + "\n" 
    out = out + "                             " + str(c['cell_parameters'][1,0]) + " " + str(c['cell_parameters'][1,1]) + " " + str(c['cell_parameters'][1,2]) + "\n"
    out = out + "                             " + str(c['cell_parameters'][2,0]) + " " + str(c['cell_parameters'][2,1]) + " " + str(c['cell_parameters'][2,2]) + "\n" 
    out = out + "alat adjusted:               " + str(c['alat_adj']) + "\n"  
    out = out + "cp adjusted:                 " + str(c['cell_parameters_adj'][0,0]) + " " + str(c['cell_parameters_adj'][0,1]) + " " + str(c['cell_parameters_adj'][0,2]) + "\n" 
    out = out + "                             " + str(c['cell_parameters_adj'][1,0]) + " " + str(c['cell_parameters_adj'][1,1]) + " " + str(c['cell_parameters_adj'][1,2]) + "\n"
    out = out + "                             " + str(c['cell_parameters_adj'][2,0]) + " " + str(c['cell_parameters_adj'][2,1]) + " " + str(c['cell_parameters_adj'][2,2]) + "\n" 

    out = out + "energy:                      " + str(s['energy']) + "\n"  
    out = out + "total_force:                 " + str(s['total_force']) + "\n"  
    out = out + "stress:                      " + str(s['stress'][0,0]) + " " + str(s['stress'][0,1]) + " " + str(s['stress'][0,2]) + "\n"  
    out = out + "stress:                      " + str(s['stress'][1,0]) + " " + str(s['stress'][1,1]) + " " + str(s['stress'][1,2]) + "\n"  
    out = out + "stress:                      " + str(s['stress'][2,0]) + " " + str(s['stress'][2,1]) + " " + str(s['stress'][2,2]) + "\n"   
    out = out + "alat positions,  crystal positions,  forces:   " + "\n" 
    for m in range(self.atom_count):
      out = out + "                             " + str(c['atomic_positions'][m,0]) + " " + str(c['atomic_positions'][m,1]) + " " + str(c['atomic_positions'][m,2]) + "    " + str(c['crystal_positions'][m,0]) + " " + str(c['crystal_positions'][m,1]) + " " + str(c['crystal_positions'][m,2]) + "    " + str(s['forces'][m,0]) + " " + str(s['forces'][m,1]) + " " + str(s['forces'][m,2]) + "\n"
         
    out = out +     "cell volume:                 " + str(c['cell_volume']) + "\n"  
    out = out +     "density:                     " + str(c['cell_density']) + "\n" 
    
    fh = open(file,'w')
    fh.write(out)
    fh.close()
    
    return out
    
#################################
# Static Methods
#################################

  @staticmethod
  def remove_spaces(input_string):
    return input_string.replace(" ", "")
    
  @staticmethod
  def extract(input_string, start=None, end=None, type=None, split=None, trim=False):
    if(start == ""):
      start = None
    if(end == ""):
      end = None
    if(trim):
      input_string = input_string.strip()
    
# Start/End
    start_n = None
    end_n = None
      
    if(start == None and end == None):   
      start_n = 0
      end_n = len(input_string)
    elif(start == None and end != None):  
      end_l = len(end)   
      start_n = 0
      for n in range(len(input_string)):
        if(input_string[n:n+end_l] == end[0:end_l]):
          end_n = n
          break
    elif(start != None and end == None):  
      start_l = len(start)
      end_n = len(input_string)
      for n in range(len(input_string)):
        if(input_string[n:n+start_l] == start[0:start_l]):
          start_n = n + start_l
    else:  
      start_l = len(start)
      end_l = len(end)  
    
      for n in range(len(input_string)):
        if(input_string[n:n+start_l] == start[0:start_l]):
          start_n = n + start_l
        if(start_n != None and input_string[n:n+end_l] == end[0:end_l]):
          end_n = n
          break
        
# Read
    result = input_string[start_n:end_n].strip()       

# Split
    if(split != None):
      if(split == " "):
#result = re.sub(r'\s\s+', ' ', result)
        result = pwscf_output.single_spaces(result)
      result = result.split(split)
      for i in range(len(result)):
        if(type.lower() == "f"):
          result[i] = float(result[i])
        elif(type.lower() == "i"):
          result[i] = int(result[i])
        
    else:  
      if(type.lower() == "f"):
        result = float(result)
      elif(type.lower() == "i"):
        result = int(result)
        
# Return
    return result
      
  @staticmethod
  def compare(line, field):
    line = line.strip()
    line = line.upper() 
    
    field = field.strip()
    field = field.upper()
    
    f_len = len(field)
    if(len(line) >= f_len and line[0:f_len] == field[0:f_len]):
      return True
    return False
    
  @staticmethod
  def read_line(line, field):
    line = line.strip()
#line = re.sub(r'\s\s+', ' ', line)
#line = re.sub(r'\s=\s', '=', line)
    line = pwscf_output.clean(line)
    line_uc = line.upper() 
    
    field = field.strip()
#field = re.sub(r'\s\s+', ' ', field)
#field = re.sub(r'\s=\s', '=', field)
    field = pwscf_output.clean(field)
    field = field.upper()
    
    f_len = len(field)
    if(len(line_uc) >= f_len and line_uc[0:f_len] == field[0:f_len]):
      output = line[f_len:].strip()
      fields = output.split(" ")
      return output, fields      
    return False, False
    
  @staticmethod
  def fields(input_string):
    input_string = input_string.strip()
    output_string = ""
    last = None
    for character in input_string:
      if(character != " " or (character == " " and last != " ")):
        output_string += character
    return output_string.split(" ")
    
  @staticmethod
  def check_keyword(line, keyword):
    if(line.upper()[0:len(keyword)] == keyword.upper()):
      return True
    return False

  @staticmethod
  def clear_screen():
    os.system('cls' if os.name == 'nt' else 'clear')

  @staticmethod
  def header(sub_title=""):
    pwscf_output.clear_screen()
    print("==========================================================")
    print("                    PWscf Input Editor                    ")
    print("==========================================================")
    print()
    print(sub_title)
    print()
    print()
    print()

  @staticmethod
  def process_keyword(str_in):
    str_in = str_in.lower().strip()
    str_in = pwscf_output.remove_spaces(str_in)
    id = None
    keyword = ""
    flag = 0
    for character in str_in:
      if(character == "("):
        id = ""
        flag = 1
      elif(character == ")"):
        flag = 2
      elif(flag == 0):
        keyword += character
      elif(flag == 1):
        id = id + character
    if(id != None):
      try:
        id = int(id)
      except:
        id = None
    return keyword, id  

  @staticmethod
  def add_keyword(keywords, keyword, id, value):
    if(id == None):
      added = False
      for i in range(len(keywords)):
        if(keywords[i][0] == keyword):
          added = True
          keywords[i][1] = keyword
      if(added == False):
        keywords.append([keyword, value])
    else:   
      n = None
      for i in range(len(keywords)):
        if(keywords[i][0] == keyword):
          n = i
          break
      if(n == None):    
        keywords.append([keyword,[None]])
        n = len(keywords) - 1
        
      while(len(keywords[n][1]) < id):
        keywords[n][1].append(None)

      keywords[n][1][id-1] = value  

  @staticmethod
  def make_line(key, value):
    output = ""
    if(value != None):
       if(isinstance(value, (list,))):
         for i in range(len(value)):
           if(value[i] != None):
             output += key + "(" + str(i+1) + ") = " + value[i] + ", \n"                
       else:
         output += key + " = " + value + ", \n"   
    return output    

  @staticmethod
  def coord_format(float_in):
    pad = "              "
    value = str(round(float_in, 6)).strip()
    return value
    
  @staticmethod
  def label_format(label):  
    pad = "              "
    label = label.strip()
    return label
    
  @staticmethod
  def is_zero(arr):
    for i in range(arr.shape[0]):
      for j in range(arr.shape[1]):
        if(arr[i, j] != 0.0):
          return False
    return True
    
  @staticmethod
  def clean(str_in):  
    str_out = ""
    l = len(str_in)
    for i in range(l):
# Last, Next, This
      if(i == 0):
        last = None
      else:
        last = str_in[i-1]
      if(i < (l-1)):
        next = str_in[i+1]
      else:  
        next = None
      char = str_in[i]
    
# Check
      ok = True
      if(last == " " and char == " "):
        ok = False
      elif(last == "\n" and char == "\n"):
        ok = False
      elif(last == "\n" and char == " "):
        ok = False
      elif(char == " " and next == "\n"):
        ok = False
      elif(last == "=" and char == " "):
        ok = False
      elif(char == " " and next == "="):
        ok = False
        
# Add to string
      if(ok):
        str_out += char
    return str_out    
    
  @staticmethod
  def electron_string(str_in):
    arr = str_in.split("(up:")
    e = arr[0]
    if(len(arr) == 1):
      return e.strip(), None, None
    if(len(arr)==2):
      arr_b = arr[1].split(", down:")
      eu = arr_b[0]
      arr_c = arr_b[1].split(")")
      ed = arr_c[0]
      return e.strip(), eu.strip(), ed.strip()
  
    print("TEST")
    return "","",""
    
  @staticmethod
  def single_spaces(str_in):
    str_out = ""
    last = None
    for char in str_in:
      if(char != " " or (char == " " and last != " ")):
        str_out = str_out + char
      last = char
    return str_out
    
  @staticmethod
  def cell_volume(alat, cp):
    cp_alat = numpy.zeros((3,3,),)
    cp_alat[:,:] = alat * cp[:,:]
    v = numpy.dot(cp_alat[0,:],numpy.cross(cp_alat[1,:], cp_alat[2,:]))
    return v
    
  @staticmethod
  def cell_density(v, mpc):
    v_m3 = v * 1.48036e-31
    m = mpc * 1.66054E-027
    rho = m / v_m3
    return rho
    
"""
        
  def aaa():
    
# Load
###################################
    n = 0
    counter = 0
    while(n < len(data)):
      n, line, line_uc = self.next_line(n, data)
      if(line != ""):
        counter += 1
        if(counter == 1):
          self.data['summary'] = line
        else:
          if(pwscf_output.compare(line, "Number of MPI processes:")):
            self.data['mpi_processes'] = pwscf_output.extract(line, ":", "", "i") 
        
          if(pwscf_output.compare(line, "bravais-lattice index     =")):
            self.data['bravais_lattice_index'] = pwscf_output.extract(line, "=", "", "i") 
            
          if(pwscf_output.compare(line, "lattice parameter (alat)  =")):
            self.data['alat'] = pwscf_output.extract(line, "=", "a.u.", "f")  
            
          if(pwscf_output.compare(line, "unit-cell volume          =")):
            self.data['volume'] = pwscf_output.extract(line, "=", "(a.u.)^3", "f")     
            
          if(pwscf_output.compare(line, "number of atoms/cell      =")):
            self.data['nat'] = pwscf_output.extract(line, "=", "", "i")  
            
          if(pwscf_output.compare(line, "number of atomic types    =")):
            self.data['types'] = pwscf_output.extract(line, "=", "", "i")  
            
          if(pwscf_output.compare(line, "number of electrons       =")):
            str_e = pwscf_output.extract(line, "=", "", "s")
            e, eu, ed = pwscf_output.electron_string(str_e)
            self.data['electrons'] = e
            self.data['electrons_up'] = eu
            self.data['electrons_down'] = ed
          
          if(pwscf_output.compare(line, "number of Kohn-Sham states=")):
            self.data['ks_states'] = pwscf_output.extract(line, "=", "", "i")   
            
          if(pwscf_output.compare(line, "kinetic-energy cutoff     =")):
            self.data['ecutwfc'] = pwscf_output.extract(line, "=", "Ry", "f")  
            
          if(pwscf_output.compare(line, "charge density cutoff     =")):
            self.data['ecutrho'] = pwscf_output.extract(line, "=", "Ry", "f")   
        
          if(pwscf_output.compare(line.strip(), "crystal axes:") and pwscf_output.is_zero(self.data['crystal_in'])):            
            for j in range(3):              
              n, line, line_uc = self.next_line(n, data)
              fields = pwscf_output.extract(line, "= (", ")", "s", " ")
              self.data['crystal_in'][j,:] = fields  
              self.data['crystal_calc'][j,:] = fields  
          
          if(pwscf_output.compare(line.strip(), "crystal axes:")):            
            for j in range(3):              
              n, line, line_uc = self.next_line(n, data)
              fields = pwscf_output.extract(line, "= (", ")", "s", " ")
              self.data['crystal_calc'][j,:] = fields            
        
          if(pwscf_output.compare(line, "!    total energy")):
            self.data['total_energy'] = pwscf_output.extract(line, "=", "Ry", "f")
            
          if(pwscf_output.compare(line, "Total force =")):
            self.data['total_force'] = pwscf_output.extract(line, "=", "T", "f")
            
          if(pwscf_output.compare(line, "total   stress  (Ry/bohr**3)")):        
            self.data['stress_sum'] = 0.0
            for j in range(3):              
              n, line, line_uc = self.next_line(n, data)   
              fields = pwscf_output.extract(line, "", "", "f", " ", True)  
              self.data['stress'][j,0] = fields[3] 
              self.data['stress'][j,1] = fields[4] 
              self.data['stress'][j,2] = fields[5]
              self.data['stress_sum'] = self.data['stress_sum'] + abs(fields[0]) + abs(fields[1]) + abs(fields[2])
            
#                  "stress": numpy.zeros((3,3)),
#      "stress_sum": None,
            
          if(pwscf_output.compare(line, "density = ")):
            self.data['density_full'] = pwscf_output.extract(line, "=", "", "s")
            self.data['density'] = pwscf_output.extract(line, "=", "g/cm^3", "f")
          
          if(pwscf_output.compare(line, "PWSCF        :")):
            self.data['cpu_time'] = pwscf_output.extract(line, ":", "CPU", "s")
            
          if(pwscf_output.compare(line, "PWSCF        :")):
            self.data['wall_time'] = pwscf_output.extract(line, "CPU", "WALL", "s")
          
          if(pwscf_output.compare(line, "JOB DONE.")):
            self.data['job_done'] = True
            
          if(pwscf_output.compare(line, "Exit code:")):
            self.data['error'] = True  
"""

###########################################
#  CLASS e_adjus
###########################################
class e_adjust:

  @staticmethod
  def load():
    if('dft' not in g.inp.keys()):
      return None
    if('e_adjust' not in g.inp['dft'].keys()):
      return None
      
    dft_file = g.inp['dft']['e_adjust']
    if(os.path.isfile(dft_file)): 
# Load csv to list
      fd = std.read_csv(dft_file)
      
      for row in fd:
        label_text = row[0].strip().upper()
        if(label_text != ''):
          label_str, label_id = labels.add(label_text)
          try:        
            relaxed_dft_ev = units.convert(row[3], "EV", float(row[2]) / int(row[1])) # relaxed per atom in eV
            coh_ev = units.convert(row[5], "EV", float(row[4]))
            apaev = coh_ev - relaxed_dft_ev # Adjustment per atom ev
            g.dft_energy_adjustments[label_id] = {
                                                  'label_id': label_id,
                                                  'label_text': label_str,
                                                  'atom_count': int(row[1]),
                                                  'relaxed_energy': float(row[2]),
                                                  'relaxed_energy_unit': row[3],
                                                  'coh_energy': float(row[4]),
                                                  'coh_energy_unit': row[5],
                                                  'calc_relaxed_dft_ev': relaxed_dft_ev,
                                                  'calc_coh_ev': coh_ev,
                                                  'calc_apaev': apaev,
                                                 }
          except:
            pass  
          
      """
      g.rd['dft_energy_adjustments'] = {}
      
      dft_energy_adjustments
      for row in fd:
        label_text = row[0].strip().upper()
        if(label_text != '' and label_text not in g.rd['dft_energy_adjustments'].keys()):
          try:            
            relaxed_dft_ev = units.convert(row[3], "EV", float(row[2]) / int(row[1])) # relaxed per atom in eV
            coh_ev = units.convert(row[5], "EV", float(row[4]))
            apaev = coh_ev - relaxed_dft_ev # Adjustment per atom ev
            g.rd['dft_energy_adjustments'][label_text] = {
                                                      'label': 0,
                                                      'label_text': label_text,
                                                      'atom_count': int(row[1]),
                                                      'relaxed_energy': float(row[2]),
                                                      'relaxed_energy_unit': row[3],
                                                      'coh_energy': float(row[4]),
                                                      'coh_energy_unit': row[5],
                                                      'calc_relaxed_dft_ev': relaxed_dft_ev,
                                                      'calc_coh_ev': coh_ev,
                                                      'calc_apaev': apaev,
                                                     }
          except:
            pass
      """

###########################################
#  CLASS unit
###########################################
class units:
  
  @staticmethod
  def convert(conv_from, conv_to, value_in):
    try:
      value_in = float(value_in)
    except:
      return None
    conv_from = conv_from.upper()
    conv_to = conv_to.upper()

# LENGTH METERS
    length = {
    'M': 1.0,
    'CM': 100,
    'MM': 1E3,
    'UM': 1E6,
    'NM': 1E9,
    'ANG': 1E10,
    'BOHR': 1.89E10,
    }

# ENERGY J
    energy = {
    'J': 1.0,
    'EV': 6.2415E18,
    'RY': 4.5874E17,
    }

# FORCE N
    force = {
    'N': 1.0,
    'RY/BOHR': 2.4276e7,
    'EV/ANG':6.2414E8,    
    }
    
# VELOCITY
    velocity = {
    'M/S': 1.0,
    'MPH': 2.25,    
    }
    
# PRESSURE
    pressure = {
    'PA': 1.0,
    'GPA': 1.0E-9,    
    'BAR': 1.0E-5,    
    'ATMOSPHERE': 9.8692E-6,    
    'PSI': 1.45038E-4, 
    'KBAR': 1.0E-8,   
    'RY/BOHR3': 6.857E-14,   
    'EV/ANG3': 6.241E-12
    }
    
# CHARGE DENSITY (UNIT CHARGE PER VOLUME - ANG^3)
    charge_density = {
    'ANG-3': 1.0,
    'BOHR-3': 0.14812,    
    }
    
# TEMPERATURE
    
    unit_list = [length, energy, force, velocity, pressure, charge_density]
    
    for l in unit_list:
      if(conv_from in l.keys() and conv_to in l.keys()):
        return round((l[conv_to] / l[conv_from]) * float(value_in),9)
  
"""  
  @staticmethod
  def convert(conv_from, conv_to, value_in):

    conv_from = conv_from.upper()
    conv_to = conv_to.upper()

# METERS
    length = {
    'M': 1.0,
    'CM': 100,
    'MM': 1E3,
    'UM': 1E6,
    'NM': 1E9,
    'ANG': 1E10,
    'BOHR': 1.89E10,
    }

# J
    energy = {
    'J': 1.0,
    'EV': 6.242E18,
    'RY': 4.5874E17,
    }

    if(conv_from in length.keys() and conv_to in length.keys()):
      return round((length[conv_to] / length[conv_from]) * float(value_in),9)

    if(conv_from in energy.keys() and conv_to in energy.keys()):
      return round((energy[conv_to] / energy[conv_from]) * float(value_in),9)
"""

###########################################
#  CLASS b_prop
###########################################
class b_props:

  @staticmethod
  def load():
    if('bp' not in g.inp.keys()):
      return None
    if('bp_file' not in g.inp['bp'].keys()):
      return None
    
# Read BP data file
    bp_inp = read_config.read_file(g.inp['bp']['bp_file'])
    
# READ IN UNITS
    try:
      bp_pressure = bp_inp['units']['pressure']
    except:
      bp_pressure = 'GPA'
    try:
      bp_length = bp_inp['units']['length']
    except:
      bp_length = 'ang'
    try:
      bp_energy = bp_inp['units']['energy']
    except:
      bp_energy = 'ev'
      
    for k in bp_inp.keys():
      if('alat' in bp_inp[k].keys()):
        label_str, label_id = labels.add(k)        
        g.bulk_properties[label_id] = b_props.make(label_id, label_str)

# There must be an alat value set
        g.bulk_properties[label_id]['alat'] = float(bp_inp[k]['alat'])
                
        try:
          if(bp_inp[k]['type'].lower() == 'sc'):
            g.bulk_properties[label_id]['type'] = 1
          elif(bp_inp[k]['type'].lower() == 'bcc'):
            g.bulk_properties[label_id]['type'] = 2
          elif(bp_inp[k]['type'].lower() == 'fcc'):
            g.bulk_properties[label_id]['type'] = 3
          elif(bp_inp[k]['type'].lower() == 'zb'):
            g.bulk_properties[label_id]['type'] = 4
        except:
          pass 
        
        try: 
          g.bulk_properties[label_id]['expansion'] = float(bp_inp[k]['expansion'])
        except:
          pass         
        try: 
          g.bulk_properties[label_id]['rcut'] = float(bp_inp[k]['rcut'])
        except:
          pass      
        
        try:
          g.bulk_properties[label_id]['amu_per_crystal'] = float(bp_inp[k]['amu_per_crystal'])
        except:
          pass  
        
        try:
          g.bulk_properties[label_id]['b0'] = float(bp_inp[k]['b0'])
        except:
          pass  
        try:
          g.bulk_properties[label_id]['e0'] = float(bp_inp[k]['e0'])
        except:
          pass        
        try:
          if('ec' in bp_inp[k].keys()):
            g.bulk_properties[label_id]['ec'] = numpy.zeros((6,6,),)
            
# CUBIC  C11 C12 C44
            if(len(bp_inp[k]['ec']) == 3):
              c11 = float(bp_inp[k]['ec'][0])
              c12 = float(bp_inp[k]['ec'][1])
              c44 = float(bp_inp[k]['ec'][2])
              
              g.bulk_properties[label_id]['ec'][0,0] = c11
              g.bulk_properties[label_id]['ec'][1,1] = c11
              g.bulk_properties[label_id]['ec'][2,2] = c11
              
              g.bulk_properties[label_id]['ec'][0,1] = c12
              g.bulk_properties[label_id]['ec'][0,2] = c12
              g.bulk_properties[label_id]['ec'][1,2] = c12
              g.bulk_properties[label_id]['ec'][1,0] = c12
              g.bulk_properties[label_id]['ec'][2,0] = c12
              g.bulk_properties[label_id]['ec'][2,1] = c12
              
              g.bulk_properties[label_id]['ec'][3,3] = c44
              g.bulk_properties[label_id]['ec'][4,4] = c44
              g.bulk_properties[label_id]['ec'][5,5] = c44
              
# ORTHORHOMBIC  C11 C22 C33 C44 C55 C66 C12 C13 C23
            if(len(bp_inp[k]['ec']) == 9):
              c11 = float(bp_inp[k]['ec'][0])
              c22 = float(bp_inp[k]['ec'][1])
              c33 = float(bp_inp[k]['ec'][2])
              c44 = float(bp_inp[k]['ec'][3])
              c55 = float(bp_inp[k]['ec'][4])
              c66 = float(bp_inp[k]['ec'][5])
              c12 = float(bp_inp[k]['ec'][6])
              c13 = float(bp_inp[k]['ec'][7])
              c23 = float(bp_inp[k]['ec'][8])
              
              g.bulk_properties[label_id]['ec'][0,0] = c11
              g.bulk_properties[label_id]['ec'][1,1] = c22
              g.bulk_properties[label_id]['ec'][2,2] = c33
              
              g.bulk_properties[label_id]['ec'][3,3] = c44
              g.bulk_properties[label_id]['ec'][4,4] = c55
              g.bulk_properties[label_id]['ec'][5,5] = c66
              
              g.bulk_properties[label_id]['ec'][0,1] = c12
              g.bulk_properties[label_id]['ec'][0,2] = c13
              g.bulk_properties[label_id]['ec'][1,2] = c23
              g.bulk_properties[label_id]['ec'][1,0] = c12
              g.bulk_properties[label_id]['ec'][2,0] = c13
              g.bulk_properties[label_id]['ec'][2,1] = c23
        except:
          pass
          
# CONVERT TO ev/angs
    for k in g.bulk_properties.keys():
      try: 
        g.bulk_properties[k]['alat'] = units.convert(bp_length, 'ang', g.bulk_properties[k]['alat'])
      except:
        pass  
      try: 
        g.bulk_properties[k]['e0'] = units.convert(bp_energy, 'EV', g.bulk_properties[k]['e0'])
      except:
        pass  
      try: 
        g.bulk_properties[k]['b0'] = units.convert(bp_pressure, 'EV/ANG3', g.bulk_properties[k]['b0'])
      except:
        pass  
      try: 
        g.bulk_properties[k]['g'] = units.convert(bp_pressure, 'EV/ANG3', g.bulk_properties[k]['g'])
      except:
        pass  
      try: 
        g.bulk_properties[k]['e'] = units.convert(bp_pressure, 'EV/ANG3', g.bulk_properties[k]['e'])
      except:
        pass  
      try: 
        for i in range(6):
          for j in range(6):
            g.bulk_properties[k]['ec'][i,j] = units.convert(bp_pressure, 'EV/ANG3', g.bulk_properties[k]['ec'][i,j])
      except:
        pass  

  def make(label_id, label_str):
    bp_d ={
           'label_id': label_id,
           'label_str': label_str,
           'alat': None,
           'uv': numpy.zeros((3,3,),),
           'type': 1,
           'expansion': 4,
           'rcut': 6.5,
           'b0': None,
           'e0': None,
           'ec': None,
           'g': None,
           'e': None,
           'amu_per_crystal': None,
          }
    bp_d['uv'][:,:] = 0.0
    bp_d['uv'][0,0] = 1.0
    bp_d['uv'][1,1] = 1.0
    bp_d['uv'][2,2] = 1.0
    return bp_d
           
###########################################################
# F2PY functions
###########################################################
  
  @staticmethod
  def bp_add():  
    for bp_n in g.bulk_properties:
#bp_id = bp.add_fcc(6.5, g.bulk_properties[bp_n]['alat'], 1)
#add_bp_config(rcut_in, alat_in, uv_in, label_in, crystal_type_in, expansion_in, bp_id)
      rcut = g.bulk_properties[bp_n]['rcut']
      alat = g.bulk_properties[bp_n]['alat']
      uv = g.bulk_properties[bp_n]['uv']
      label = g.bulk_properties[bp_n]['label_id']
      type = g.bulk_properties[bp_n]['type']
      expansion = g.bulk_properties[bp_n]['expansion']
      
# Add Config
      bp_id = bp.add_bp_config(rcut, alat, uv, label, type, expansion)
    
# Add known data
      bp.add_alat(bp_id, g.bulk_properties[bp_n]['alat'])
      bp.add_e0(bp_id, g.bulk_properties[bp_n]['e0'])
      bp.add_b0(bp_id, g.bulk_properties[bp_n]['b0'])
      bp.add_ec(bp_id, g.bulk_properties[bp_n]['ec'])
      bp.add_amu_per_crystal(bp_id, g.bulk_properties[bp_n]['amu_per_crystal'])
    
# Add rss
    try:
      bp.set_rss_alat(g.inp['rss']['alat'])
    except:
      pass
    try:
      bp.set_rss_e0(g.inp['rss']['e0'])
    except:
      pass
    try:
      bp.set_rss_b0(g.inp['rss']['b0'])
    except:
      pass
    try:
      bp.set_rss_ec(g.inp['rss']['ec'])
    except:
      pass
    try:
      bp.set_rss_g(g.inp['rss']['g'])
    except:
      pass
    try:
      bp.set_rss_e(g.inp['rss']['e'])
    except:
      pass
    try:
      bp.set_rss_v(g.inp['rss']['v'])
    except:
      pass

  @staticmethod
  def bp_output():  
  
    for bp_id in range(bp.bp_configs_count):
      fh = open(g.dirs['results'] + '/' + 'bp_' + str(bp_id + 1) + '.dat', 'w')

      t_pad = 30
      f_pad = 18
      
      std.write_file_line(fh, '######################################################', 1, '', 1)
      std.write_file_line(fh, 'Known Properties', 1, '', 1)
      std.write_file_line(fh, '######################################################', 1, '', 1)
      std.write_file_line(fh, '', 1, '', 1)
      std.write_file_line(fh, 'All units are in ev/Ang unless specified', 1, '', 1)
      std.write_file_line(fh, 'Energy: eV', 1, '', 1)
      std.write_file_line(fh, 'Length: ang', 1, '', 1)
      std.write_file_line(fh, 'Force: eV/ang', 1, '', 1)
      std.write_file_line(fh, 'Pressure: eV/ang3', 1, '', 1)
      std.write_file_line(fh, '', 1, '', 1)
      
      std.write_file_line(fh, 'Atoms per crystal:', t_pad, bp.known_atoms_per_crystal[bp_id], f_pad)
      std.write_file_line(fh, 'Expansion:', t_pad, bp.known_expansion[bp_id], f_pad)
      std.write_file_line(fh, '', 1, '', 1)
      std.write_file_line(fh, 'Equation of State', 1, '', 1)
      std.write_file_line(fh, '#################', 1, '', 1)
      std.write_file_line(fh, 'alat:', t_pad, bp.known_alat[bp_id], f_pad)
      std.write_file_line(fh, 'e0:', t_pad, bp.known_e0[bp_id], f_pad)
      std.write_file_line(fh, 'b0:', t_pad, bp.known_b0[bp_id], f_pad)
      std.write_file_line(fh, 'b0/GPA:', t_pad, 160.230732254e0 * bp.known_b0[bp_id], f_pad)
      std.write_file_line(fh, '', 1, '', 1)
      std.write_file_line(fh, 'Stiffness Matrix', 1, '', 1)
      std.write_file_line(fh, '#################', 1, '', 1)
      std.write_file_line(fh, 'Stiffness:', t_pad, bp.known_ec[bp_id,0,:], f_pad)
      for i in range(1,6):   
        std.write_file_line(fh, '', t_pad, bp.known_ec[bp_id,i,:], f_pad)
      std.write_file_line(fh, '', 1, '', 1)
      
      std.write_file_line(fh, 'Stiffness (GPA):', t_pad, 160.230732254e0 * bp.known_ec[bp_id,0,:], f_pad)
      for i in range(1,6):   
        std.write_file_line(fh, '', t_pad, 160.230732254e0 * bp.known_ec[bp_id,i,:], f_pad)
        
      std.write_file_line(fh, '', 1, '', 1)
      std.write_file_line(fh, '', 1, '', 1)
      
      std.write_file_line(fh, '######################################################', 1, '', 1)
      std.write_file_line(fh, 'Calculated Properties', 1, '', 1)
      std.write_file_line(fh, '######################################################', 1, '', 1)
      std.write_file_line(fh, '', 1, '', 1)
      
      std.write_file_line(fh, 'Equation of State', 1, '', 1)
      std.write_file_line(fh, '#################', 1, '', 1)
      std.write_file_line(fh, 'alat:', t_pad, bp.calc_alat[bp_id], f_pad)
      std.write_file_line(fh, 'v0:', t_pad, bp.calc_v0[bp_id], f_pad)
      std.write_file_line(fh, 'e0:', t_pad, bp.calc_e0[bp_id], f_pad)
      std.write_file_line(fh, 'b0:', t_pad, bp.calc_b0[bp_id], f_pad)
      std.write_file_line(fh, 'b0/GPA:', t_pad, 160.230732254e0 * bp.calc_b0[bp_id], f_pad)
      std.write_file_line(fh, '', 1, '', 1)  
      
      std.write_file_line(fh, 'Stiffness Matrix', 1, '', 1)
      std.write_file_line(fh, '#################', 1, '', 1)
      std.write_file_line(fh, 'Stiffness:', t_pad, bp.calc_ec[bp_id,0,:], f_pad)
      for i in range(1,6):   
        std.write_file_line(fh, '', t_pad, bp.calc_ec[bp_id,i,:], f_pad)
      std.write_file_line(fh, '', 1, '', 1)
      
      std.write_file_line(fh, 'Stiffness (GPA):', t_pad, 160.230732254e0 * bp.calc_ec[bp_id,0,:], f_pad)
      for i in range(1,6):   
        std.write_file_line(fh, '', t_pad, 160.230732254e0 * bp.calc_ec[bp_id,i,:], f_pad)   
      std.write_file_line(fh, '', 1, '', 1)
      
      std.write_file_line(fh, 'Compliance Matrix', 1, '', 1)
      std.write_file_line(fh, '#################', 1, '', 1)
      std.write_file_line(fh, 'Compliance:', t_pad, bp.calc_sc[bp_id,0,:], f_pad)
      for i in range(1,6):   
        std.write_file_line(fh, '', t_pad, bp.calc_sc[bp_id,i,:], f_pad)
      std.write_file_line(fh, '', 1, '', 1)
      
      std.write_file_line(fh, 'Compliance (1/GPA):', t_pad, bp.calc_sc_gpa[bp_id,0,:], f_pad)
      for i in range(1,6):   
        std.write_file_line(fh, '', t_pad, bp.calc_sc_gpa[bp_id,i,:], f_pad)
      std.write_file_line(fh, '', 1, '', 1)
          
      std.write_file_line(fh, 'Bulk Modulus', 1, '', 1)
      std.write_file_line(fh, '#################', 1, '', 1)
      std.write_file_line(fh, 'b0 reuss:', t_pad, bp.calc_b0_r[bp_id], f_pad)
      std.write_file_line(fh, 'b0 voight:', t_pad, bp.calc_b0_v[bp_id], f_pad)
      std.write_file_line(fh, 'b0 avg:', t_pad, bp.calc_b0_avg[bp_id], f_pad)
      std.write_file_line(fh, 'b0 reuss (GPA):', t_pad, 160.230732254e0 * bp.calc_b0_r[bp_id], f_pad)
      std.write_file_line(fh, 'b0 voight (GPA):', t_pad, 160.230732254e0 * bp.calc_b0_v[bp_id], f_pad)
      std.write_file_line(fh, 'b0 avg (GPA):', t_pad, 160.230732254e0 * bp.calc_b0_avg[bp_id], f_pad)
      std.write_file_line(fh, '', 1, '', 1)
      
      std.write_file_line(fh, 'Shear Modulus', 1, '', 1)
      std.write_file_line(fh, '#################', 1, '', 1)
      std.write_file_line(fh, 'G reuss:', t_pad, bp.calc_g_r[bp_id], f_pad)
      std.write_file_line(fh, 'G voight:', t_pad, bp.calc_g_v[bp_id], f_pad)
      std.write_file_line(fh, 'G avg:', t_pad, bp.calc_g_avg[bp_id], f_pad)
      std.write_file_line(fh, 'G reuss (GPA):', t_pad, 160.230732254e0 * bp.calc_g_r[bp_id], f_pad)
      std.write_file_line(fh, 'G voight (GPA):', t_pad, 160.230732254e0 * bp.calc_g_v[bp_id], f_pad)
      std.write_file_line(fh, 'G avg (GPA):', t_pad, 160.230732254e0 * bp.calc_g_avg[bp_id], f_pad)
      std.write_file_line(fh, '', 1, '', 1)  
      
      std.write_file_line(fh, 'Young Modulus', 1, '', 1)
      std.write_file_line(fh, '#################', 1, '', 1)
      std.write_file_line(fh, 'E:', t_pad, bp.calc_e[bp_id], f_pad)
      std.write_file_line(fh, 'E vec:', t_pad, bp.calc_e_vec[bp_id, :], f_pad)
      std.write_file_line(fh, '', 1, '', 1) 
      
      std.write_file_line(fh, 'Poisson Ratio', 1, '', 1)
      std.write_file_line(fh, '#################', 1, '', 1)
      std.write_file_line(fh, 'v:', t_pad, bp.calc_v[bp_id], f_pad)
      std.write_file_line(fh, '', 1, '', 1) 
      
      std.write_file_line(fh, 'Temperatures', 1, '', 1)
      std.write_file_line(fh, '#################', 1, '', 1)
      std.write_file_line(fh, 'Melting (K):', t_pad, bp.calc_melting[bp_id], f_pad)
#std.write_file_line(fh, 'Debye:', t_pad, bp.calc_debye[bp_id], f_pad)   # Check, calc might be wrong
      std.write_file_line(fh, '', 1, '', 1) 
             
      std.write_file_line(fh, '', 1, '', 1)
      std.write_file_line(fh, '', 1, '', 1)
      
      std.write_file_line(fh, '######################################################', 1, '', 1)
      std.write_file_line(fh, 'RSS', 1, '', 1)
      std.write_file_line(fh, '######################################################', 1, '', 1)
      std.write_file_line(fh, '', 1, '', 1) 
      if(bp.known_set[bp_id,0]):
        std.write_file_line(fh, 'Alat rss:', t_pad, bp.rss[bp_id,0], f_pad)
      if(bp.known_set[bp_id,1]):
        std.write_file_line(fh, 'e0 rss:', t_pad, bp.rss[bp_id,1], f_pad)
      if(bp.known_set[bp_id,2]):
        std.write_file_line(fh, 'b0 rss:', t_pad, bp.rss[bp_id,2], f_pad)
      if(bp.known_set[bp_id,3]):
        std.write_file_line(fh, 'Stiffness rss:', t_pad, bp.rss[bp_id,3], f_pad)
      if(bp.known_set[bp_id,4]):
        std.write_file_line(fh, 'G rss:', t_pad, bp.rss[bp_id,4], f_pad)
      if(bp.known_set[bp_id,5]):
        std.write_file_line(fh, 'E rss:', t_pad, bp.rss[bp_id,5], f_pad)
      if(bp.known_set[bp_id,6]):
        std.write_file_line(fh, 'Poisson rss:', t_pad, bp.rss[bp_id,6], f_pad)
      std.write_file_line(fh, 'Total rss for ' +str(bp_id+1) + ':', t_pad, bp.rss_total[bp_id], f_pad)
      std.write_file_line(fh, 'Total rss for All:', t_pad, [bp.rss_total_rss], f_pad)
        
      std.write_file_line(fh, '', 1, '', 1)
      std.write_file_line(fh, '', 1, '', 1)
      
      std.write_file_line(fh, '######################################################', 1, '', 1)
      std.write_file_line(fh, 'RSS Weighted', 1, '', 1)
      std.write_file_line(fh, '######################################################', 1, '', 1)
      std.write_file_line(fh, '', 1, '', 1) 
      if(bp.known_set[bp_id,0]):
        std.write_file_line(fh, 'Alat rss:', t_pad, bp.rss_w[bp_id,0], f_pad)
      if(bp.known_set[bp_id,1]):
        std.write_file_line(fh, 'e0 rss:', t_pad, bp.rss_w[bp_id,1], f_pad)
      if(bp.known_set[bp_id,2]):
        std.write_file_line(fh, 'b0 rss:', t_pad, bp.rss_w[bp_id,2], f_pad)
      if(bp.known_set[bp_id,3]):
        std.write_file_line(fh, 'Stiffness rss:', t_pad, bp.rss_w[bp_id,3], f_pad)
      if(bp.known_set[bp_id,4]):
        std.write_file_line(fh, 'G rss:', t_pad, bp.rss_w[bp_id,4], f_pad)
      if(bp.known_set[bp_id,5]):
        std.write_file_line(fh, 'E rss:', t_pad, bp.rss_w[bp_id,5], f_pad)
      if(bp.known_set[bp_id,6]):
        std.write_file_line(fh, 'Poisson rss:', t_pad, bp.rss_w[bp_id,6], f_pad)
      std.write_file_line(fh, 'Total rss for ' +str(bp_id+1) + ':', t_pad, bp.rss_total_w[bp_id], f_pad)
      std.write_file_line(fh, 'Total rss for All:', t_pad, [bp.rss_total_rss_w], f_pad)
      
      std.write_file_line(fh, '', 1, '', 1)
      std.write_file_line(fh, '', 1, '', 1)
      
      fh.close()

  @staticmethod
  def bp_output_terminal():  
  
    for bp_id in range(bp.bp_configs_count):
      fh = open(g.dirs['results'] + '/' + 'bp_' + str(bp_id + 1) + '.dat', 'w')

      t_pad = 30
      f_pad = 18
      
      std.print_file_line('######################################################', 1, '', 1)
      std.print_file_line('Known Properties', 1, '', 1)
      std.print_file_line('######################################################', 1, '', 1)
      std.print_file_line('', 1, '', 1)
      std.print_file_line('All units are in ev/Ang unless specified', 1, '', 1)
      std.print_file_line('Energy: eV', 1, '', 1)
      std.print_file_line('Length: ang', 1, '', 1)
      std.print_file_line('Force: eV/ang', 1, '', 1)
      std.print_file_line('Pressure: eV/ang3', 1, '', 1)
      std.print_file_line('', 1, '', 1)
      
      std.print_file_line('Atoms per crystal:', t_pad, bp.known_atoms_per_crystal[bp_id], f_pad)
      std.print_file_line('Expansion:', t_pad, bp.known_expansion[bp_id], f_pad)
      std.print_file_line('', 1, '', 1)
      std.print_file_line('Equation of State', 1, '', 1)
      std.print_file_line('#################', 1, '', 1)
      std.print_file_line('alat:', t_pad, bp.known_alat[bp_id], f_pad)
      std.print_file_line('e0:', t_pad, bp.known_e0[bp_id], f_pad)
      std.print_file_line('b0:', t_pad, bp.known_b0[bp_id], f_pad)
      std.print_file_line('b0/GPA:', t_pad, 160.230732254e0 * bp.known_b0[bp_id], f_pad)
      std.print_file_line('', 1, '', 1)
      std.print_file_line('Stiffness Matrix', 1, '', 1)
      std.print_file_line('#################', 1, '', 1)
      std.print_file_line('Stiffness:', t_pad, bp.known_ec[bp_id,0,:], f_pad)
      for i in range(1,6):   
        std.print_file_line('', t_pad, bp.known_ec[bp_id,i,:], f_pad)
      std.print_file_line('', 1, '', 1)
      
      std.print_file_line('Stiffness (GPA):', t_pad, 160.230732254e0 * bp.known_ec[bp_id,0,:], f_pad)
      for i in range(1,6):   
        std.print_file_line('', t_pad, 160.230732254e0 * bp.known_ec[bp_id,i,:], f_pad)
        
      print()
      print()

      std.print_file_line('######################################################', 1, '', 1)
      std.print_file_line('Calculated Properties', 1, '', 1)
      std.print_file_line('######################################################', 1, '', 1)
      std.print_file_line('', 1, '', 1)
      
      std.print_file_line('Equation of State', 1, '', 1)
      std.print_file_line('#################', 1, '', 1)
      std.print_file_line('alat:', t_pad, bp.calc_alat[bp_id], f_pad)
      std.print_file_line('v0:', t_pad, bp.calc_v0[bp_id], f_pad)
      std.print_file_line('e0:', t_pad, bp.calc_e0[bp_id], f_pad)
      std.print_file_line('b0:', t_pad, bp.calc_b0[bp_id], f_pad)
      std.print_file_line('b0/GPA:', t_pad, 160.230732254e0 * bp.calc_b0[bp_id], f_pad)
      std.print_file_line('', 1, '', 1)  
      
      std.print_file_line('Stiffness Matrix', 1, '', 1)
      std.print_file_line('#################', 1, '', 1)
      std.print_file_line('Stiffness:', t_pad, bp.calc_ec[bp_id,0,:], f_pad)
      for i in range(1,6):   
        std.print_file_line('', t_pad, bp.calc_ec[bp_id,i,:], f_pad)
      std.print_file_line('', 1, '', 1)
        
  @staticmethod
  def bp_eos_plot():       
  
    for bp_id in range(bp.bp_configs_count):
    
# EQUATION OF STATE
    
      s = bp.calc_sizes[bp_id, 0]
    
      plt.clf()
    
      plt.rc('font', family='serif')
      plt.rc('xtick', labelsize='x-small')
      plt.rc('ytick', labelsize='x-small')

      fig, axs = plt.subplots(1, 1, figsize=(12,9))
      fig.tight_layout(pad=5.0)
      fig.suptitle('Equation of State')  
      
      plt.xlabel('Volume (ang3)')    
      plt.ylabel('Energy (eV)')
          
      plt.plot(bp.calc_volumes[bp_id, 0, 0:s], bp.calc_energies[bp_id, 0, 0:s], color='k',  marker="x", ls='')
      plt.plot(bp.calc_volumes[bp_id, 0, 0:s], bp.calc_energies_fit[bp_id, 0, 0:s], color='k', ls='solid')

      plt.savefig(g.dirs['eos'] + '/' + 'eos.svg')
      plt.savefig(g.dirs['eos'] + '/' + 'eos.eps')
      
# ELASTIC CONSTANTS
      
      plt.clf()
    
      plt.rc('font', family='serif')
      plt.rc('xtick', labelsize='x-small')
      plt.rc('ytick', labelsize='x-small')

      fig, axs = plt.subplots(3, 3, figsize=(12,9))
      fig.tight_layout(pad=5.0)
      fig.suptitle('Elastic Constant Curves')    
    
      for dn in range(9):
        s = bp.calc_sizes[bp_id, dn + 1]
    
        axs[int(numpy.floor(dn/3)), dn % 3].plot(bp.calc_strains[bp_id, dn + 1, 0:s],
                                                 bp.calc_energies[bp_id, dn + 1, 0:s],
                                                 color='k',  marker="x", ls='')
        axs[int(numpy.floor(dn/3)), dn % 3].plot(bp.calc_strains[bp_id, dn + 1, 0:s],
                                                 bp.calc_energies_fit[bp_id, dn + 1, 0:s],
                                                 color='k', ls='solid')
        axs[int(numpy.floor(dn/3)), dn % 3].set_title('Distortion D' + str(dn + 1))
        axs[int(numpy.floor(dn/3)), dn % 3].set_xlabel('Strain (Expanded Alat)')
        axs[int(numpy.floor(dn/3)), dn % 3].set_ylabel('Energy (eV)')
               
      plt.savefig(g.dirs['ec'] + '/' + 'ec.svg')
      plt.savefig(g.dirs['ec'] + '/' + 'ec.eps')
      
    """
    n = globals.d['eos_data_size']
    x = np.linspace(globals.d['eos_data'][1,0], globals.d['eos_data'][1,n-1], 101)    
    y = np.zeros((101,),)
    p = numpy.zeros((4,),)
    p[:] = globals.d['eos_fitting'][:]
    y[:] = run_eos.bm_calc(x[:], p)
    
    n = globals.d['eos_data_size'] 
    plt.plot(globals.d['eos_data'][1,:n], globals.d['eos_data'][2,:n], color='k',  marker="x", ls='')
    plt.plot(x[:], y[:], color='k', ls='solid')
    """
  
###########################################
#  CLASS efs_cal
###########################################
class efs_calc:

  def run_energy():
  
    print("Calc Energy") 
  
# Setup EFS
    efs.init()                         # Initialise (allocate arrays)
    potential.efs_add_potentials()     # Load potentials
    configs.efs_add_config()           # Add configs
    efs.energy()
    efs_calc.output_energy()
    potential.plot_fortran_potentials()
    potential.plot_python_potentials()
    
  def run_energy_force():
  
    print("Calc Energy and Forces") 
    
# Setup EFS
    efs.init()                         # Initialise (allocate arrays)
    potential.efs_add_potentials()     # Load potentials
    configs.efs_add_config()           # Add configs
    efs.energy_force() 
    efs_calc.output_energy()
    efs_calc.output_forces()
    potential.plot_fortran_potentials()
    potential.plot_python_potentials()
  
  def run_energy_force_stress():
  
    print("Calc Energy, Forces and Stress") 
    
# Setup EFS
    efs.init()                         # Initialise (allocate arrays)
    potential.efs_add_potentials()     # Load potentials
    configs.efs_add_config()           # Add configs
    efs.energy_force_stress() 
    efs_calc.output_energy()
    efs_calc.output_forces()
    efs_calc.output_stress()
    potential.plot_fortran_potentials()
    potential.plot_python_potentials()
  
#efs.max_density_calc()

  def output_energy():
  
    t_pad = 12
    f_pad = 18
  
    fh = open(g.dirs['results'] + '/' + 'config_energies.txt', 'w')
    fh.write('###############################################################')
    fh.write('###############################################################\n')
    fh.write('ENERGY RESULTS\n')
    fh.write('Config Count: ' + str(efs.cc) + '\n')
    fh.write('###############################################################')
    fh.write('###############################################################\n')
    for n in range(efs.cc):    
      std.write_file_line(fh, 'Config ' + str(n+1) + ':', t_pad, efs.config_energy[n,:], f_pad)
    fh.write('###############################################################')
    fh.write('###############################################################\n')
    fh.close()
  
  def output_forces():
  
    t_pad = 12
    f_pad = 18
  
    fh = open(g.dirs['results'] + '/' + 'config_forces.txt', 'w')
    fh.write('###############################################################')
    fh.write('###############################################################\n')
    fh.write('FORCE RESULTS\n')
    fh.write('Config Count: ' + str(efs.cc) + '\n')
    fh.write('###############################################################')
    fh.write('###############################################################\n')
    fh.write('\n')
    for n in range(efs.cc): 
      fh.write('##################\n')
      fh.write('Config ' + str(n) + '\n')
      fh.write('##################\n')
      a = efs.key[n, 0] - 1
      b = efs.key[n, 1]
      for l in range(a, b):
        std.write_file_line(fh, str(efs.labels[l]) + ':', t_pad, efs.config_forces[l,:], f_pad)
      fh.write('\n')
     
    fh.write('\n')
    fh.write('###############################################################')
    fh.write('###############################################################\n')
    fh.close()

  def output_stress():
  
    t_pad = 12
    f_pad = 18
  
    fh = open(g.dirs['results'] + '/' + 'config_stresses.txt', 'w')
    fh.write('###############################################################')
    fh.write('###############################################################\n')
    fh.write('STRESS RESULTS\n')
    fh.write('Config Count: ' + str(efs.cc) + '\n')
    fh.write('###############################################################')
    fh.write('###############################################################\n')
    for n in range(efs.cc):    
      a = n * 3
      std.write_file_line(fh, 'Config ' + str(n+1) + ':', t_pad, efs.config_stresses[a,:], f_pad)
      std.write_file_line(fh, '', t_pad, efs.config_stresses[a+1,:], f_pad)
      std.write_file_line(fh, '', t_pad, efs.config_stresses[a+2,:], f_pad)
      fh.write('\n')
    fh.write('###############################################################')
    fh.write('###############################################################\n')
    fh.close()
  
  def output_rss():
    fh = open(g.dirs['results'] + '/' + 'rss_configs.txt', 'w')
    fh.write('###############################################################\n')
    fh.write('ENERGY RESULTS\n')
    fh.write('Config Count: ' + str(efs.cc) + '\n')
    fh.write('###############################################################\n')
    for n in range(efs.cc):
      fh.write(str(efs.config_energy[n,0]) + ' ')
      fh.write(str(efs.config_energy[n,1]) + ' ')
      fh.write(str(efs.config_energy[n,2]) + ' ')
      fh.write(str(efs.config_energy[n,3]) + ' ')
      fh.write(str(efs.config_energy[n,4]) + ' ')
      fh.write(str(efs.config_energy[n,5]) + ' ')
      fh.write('\n')
    fh.write('###############################################################\n')
    fh.close()

###########################################
#  CLASS bp_cal
###########################################
class bp_calc:

  def run():  
  
    print("Calc Bulk Properties") 
    
# Setup BP
    bp.init()
    potential.bp_add_potentials()
    b_props.bp_add()
    
# Calculate
    bp.energy()
    bp.calculate_bp()    
    
# Output to Terminal
    b_props.bp_output_terminal()
    
# Output to File
    b_props.bp_output()
    
# Plots
    b_props.bp_eos_plot()
    potential.plot_fortran_potentials()
    potential.plot_python_potentials()
    
###########################################
#  CLASS es_cal
###########################################
class es_calc:

  def run():  
  
    print("Energies: Surface, Vacany etc") 

# Setup ES
    es.init()
    potential.es_add_potentials()
    es_calc.es_add()
    es.energy()
    for i in range(es.cc):
      print(es.config_energy[i,:])

  @staticmethod
  def es_add():  
# Test values
    
    rcut = 6.5
    alat = 4.04
    label = 1
    type = 3
      
# Add Config
    es_id = es.add_es_config(rcut, alat, label, type)

###########################################
#  CLASS rss_cal
###########################################
class rss_calc:

  def run():
    print("Calc RSS") 
    
# Setup EFS
    efs.init()                         # Initialise (allocate arrays)
    potential.efs_add_potentials()     # Load potentials
    configs.efs_add_config()           # Add configs
    
# Setup BP
    bp.init()
    potential.bp_add_potentials()
    b_props.bp_add()

    rss = rss_calc.get_rss()
    
# Output to File
    efs_calc.output_energy()
    efs_calc.output_forces()
    efs_calc.output_stress()
    
# Plots
    potential.plot_fortran_potentials()
    potential.plot_python_potentials()
    b_props.bp_output()
    b_props.bp_eos_plot()
    
    print('')     
    print('CONFIGS')    
    for n in range(efs.cc):    
      print('Config ' + str(n+1) + ':', efs.config_energy[n,2], efs.energies[n], (efs.config_energy[n,2]-efs.energies[n])**2)
    
    print('')   
    for bp_id in range(bp.bp_configs_count):  
      print('BP') 
      print('alat:', bp.calc_alat[bp_id], bp.known_alat[bp_id], (bp.calc_alat[bp_id] - bp.known_alat[bp_id])**2)
      print('v0:', bp.calc_v0[bp_id])
      print('e0:', bp.calc_e0[bp_id], bp.known_e0[bp_id], (bp.calc_e0[bp_id] - bp.known_e0[bp_id])**2)
      print('b0:', bp.calc_b0[bp_id], bp.known_b0[bp_id], (bp.calc_b0[bp_id] - bp.known_b0[bp_id])**2)
      print("Calculated Stiffness Matrix (GPA)")
      for i in range(6):
        print(160.230732254e0 * bp.calc_ec[bp_id,i,:])
      print("Known Stiffness Matrix (GPA)")
      for i in range(6):
        print(160.230732254e0 * bp.known_ec[bp_id,i,:])

    print('')
    print('RSS: ' + str(rss))
    
  def get_rss():
    
    rss = 0.0
    
# Run efs calc
    if(efs.cc > 0):
      efs.rss_calc()
      rss = efs.total_rss_weighted / efs.cc
    else:
      rss = 0.0
    
    bp.energy()
    bp.calculate_bp()    
    
    try:
      rss = rss + g.inp['rss']['bp'] * bp.rss_total_rss
    except: 
      rss = rss + bp.rss_total_rss
    
    return rss
    
###########################################
#  CLASS p
###########################################
class pf:

  """
  pv_in = 0.1
  pop_half_size = 8
  pop_size = None
  tab_nodes = 10
  generations = 2
  pv = []             # VARIANCE
  ps = []             # PARAMETERS
  pt = []             # 2 CHILD PARAMETERS
  rss = None          # RSS
  rss_child = None    # RSS
  """
  
# Array Sizes
  pop_size = 50  
  fresh_size = 10  
  pop_overflow = 10
  pop_size_d = None
  fresh_size_d = None 
  arr_len = None 
  
# Parameter arrays
  ps = None
  p_var = None
  p_width = 0

# Rescale Density  0.0 to 1.0 (max is approximated)
  rescale_on = True
  
# Variance Coeffs
  vc_extinction = 0.2

  ncv = 0.05
  tab_nodes = 10
  generations = 200
  spline_generations = 20
  breed_switch = 0.5

  extinction_counter = 0
  extinction_percentage = 50
  extinction_frequency = 3
  extinction_threshold_t = 0.0
  extinction_threshold_k = 0.0

# Record alat, e0, b0
  bulk_properties = None

# PROGRESS COUNTERS
  stage = 'INITIALISING'
  start_time = 0.0
  time_per_calc = 0.0
  time_estimate = 0.0
  time_remaining = 0.0
  since_improvement = 0
  this_gen = 0
  this_cycle = 0
  total_generations = 0
  last_rss = 0.0
  best_rss = 0.0
  calc_counter = 0
  calc_counter_expected = 0
  progress = 0

# Potential Output
  output_counter = 0
  
  def run():    
    
# Read parameter variance
    try:
      pf.pv_in = g.inp['fit']['variance']
    except:
      pass
      
# Read from input into globals
    try:
      g.fit['generations'] = g.inp['fit']['gens']
    except:
      pass
    try:
      g.fit['spline_generations'] = g.inp['fit']['sgens']
    except:
      pass
    try:
      g.fit['spline_cycles'] = g.inp['fit']['scycles']
    except:
      pass
    try:
      g.fit['pop_size'] = g.inp['fit']['ps']
    except:
      pass
    try:
      g.fit['fresh_size'] = g.inp['fit']['fs']
    except:
      pass
    try:
      g.fit['extinction_percentage'] = g.inp['fit']['ep']
    except:
      pass
    try:
      g.fit['extinction_frequency'] = g.inp['fit']['ef']
    except:
      pass
    try:
      g.fit['enhance_frequency'] = g.inp['fit']['ehf']
    except:
      pass
      
# Set pf globals
    pf.generations = g.fit['generations']
    pf.spline_generations = g.fit['spline_generations']
    pf.spline_cycles = g.fit['spline_cycles']
    pf.pop_size = g.fit['pop_size']
    pf.fresh_size = g.fit['fresh_size']
    pf.extinction_percentage = g.fit['extinction_percentage']
    pf.extinction_frequency = g.fit['extinction_frequency']
    pf.enhance_frequency = g.fit['enhance_frequency']
    
# Set Up EFS and BP
    pf.set_up()
    
# Run start up screen
    pf.startup()
    
# Start fit
    pf.fit()
    
#for i in range(10):
#  rss = pot_fit.get_rss()
#  print(i+1, rss)
      
  def startup():
  
    os.system('cls' if os.name == 'nt' else 'clear')    
    start_rss = pf.get_rss_startup()
    
    print("###################################################")      
    print("Starting RSS: ", start_rss)   
    print("###################################################") 
    print("Potential Parameters")
    print("###################################################")   
    pf.print_parameters() 
    
    time.sleep(2.5)
    
  def print_parameters():  
    for fn in range(len(g.pot_functions['functions'])):          
      if(g.pot_functions['functions'][fn]['fit_type'] == 1):     # SPLINE
        print("fn: " + str(fn))
        print("type: spline")
        for i in range(g.pot_functions['functions'][fn]['fit_size']):
          print("Parameter " + str(i) + ": ", end='')
          print(g.pot_functions['functions'][fn]['fit_parameters'][0,i], end='')
          print(", ", end='')
          print(g.pot_functions['functions'][fn]['fit_parameters'][1,i], end='')
          print()
#pf.p_width = pf.p_width + g.pot_functions['functions'][fn]['fit_size']
      elif(g.pot_functions['functions'][fn]['fit_type'] == 2):   # ANALYTIC
        print("fn: " + str(fn))
        print("type: analytic")
        for i in range(g.pot_functions['functions'][fn]['fit_size']):
          print("Parameter " + str(i) + ": ", end='')
          print(g.pot_functions['functions'][fn]['a_params'][i], end='')
          print(", ", end='')
          print(g.pot_functions['functions'][fn]['fit_parameters'][0,i], end='')
          print(", ", end='')
          print(g.pot_functions['functions'][fn]['fit_parameters'][1,i], end='')
          print()
      else:
        print("fn: " + str(fn))
        print("type: no variance")
    
  def set_up():
  
# If a spline fit, convert into a spline
    for fn in range(len(g.pot_functions['functions'])):       
      if(g.pot_functions['functions'][fn]['fit_type'] == 1): 
        potential.vary_tabulated_points(fn)
  
# Setup EFS
    efs.init()                         # Initialise (allocate arrays)
    potential.efs_add_potentials()     # Load potentials
    configs.efs_add_config()           # Add configs
    
# Setup BP
    bp.init()
    potential.bp_add_potentials()
    b_props.bp_add()
    
    rescale_density.run()
    
  def get_rss():
    pf.calc_counter = pf.calc_counter  + 1
    pf.progress = min(int(numpy.floor(100 * (pf.calc_counter / pf.calc_counter_expected))), 100)
    pf.time_per_calc = (time.time() - pf.start_time) / pf.calc_counter
    remaining = max(0, pf.calc_counter_expected - pf.calc_counter)
    pf.time_remaining = pf.time_per_calc * remaining
    pf.time_estimate = pf.time_per_calc * pf.calc_counter_expected
    
# Run efs calc
    if(efs.cc > 0):
      efs.rss_calc()
      rss = efs.total_rss_weighted / efs.cc
    else:
      rss = 0.0    
      
# Bulk property rss
    bp.energy()
    bp.calculate_bp() 
    try:
      rss = rss + g.inp['rss']['bp'] * bp.rss_total_rss
    except: 
      rss = rss + bp.rss_total_rss  
      
# Store last bp
    try:
      if(pf.bulk_properties is None):
        bp_id_max = 0
        for bp_id in range(bp.bp_configs_count):
          if(bp_id + 1 > bp_id_max):
            bp_id_max = bp_id + 1 
        pf.bulk_properties = numpy.zeros((bp_id_max, 20,),)
      for bp_id in range(bp.bp_configs_count):
        pf.bulk_properties[bp_id,0] = bp.known_alat[bp_id]
        pf.bulk_properties[bp_id,1] = bp.known_e0[bp_id]
        pf.bulk_properties[bp_id,2] = 160.230732254e0 * bp.known_b0[bp_id]
        pf.bulk_properties[bp_id,3] = bp.calc_alat[bp_id]
        pf.bulk_properties[bp_id,4] = bp.calc_e0[bp_id]
        pf.bulk_properties[bp_id,5] = 160.230732254e0 * bp.calc_b0[bp_id]
    except:
      pass
# print(pf.bulk_properties)

    return rss
   
  def get_rss_startup():    
# Run efs calc
    if(efs.cc > 0):
      efs.rss_calc()
      rss = efs.total_rss_weighted / efs.cc
    else:
      rss = 0.0    
# Bulk property rss
    bp.energy()
    bp.calculate_bp() 
    try:
      rss = rss + g.inp['rss']['bp'] * bp.rss_total_rss
    except: 
      rss = rss + bp.rss_total_rss          
    return rss
    
  def update_potential(p):    
# Update potential
    a = 0
    for fn in range(len(g.pot_functions['functions'])): 
      if(g.pot_functions['functions'][fn]['fit_type'] == 1):     # SPLINE    
###NEED TO REDO MAYBE??###
      
# Calc b
        b = a + g.pot_functions['functions'][fn]['fit_size']          
# LOAD ORIGINAL
        g.pot_functions['functions'][fn]['points'] = numpy.copy(g.pot_functions['functions'][fn]['points_original'])        
# VARY SPLINE
        potential.vary_tabulated_points(fn, p[a:b])
# Update a
        a = b        
      elif(g.pot_functions['functions'][fn]['fit_type'] == 2):   # ANALYTIC  
        b = a + g.pot_functions['functions'][fn]['fit_size'] 
# Make Analytic Points
        g.pot_functions['functions'][fn]['a_params'][:] = p[a:b]
        potential.make_analytic_points_inner(fn)
        a = b    
      
# Rescale density functions 0.0 to 1.0
    if(pf.rescale_on):
      rescale_density.run()

# Update efs and bp modules
    potential.efs_add_potentials()     # Load potentials
    potential.bp_add_potentials()      # Load potentials
    
  def spline_cycles_prep(cycle):
# Update with best
    pf.update_potential(pf.ps[-5,:])
    pw = 0
    pw_per_f = 15
   
# Update potential fit settings
    for fn in range(len(g.pot_functions['functions'])): 
# Replace original points with new points
      g.pot_functions['functions'][fn]['points_original'] = numpy.copy(g.pot_functions['functions'][fn]['points'])  
# fit type to spline
      g.pot_functions['functions'][fn]['fit_type'] = 1  
      g.pot_functions['functions'][fn]['fit_size'] = pw_per_f    
      g.pot_functions['functions'][fn]['fit_parameters'] = numpy.zeros((pw_per_f, ),)
      g.pot_functions['functions'][fn]['fit_parameters'][:] = 10 * 0.1**cycle
      g.pot_functions['functions'][fn]['fit_parameters'][-1] = 0.0  # Don't move last node
      g.pot_functions['functions'][fn]['fit_mult'] = 1.0
      pw = pw + pw_per_f
      
# Update
    z = numpy.zeros((pw, ),)
    pf.update_potential(z)
    
  def fit():
# Init counters
    pf.start_time = time.time()
    pf.since_improvement = 0
    pf.improvement_time = time.time()
    pf.this_gen = 0
    pf.extinction_counter = 0
    pf.extinction_threshold_k = 0.0
    pf.extinction_threshold_t = 0.0
    pf.this_cycle = 0
    pf.total_generations = 0
    pf.last_rss = 0.0
    pf.best_rss = 0.0
    pf.calc_counter = 0    
    pf.estimate_calc_counter()

    pf.cycle(1)
    
    """
# OUTPUT FIT
    pf.update_potential(pf.ps[-5,:])
    
    if(pf.spline_cycles > 0):
      for i in range(pf.spline_cycles):
        pf.update_potential(pf.ps[-5,:])
        pf.spline_cycles_prep(i)
        pf.cycle(i+2)
        pf.update_potential(pf.ps[-5,:])
    
    potential.pf_output()
    potential.plot_python_potentials()
    potential.plot_fortran_potentials()
    """
    
  def start_cycle(): 
# Pop size
    pf.pop_size_d = 2 * pf.pop_size
    pf.fresh_size_d = 2 * pf.fresh_size    
    pf.arr_len = pf.pop_size_d + pf.fresh_size_d + pf.pop_overflow
    pf.child_len = pf.pop_size_d + pf.fresh_size_d + pf.fresh_size_d
    
# Find width
    pf.p_width = 0
    for fn in range(len(g.pot_functions['functions'])):          
      if(g.pot_functions['functions'][fn]['fit_type'] == 1):     # SPLINE
        pf.p_width = pf.p_width + g.pot_functions['functions'][fn]['fit_size'] 
      elif(g.pot_functions['functions'][fn]['fit_type'] == 2):   # ANALYTIC
        pf.p_width = pf.p_width + g.pot_functions['functions'][fn]['fit_size']
    
# MAKE ARRAYS
    pf.p_var = numpy.zeros((2, pf.p_width,),)
    pf.ps = numpy.zeros((pf.arr_len, pf.p_width,),)
    pf.rss = numpy.zeros((pf.arr_len,),)   
    pf.p_children = numpy.zeros((pf.child_len, pf.p_width,),)
    pf.rss_children = numpy.zeros((pf.child_len,),)
    pf.rss_list = numpy.zeros((pf.pop_size_d + pf.child_len,),)
    pf.p_temp = numpy.zeros((pf.pop_size_d, pf.p_width,),)
    pf.rss_temp = numpy.zeros((pf.pop_size_d,),)
    pf.p_top = numpy.zeros((int(pf.pop_size_d / 4), pf.p_width,),)
    pf.rss_sort = numpy.zeros((pf.pop_size_d,),)
    
# SET Variation
    a = 0
    for fn in range(len(g.pot_functions['functions'])):        
      if(g.pot_functions['functions'][fn]['fit_type'] == 1):     # TABULATED      
        b = a + g.pot_functions['functions'][fn]['fit_size']     
        pf.p_var[0,a:b] = g.pot_functions['functions'][fn]['fit_parameters'][0,:]  # Lower
        pf.p_var[1,a:b] = g.pot_functions['functions'][fn]['fit_parameters'][1,:]  # Upper
        a = b 
      elif(g.pot_functions['functions'][fn]['fit_type'] == 2):   # ANALYTIC  
        b = a + g.pot_functions['functions'][fn]['fit_size'] 
        pf.p_var[0,a:b] = g.pot_functions['functions'][fn]['fit_parameters'][0,:]  # Lower
        pf.p_var[1,a:b] = g.pot_functions['functions'][fn]['fit_parameters'][1,:]  # Upper
        a = b  
        
#    yvar = numpy.zeros((len(g.pot_functions['functions'][fn]['fit_parameters']),),)
#    potential.vary_tabulated_points(fn, yvar)
        
# SET -1.0 for best rss
    pf.rss[-5] = -1.0
        
  def cycle(c):
    """
    pv_in = 0.1
    pop_half_size = 8
    pop_size = None
    tab_nodes = 10
    generations = 2
    """
    
    """
    -1  Initial input parameters
    -2  Mutation
    -3  child a
    -4  child b
    -5  best
    """
    time_start = time.time()
    
# Counter
    pf.this_cycle = pf.this_cycle + 1
    pf.this_gen = 0
    
# Start Cycle   (calculate width,
    pf.start_cycle()
    
# SAVE Starting Parameters
    a = 0
    for fn in range(len(g.pot_functions['functions'])):        
      if(g.pot_functions['functions'][fn]['fit_type'] == 1):     # TABULATED      
        b = a + g.pot_functions['functions'][fn]['fit_size']     
        pf.ps[-1,a:b] = numpy.zeros((g.pot_functions['functions'][fn]['fit_size'],),)
        a = b
      elif(g.pot_functions['functions'][fn]['fit_type'] == 2):   # ANALYTIC  
        b = a + g.pot_functions['functions'][fn]['fit_size'] 
        pf.ps[-1,a:b] = g.pot_functions['functions'][fn]['a_params'][:]
        a = b

# Update and calculate starting rss
    pf.update_potential(pf.ps[-1,:])
    pf.rss[-1] = pf.get_rss() 
    
# Try within the range provided
    
# Create initial population - First half
    pf.stage = 'Initialising population - First Half'
    for p in range(pf.pop_size):
      if(p == 0):
        pf.ps[p,:] = pf.ps[-1,:]
      else:
        pf.ps[p,:] = pf.ps[-1,:]
                
# Calculate and save RSS
      pf.ps[p,:] = pf.random_p(0.0, 1.0)
      pf.update_potential(pf.ps[p,:])
      pf.rss[p] = pf.get_rss() 
      pf.check_improvement(pf.ps[p,:], pf.rss[p])

# Look in a wider range

# Create initial population - Second half
    pf.stage = 'Initialising population - Second Half'
    m = 0.5
    m_inc = (10.0-0.5) / (pf.pop_size - 1)
    for p in range(pf.pop_size, pf.pop_size_d):
      loop = True
      while(loop):  
        pf.ps[p,:] = pf.random_p(0.0, m)
        pf.update_potential(pf.ps[p,:])
        pf.rss[p] = pf.get_rss() 
        if(not numpy.isnan(pf.rss[p])):
          loop = False
          pf.check_improvement(pf.ps[p,:], pf.rss[p])
          m = m + m_inc

###################################
# LOOP THROUGH GENERATIONS
###################################
    
    gens = pf.generations
    if(c > 1):    
      gens = pf.spline_generations
    pf.since_improvement = 0
    for gen in range(gens):
      pf.this_gen = gen + 1
    
      parents = numpy.arange(pf.pop_size_d)
      numpy.random.shuffle(parents)

      pf.stage = 'Looping population'
      
#######################
# Parents + Parents
#######################
      
      ca = 0
# Loop through population
      for p in range(pf.pop_size):   
        pa = parents[p]
        pb = parents[p + pf.pop_size]
        
# Breed
        pf.breed(pa, pb, ca, ca + 1)
        
# Run calculations
        pf.update_potential(pf.p_children[ca, :])        
        pf.rss_children[ca] = pf.get_rss()  
        pf.check_improvement(pf.p_children[ca,:], pf.rss_children[ca])      
          
        pf.update_potential(pf.p_children[ca+1, :])        
        pf.rss_children[ca+1] = pf.get_rss()  
        pf.check_improvement(pf.p_children[ca+1,:], pf.rss_children[ca+1])
                
# Increment
        ca = ca + 2
 
# MAKE FRESH PARAMETERS
      m = 0.5
      m_inc = (10.0-0.5) / (pf.fresh_size_d - 1)
      for p_fresh in range(pf.pop_size_d, pf.pop_size_d + pf.fresh_size_d):          
        pf.ps[p_fresh,:] = pf.random_p(pf.ps[-5,:], m)
        m = m + m_inc
      
#######################
# Parents + Fresh
#######################
      
# Pick parents
      parents = numpy.arange(pf.pop_size_d)
      numpy.random.shuffle(parents)
      
      pf.stage = 'Looping population + fresh pool'
# Breed random parents with fresh pool
      pa_n = 0
      for p_fresh in range(pf.pop_size_d, pf.pop_size_d + pf.fresh_size_d):
        pa = parents[pa_n]
        
# Breed
        pf.breed(pa, p_fresh, ca, ca + 1)
        
# Run calculations
        pf.update_potential(pf.p_children[ca, :])        
        pf.rss_children[ca] = pf.get_rss()  
        pf.check_improvement(pf.p_children[ca, :], pf.rss_children[ca])     
          
        pf.update_potential(pf.p_children[ca+1, :])        
        pf.rss_children[ca + 1] = pf.get_rss()  
        pf.check_improvement(pf.p_children[ca+1, :], pf.rss_children[ca+1])
         
# Increment
        ca = ca + 2
        
#######################
# Select Best
#######################
        
      pf.rss_list[:pf.pop_size_d] = pf.rss[:pf.pop_size_d]
      pf.rss_list[pf.pop_size_d:] = pf.rss_children[:]

      rss_sorted = sort.sort_1d_dp_asc(pf.rss_list)
      rt = rss_sorted[pf.pop_size_d-1]

# SELECT BEST FROM PARENTS AND CHILDREN FOR NEXT GEN
      n = 0
      for i in range(pf.pop_size_d):
        if(pf.rss[i] < rt):
          pf.p_temp[n,:] = pf.ps[i,:]
          pf.rss_temp[n] = pf.rss[i]
          n = n + 1          
      for i in range(len(pf.rss_children)):
        if(pf.rss_children[i] < rt):
          pf.p_temp[n,:] = pf.p_children[i,:]
          pf.rss_temp[n] = pf.rss_children[i]
          n = n + 1
      for i in range(pf.pop_size_d):
        if(pf.rss[i] == rt and n < pf.pop_size_d):
          pf.p_temp[n,:] = pf.ps[i,:]
          pf.rss_temp[n] = pf.rss[i]
          n = n + 1      
      for i in range(len(pf.rss_children)):
        if(pf.rss_children[i] == rt and n < pf.pop_size_d):
          pf.p_temp[n,:] = pf.p_children[i,:]
          pf.rss_temp[n] = pf.rss_children[i]
          n = n + 1
          
# UPDATE RSS AND PS
      pf.rss[1:pf.pop_size_d] = numpy.copy(pf.rss_temp[1:pf.pop_size_d])  
      pf.ps[1:pf.pop_size_d] = numpy.copy(pf.p_temp[1:pf.pop_size_d])  

# RUN EXTINCTION EVENT
      if((gen % pf.extinction_frequency) == 0 and gen > 0 and gen < pf.generations-1):
        pf.stage = 'Extinction event'
        pf.extinction()

# RUN ENHANCE EVENT
      if(((gen + 1) % pf.enhance_frequency) == 0):
        pf.stage = 'Enhancing Top 10% With Steepest Descent'
        pf.enhance()
        
# SAVE POTENTIAL
      pf.output()
        
# Run end function
    pf.stage = 'End' 
    results = pf.end()   
    pf.progress = 100
    pf.time_remaining = 0.0
    progress.display(results) 
    
  def rss_threshold(f=50):
    rss_sorted = sort.sort_1d_dp_asc(pf.rss[0:pf.pop_size_d])
    t = int(numpy.floor(len(rss_sorted) * (f / 100)))
    return rss_sorted[t]
    
  def check_improvement(params, rss): 
# Increment Improvement Counter
    pf.since_improvement = pf.since_improvement + 1 
   
# Update best (if better)
    if(pf.rss[-5] < 0 or rss < pf.rss[-5]):
      pf.ps[-5,:] = params
      pf.rss[-5] = rss  
      pf.improvement_time = time.time()  
      pf.since_improvement = 0
      pf.best_rss = pf.rss[-5]
      
      pf.bulk_properties[:,6:9] = pf.bulk_properties[:,3:6]
    pf.last_rss = rss 
    
# Save to print in display
    pf.last_ps = params
    pf.best_ps = pf.ps[-5,:]
    pf.density_list = rescale_density.max_densities()
    progress.display()
    
  def breed(pa, pb, ca, cb):    
    state = True    
    for i in range(pf.p_width):
      rn = random.random()
      if(rn <= pf.breed_switch):
        if(state):
          state = False
        else:
          state = True
      if(state):
        pf.p_children[ca, i] = pf.ps[pa,i]
        pf.p_children[cb, i] = pf.ps[pb,i]
      else:      
        pf.p_children[cb, i] = pf.ps[pa,i]
        pf.p_children[ca, i] = pf.ps[pb,i]
        
# Mutation (possibly)
    pf.p_children[ca, :] = pf.mutate(pf.p_children[ca, :], chance=5) 
    pf.p_children[cb, :] = pf.mutate(pf.p_children[cb, :], chance=5)    

# No clones - Child A
    vary = False
    for p in range(pf.pop_size_d):
      if((pf.p_children[ca, :] == pf.ps[p,:]).all()):
        vary = True
        break
    if(vary):
      r = numpy.random.rand(pf.p_width)
      pf.p_children[ca, :] = pf.p_children[ca, :] * pf.ncv * (0.5 - r)
    
# No clones - Child B
    vary = False
    for p in range(pf.pop_size_d):
      if((pf.p_children[cb, :] == pf.ps[p,:]).all()):
        vary = True
        break
    if(vary):
      r = numpy.random.rand(pf.p_width)
      pf.p_children[cb, :] = pf.p_children[cb, :] * pf.ncv * (0.5 - r)

  def extinction():
# Extinction for those below rss_thresh
# Replace with variations of the best 25% parameters
  
    pf.extinction_counter = pf.extinction_counter + 1
    rss_thresh = pf.rss_threshold(pf.extinction_percentage)
    pf.extinction_threshold_k = rss_thresh
    
    top_size = int(pf.pop_size_d / 4)    
    pf.rss_sort = sort.sort_1d_dp_asc(pf.rss[0:pf.pop_size_d])
    top_rss_min = pf.rss_sort[top_size-1]
    pf.extinction_threshold_t = top_rss_min
    
    n = 0
    for i in range(pf.pop_size_d):
      if(pf.rss[i] <= top_rss_min and n < top_size):
        pf.p_top[n,:] = pf.ps[i,:]
        n = n + 1
      if(n == top_size):
        break
    
#pf.rss_sort = numpy.zeros((pf.pop_size_d,),)

    for p in range(pf.pop_size_d):
      if(pf.rss[p]>rss_thresh):   
        r_top = numpy.random.randint(top_size)
        
        pf.ps[p,:] = pf.random_p(pf.p_top[r_top], 0.01)

        pf.update_potential(pf.ps[p,:])
        pf.rss[p] = pf.get_rss() 
        pf.check_improvement(pf.ps[p,:], pf.rss[p])
        
  def enhance():
# Enhance best 10% with gradient descent
    top_size = int(pf.pop_size_d / 10)    
    pf.rss_sort = sort.sort_1d_dp_asc(pf.rss[0:pf.pop_size_d])
    top_rss = pf.rss_sort[top_size-1]
    
    for p in range(pf.pop_size_d):
      if(pf.rss[p] <= top_rss):
        params = gd.opt(pf.gd_rss, pf.ps[p])
        if(gd.rss_out < pf.rss[p]):
          pf.rss[p] = gd.rss_out
          pf.ps[p, :] = numpy.copy(params)
          
  def gd_rss(params):
    pf.update_potential(params)
    return pf.get_rss() 

  def estimate_calc_counter():
  
    e = 0
    
    for c in range(pf.spline_cycles + 1):
      e = e + pf.pop_size
      e = e + pf.pop_size
            
      gens = pf.generations
      if(c > 1):    
        gens = pf.spline_generations
      for gen in range(gens):
        e = e + 2 * pf.pop_size
        e = e + 2 * pf.fresh_size
        if((gen % pf.extinction_frequency) == 0 and gen > 0 and gen < pf.generations-1):
          e = e + pf.pop_size
        if((gen % pf.enhance_frequency) == 0):
          e = e + (300 * (pf.pop_size // 4))
    
    pf.calc_counter_expected = e
  
  def random_p(c=0.0, m=1.0):
#
    lower = pf.p_var[0,:]
    upper = pf.p_var[1,:]
    range = upper - lower
    
# If there's no center, take midpoint of upper/lower - else center it on the parameters c
    if(type(c) != numpy.ndarray and c == 0.0):
      c = lower + 0.5 * range 
    
# Multiply range
    m_range = m * range
    
# Get random parameters 0 to 1
    r = numpy.random.rand(pf.p_width)
    
# New parameters
    p_new = c + (r - 0.5) * m_range
    
# Return
    return p_new
    
  def mutate(params, chance=5):
    mutant = pf.random_p(params, 1.0)
    for i in range(len(params)):
      if(random.randrange(0, 100) <= chance):
        params[i] = mutant[i]
    return params
    
  def end():
    pf.update_potential(pf.ps[-5,:])
    rss = pf.get_rss() 
    
    results = []
    results.append("RSS: " + str(rss))
    
    for bp_id in range(bp.bp_configs_count):
      results.append("alat: " + str(bp.calc_alat[bp_id]))
      results.append("v0: " + str(bp.calc_v0[bp_id]))
      results.append("e0: " + str(bp.calc_e0[bp_id]))
      results.append("b0: " + str(bp.calc_b0[bp_id]))
      results.append("b0/GPA: " + str(160.230732254e0 * bp.calc_b0[bp_id]))
      return results
  
# To Do
  def sort_parameters():
    pf.rss[0:pf.pop_size_d] = sort.sort_1d_dp_asc(pf.rss[0:pf.pop_size_d])
    
#kt = sort.key_table
#pf.ps[0:pf.pop_size_d] = sort.match_key_table_2d_dp_kt(pf.ps[0:pf.pop_size_d], kt)
    
  def output():
    pot_dir = g.dirs['fitting'] + "/fits"
    std.make_dir(pot_dir)   
  
    pf.output_counter = pf.output_counter + 1
    pf.update_potential(pf.ps[-5,:])
  
    c = str(pf.output_counter)
    
    while(len(c)<4):
      c = "0" + c
    
    fit_dir = pot_dir + "/" + c
    std.make_dir(fit_dir)   
    
    potential.plot_fortran_potentials(fit_dir)
    
###########################################
#  CLASS progres
###########################################
class progress:

  def display(results=None):
  
# CLEAR
    os.system('cls' if os.name == 'nt' else 'clear')

    print("###################################################################")
    print("#  Pop/Fresh size  :   " + progress.pad_l(pf.pop_size)  + " /  " +  progress.pad_l(pf.fresh_size))
    print("#  Generations:        " + progress.pad_l(pf.generations) )
    print("#  Spline Cycles:      " + progress.pad_l(pf.spline_cycles) + " /  " + progress.pad_l(pf.spline_generations) + "")
    print("#  Ext/Enhance Freq:   " + progress.pad_l(pf.extinction_frequency) + " /  " + progress.pad_l(pf.enhance_frequency) + "")    
    print("###################################################################")
    print("# Timer:              " + progress.pad_l(time.time() - pf.start_time, 10) + " (Estimate: " + progress.pad_l(pf.time_estimate, 10) + ")")
    print("# Stage:              " + pf.stage) 
    print("# Cycle:              " + progress.pad_l(pf.this_cycle, 10))  
    print("###################################################################")
    print("# Gen:                " + progress.pad_l(pf.this_gen, 10))  
    print("# Exctinctions:       " + progress.pad_l(pf.extinction_counter, 10) 
          + " (" + progress.pad_l(pf.extinction_threshold_t, 10) + "/" 
          + "" + progress.pad_l(pf.extinction_threshold_k, 10) + ")") 
    print("# Since Improvement:  " + progress.pad_l(pf.since_improvement, 10))  
    print("#    ")  
    print("# Best RSS:           ", pf.best_rss) 
    print("#    ")   
    print("# Last RSS:           ", pf.last_rss)  
    for bp_id in range(len(pf.bulk_properties)):
      print("#                     a0  " + progress.pad_l(pf.bulk_properties[bp_id,3], 10) + "   " + progress.pad_l(pf.bulk_properties[bp_id,6], 10) + "   " + progress.pad_l(pf.bulk_properties[bp_id,0], 10) + " ") 
      print("#                     e0  " + progress.pad_l(pf.bulk_properties[bp_id,4], 10) + "   " + progress.pad_l(pf.bulk_properties[bp_id,7], 10) + "   " + progress.pad_l(pf.bulk_properties[bp_id,1], 10) + " ") 
      print("#                     b0  " + progress.pad_l(pf.bulk_properties[bp_id,5], 10) + "   " + progress.pad_l(pf.bulk_properties[bp_id,8], 10) + "   " + progress.pad_l(pf.bulk_properties[bp_id,2], 10) + " ") 
      
    print("#    ")  
    print("# PROGRESS:   ", progress.bar(pf.progress, 25))  
    print("# TIME LEFT:  ", round(pf.time_remaining,2))  
    print("#   ")  
    print("# Parameters (this): ", end="")
    if(type(pf.last_ps) == numpy.ndarray):
      for pn in range(len(pf.last_ps)):
        print(pf.last_ps[pn], end=" ")
        if((pn+1) % 5 == 0):
          print()
          print("#                   ", end="")
    print("")
    print("# Parameters (best): ", end="")
    if(type(pf.best_ps) == numpy.ndarray):
      for pn in range(len(pf.best_ps)):
        print(pf.best_ps[pn], end=" ")
        if((pn+1) % 5 == 0):
          print()
          print("#                   ", end="")
    print("")
    print("#   ")  
    print("# Density List: ", pf.density_list)  
    print("#   ")  
    
#    pf.last_ps = params
#pf.best_ps = params
    
    print("###################################################################")
    
    if(type(results) == list):
      for r in results:
        print("# " + str(r))  
      print("###################################################################")
        
#for bp_id in range(bp.bp_configs_count):
    
#pf.since_improvement
    """ 
    print("Running Fit")
    print("Pop size: ", str(pf.pop_size), " (", str(pf.pop_size_d), ")", sep="")
    print("Fresh size: ", str(pf.fresh_size), " (", str(pf.fresh_size_d), ")", sep="")
    print("Generations: ", str(pf.generations), sep="")
 pf.spline_cycles
    """
 
  @staticmethod
  def pad_r(inp, p=7):
    if(inp == None):
      return ""      
    out = str(inp).strip()  
    while(len(out)<p):
      out = out + " "      
    return out[0:p]
    
  @staticmethod
  def pad_l(inp, p=7):
    if(inp == None):
      return ""      
    out = str(inp).strip()  
    while(len(out)<p):
      out = " " + out     
    return out[0:p]
    
  @staticmethod
  def bar(p, w=25):
    out = ''
    p_num = p
    p = p * (25 / 100)
    for i in range(w):
      if(i <= p):
        out = out + '#'
      else:
        out = out + '_'
    out = out + '   ' + str(p_num) + '%'
    return out
    
  """
  start_time = 0
  since_improvement = 0
  this_gen = 0
  this_cycle = 0
  total_generations = 0
  last_rss = 0.0
  best_rss = 0.0
  """

###########################################
#  CLASS g
###########################################
class gd:

  precision = 1.0e-6
  max_iterations = 10
  h = 1.0e-8
  momentum = 0.1
  p_in = None
  rss_in = None
  p_out = None
  rss_out = None
  
# Central Difference
  def df(p):
    d = numpy.zeros((len(p),),)
    for i in range(len(p)):
      p_f = numpy.copy(p)
      p_f[i] = p_f[i] + gd.h
      p_b = numpy.copy(p)
      p_b[i] = p_b[i] - gd.h
      d[i] = (gd.rss(p_f) - gd.rss(p_b)) / (2 * gd.h)
    return d

  def line_search(p, dp):
# Back track
    best_rss = gd.rss_in
    best_gamma = 0.0
    gamma = 20.0 * gd.last_gamma 
    if(gamma <1.0e-8):
      gamma <1.0e-6
    loop = True
    n = 0
    while(gamma > 1.0e-12):
      n = n + 1
      p_test = p - gamma * dp
      rss = gd.rss(p_test)
      pf.check_improvement(p + gamma * dp, rss)
      gamma = 0.5 * gamma
      if(best_rss is None or rss < best_rss):
        p_best = p_test
        best_rss = rss
        best_gamma = gamma
    gd.last_gamma = best_gamma
    return best_rss, best_gamma * dp
    
# Gradient Descent
  def opt(f_rss, p0):
    gd.rss = f_rss   
    gd.last_gamma = 1.0
    p = p0
    
    gd.p_in = numpy.copy(p)
    gd.rss_in = gd.rss(p)
    
    best_p = numpy.copy(gd.p_in)
    best_rss = gd.rss_in

    n = 0
    last_dp = 0.0
    while(n < gd.max_iterations):
      df = gd.df(p)
      rss, dp = gd.line_search(p, df)
      if(rss > best_rss):
        n = gd.max_iterations
      else:
        p = p - (gd.momentum * last_dp + dp)
        last_dp = dp
        rss = gd.rss(p)
        best_p = numpy.copy(p)
        best_rss = rss
        n = n + 1

    pf.check_improvement(best_p, best_rss)
    gd.p_out = numpy.copy(best_p)
    gd.rss_out = best_rss
    
    return p

###########################################
###########################################
#  MAIN
###########################################
###########################################

def main():

  

# RECORD START TIME
  g.times['start'] = time.time()

  

# MAKE DIRS
  for d in g.dirs.keys():

    dir = g.dirs[d]

    std.make_dir(dir)

  

# OPEN LOG
  g.log_fh = open(g.dirs['log'] + '/main.log', 'w')

  g.log_fh.write('###########################################################################\n')

  g.log_fh.write(time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime()) + '\n')

  g.log_fh.write('###########################################################################\n')

  g.log_fh.write('\n')

  g.log_fh.write('Script: ' + str(sys.argv[0]) + '\n')

  if(len(sys.argv)>1):

    shutil.copyfile(sys.argv[1], g.dirs['log'] + '/input.log')

    run_program = False

    try:

      g.inp = read_config.read_file(sys.argv[1])

      g.log_fh.write('Loaded: ' + str(sys.argv[1]) + '\n')

      run_program = True

    except:

      g.log_fh.write('Unable to load, exiting\n')

      

# RUN
    if(run_program):

      eampa.run()

    

# CLOSE LOG
  g.times['end'] = time.time()

  g.log_fh.write('\n')

  g.log_fh.write('###########################################################################\n')

  g.log_fh.write('Duration: ' + str(g.times['end'] - g.times['start']) + '\n')

  g.log_fh.write('###########################################################################\n')

  g.log_fh.close()



# Run
main()

