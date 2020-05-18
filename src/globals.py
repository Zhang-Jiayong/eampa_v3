##########
# GLOBALS

import numpy

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
         
         