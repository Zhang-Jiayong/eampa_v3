######################################################
#  Ben Palmer University of Birmingham 2020
#  Free to use
######################################################

import numpy
#from tendl import tendl
#from isotopes import isotopes
#import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from labels import labels
from potential import potential
from potential_vary import potential_vary
from configs import configs
from e_adjust import e_adjust
from b_props import b_props
from efs_calc import efs_calc
from bp_calc import bp_calc
from es_calc import es_calc
from rss_calc import rss_calc
from pot_fit import pf

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
  
  
  
  