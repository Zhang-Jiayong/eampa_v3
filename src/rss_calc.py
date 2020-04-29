######################################################
#  Ben Palmer University of Birmingham 2020
#  Free to use
######################################################

import numpy
from f2py_lib.f_efs import efs
from f2py_lib.f_bp import bp
from potential import potential
import matplotlib.pyplot as plt
import time


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
    print ('RSS: ' + str(rss))
    
    # Output to File
    efs_calc.output_energy()
    efs_calc.output_forces()
    efs_calc.output_stress()
    
    # Plots
    potential.plot_fortran_potentials()
    potential.plot_python_potentials()
    b_props.bp_output()
    b_props.bp_eos_plot()
    
    
    
    
    
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
    
    