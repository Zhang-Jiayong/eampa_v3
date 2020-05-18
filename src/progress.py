import numpy
import os
from pot_fit import pf

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

