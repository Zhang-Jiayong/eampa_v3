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
import random





"""
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
"""


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
  ps_key = []
  p_width = 0
  since_improvement = 0

  
  # Variance Coeffs
  vc_initial_mult = 2.0
  vc_initial_add = 2.0
  vc_mutation = 0.2
  vc_extinction = 0.2


  pv_initial = 2.0
  pv_mult = 10.0
  pv_add = 10.0
  pv_fresh = 5.0
  
  ext_thresh = 0.5
  gen_per_ext = 10
  
  ncv = 0.05
  tab_nodes = 10
  generations = 200
  breed_switch = 0.5



  def run():
    print("POT FITTING")
    
    # Read parameter variance
    try:
      pf.pv_in = g.inp['fit']['variance']
    except:
      pass
    
    # Set Up EFS and BP
    pf.set_up()
    
    # Start fit
    pf.fit()
    
    
    #for i in range(10):
    #  rss = pot_fit.get_rss()
    #  print(i+1, rss)
      
      
    
  def set_up():
    # Setup EFS
    efs.init()                         # Initialise (allocate arrays)
    potential.efs_add_potentials()     # Load potentials
    configs.efs_add_config()           # Add configs
    
    # Setup BP
    bp.init()
    potential.bp_add_potentials()
    b_props.bp_add()
    
   
    
  def get_rss():
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
    for fn in range(len(g.pot_functions['functions'])):  
      a = pf.ps_key[fn][1]
      b = pf.ps_key[fn][2]        
      if(g.pot_functions['functions'][fn]['function_type'] == 1):     # TABULATED
        pass
      elif(g.pot_functions['functions'][fn]['function_type'] == 2):   # ANALYTIC
        g.pot_functions['functions'][fn]['a_params'][:] = p[a:b]
        
    # Make Tabulated Points
    potential.make_tabulated_points()
    
    # Make Analytic Points
    potential.make_analytic_points()
    
    # Update efs and bp modules
    potential.efs_add_potentials()     # Load potentials
    potential.bp_add_potentials()      # Load potentials
    
    
    
    
  def fit():
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
    
    
    spline_temp = numpy.linspace(1.0,10.0, pf.tab_nodes)
    
    # Pop size
    pf.pop_size_d = 2 * pf.pop_size
    pf.fresh_size_d = 2 * pf.fresh_size
    
    pf.arr_len = pf.pop_size_d + pf.fresh_size_d + pf.pop_overflow
    
    # Fill ps_key and find width of parameter table
    pf.ps_key = []
    for fn in range(len(g.pot_functions['functions'])):   
      pf.ps_key.append([None,None,None,None])
    
    pf.p_width = 0
    a = 0
    for fn in range(len(g.pot_functions['functions'])):          
      if(g.pot_functions['functions'][fn]['function_type'] == 1):     # TABULATED
        pf.p_width = pf.p_width + pf.tab_nodes   
        b = a + pf.tab_nodes
        pf.ps_key[fn][0] = 1
        pf.ps_key[fn][1] = a
        pf.ps_key[fn][2] = b
        a = a + pf.tab_nodes
      elif(g.pot_functions['functions'][fn]['function_type'] == 2):   # ANALYTIC
        pf.p_width = pf.p_width + len(g.pot_functions['functions'][fn]['a_params'])
        b = a + len(g.pot_functions['functions'][fn]['a_params'])
        pf.ps_key[fn][0] = 2
        pf.ps_key[fn][1] = a
        pf.ps_key[fn][2] = b
        a = a + len(g.pot_functions['functions'][fn]['a_params'])
    #ps
    
    # Make Array
    pf.ps = numpy.zeros((pf.arr_len, pf.p_width,),)
    pf.rss = numpy.zeros((pf.arr_len,),)
    
    # Load Starting Parameters
    for fn in range(len(g.pot_functions['functions'])):          
      if(g.pot_functions['functions'][fn]['function_type'] == 1):     # TABULATED
        a = pf.ps_key[fn][1]
        b = pf.ps_key[fn][2]
        pf.ps[-1,a:b] = spline_temp[:] # Temporary 
      elif(g.pot_functions['functions'][fn]['function_type'] == 2):   # ANALYTIC
        a = pf.ps_key[fn][1]
        b = pf.ps_key[fn][2]
        pf.ps[-1,a:b] = g.pot_functions['functions'][fn]['a_params'][:]
   
    
    # Create initial population
    for p in range(pf.pop_size_d):
      if(p == 0):
        pf.ps[p,:] = pf.ps[-1,:]
      else:
        r = numpy.random.rand(pf.p_width)
        if(p % 2 == 0):
          pf.ps[p,:] = (0.5 - r) * pf.vc_initial_mult * pf.ps[-1,:]
        else:
          pf.ps[p,:] = (1.0 + (0.5 - r) * pf.vc_initial_add) * pf.ps[-1,:]

    
    # Calc rss for initial population
    for p in range(pf.pop_size_d):
      pf.update_potential(pf.ps[p,:])
      pf.rss[p] = pf.get_rss() 
      if(p==0 or pf.rss[p] < pf.rss[-5]):
        pf.ps[-5,:] = pf.ps[p,:]
        pf.rss[-5] = pf.rss[p]
      print(pf.rss[p], pf.rss[-5])
    pf.rss[-1] = pf.rss[0]
    
    # Make Mutation Array
    pf.make_mutation()
    
    # LOOP THROUGH GENERATIONS  
    pf.since_improvement = 0
    for gen in range(pf.generations):
    
      parents = numpy.arange(pf.pop_size_d)
      numpy.random.shuffle(parents)

      # Loop through population
      for p in range(pf.pop_size):   
        pa = parents[p]
        pb = parents[p + pf.pop_size]
        
        # Breed
        pf.since_improvement = pf.since_improvement + 1
        pf.breed(pa, pb)
        
        # Run with children and save if better
        pf.pp_children(pa, pb)
        
        print(pf.rss[-3], pf.rss[-4], pf.rss[-5], pf.since_improvement)
        
      # MAKE FRESH PARAMETERS
      for p_fresh in range(pf.pop_size_d, pf.pop_size_d + pf.fresh_size_d):  
        if(p_fresh % 2 == 0):
          r = numpy.random.rand(pf.p_width)
          pf.ps[p_fresh,:] = (0.5 - r) * pf.pv_mult * pf.ps[-5,:]      
        else:
          pf.ps[p_fresh,:] = (0.5 - r) * pf.pv_add + pf.ps[-5,:] 
          
      # Pick parents
      parents = numpy.arange(pf.pop_size_d)
      numpy.random.shuffle(parents)
      
      # Breed random parents with fresh pool
      pa_n = 0
      for p_fresh in range(pf.pop_size_d, pf.pop_size_d + pf.fresh_size_d):
        pa = parents[pa_n]
        
        # Breed
        pf.since_improvement = pf.since_improvement + 1
        pf.breed(pa, p_fresh)
        
        # Run with children and save if better
        pf.pf_children(pa)
        print(pf.rss[-3], pf.rss[-4], pf.rss[-5], pf.since_improvement)
        
        
      if((gen % pf.gen_per_ext) == 0 and gen > 0 and gen < pf.generations-1):
        pf.extinction()
        
    # Run end function
    pf.end()    
        
        
        
        

      
          
  def breed(pa, pb):    
    state = True    
    for i in range(pf.p_width):
      rn = random.random()
      if(rn <= pf.breed_switch):
        if(state):
          state = False
        else:
          state = True
      if(state):
        pf.ps[-3,i] = pf.ps[pa,i]
        pf.ps[-4,i] = pf.ps[pb,i]
      else:      
        pf.ps[-4,i] = pf.ps[pa,i]
        pf.ps[-3,i] = pf.ps[pb,i]
    
    # Mutation
    ri = random.randrange(0, 100)     
    if(ri > 3):
      ri = 0
    else:      
      wa = numpy.arange(pf.p_width)
      numpy.random.shuffle(wa)
      wb = numpy.arange(pf.p_width)
      numpy.random.shuffle(wb)
      ra = numpy.random.rand(pf.p_width)
      rb = numpy.random.rand(pf.p_width)
    for m in range(ri):
      ri = random.randrange(0, 100) 
      if(ri < 90):
        pf.ps[-3,wa[m]] = pf.ps[-3,wa[m]] + (0.5 - ra[m]) * pf.ps[-2,wa[m]]
      else:  
        pf.ps[-3,wa[m]] = pf.ps[-3,wa[m]] * (0.5 - ra[m]) * pf.ps[-2,wa[m]]
      if(ri < 90):        
        pf.ps[-4,wb[m]] = pf.ps[-4,wb[m]] + (0.5 - rb[m]) * pf.ps[-2,wb[m]]
      else:  
        pf.ps[-4,wb[m]] = pf.ps[-4,wb[m]] * (0.5 - rb[m]) * pf.ps[-2,wb[m]]
    
    # No clones
    vary = False
    for p in range(pf.pop_size_d):
      if((pf.ps[-3,:] == pf.ps[p,:]).all()):
        vary = True
        break
    if(vary):
      r = numpy.random.rand(pf.p_width)
      pf.ps[-3,:] = pf.ps[-3,:] * pf.ncv * (0.5 - r)
    vary = False
    for p in range(pf.pop_size_d):
      if((pf.ps[-4,:] == pf.ps[p,:]).all()):
        vary = True
        break
    if(vary):
      r = numpy.random.rand(pf.p_width)
      pf.ps[-4,:] = pf.ps[-4,:] * pf.ncv * (0.5 - r)
    
    
    
    
 
  def pp_children(pa, pb): 
    # Run child A
    pf.update_potential(pf.ps[-3,:])
    pf.rss[-3] = pf.get_rss()  
    pr = pa
    if(pf.rss[pa]<pf.rss[pb]):
      pr = pb        
    if(pf.rss[-3]<pf.rss[pr]):
      pf.ps[pr,:] = pf.ps[-3,:]
      pf.rss[pr] = pf.rss[-3]          
      if(pf.rss[pr] < pf.rss[-5]):
        pf.since_improvement = 0
        pf.ps[-5,:] = pf.ps[pr,:]
        pf.rss[-5] = pf.rss[pr]
        
    # Run child B
    pf.update_potential(pf.ps[-4,:])        
    pf.rss[-4] = pf.get_rss()  
    pr = pa
    if(pf.rss[pa]<pf.rss[pb]):
      pr = pb        
    if(pf.rss[-4]<pf.rss[pr]):
      pf.ps[pr,:] = pf.ps[-4,:]
      pf.rss[pr] = pf.rss[-4]        
      if(pf.rss[pr] < pf.rss[-5]):
        pf.since_improvement = 0
        pf.ps[-5,:] = pf.ps[pr,:]
        pf.rss[-5] = pf.rss[pr]
    
  def pf_children(pa): 
    # Run child A
    pf.update_potential(pf.ps[-3,:])
    pf.rss[-3] = pf.get_rss()  
    if(pf.rss[-3]<pf.rss[pa]):
      pf.ps[pa,:] = pf.ps[-3,:]
      pf.rss[pa] = pf.rss[-3]          
      if(pf.rss[pa] < pf.rss[-5]):
        pf.since_improvement = 0
        pf.ps[-5,:] = pf.ps[pa,:]
        pf.rss[-5] = pf.rss[pa]
    # Run child B
    pf.update_potential(pf.ps[-4,:])
    pf.rss[-4] = pf.get_rss()  
    if(pf.rss[-4]<pf.rss[pa]):
      pf.ps[pa,:] = pf.ps[-4,:]
      pf.rss[pa] = pf.rss[-4]          
      if(pf.rss[pa] < pf.rss[-5]):
        pf.since_improvement = 0
        pf.ps[-5,:] = pf.ps[pa,:]
        pf.rss[-5] = pf.rss[pa]
            
    
  def make_mutation():   
    # Mutation Array
    pf.ps[-2,:] = pf.vc_mutation * pf.ps[-5,:]   # Use best fit parameters
    for n in range(pf.p_width):
      if(pf.ps[-2,n] == 0.0):
        pf.ps[-2,n] = pf.vc_mutation
    
    
  def extinction():
    print("Extinction Event")
    rss_thresh = pf.ext_thresh * numpy.average(pf.rss[:pf.pop_size_d])
    init = True
    for p in range(pf.pop_size_d):
      if(pf.rss[p]<=rss_thresh):    
        if(init):
          pf.ps[-6,:] = pf.ps[p,:]
          pf.ps[-7,:] = pf.ps[p,:]
          init = False
        else:
          for n in range(pf.p_width):
            if(pf.ps[p,n]<pf.ps[-6,n]):
              pf.ps[-6,n] = pf.ps[p,n]
            if(pf.ps[p,n]>pf.ps[-7,n]):
              pf.ps[-7,n] = pf.ps[p,n]
    pf.ps[-8,:] = pf.ps[-7,:] - pf.ps[-6,:]
    
    for p in range(pf.pop_size_d):
      if(pf.rss[p]>rss_thresh):    
        r = numpy.random.rand(pf.p_width)
        pf.ps[p,:] = pf.ps[-5,:] + pf.vc_extinction * (0.5 - r) * pf.ps[-8,:]
        pf.update_potential(pf.ps[p,:])
        pf.rss[p] = pf.get_rss()  
        
        if(pf.rss[p] < pf.rss[-5]):
          pf.since_improvement = 0
          pf.ps[-5,:] = pf.ps[p,:]
          pf.rss[-5] = pf.rss[p]
        print(pf.rss[-3], pf.rss[-4], pf.rss[-5], pf.since_improvement)
    
    
    
    
  def end():
    pf.update_potential(pf.ps[-5,:])
    rss = pf.get_rss() 
    
    print("RSS: ", rss)
    
  
    for bp_id in range(bp.bp_configs_count):
      print('alat:', bp.calc_alat[bp_id])
      print('v0:', bp.calc_v0[bp_id])
      print('e0:', bp.calc_e0[bp_id])
      print('b0:', bp.calc_b0[bp_id])
      print('b0/GPA:', 160.230732254e0 * bp.calc_b0[bp_id])
  
  