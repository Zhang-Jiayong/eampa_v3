######################################################
#  Ben Palmer University of Birmingham 2020
#  Free to use
######################################################

import numpy
import os
from f2py_lib.f_fnc import fnc


################################################
#  Define functions here, or define in fnc
#  better for vectorising?
#  Heaviside step etc are in fnc
################################################
 
 
 
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
    
    
