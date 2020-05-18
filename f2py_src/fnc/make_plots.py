from f_fnc import fnc
import numpy
import matplotlib.pyplot as plt




class make_plots:


  
  @staticmethod
  def run():
    x = numpy.linspace(1.0, 7, 100)
    
    # LJ
    y = fnc.lennard_jones_v(x, numpy.asarray([2.3,3.5], dtype=numpy.float32))
    make_plots.plot('lennard_jones', 'Lennard-Jones', x, y) 
    
    # Summed spline (Olsson et al)
    p = [0.976,-165.0,1.15,-78.499908, 1.216,-78.15495,1.650,1.8679553]
    y = fnc.summed_spline_v(x, numpy.asarray(p, dtype=numpy.float32))
    make_plots.plot('summed_spline_pair', 'Summed Spline', x, y) 
    
    y = fnc.summed_spline_v(x, numpy.asarray([0.963,-11.0828,1.284,0.013905,1.685,-0.447541], dtype=numpy.float32))
    make_plots.plot('summed_spline_density', 'Summed Spline - Density', x, y) 
    
    
    x = numpy.linspace(0.0, 7, 100)
    y = fnc.simple_spline_v(x, numpy.asarray([50.0,-2.0,-0.5,-0.011,0.0,0.0], dtype=numpy.float32))
    make_plots.plot('simple_spline_density', 'Simple Spline - Density', x, y) 
    
    
    
    x = numpy.linspace(0.0, 1.0, 100)
    y = fnc.fs_embedding_v(x, numpy.asarray([50.0,-2.0,-0.5,-0.011,0.0,0.0], dtype=numpy.float32))
    make_plots.plot('fs_embedding', 'Finnis-Sinclair Embedding', x, y) 

  @staticmethod
  def plot(plot_file, plot_name, x, y, y_a = None, y_b = None):
    plt.clf()   
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig, axs = plt.subplots(1, 1, figsize=(7,5))
    fig.tight_layout(pad=5.0)    
    fig.suptitle(plot_name)
    axs.plot(x, y, color='k', ls='solid', label='potential')
    if(type(y_a) == numpy.ndarray):
      axs.plot(x, y_a, color='k', ls='dashed', label='repulsion')
    if(type(y_b) == numpy.ndarray):
      axs.plot(x, y_b, color='k', ls='dotted', label='attraction')
    axs.set_ylim(-5.0,5.0)
    axs.legend()
    plt.savefig('plots/' + plot_file + '.eps', format='eps')




make_plots.run()