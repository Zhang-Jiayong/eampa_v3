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
        