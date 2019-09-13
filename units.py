# some basic handy conversions. 

def waterPpmvDry2GmKgDry(x):
    """
    Convert water vapor parts per million dry air to Grams H2o/ Kg dry air.
    """
    Mair = 28.9648
    Mh2o = 18.01528
    return (x*Mh2o)/(1e3*Mair)    
def waterGmKgDryToPpmvDry(q):
    """
    Convert water vapor Grams H2o / Kg dry air to ppmv dry air.
    """
    Mair = 28.9648
    Mh2o = 18.01528
    return (q*1e3*Mair)/Mh2o
          
def gasKgKgMoistToDry(q,qh2o):
    """
    Take Kg/Kg moist air to Kg/Kg dry air  
    """
    r = q/(1-qh2o)
    return r

def gasPpmvMoistToDry(x,xh2o):
    """
    Take ppmv moist air to ppmv dry air  
    """
    xd = x/(1-(1e-6*xh2o))
    return xd

def gasKgKgDryToMoist(q,qh2o):
    """
    Take Kg/Kg dry air to Kg/Kg moist air.
    """
    r = q/(1+qh2o)
    return r

def gasPpmvDryToMoist(x,xh2o):
    """
    Take ppmv dry air to ppmv moist air.
    """
    xd = x/(1+(1e-6*xh2o))
    return xd
