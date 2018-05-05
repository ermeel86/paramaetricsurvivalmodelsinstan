"""
c.f. https://commons.wikimedia.org/wiki/File:Mspline_order3.svg
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

def mspline(x, k, i, T):
    """Evaluate the kth order m-spline at x.                                    
                                                                                
    Arguments:                                                                  
    ----------                                                                  
    x : the point at which the spline is evaluated                              
    k : the order of the spline (k>=1)                                          
    i : the interval on which the spline is defined (i<len(T)-k)                
    T : the knots                                                               
                                                                                
    Returns:                                                                    
    --------                                                                    
    The value of the spline at x                                                
                                                                                
    Note: Use mspline_knots(...) to create T.                                   
    """

    if x<T[i]: return 0.
    if x>=T[i+k]: return 0.

    if k==1: return 1/(T[i+1]-T[i])

    d1 = x-T[i]
    d2 = T[i+k]-x

    v = d1*mspline(x, k-1, i, T) + d2*mspline(x, k-1, i+1, T)
    v *= k
    v /= (k-1)*(T[i+k]-T[i])

    return v

def mspline_knots(xmin, xmax, k, m):
    """Create an order list of points in the interval (xmin,xmax)                
    suitable for use as knots in m-splines and i-splines.                       
                                                                                
    Arguments:                                                                  
    ----------                                                                  
    xmin : the lower bound of the interval                                      
    xmax : the upper bount of the interval                                      
    k : the order of the m-spline for which the knots will be used              
    m : the number of intervals into which (xmin,xmax) will be partitioned      
    """

    T = [xmin]*(k-1)
    h = (xmax-xmin)/m
    T.extend(np.arange(xmin, xmax, h))
    T.extend([xmax]*k)

    return T

