#input between 0 and 1 converted to normalised linearized output, with respect to polynomial factors

import numpy as np
import cmath
from scipy import interpolate

def linz(input):
    
    red=    [1.127339114,-0.163569181,0.017654147]
    green=  [1.107403255,-0.138534706,0.009456395]
    blue=   [1.095060654,-0.140152798,0.022707094]

    a=red[0]
    b=red[1]
    c=red[2]
    red_min=(-b+cmath.sqrt(b**2-4*a*c+4*a*0))/(2*a)
    red_max=(-b+cmath.sqrt(b**2-4*a*c+4*a*1))/(2*a)
    xred=(-b+cmath.sqrt(b**2-4*a*c+4*a*input[0]))/(2*a)
    red_output=xred.real
    #(xred.real-red_min.real)*(1/(red_max.real-red_min.real))

    a=green[0]
    b=green[1]
    c=green[2]
    green_min=(-b+cmath.sqrt(b**2-4*a*c+4*a*0))/(2*a)
    green_max=(-b+cmath.sqrt(b**2-4*a*c+4*a*1))/(2*a)
    xgreen=(-b+cmath.sqrt(b**2-4*a*c+4*a*input[1]))/(2*a)
    green_output=xgreen.real
    (xgreen.real-green_min.real)*(1/(green_max.real-green_min.real))
    
    a=blue[0]
    b=blue[1]
    c=blue[2]
    blue_min=(-b+cmath.sqrt(b**2-4*a*c+4*a*0))/(2*a)
    blue_max=(-b+cmath.sqrt(b**2-4*a*c+4*a*1))/(2*a)
    xblue=(-b+cmath.sqrt(b**2-4*a*c+4*a*input[2]))/(2*a)
    blue_output=xblue.real
    (xblue.real-blue_min.real)*(1/(blue_max.real-blue_min.real))
    
    return [red_output,green_output,blue_output]

#def linzs(input):
#red=np.array([0.81,0.98,1.27,2.08,3.25,5.09,7.54,10.67,14.36,18.81,23.65,29.32,35.15,41.54,48.32,56.06,64.44,75.56])
#red=red/max(red)
#green=np.array([0.98,1.38,2.79,5.52,9.99,16.44,25.33,36.59,50,65.11,81.72,100.13,120,141.64,163.34,189.38,222.15,258.75])
#green=green/max(green)
#blue=np.array([1.1,1.09,1.39,1.99,3.01,4.35,6.17,8.68,11.67,15.19,18.87,22.87,27.6,32.41,37.62,42.86,49.39,58.29])
#blue=blue/max(blue)
#
#print red
#print green
#print blue
#
#x = np.linspace(0,1,18)
#xx = np.linspace(0,1,100)
#yyr = spline(x,red,xx);
#yyg = spline(x,green,xx);
#yyb = spline(x,blue,xx);
#
#yy = interpolate.interp1d(x, red, kind='cubic')(xx)
#print yy
#print len(yy)


######
#print linz([0,0,0])
#print linz([ 0.5,0.5,0.5 ])
#print linz([1,1,1])



