from scipy import linspace
from scipy import pi,sqrt,exp
from scipy.special import erf

from pylab import plot,show,legend

def pdf(x):
    return 1/sqrt(2*pi) * exp(-x**2/2)

def cdf(x):
    return (1 + erf(x/sqrt(2))) / 2

def skew(x,e=0,w=1,a=0):
    t = (x-e) / w
    return 2 / w * pdf(t) * cdf(a*t)
    # You can of course use the scipy.stats.norm versions
    # return 2 * norm.pdf(t) * norm.cdf(a*t)


n = 2**10

mu=3000.0
sigma=300.0

x = linspace(2000,4000,n) 

pl=[]
la=[]

i=0
for a in range(-20,21,4):
    d = a/sqrt(1.0 + a**2)
    #a = d/sqrt(1 - d**2)
    w = sigma/sqrt(1.0 - (2.0/pi)*d**2)
    e = mu - w*d*sqrt(2.0/pi)
    p = skew(x,e,w,a)
    pl.append(plot(x,p))
    la.append('a='+str(a))
    i=i+1

legend(la)
show()

