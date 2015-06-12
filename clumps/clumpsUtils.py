


def to_gauss((a,b,x0,y0,v0,phi,sx,sy,sv,dvx,dvy)):
   sphi2=np.square(np.sin(phi))
   cphi2=np.square(np.cos(phi))
   s2phi=np.sin(2*phi)
   sx2=np.square(sx)
   sy2=np.square(sy)
   sv2=np.square(sv)
   La=cphi2/sx2 + sphi2/sy2 + np.square(dvx)/sv2
   Ld=sphi2/sx2 + cphi2/sy2 + np.square(dvy)/sv2
   Lb=-s2phi/(2*sx2) + s2phi/(2*sy2) + dvx*dvy/sv2
   Lc=-dvx/sv2
   Le=-dvy/sv2
   Lf=1.0/sv2
   L=np.array([[La,Lb,Lc],[Lb,Ld,Le],[Lc,Le,Lf]])
   mu=[x0,y0,v0]
   return (a,b,mu,L)



