\\ Pari/GP script to generate one large planar difference set

p=2096993;g=Mod(x,Mod(x^3+3*x-3,p));
m=p^2+p+1;r=g^0;print1(0);
for(j=1,m-1,r*=g;if(poldegree(r.pol)<2,print1(" ",j)));
print("");
