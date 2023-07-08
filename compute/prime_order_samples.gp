\\ Pari/GP script to generate prime order sample sets.

forprime(p=2,4096,\
    ge=Mod(x,minpoly(ffprimroot(ffgen([p,3]))));\
    r=ge^0;print1(0);\
    for(j=1,p^2+p,\
        r*=ge;\
        if(poldegree(r.pol)<2,print1(" ",j)));\
    print("")\
)
