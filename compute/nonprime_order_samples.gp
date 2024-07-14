is_max_o = (x1, m, e)-> {
  x0 = x1 ^ 0;
  for(i = 1, #e, if(x1^e[i] == x0, return(0))); x1^m == x0
};
mk_gen = (q)-> {
  z = ffprimroot(ffgen(q, 'y));
  m = q^3 - 1; f = factor(m); d = #f~;
  e = vector(d, i, m/f[d + 1 - i, 1]);
  co = vector(q - 1, i, z^(i - 1));
  r = 0;
  for(b = 1, q - 1,
    for(c = 1, q - 1,
      p = co[1]*x^3 + co[b]*x + co[c];
      x1 = Mod(x, p);
      if(is_max_o(x1, m, e) && polisirreducible(p), return(x1))
    )
  );
  for(a = 1, q -1,
    for(b = 1, q - 1,
      for(c = 1, q - 1,
        p = co[1]*x^3 + co[a]*x^2 + co[b]*x + co[c];
        x1 = Mod(x, p);
        if(is_max_o(x1, m, e) && polisirreducible(p), return(x1))
      )
    )
  )
};
pr_diffset = (q)-> {
  ge = mk_gen(q);
  r=ge^0; print1(0);
  for(j=1,q^2+q,
    r*=ge;
    if(poldegree(r.pol) < 2, print1(" ", j))
  );
  print("")
};
for(q=4, 1024, if(isprimepower(q) > 1, pr_diffset(q)))
