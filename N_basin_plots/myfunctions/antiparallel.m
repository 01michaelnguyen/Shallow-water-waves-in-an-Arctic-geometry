%relation to get antiparallel (trivial) eigenvalues (if inner basins are the same size)
function sym=antiparallel(w,b,s_0,s_1,i)
d_i_plus=(-b^2+4*(b*(i/w)-(i^2)))^(1/2);
sym=tan(d_i_plus*(s_1-s_0)/2)-d_i_plus/(2*(i+b/2));
end