#path(LOADPATH,"..");
function orbita(x0,alpha)
#x0=[  1.1154398795E+00   0.0000000000E+00  -5.5771993973E-01   3.8390157495E-01  -5.5771993973E-01  -3.8390157495E-01 0.0000000000E+00  -5.4749389715E-01  -1.1930727283E+00   2.7374694857E-01   1.1930727283E+00   2.7374694857E-01 ];
global m1 m2 m3 alpha;
m1=m2=m3=1;
#alpha=1.3424242750;
lsode_options('absolute tolerance',1.e-14);
lsode_options('relative tolerance',1.e-14);
tspan=linspace(0.,2*pi,400);
x=lsode('tbp2d',x0, tspan);
for i=1:length(tspan)
	printf("%17.14e ",tspan(i),x(i,:));
	printf("\n");
end

return