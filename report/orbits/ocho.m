path(LOADPATH,"..");
x0=[ 0.1076143733510925E+01  0.0000000000000000E+00  -0.5380718667554619E+00  0.3437068277582447E+00  -0.5380718667554627E+00 -0.3437068277582447E+00  0.0000000000000000E+00 -0.4682662184090647E+00   -0.1099603752075198E+01 0.2341331092045306E+00   0.1099603752075198E+01  0.2341331092045340E+00 ];
global m1 m2 m3;
m1=m2=m3=1;

lsode_options('absolute tolerance',1.e-14);
lsode_options('relative tolerance',1.e-14);
tspan=linspace(0.,2*pi,400);
x=lsode('tbp2d',x0, tspan);
for i=1:length(tspan)
	printf("%17.14e ",tspan(i),x(i,:));
	printf("\n");
end
