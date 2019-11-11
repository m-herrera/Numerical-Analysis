hold on;
x = 1:0.01:6;
y=sin(6-x)./ (sin(5) * sqrt(x));
plot(x, y );
pause (1)
[j,f]=edo2("@(x) -1/x","@(x) 1/(4*x**2)-1","@(x) 0",0.1,1,6,1,0);
plot(j,f)
legend ({"Verdadera", "Con h=0.1"}, "location", "north");
legend hide
legend show
pause (1)
[j,f]=edo2("@(x) -1/x","@(x) 1/(4*x**2)-1","@(x) 0",0.01,1,6,1,0);
plot(j,f)
legend ({"Verdadera", "Con h=0.1","Con h=0.01"}, "location", "north");
legend hide
legend show
pause (1)
[j,f]=edo2("@(x) -1/x","@(x) 1/(4*x**2)-1","@(x) 0",0.001,1,6,1,0);
plot(j,f)
legend ({"Verdadera", "Con h=0.1","Con h=0.01","Con h=0.001"}, "location", "north");
legend hide
legend show
pause (0.1)
[j,f]=edo2("@(x) -1/x","@(x) 1/(4*x**2)-1","@(x) 0",0.0001,1,6,1,0);
plot(j,f)
legend ({"Verdadera", "Con h=0.1","Con h=0.01","Con h=0.001","Con h=0.0001"}, "location", "north");
legend hide
legend show