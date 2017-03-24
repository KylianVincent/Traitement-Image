clf();
figure(0);
y1=read("timelog1.txt", -1, 1);
y2=read("timelog2.txt", -1, 1);
y3=read("timelog3.txt", -1, 1);
y4=read("timelog4.txt", -1, 1);
y5=read("timelog5.txt", -1, 1);
y=(y1'+y2'+y3'+y4'+y5')/5;

x=linspace(0, 10, 100);
plot(x, y);
[a, b, r] = reglin(x, y);
plot(x, a*x+b, "r");
drawaxis(x=0:10, y=0, dir='u', tics='v')

ag=gca();
ag.x_label.text="sigma";
ag.x_label.font_size = 3;
ag.y_label.text="Ecart relatif en %";
ag.y_label.font_size = 3;
ag.title.text="Ecart de temps relatif entre filtrage fréquentiel et spatial en fonction de sigma (W = 3sigma)";
ag.title.font_size = 4;
legend('(Tfreq-Tspat)/Tfreq (%)', 'Droite de régression : y='+string(a(1))+'*x+'+string(b), pos=2, boxed=%f);
