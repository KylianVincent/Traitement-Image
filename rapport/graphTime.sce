clf();
figure(0);
colordef('white');
y=read("timelog.txt", -1, 1);
plot(y05);
a=gca();
a.x_label.text="sigma";
a.y_label.text="Temps en s";
a.title.text="Ecart de temps entre filtrage fréquentiel et spatial en fonction de sigma (W = 3sigma)";
a.title.font_size = 4;
