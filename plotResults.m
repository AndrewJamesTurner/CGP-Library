data = csvread("tmp.csv",1,0);

linewidth = 3;
fontsize = 16;

h = figure;
hold on;
plot(data(:,2), "r-", "linewidth", linewidth);
plot(data(:,3), "b--", "linewidth", linewidth);
hold off;

xlim([0,5000]);
ylim([-1,1]);

set(gca, 'fontsize', fontsize);

legend ("Desired Output", "Actual Output", "location", "southeast");

set(gca, 'fontsize', fontsize);
copied_legend = findobj(gcf(),"type","axes","Tag","legend");
set(copied_legend, "fontsize", fontsize);

saveas(h,"tmp",'pdf');
system("pdfcrop tmp.pdf tmp.pdf");
rename("tmp.pdf", "results.pdf");

h = figure;
hold on;
plot(data(:,1), "r-", "linewidth", linewidth);
plot(data(:,2), "b--", "linewidth", linewidth);
hold off;

xlim([0,5000]);
ylim([-1,1]);

legend ("Input", "Desired Output", "location", "southeast");

set(gca, 'fontsize', fontsize);
copied_legend = findobj(gcf(),"type","axes","Tag","legend");
set(copied_legend, "fontsize", fontsize);

saveas(h,"tmp",'pdf');
system("pdfcrop tmp.pdf tmp.pdf");
rename("tmp.pdf", "sin2saw.pdf");
