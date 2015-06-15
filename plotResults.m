data = csvread("tmp.csv",1,0);

h = figure;
hold on;
plot(data(:,1), "r");
plot(data(:,2), "g");
plot(data(:,3), "b");
hold off;

saveas(h,"tmp",'pdf');
%system("pdfcrop tmp.pdf tmp.pdf");
%rename("tmp.pdf", strcat(type,'-',ff,'-boxBlots.pdf'));
