
clear all;
close all;
clc;

pkg load signal

i = 1;
saw = -1;

for x = 0:0.01:50

	inputs(i) = i;
	sinWave(i) = sin(x);
	sawWave(i) = sawtooth(x);

	i = i + 1;
	saw = saw + 0.01;

	if saw > 1
		saw = -1;
	end
	
	

end

h = figure;
hold on;
plot(inputs, sinWave, 'r');
plot(inputs, sawWave);
hold off;

file = "sin2saw.csv";
delete(file);

csvwrite(file, [1,1,length(sinWave)]);

for i=1:length(sinWave)
	csvwrite(file, [sinWave(i), sawWave(i)], "-append");
end


