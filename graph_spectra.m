
groups = {'day2', 'day5', 'day10'};
graph_names = {'2-day', '5-day', '10-day'};
% blocks = {'B1_baseline', 'B1_baseline', 'B1_baseline'};
% blocks = {'B3', 'B3', 'B3'};
% blocks = {'B9', 'B24', 'B49'};
blocks = {'B11_habit', 'B26_habit', 'B51_habit'};
allSubj = [13 14 4];

for i = 1:3
    for j = 1:allSubj(i)
        amplitudeX{i}(:,j) = mean(abs(data.(groups{i}){j}.(blocks{i}).cursor.x_fft), 2);
        amplitudeY{i}(:,j) = mean(abs(data.(groups{i}){j}.(blocks{i}).cursor.y_fft), 2);
    end
end

a = data.day2{1}.B1_baseline;
fs = 130;
Nsamples = length(a.time);
x_axis = fs*(0:Nsamples/2)/Nsamples;
freqX = a.freqX;
ampX = a.ampX;
freqY = a.freqY;
ampY = a.ampY;

figure(5); clf
for i = 1:3
    subplot(1,3,i); hold on
    plot(freqX, ampX, '.k', 'MarkerSize', 20)
    plot(x_axis, mean(amplitudeX{i},2))
    axis([0 1.7 0 0.024])
    title(graph_names{i})
    if i == 1
        ylabel('Amplitude (m)')
    elseif i == 2
        xlabel('Frequency (Hz)')
    end
end

figure(6); clf
for i = 1:3
    subplot(1,3,i); hold on
    plot(freqY, ampY, '.k', 'MarkerSize', 20)
    plot(x_axis, mean(amplitudeY{i},2))
    axis([0 1.7 0 0.024])
    title(graph_names{i})
    if i == 1
        ylabel('Amplitude (m)')
    elseif i == 2
        xlabel('Frequency (Hz)')
    end
end
