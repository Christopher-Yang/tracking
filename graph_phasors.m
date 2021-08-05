function graph_phasors(data)

group = {'day2','day5','day10'};
names = {'2-day','5-day','10-day'};
allSubj = [13 14 5];
Ngroup = 3;
block{1} = {'B9', 'B11_habit'};
block{2} = {'B24', 'B26_habit'};
block{3} = {'B49', 'B51_habit'};
Nfreq = 6;
col = copper;
col = col(floor((size(col,1)/(Nfreq))*(1:Nfreq)),:);

for i = 1:3
    for j = 1:allSubj(i)
        for k = 1:2
            phasor{i}(:,j,k) = mean(data.(group{i}){j}.(block{i}{k}).cursor.phasors.x_x.ratio,2);
        end
    end
end

figure(1); clf
for k = 1:2
    for j = 1:Ngroup
        subplot(2,3,3*(k-1) + j); hold on
        plot([0 0], [-1 1], 'k')
        plot([-1 1], [0 0], 'k')
        for i = 1:6
            plot(real(phasor{j}(i,:,k)), imag(phasor{j}(i,:,k)), '.', 'MarkerSize', 10, 'Color', col(i,:))
            plot(real(mean(phasor{j}(i,:,k), 2)), imag(mean(phasor{j}(i,:,k), 2)), '.', 'MarkerSize', 30, 'Color', col(i,:))
        end
        axis([-1 1 -1 1])
        axis square
        if k == 1
            title(names{j})
            if j == 1
                ylabel('Late')
            end
        elseif k == 2 && j == 1
            ylabel('Flip')
        end
    end
end

end