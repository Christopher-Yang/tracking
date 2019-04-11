function graph_phasor(data, block_name, groups, graph_name, subj_rot, subj_rot_i)

    Nfreqs = length(data.(groups{2}).(subj_rot_i{1}).(block_name{1}).freqs_x);
    Nreps = size(data.(groups{2}).(subj_rot_i{1}).(block_name{1}).x_x_all.ratio,2);
    leg.rot = subj_rot;
    leg.rot_i = subj_rot_i;
    names = {'x_x','y_y','x_y','y_x'};
    names2 = {'x_x_all','y_y_all','x_y_all','y_x_all'};
    for i = 1:length(groups)
        Nsubj = length(fieldnames(data.(groups{i})))-1;
        g1 = NaN(Nfreqs,length(block_name),Nsubj);
        g2 = NaN(Nfreqs,Nreps,length(block_name),Nsubj);        
        for j = 1:length(names)
            d.(groups{i}).(names{j}).amp = g1;
            d.(groups{i}).(names{j}).phase = g1;
            d.(groups{i}).(names{j}).amp_all = g2;
            d.(groups{i}).(names{j}).phase_all = g2;
        end
    end
    for i = 1:length(groups)
        for j = 1:Nsubj
            subj_name = fieldnames(data.(groups{i}));
            subj_name = subj_name(1:Nsubj);
            for k = 1:length(block_name)
                for l = 1:length(names)
                    d.(groups{i}).(names{l}).amp(:,k,j) = data.(groups{i}).(subj_name{j}).(block_name{k}).(names{l}).amplitude;
                    d.(groups{i}).(names{l}).phase(:,k,j) = data.(groups{i}).(subj_name{j}).(block_name{k}).(names{l}).phase;
                    d.(groups{i}).(names{l}).amp_all(:,:,k,j) = data.(groups{i}).(subj_name{j}).(block_name{k}).(names2{l}).amplitude;
                    d.(groups{i}).(names{l}).phase_all(:,:,k,j) = data.(groups{i}).(subj_name{j}).(block_name{k}).(names2{l}).phase;
                end
            end
        end
    end

    group_names = {'Rotation','Mirror Reversal'};
    freqs_x = {'0.10 Hz','0.25 Hz','0.55 Hz','0.85 Hz','1.15 Hz','1.55 Hz','2.05 Hz'};
    freqs_y = {'0.15 Hz','0.35 Hz','0.65 Hz','0.95 Hz','1.45 Hz','1.85 Hz','2.15 Hz'};

%% across all frequencies and blocks; subjects averaged together

    figure('Name','Average Response (Rotation)','NumberTitle','off');
    subplot(1,2,1);
    polarplot(data.rot.avg.x_x.d.phase',data.rot.avg.x_x.d.amplitude','-o');
    hold on;
    title('X movements');
    rlim([0 1]);
    
    subplot(1,2,2);
    polarplot(data.rot.avg.y_y.d.phase',data.rot.avg.y_y.d.amplitude','-o');
    hold on;
    title('Y movements');
    rlim([0 1]);
    legend(graph_name, 'Position',[0.9 0.4 0.1 0.2]);
    
    figure('Name','Average Response (Mirror Reversal)','NumberTitle','off');
    subplot(1,2,1);
    polarplot(data.rot_i.avg.x_x.d.phase',data.rot_i.avg.x_x.d.amplitude','-o');
    hold on;
    title('X movements');
    rlim([0 1]);
    
    subplot(1,2,2);
    polarplot(data.rot_i.avg.y_y.d.phase',data.rot_i.avg.y_y.d.amplitude','-o');
    hold on;
    title('Y movements');
    rlim([0 1]);
    legend(graph_name, 'Position',[0.9 0.4 0.1 0.2]);

%% across all frequencies and blocks; keep subject and rotation constant    

%     for i = 1:length(groups)
%         polarplot(d.x_x.phase(:,:,1,1),d.x_x.amp(:,:,1,1),'-o');
    
%% graph across all frequencies and subjects; keep training block and rotation constant
%     for i = 1:length(groups)
%         figure('Name',[group_names{i},' (X)'],'NumberTitle','off');
%         for j = 1:length(block_name)
%             subplot(2,3,j); 
%             polarplot(squeeze(d.(groups{i}).x_x.phase(:,j,:))',squeeze(d.(groups{i}).x_x.amp(:,j,:))','.','MarkerSize',15);
% %             polarplot(squeeze(d.(groups{i}).x_x.phase(4,j,:))',squeeze(d.(groups{i}).x_x.amp(4,j,:))','.');
%             hold on;
%             title(graph_name{j});
%             rlim([0 1.75]);
% %             rlim([0 0.5]);
%         end
%         legend(freqs_x,'Position',[0.9 0.4 0.1 0.2]);
%         figure('Name',[group_names{i},' (Y)'],'NumberTitle','off');
%         for j = 1:length(block_name)
%             subplot(2,3,j);
%             polarplot(squeeze(d.(groups{i}).y_y.phase(:,j,:))',squeeze(d.(groups{i}).y_y.amp(:,j,:))','.','MarkerSize',15);
% %             polarplot(squeeze(d.(groups{i}).y_y.phase(4,j,:))',squeeze(d.(groups{i}).y_y.amp(4,j,:))','.');
%             hold on;
%             title(graph_name{j});
%             rlim([0 1.75]);
% %             rlim([0 0.5]);
%         end
%         legend(freqs_y,'Position',[0.9 0.4 0.1 0.2]);
%     end
    
%% graph all trials across all frequencies and subjects; keep trianing block and rotation constant
%     freq = 4;
%     for i = 1:length(groups)
%         subjs = length(fieldnames(data.(groups{i})))-1;
%         figure('Name',[group_names{i},' (X)'],'NumberTitle','off');
%         for j = 1:length(block_name)
%             subplot(2,3,j);
%             polarplot(squeeze(d.(groups{i}).x_x.phase_all(freq,:,j,subjs)),squeeze(d.(groups{i}).x_x.amp_all(freq,:,j,subjs)),'.','MarkerSize',15);
%             hold on;
%             title(graph_name{j});
%             rlim([0 2]);
% %             rlim([0 1]);
%         end
%         legend(leg.(groups{i}),'Position',[0.9 0.4 0.1 0.2]);
%         figure('Name',[group_names{i},' (Y)'],'NumberTitle','off');
%         for j = 1:length(block_name)
%             subplot(2,3,j);
%             polarplot(squeeze(d.(groups{i}).y_y.phase_all(freq,:,j,subjs)),squeeze(d.(groups{i}).y_y.amp_all(freq,:,j,subjs)),'.','MarkerSize',15);
%             hold on;
%             title(graph_name{j});
%             rlim([0 2]);
% %             rlim([0 1]);
%         end
%         legend(leg.(groups{i}),'Position',[0.9 0.4 0.1 0.2]);
%     end
end