% This script calculate individual levle inter-session CCC 
%% 
% ----load data using load_grandAvg_tep.m
% #### assign data ####
% S1 = GA_TEPs_clean_S1_PRE.individual;
% S2 = GA_TEPs_clean_S2_PRE.individual;
% S3 = GA_TEPs_clean_S3_PRE.individual;

% S1 = GA_TEPs_SHAM_S1_PRE.individual;
% S2 = GA_TEPs_SHAM_S2_PRE.individual;
% S3 = GA_TEPs_SHAM_S3_PRE.individual;
 
S1 = GA_TEPs_ACTIVE_S1_PRE.individual;
S2 = GA_TEPs_ACTIVE_S2_PRE.individual;
S3 = GA_TEPs_ACTIVE_S3_PRE.individual;

alldata = cat(4,S1,S2,S3); % subject x channel x time x session

times = GA_TEPs_clean_S1_PRE.time;   
% post stimulation
[~,index_t1] = min(abs(times*1000 - 20)); % 20ms post pulse
[~,index_t2] = min(abs(times*1000 - 300)); % 300ms post pulse

% pre stimulation
% [~,index_t1] = min(abs(times*1000 - (-300))); 
% [~,index_t2] = min(abs(times*1000 - (-20))); 
time_window = index_t1:index_t2;

nSubjects = size(S1, 1);
resultMatrix = nan(3*nSubjects, 3*nSubjects);
labels = cell(1,nSubjects*3);

for iSubject = 1:nSubjects
    for iSession = 1:3
        labels{3*(iSubject-1) + iSession} = sprintf('%d/%d', iSubject, iSession);
        for jSubject = 1:nSubjects
            for jSession = 1:3
                row = 3*(iSubject-1)+iSession;
                col = 3*(jSubject-1)+jSession;
                resultMatrix(row, col) = concordanceSimilarity(squeeze(alldata(iSubject,:,time_window, iSession)), ...
                    squeeze(alldata(jSubject,:,time_window, jSession)));
            end
        end
    end
end
%% PLOT
fig = figure('color','w');
set(gcf,'position',[0.01 0.01 2 1.8]*1000);
tiledlayout(1,3,'TileSpacing','compact');

cond ='(Active)';
map = colorcet('D1A');

% #### 1st plot: s1 s2 s3 session ####
nexttile;
imagesc(resultMatrix)
% axis equal;
hold on
for i = 1:nSubjects
    xline(3*i+0.5)
    yline(3*i+0.5)
end
xticks(1:size(resultMatrix))
xticklabels(labels)
xtickangle(90);
yticks(1:size(resultMatrix))
yticklabels(' ')
set(gca,'box','off','tickdir','out','linewidth',1,'fontsize',7);
title(['PairWise Similarity ' cond],FontSize=11)

clim([-1 1])
colormap(map)

% #### 2nd plot: average across three sessions without autosimilarity (self-similarity) ####
resultAveragedCross = nan(nSubjects, nSubjects);
for iSubject = 1:nSubjects
    for jSubject = 1:nSubjects
        rows = 3*(iSubject-1) + (1:3);
        cols = 3*(jSubject-1) + (1:3);
        if iSubject == jSubject
            withinSubjectMatrix = resultMatrix(rows, cols); 
            withinSubjectMatrix = withinSubjectMatrix(~eye(3)); % remove the values on the diag (autosimilarity)
            resultAveragedCross(iSubject, jSubject) = mean(withinSubjectMatrix, 'all');
        else
            resultAveragedCross(iSubject, jSubject) = mean(resultMatrix(rows, cols), 'all');
        end
    end
end

nexttile;
imagesc(resultAveragedCross)
hold on
for i = 1:nSubjects
    xline(i+0.5)
    yline(i+0.5)
end
xticks(1:size(resultAveragedCross))
xticklabels(1:nSubjects)
yticks(1:size(resultAveragedCross))
yticklabels(1:nSubjects)
set(gca,'box','off','tickdir','out','linewidth',1,'fontsize',8);
title(['PairWise Similarity Averaged Across Sessions ' cond],FontSize=11)

clim([-1 1])
colormap(map)

% #### 3rd plot: histogram distribution ####
nexttile
histogram(resultAveragedCross(tril(true(nSubjects) & ~eye(nSubjects))),...
    BinEdges=linspace(-1,1,21))% only lower diag values without diag
hold on
histogram(diag(resultAveragedCross),...
    BinEdges=linspace(-1,1,21))
legend('betweenSubjects','withinSubject','Location','northwest')
xlabel('Cosine Similarity')
ylabel('Count')
set(gca,'box','off','tickdir','out','linewidth',1,'fontsize',8);
