%% Project: Plearning melodies
% MelodyAnalysisPupil_preprocessing
% August 2017 Montreal, Roberta Bianco
%
% input: preprocessed data structure: subject x trial x (data, liking, '')
% output: bin data (final matrix: (27 subj*16 trials, 13 bins)
%%-----------------------------------

clear all
%% ================== SETTINGs
checksubdata =1;
selectnote = 0;  %1 if want to look at last N notes

load('data.mat');
subj=27;
noteDur= 428;
ntrial= 16;
seqend = 400+(428*13)+428;

howmanynotes= 13;
baseline_sample = 400; % 400 millisecond baseline @ 1000Hz
seqlength = baseline_sample +(noteDur*howmanynotes);  %baseline, 13 notes
seqstart= baseline_sample;  %start after baseline


%% ================== slect subjects
dataor=data;
% data= data(:,:,3);
% excludedSub = [3, 10, 27];
% data(excludedSub,:) = [];

%% ================= start

%% check single subject data
avg_data_trial = {};
if checksubdata
    for s=1:subj
        data_sub=reshape(dataor(s,:,:),ntrial,3);
        
        for t = 1:ntrial
            
            Cell = {cat(3,data_sub{t,3})};     
            Cell=cell2mat(Cell);
            Cell = Cell(1:seqend);
            
            baselineCell=median(Cell(1:baseline_sample));
            Cell=Cell-baselineCell;
            
            avg_data_trial{s,t}=Cell;
%             figure(s);
%             xval=(1:(size(Cell,2)))-baseline_sample;
%             subplot(4,4,t);
%             plot(xval,Cell,'k-');
%             ax = gca;
%             ax.XAxisLocation = 'origin';
%             ax.YAxisLocation = 'origin';
        end
%         path = 'pic/'; mkdir(path);
%         filename = [path 'trials_' num2str(s) '.jpg'];
%         saveas(figure(s),filename)
%         close all
    end

end

%% extract data for ANOVA with all trials
rob=avg_data_trial;

slots = baseline_sample:noteDur:((howmanynotes+1)*noteDur);
for i = 1:subj
    for t = 1:ntrial
        avgtrial(i,t)= nanmean(rob{i,t});
        for b = 2:1:size(slots,2)
            bintrial{i,b-1,t}=nanmean(rob{i,t}(slots(b-1):slots(b)));
        end
    end
end

bintrial_save = [];
da=1;
la=subj;
for t = 1:ntrial
    bintrial_save(da:la,:)=cell2mat(bintrial(:,:,t));
    da=da+subj;
    la=la+subj;
end
s=(1:subj)';
subid=repmat(s,ntrial,1)';
trialid = repelem(1:ntrial,subj);
bintrial_final=[bintrial_save, subid', trialid'];


%csvwrite(['avgtrial_' name '.txt'], avgtrial)
csvwrite(['pupil_bin_singletrial.txt'],bintrial_final) %(final matrix: (27SUbj*16Trials, 13bins)