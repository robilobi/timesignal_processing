%% Project: Plearning melodies
% MelodyAnalysisPupil_preprocessing
% August 2017 Montreal, Roberta Bianco, Aaron Johonson
%
% input: raw data from Eyelink
% output: preprocessed data structure: subject x trial x (data, ratings, '')
%%-----------------------------------

baseline_sample = 400; % 400 millisecond baseline @ 1000Hz
em_surpress_window=200;% region of window to surpress around a blink (in ms) 


fid = fopen('pupil_raw_400ms_27subjects.txt');
pause(.1);
C = textscan(fid, '%s	%s	%u8	%u8	%s	%s	%u8	%u8	%u8	%u8	%n	%n	%n', 'delimiter','\t','TreatAsEmpty',{'.'},'headerLines', 1);  

% C = textscan(fid, '%s %s %u8 %u8 %s %s %u8 %u8 %n %n %n %n %n', 'delimiter','\t','TreatAsEmpty',{'.'},'headerLines', 1);  
fclose(fid);   
    
% now convert cell array into numeric array.  
data_matrix=zeros(size(C{4},1),8);
data_matrix(:,1)=C{4}; % trial - 1 to 16 trials
data_matrix(:,2)=C{7}; % PLEASANT_RATING    
data_matrix(:,3)=C{9}; % RT
data_matrix(:,4)=C{10}; % time stamp counter(in ms)
data_matrix(:,5)=C{11}; % AVERAGE_PUPIL_SIZE
data_matrix(:,6)=C{12}; % LEFT_PUPIL_SIZE
data_matrix(:,7)=C{13}; % RIGHT_PUPIL_SIZE

% find subject and replace string with number
subj=C{1};
subj_unique=unique(subj);
for s=1:size(subj_unique,1)
    subj_name=subj_unique{s}; 
    IndexC = strfind(subj,subj_name);
    IndexC = find(not(cellfun('isempty', IndexC)));
    data_matrix(IndexC,8)=s;
end   

% subject data 
data = cell(size(subj_unique,1),16,3);

for s=1:size(subj_unique,1)
    sT= data_matrix(:,8)==s;
    subj_data_matrix=data_matrix(sT,:);
    
    %% clean pupil data
    for nTrial=1:16
        nT = subj_data_matrix(:,1)==nTrial;
        trial_data_matrix=subj_data_matrix(nT,:);
        nLength=size(trial_data_matrix,1);

        %% psychophysics data
        data{s,nTrial,1}=trial_data_matrix(1,2); % pleasant rating
        data{s,nTrial,2}=trial_data_matrix(1,3); % reaction time
        if nTrial==1
            disp(['Subject ',num2str(s)]);
        end
        disp(['      Cleaning data...Trial ',num2str(nTrial)]);
        
        curveL=trial_data_matrix(:,6)';
        curveR=trial_data_matrix(:,7)';

        %% 1) Interpolate blinks from the signal. 
        % These are characterised by rapid declines towards 0 at blink onset, and rapid rises from 0 back to a regular value at blink offset.
        % Note that the eyelink will automatically remove these zeros and replace with NaN.
        % Need to find blinks, and remove 100ms either side: removes artifacts in eye tracking data 
        % see http://jov.arvojournals.org/article.aspx?articleid=2193211 for example of doing this with 200ms (given our data recored on
        % Eyelink1000 - has higher sample rate so blink algorithm more accurate)
        if nansum(curveL>0) % so long as there is data from the left eye....
            test2 = diff(~isnan([ NaN curveL NaN ]) );
            NumBlockStart = find(test2>0)-0;
            NumBlockEnd = find(test2<0)-1;
            for blinkremove=1:size(NumBlockEnd,2)-1
                if NumBlockEnd(blinkremove)-em_surpress_window>=2
                    curveL(:,NumBlockEnd(blinkremove)-em_surpress_window:NumBlockStart(blinkremove+1)+em_surpress_window)=NaN;
                else
                    curveL(:,2:NumBlockStart(blinkremove+1)+em_surpress_window)=NaN;
                end
            end

            % now for those conditions where the blink occurs at the end of the trial....
            if max(NumBlockEnd)<nLength
                curveL(:,max(NumBlockEnd)-em_surpress_window:nLength)=NaN;
            end
        end

        if nansum(curveR>0) % so long as there is data from the left eye....
            test2 = diff(~isnan([ NaN curveR NaN ]) );
            NumBlockStart = find(test2>0)-0;
            NumBlockEnd = find(test2<0)-1;
            for blinkremove=1:size(NumBlockEnd,2)-1
                if NumBlockEnd(blinkremove)-em_surpress_window>=2
                    curveR(:,NumBlockEnd(blinkremove)-em_surpress_window:NumBlockStart(blinkremove+1)+em_surpress_window)=NaN;
                else
                    curveR(:,2:NumBlockStart(blinkremove+1)+em_surpress_window)=NaN;
                end
            end

            % now for those conditions where the blink occurs at the end of the trial....
            if max(NumBlockEnd)<nLength
                curveR(:,max(NumBlockEnd)-em_surpress_window:nLength)=NaN;
            end
        end


        %% now interpolate to fill in the missing information
        % For each blink, four equally spaced time points are required. t2 is the blink onset; t3 is
        % the blink offset; t1=t2-t3+t2; t4=t3-t2+t3. Based on these four time points and the
        % associated pupil sizes (from the original, unsmoothed signal), a cubic-spline fit is
        % generated. The original signal between t2 and t3 is replaced by the cubic spline. Thus,
        % the signal is left unchanged, except for the blink period.
        % based on http://dx.doi.org/10.6084/m9.figshare.688001
        if nansum(curveL>0) % so long as there is data from the left eye....
            test2 = diff(~isnan([ NaN curveL NaN ]) );
            NumBlockStart = find(test2>0)-0;
            NumBlockEnd = find(test2<0)-1;
            for blinkremove=1:size(NumBlockEnd,2)-1
                t2=NumBlockEnd(blinkremove);
                t3=NumBlockStart(blinkremove+1);
                t1=t2-t3+t2;
                t4=t3-t2+t3;
                if t4<nLength && t1>1
                    eye_i=[curveL(:,t1) curveL(:,t2) curveL(:,t3) curveL(:,t4)];
                elseif t4>=nLength
                    eye_i=[curveL(:,t1) curveL(:,t2) curveL(:,t3) curveL(:,nLength)];
                elseif t1<=1
                    eye_i=[curveL(:,1) curveL(:,t2) curveL(:,t3) curveL(:,t4)];
                end
                eye_spline=interp1([t1 t2 t3 t4],eye_i,t1:1:t4,'spline');
                tmin=t2-t1;
                tmax=t3-t1;
                curveL(:,t2:t3)=eye_spline(:,tmin:tmax);
            end
        end

        if nansum(curveR>0) % so long as there is data from the left eye....
            test2 = diff(~isnan([ NaN curveR NaN ]) );
            NumBlockStart = find(test2>0)-0;
            NumBlockEnd = find(test2<0)-1;
            for blinkremove=1:size(NumBlockEnd,2)-1
                t2=NumBlockEnd(blinkremove);
                t3=NumBlockStart(blinkremove+1);
                t1=t2-t3+t2;
                t4=t3-t2+t3;
                if t4<nLength && t1>1
                    eye_i=[curveR(:,t1) curveR(:,t2) curveR(:,t3) curveR(:,t4)];
                elseif t4>=nLength
                    eye_i=[curveR(:,t1) curveR(:,t2) curveR(:,t3) curveR(:,nLength)];
                elseif t1<=1
                    eye_i=[curveR(:,1) curveR(:,t2) curveR(:,t3) curveR(:,t4)];
                end
                eye_spline=interp1([t1 t2 t3 t4],eye_i,t1:1:t4,'spline');
                tmin=t2-t1;
                tmax=t3-t1;
                curveR(:,t2:t3)=eye_spline(:,tmin:tmax);
            end
        end

        %% 2) Reject artifacts, e.g. by Hampel filtering. 
        DX = 3; % Filter Half Size 
        T = 0; % Threshold 
        X = 1:nLength; % Pseudo Time 
        if size(curveL,2)>nLength
            curveL=curveL(:,1:nLength);
        end
        if size(curveR,2)>nLength
            curveR=curveR(:,1:nLength);
        end
        [~,~,L0] = hampel(X,curveL,DX,T);
        curveL=L0';% Median Filtered Data
        [~,~,R0] = hampel(X,curveR,DX,T);
        curveR=R0';% Median Filtered Data

        %% 3) Optionally smooth the data (depending on your parameters, the Hampel filter might actually act as a smoothing function too). 
        % A popular approach is to use a moving window, e.g. a Savitzky-Golay Filters over an 11ms timeframe. 
        curveL = sgolayfilt(curveL,1,11);
        curveR = sgolayfilt(curveR,1,11);


        %% 4) Baseline correction: 
        % Divide the signal (e.g. timepoint 0 ms to timepoint 3000 ms) by the median pupil size during a baseline period 
        % (e.g. timepoint -200 ms to timepoint 0 ms). This is an important step, as most trackers tend to work with arbitrary numbers, 
        % whereas most papers report changes as a proportional change. 
        baselineL=nanmedian((curveL(:,1:baseline_sample)),2);
        curveL=curveL-baselineL;
        baselineR=nanmedian((curveR(:,1:baseline_sample)),2);
        curveR=curveL-baselineR;

        if size(curveL,2)>7446
           curveL=curveL(:,1:7446);
           nLength=7446;
        end

        if size(curveR,2)>7446
           curveR=curveR(:,1:7446);
           nLength=7446;
        end

        if size(curveL,2)<7446
            curveL = padarray(curveL,[0 (7446-size(curveL,2))],'replicate','post');
            nLength=7446;
        end

        if size(curveR,2)<7446
            curveR = padarray(curveR,[0 (7446-size(curveR,2))],'replicate','post');
            nLength=7446;
        end

        %% 5) plot pupillometry response
%         xval=(1:nLength)-baseline_sample;
%         subplot(4,4,nTrial);
%         plot(xval,curveL,'r',xval,curveR,'b');
%         title(['Trial ',num2str(nTrial)]);
%         legend('Left','Right');
%         ax = gca;
%         ax.XAxisLocation = 'origin';
%         ax.YAxisLocation = 'origin';

        %% 6) save into datamatrix
        % save filtered data into a matrix - use left or right eye depending on
        % data (right eye if both eyes are tracked)
        if nansum(curveL)>0
            data{s,nTrial,3}=curveL;
        else
            data{s,nTrial,3}=curveR;
        end
    end
end

save('data.mat', 'data');  %save data structure subject x trial x (data, liking, '')