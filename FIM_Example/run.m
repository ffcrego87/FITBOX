%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run.m - main file to run the FDI example
%       - data from each Monte-Carlo simulation is saved onto a cell
%       variable and a html file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);
Ts            = 0.01;                     %sample time [s]
mc_runs       = 1;                        %number of Monte-Carlo runs
nofaults      = 3;                       
%auxiliary variables to count the detections/isolations
dtc           = zeros(nofaults,mc_runs);
iso           = dtc;
max_dtc       = zeros(nofaults,1);
min_dtc       = max_dtc+inf;
max_iso       = zeros(nofaults,1);
min_iso       = max_iso+inf;
time_out_dtc  = dtc;
time_out_iso  = dtc;
falseDetect   = dtc;
overlimit_dtc = dtc;
overlimit_iso = dtc;
idx           = 1:nofaults;
tFault        = ones(nofaults,1);         %fault occurs at t=tFault s 
for J=1:mc_runs
    for I=1:nofaults
        T_init = 0; 
        T_limit = tFault(I)+101;
        time  = (T_init:Ts:T_limit)';
        fsgn1 = ones(size(time));
        fsgn2 = ones(size(time));
        fsgn3 = ones(size(time));
        if I < 3
            eval(['fsgn' num2str(I) '=[ones(tFault(I)/Ts,1);zeros((T_limit-tFault(I))/Ts+1,1)];']);
        else
            rampLength = 30; %fault level increases for 30s
            fsgn3=[ones(1,tFault(I)/Ts),... %
                1:-(1/(rampLength/Ts)):max([0,1-(T_limit-tFault(I))/rampLength]), ...
                zeros(1,((T_limit-tFault(I)-rampLength)/Ts)*((T_limit-tFault(I))>rampLength))]';
        end
        bladeData; % reset data
        sim('blade.mdl');
        sij = strcat('(',num2str(idx(I)),',',num2str(J),')');
        save(sij,'x_lb','x_ub','beta','T_limit','T_init','Ts','seed','mode','false_dtc','fault','u');
        if (~isempty(find(mode(:,2)==2,1)) && (mode(find(mode(:,2)==2,1),1)-tFault(idx(I),1))<0) ||...
                (~isempty(find(mode(:,2)==3,1)) && (mode(find(mode(:,2)==3,1),1)-tFault(idx(I),1))<0)
            %false detection
            falseDetect(idx(I),J) = 1;
        elseif isempty(find(mode(:,2)==3,1))
            %fault was not isolated
            time_out_iso(idx(I),J) = 1;
            if isempty(find(mode(:,2)==2,1))
                %fault was not detected
                time_out_dtc(idx(I),J) = 1;
            else 
                dtc(idx(I),J) = mode(find(mode(:,2)==2,1),1)-tFault(idx(I),1);
                if dtc(idx(I),J)>max_dtc(I)
                    max_dtc(I) = dtc(idx(I),J);
                end
                if dtc(idx(I),J)<min_dtc(I)
                    min_dtc(I) = dtc(idx(I),J);
                end
            end
        else
            %fault was isolated
            if isempty(find(mode(:,2)==2,1)) || find(mode(:,2)==2,1)>find(mode(:,2)==3,1)
                %fault was isolated at the same sampling time it was
                %detected
                dtc(idx(I),J) = mode(find(mode(:,2)==3,1),1)-tFault(idx(I),1);
                iso(idx(I),J) = mode(find(mode(:,2)==3,1),1)-tFault(idx(I),1);
            else 
                dtc(idx(I),J) = mode(find(mode(:,2)==2,1),1)-tFault(idx(I),1);
                iso(idx(I),J) = mode(find(mode(:,2)==3,1),1)-tFault(idx(I),1);
            end
            if iso(idx(I),J)>max_iso(I)
                max_iso(I) = iso(idx(I),J);
            end
            if iso(idx(I),J)<min_iso(I)
                min_iso(I) = iso(idx(I),J);
            end
            if dtc(idx(I),J)>max_dtc(I)
                max_dtc(I) = dtc(idx(I),J);
            end
            if dtc(idx(I),J)<min_dtc(I)
                min_dtc(I) = dtc(idx(I),J);
            end
            if min_dtc(I)<0
                disp('')
            end
        end
                
    end
%%
results = cell(numel(idx)+1,10);
results{1,2} = 'Detected';
results{1,3} = 'Median detection time';
results{1,4} = 'Missed detections';
results{1,5} = 'Max Detection Time';
results{1,6} = 'Min Detection Time';
results{1,7} = 'Isolated';
results{1,8} = 'Median isolation time';
results{1,9} = 'False detections';
results{1,10} = 'Missed isolations';
results{1,11} = 'Max Isolation Time';
results{1,12} = 'Min Isolation Time';
for I=1:numel(idx)
    index_dtc=find(dtc(idx(I),:)>=0);
    index_iso=find(iso(idx(I),:)>=0);
    results{idx(I)+1,1} = ['Fault no. ' num2str(idx(I))];
    results{idx(I)+1,2} = sum(dtc(idx(I),:)>=0);
    results{idx(I)+1,3} = median(dtc(idx(I),index_dtc));
    results{idx(I)+1,4} = sum(time_out_dtc(idx(I),:));
    results{idx(I)+1,5} = max_dtc(I);
    results{idx(I)+1,6} = min_dtc(I);
    results{idx(I)+1,7} = sum(iso(idx(I),:)>=0);
    results{idx(I)+1,8} = median(iso(idx(I),index_iso));
    results{idx(I)+1,9} = sum(falseDetect(idx(I),:));
    results{idx(I)+1,10} = sum(time_out_iso(idx(I),:));
    results{idx(I)+1,11} = max_iso(I);
    results{idx(I)+1,12} = min_iso(I);
end
save results results
cols = results(1,1:end);
rows = results(2:end,1);
mat = cell2mat(results(2:end,2:end));
GTHTMLtable(mat,cols,'%s',rows,'%s','save')
end