N=1000; %the number of individuals in the population
time=1000; %the number of generations to run the simulation
t_split = 40; %if patch='TRUE', t_split is the number of generations (40 in this case) the population stays in one environment or the other
radius=1; %distance from the hypersphere origin
nugen = 10; %number of genetic (HT) mutations for each generation
nuepi = 90; %number of epigenetic (LT) mutations for each generation
epiback = 0.0205; %epigenetic reversion rate
mlimitgen=2; %maximum effect a genetic mutation can have
mlimitepi = 0.3; %maximum effect an epigenetic mutation can have
oligo = 1; %always leave at 1 for this version of the model. It is a switch that determines number of dimensions of the hypersphere (in our case 5 dimensions)
epigenetics = 'TRUE'; %switch to toggle effects of epigenetics on or off
patch = 'TRUE'; %switch to determine if population switches between one environment and another. If TRUE, population goes in and out of selection environment at a frequency of t_split generations. If FALSE, only stays in the selection environment.
mode = 'negative'; %When patch = 'TRUE', this will run the model in two modes that determines the selection rules for environment 2: random or negative. 'negative' turns on randomly sampling individuals weighted by the reciprocal number of genetic mutations when in environment 2 (stabilizing selection). 'random' turns on randomly sampling indviduals in environment 2 with no weighting. Random sampling weighted by fitness always occurrs in environment 1 
runnumber = 1; %replicate number of model run

patch_intervals = cell(length(t_split), 1); %Create cell array of matrices to store interval runs
patch_removed_intervals = cell(length(t_split), 1); %Create cell array of matrices to store interval runs of removed individuals. These are individuals that were not sampled to go into the next generation.
results_mat = cell(length(t_split), 1); %Create cell array of matrices to store results matrix runs
results_removed_mat = cell(length(t_split), 1); %Create cell array of matrices to store results matrix runs of removed individuals
epi_mat = cell(length(t_split), 1); %Create cell array of matrices to store epigenetic distribution run
gen_mat = cell(length(t_split), 1); %Create cell array of matrices to store genetic distribution run
epi_removed_mat = cell(length(t_split), 1); %Create cell array of matrices to store epigenetic distribution run
gen_removed_mat = cell(length(t_split), 1); %Create cell array of matrices to store genetic distribution run


%Driver is a vector of 0's and 1's and has length = time in generations. The loop below will step through
%every element of Driver. When patch = 'TRUE' and Driver == 1, then select
%= 1, and the population is randomly sampled to the next generation weighted
%by fitness, i.e. is in environment 1. When patch = 'TRUE', mode = 'negative', and
%Driver == 0, then the population is randomly sampled to the next
%generation weighted by the reciprocal of the number of genetic mutations (stabilizing selection),
%i.e. is in environment 2.
Driver=ones(time,1);
for q=t_split:t_split*2:time
    Driver(q:q+t_split-1)=0; %generate time intervals to switch back and forth from
end

%Run the model
[results, population_disgen, population_disepi, results_removed, population_removed_disgen, population_removed_disepi]=EpiGen_simulation(N,time,radius,nugen,nuepi, ...
    epiback, mlimitgen, mlimitepi, epigenetics, patch, oligo, Driver, mode, runnumber, t_split);

%PULL DISTRIBUTION OF GENETIC EFFECTS FOR 'KEPT' (NON-REMOVED) POPULATIONS
dis_gen = zeros(N, time);
[mgen, ngen] = size(population_disgen);
for g=1:mgen

    %check if cell array is empty
    gtmp = {population_disgen{g,:}};
    tf = cellfun('isempty',gtmp);
    gtmp(tf) = {0};
    %gtmp = [gtmp{:}];
    for x=1:length(gtmp)
        %gtmp{x}
        dis_gen(g,x) = str2double(gtmp{x});
    end
end
%Save each genetic matrix into cell array
dis_gen(isnan(dis_gen))=0;
dlmwrite(['disgen_' num2str(runnumber) '_' mode '.txt'],dis_gen, 'delimiter', '\t');
gen_mat = dis_gen;


%PULL DISTRIBUTION OF GENETIC EFFECTS FOR 'REMOVED' POPULATIONS
dis_removed_gen = zeros(N, time);
[row, col] = size(population_removed_disgen);
for gr=1:row

    grtmp = {population_removed_disgen{gr,:}};
    tf = cellfun('isempty',grtmp);
    grtmp(tf) = {0};
    %gtmp = [gtmp{:}];
    for x=1:length(grtmp)
        %gtmp{x}
        dis_removed_gen(gr,x) = str2double(grtmp{x});
    end
end
%Save each genetic matrix into cell array
dis_removed_gen(isnan(dis_removed_gen))=0;
dlmwrite(['dis_removed_gen_' num2str(runnumber) '_' mode '.txt'],dis_removed_gen, 'delimiter', '\t');
gen_removed_mat = dis_removed_gen;


%PULL DISTRIBUTION OF EPIGENETIC EFFECTS FOR 'KEPT' POPULATIONS
dis_epi = zeros(N, time);
[mepi, nepi] = size(population_disepi);
for e=1:mepi
    etmp = {population_disepi{e,:}};
    tf = cellfun('isempty',etmp);
    etmp(tf) = {0};
    etmp = [etmp{:}];
    dis_epi(e,:) = etmp;
end
%Save each epigenetic matrix into cell array
dlmwrite(['disepi_' num2str(runnumber) '_' mode '.txt'],dis_epi, 'delimiter', '\t');
epi_mat = dis_epi;


%PULL DISTRIBUTION OF EPIGENETIC EFFECTS FOR 'REMOVED' POPULATIONS
dis_removed_epi = zeros(N, time);
[row, col] = size(population_removed_disepi);
for er=1:row
    ertmp = {population_removed_disepi{er,:}};
    tf = cellfun('isempty',ertmp);
    ertmp(tf) = {0};
    ertmp = [ertmp{:}];
    dis_removed_epi(er,:) = ertmp;

end
%Save each epigenetic matrix into cell array
dlmwrite(['dis_removed_epi_' num2str(runnumber) '_' mode '.txt'],dis_removed_epi, 'delimiter', '\t');
epi_removed_mat = dis_removed_epi;



t_hold = 0;


tmp_matrix = zeros(t_split, 5);
tmp_matrix_r = zeros(t_split, 5);


for s=t_split:t_split*2:time

    if Driver(s-1) == 1

        %disp('Driver=1')

        %SELECT KEPT 
        tmp_matrix(s,1) = results.co(s,1);
        tmp_matrix(s,2) = results.co(s,2);

        %SELECT REMOVED 
        tmp_matrix_r(s,1) = results_removed(s,1);
        tmp_matrix_r(s,2) = results_removed(s,2);

        %NONSELECT KEPT
        tmp_matrix(s,3) = results.co(s+(t_split)-1,1);
        tmp_matrix(s,4) = results.co(s+(t_split)-1,2);

        %NONSELECT REMOVED
        tmp_matrix_r(s,3) = results_removed(s+(t_split)-1,1);
        tmp_matrix_r(s,4) = results_removed(s+(t_split)-1,2);
        
        %RUNNUMBER
        tmp_matrix(s,5) = runnumber;
        tmp_matrix_r(s,5) = runnumber;

    end
end

%Save each interval run matrix into cell array
patch_intervals = tmp_matrix;
dlmwrite(['patch_intervals_' num2str(runnumber) '_' mode '.txt'], tmp_matrix, 'delimiter', '\t');
patch_removed_intervals = tmp_matrix_r;
dlmwrite(['patch_removed_intervals_' num2str(runnumber) '_' mode '.txt'], tmp_matrix_r, 'delimiter', '\t');
results_mat = results.co;
dlmwrite(['results_' num2str(runnumber) '_' mode '.txt'], results.co, 'delimiter', '\t');
results_removed_mat = results_removed;
dlmwrite(['results_removed_' num2str(runnumber) '_' mode '.txt'], results_removed, 'delimiter', '\t');


save patch_intervals
