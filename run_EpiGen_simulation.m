N = 1000;
time=15000;
%time = 10;
%t_split = [10, 15, 20, 25, 30, 40, 50, 100, 200];
t_split = mytsplit;
%t_split = 40;
radius=1;
%params = [100 0;10 90;50 50];
%Nugen = [100 0;10 90;50 50]
nugen = 10;
nuepi = 90;
epiback=0.0205; %ran cost function
%epiback = epiopt;
mlimitgen=2;
%mlimitepi=2;
%mlimitepi=epi_limit;
mlimitepi = 0.3;
%mlimitepi=[2,1,0.5,0.3,0.1];
%oligo= [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
oligo = 1;
epigenetics = 'TRUE';
patch = 'TRUE';
%mode = mymode; 
mode = 'negative'; %random or negative
%runnumber = myrun;
runnumber = 1;
tfact = 0.1;
%tfact = mytfact;

patch_intervals = cell(length(t_split), 1); %Create cell array of matrices to store interval runs
patch_removed_intervals = cell(length(t_split), 1); %Create cell array of matrices to store interval runs of removed individuals
results_mat = cell(length(t_split), 1); %Create cell array of matrices to store results matrix runs
results_removed_mat = cell(length(t_split), 1); %Create cell array of matrices to store results matrix runs of removed individuals
epi_mat = cell(length(t_split), 1); %Create cell array of matrices to store epigenetic distribution run
gen_mat = cell(length(t_split), 1); %Create cell array of matrices to store genetic distribution run
epi_removed_mat = cell(length(t_split), 1); %Create cell array of matrices to store epigenetic distribution run
gen_removed_mat = cell(length(t_split), 1); %Create cell array of matrices to store genetic distribution run


    
Driver=ones(time,1);
for q=t_split:t_split*2:time
    Driver(q+1:q+t_split-1)=0; %generate time intervals to switch back and forth from
end

[results, population_disgen, population_disepi, results_removed, population_removed_disgen, population_removed_disepi, population_numgen]= EpiGen_simulation(N,time,radius,nugen,nuepi, ...
    epiback, mlimitgen, mlimitepi, epigenetics, patch, oligo, Driver, mode, runnumber, tfact);

%PULL NUMBER OF GENETIC MUTATIONS FOR KEPT POPULATIONS
dlmwrite(['popnumgen_' num2str(runnumber) '_' mode '.txt'],population_numgen, 'delimiter', '\t');

%PULL DISTRIBUTION OF GENETIC EFFECTS FOR 'KEPT' POPULATIONS
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

    
    
    
    %PLOT DISTANCE OF GENETIC/TOTAL AND EPIGENETIC/TOTAL MUTATIONS%%%
    %walk_epi_dist = sum(abs(results.co(:,6)));
    
    %for g=1:length(results.co(:,6))
        
        %calculate epigenetic potential
    %    current_gen(g) = (walk_epi_dist - sum(abs(results.co(1:g,6))))/walk_epi_dist; 

    %end
    %plot(results.co(:,1), current_gen);
    %gen_18 = (walk_epi_dist - sum(abs(results.co(1:18,6))))/walk_epi_dist;
    %gen_50 = (walk_epi_dist - sum(abs(results.co(1:50,6))))/walk_epi_dist;
    %gen_100 = (walk_epi_dist - sum(abs(results.co(1:100,6))))/walk_epi_dist;
    %gen_432 = (walk_epi_dist - sum(abs(results.co(1:432,6))))/walk_epi_dist;
    %gen_756 = (walk_epi_dist - sum(abs(results.co(1:756,6))))/walk_epi_dist;
    %m_data = [18 gen_18; 432 gen_432; 756 gen_756];
    %m_data = [18 gen_18; 50 gen_50; 100 gen_100];
    %plot(m_data(:,1), m_data(:,2));
    %PLOT GENETIC/TOTAL AND EPIGENETIC/TOTAL MUTATIONS%%%
    
    
    %%%%%shaded error tutorial%%%%
    %https://www.mathworks.com/matlabcentral/answers/39540-continuous-error-bars
    %x = linspace(0,1,20)';
    %y = sin(x);
    %dy = .1*(1+rand(size(y))).*y;  % made-up error values
    %fill([x;flipud(x)],[y-dy;flipud(y+dy)],[.9 .9 .9],'linestyle','none');
    %line(x,y)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %{
    if oligo(i) == 1

        cond = 'n = 5';

        %Save n = 5 results to c matrix 
        n5 = results.co;
        
        dlmwrite('n5.txt',n5,'delimiter','\t','precision',3);
        
    elseif oligo(i) == 2

        cond = 'n = 10';

        %Save n = 10 results to cp matrix
        n10 = results.co;
        
        %dlmwrite('n10.txt',n10,'delimiter','\t','precision',3);
        
    elseif oligo(i) == 3

        cond = 'n = 15';

        %Save n = 15 results to ci matrix
        n15 = results.co;
        
        %dlmwrite('n15.txt',n15,'delimiter','\t','precision',3);
        
    elseif oligo(i) == 4

        cond = 'n = 20';

        %Save n = 20 results to cpi matrix
        n20 = results.co;
        
        %dlmwrite('n20.txt',n20,'delimiter','\t','precision',3);
       
    elseif oligo(i) == 5

        cond = 'n = 25';

        %Save n = 25 results to cpi matrix
        n25 = results.co;
        
        %dlmwrite('n25.txt',n25,'delimiter','\t','precision',3);
        
    elseif oligo(i) == 6

        cond = 'n = 30';

        %Save n = 30 results to cpi matrix
        n30 = results.co;
        
        %dlmwrite('n30.txt',n30,'delimiter','\t','precision',3);
        
    elseif oligo(i) == 7

        cond = 'n = 35';

        %Save n = 35 results to cpi matrix
        n35 = results.co;
        
        %dlmwrite('n35.txt',n35,'delimiter','\t','precision',3);
        
    elseif oligo(i) == 8

        cond = 'n = 40';

        %Save n = 40 results to cpi matrix
        n40 = results.co;
        
        %dlmwrite('n40.txt',n40,'delimiter','\t','precision',3);
        
    elseif oligo(i) == 9

        cond = 'n = 5 + 20';

        %Save n = 5 + 20 results to cpi matrix
        n_5_20 = results.co;
        
        %dlmwrite('n_5_20.txt',n_5_20,'delimiter','\t','precision',3);
        
    elseif oligo(i) == 10

        cond = 'n = 5 + 40';

        %Save n = 5 + 40 results to cpi matrix
        n_5_40 = results.co;
        
        %dlmwrite('n_5_40.txt',n_5_40,'delimiter','\t','precision',3);
        
    end


    legendinfo{i}= ... 
        ['Oligotrophic conditions = ' cond];
        %['EpiGenetic mutational limit = ' num2str(mlimitepi(i))];
            
    
    
    %figure
    
    
%}
    


%save results_mat 
%save results_removed_mat
save patch_intervals
%save patch_removed_intervals
%save gen_mat
%save gen_removed_mat
%save epi_mat
%save epi_removed_mat


%legend(legendinfo, 'Location', 'northwest');
%title('Fitness increases with different Epigenetic mutational effects Genetic Supply = 10 -- Epigenetic Supply = 90 ');
%xlabel('Generations')
%ylabel('Relative Number of mutations to total pool')
%hold off



%saveas(figure1,'test.pdf')

%quit
