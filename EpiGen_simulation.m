function [results, population_disgen, population_disepi, results_removed, population_removed_disgen, population_removed_disepi]=EpiGen_simulation(N, time, radius, nugen, nuepi, epiback, mlimitgen, mlimitepi, epigenetics, patch, oligo, Driver, mode, runnumber, t_split)



%{
N=1000; %the number of individuals in the population
time=1000; %the number of generations in the population

radius=1; %distance from the hypersphere origin
nugen=10; %number of genetic (HT) mutations for each generation
nuepi=90; %number of epigenetic (LT) mutations for each generation
epiback=0.0205; %epigenetic reversion rate
mlimitgen=2; %maximum effect a genetic mutation can have
mlimitepi=0.3; %maximum effect an epigenetic mutation can have
oligo = 1; % switch to determine number of dimensions of the hypersphere. 
epigenetics = 'TRUE'; %switch to toggle effects of epigenetics on or off
patch = 'TRUE'; %switch to determine if population switches between one
environment and another. If TRUE, population goes in and out of selection
environment at a frequency of t_split generations. If FALSE, only stays in
the selection environment.
t_split = 10; %if patch='TRUE', t_split is the number of generations the
population stays in one environment or the other
mode = 'negative'; When patch = 'TRUE', this will run the model in two modes that determines the selection rules for environment 2: random or negative. 'negative' turns on randomly sampling individuals weighted by the reciprocal number of genetic mutations when in environment 2 (stabilizing selection). 'random' turns on randomly sampling indviduals in environment 2 with no weighting.
runnumber = 1; %replicate number of model run


Driver is a vector of 0's and 1's and has length = time in generations. The loop below will step through
every element of Driver. When patch = 'TRUE' and Driver == 1, then select
= 1, and the population is randomly sampled to the next generation weighted
by fitness, i.e. is in environment 1. When patch = 'TRUE' and
Driver == 0, then the population is randomly sampled to the next
generation weighted by the reciprocal of the number of genetic mutations,
i.e. is in environment 2.

Driver=ones(time,1);
for q=t_split(1):t_split(1)*2:time
    Driver(q:q+t_split(1)-1)=0; %generate time intervals to switch back and forth from
end
    
%}  

%this switch was originally encoded to switch between the number of n
%dimensions of the hypersphere. For all model runs, we kept n = 5. So, below when
%rep = 1, it sets the number of dimensions to n = 5
rep = 1; 

%Population information
population.co=zeros(N,25); %%create populations
population.co(:,1) = 1:N; %Index for individuals
population.co(:,2)=1; %initial distance from the optimum is radius
population.co(:,3)=exp((-population.co(:,2).^2)/2); %Fitness - Using a Gaussian fitness function
%Column 4 is for distance traveled by genetic mutation
%Column 5 is for number of genetic mutations
%Column 6 is last time point of genetic mutations
%Column 7 is for distance traveled by epigenetic mutations
%Column 8 is for number of epigenetic mutations
%Column 9 is last time point of epigenetic mutations


%Information of individuals in the population that were removed from the main reproducing population
population_removed=zeros(N,25); %%create populations
%Column 1 %Index for individuals
%Column 2 %initial distance from the optimum is radius
%Column 3 %fitness function
%Column 4 is for distance traveled by genetic mutation
%Column 5 is for number of genetic mutations
%Column 6 is last time point of genetic mutations
%Column 7 is for distance traveled by epigenetic mutations
%Column 8 is for number of epigenetic mutations
%Column 9 is last time point of epigenetic mutations


%%%matrix to populate in loop
population_disgen=cell(N, time);  %distribution of genetic mutations of kept population (i.e. sampled to go onto next generation)
population_disepi=cell(N, time);  %distribution of epigenetic mutations of kept population
population_removed_disgen = cell(N, time); %distribution of genetic mutations of removed population (i.e. not sampled to go onto next generation)
population_removed_disepi = cell(N, time); %distribution of epigenetic mutations of removed population

%%
%Initialize geometric model
%dimensions of spherical phenotypic space 
%Z-values have been solved manually, dimensions solved for are 5, 10, 15,
%20, 25, 30, 35, 40
%If some other value of n is needed, need to solve the integral for Z and add the value to Z.mat
zmat=[5 0.75; 10 1.1641; 15 1.4663; 20 1.7162; 25 1.9342; 30 2.1299; 35 2.3092; 40 2.4755];

%Store Z-value for 5 dimension. We use 5 dimensions for all model runs
Z5=zmat(1,2);  

%the angle between the mutational vector and the vector running from the current phenotype to the origin
%Angles can take on values between -pi/2 to pi/2
phi_values=linspace(-pi/2,pi/2,50000); 

%Calculate the n = 5 probability density function for the random angle
prob_density5=Z5*cos(phi_values).^(zmat(1,1) - 2);

%Initialize Results of kept population
%Column 1 generations
%Column 2 mean population fitness
%Column 3 mean number of genetic mutations
%Column 4 mean number of epigenetic mutations
%Column 5 mean distance traveled by genetic mutations
%Column 6 mean distance traveled by epigenetic mutations
results.co=zeros(time,28);
results.co(:,1)=1:time;
%colname={'time','meanfitness','nogenmut','noepimut','distgen','distepi','net_fi_carbon_gen','net_fi_p_gen','net_fi_i_gen','net_fi_carbon_epi','net_fi_p_epi','net_fi_i_epi'};
%colindex={1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};

%Initialize Results of removed population (i.e. individuals not sampled to
%next generation
results_removed=zeros(time,28);
results_removed(:,1)=1:time;


%select = 1; %Switch to randomly resample population either not weighted (0) or weighted (1) by fitness 
check_count = 100; %print results every check_count generations 
%%%%%%%%%%%%%%%%%%%%%
%%Looping over time
%%%%%%%%%%%%%%%%%%%%%
for i=1:length(Driver)
    

    %Check if patch is TRUE or FALSE. If patch_tf == 1, run model with only
    %1 type of selection (i.e. 1 environment) where sampling is weighted by
    %fitness. If patch_ft == 0, run the model for two types of selection (i.e. two
    %environments). When Driver = 1, sample weighted by fitness. When
    %Driver = 0, sample weighted by the reciprocal number of genetic
    %mutations
    patch_tf = strcmp(patch,'FALSE'); 

    if patch_tf == 1
        
        select = 1; %only sample with weighting as fitness

    elseif (patch_tf == 0) && (Driver(i) == 1)
        
        
        select = 1; %resample population with weight by fitness
       
        
        
    elseif (patch_tf == 0) && (Driver(i) == 0)
    
        
        select = 0; %resample population with other weighting parameters (e.g., by the reciprocal number of genetic mutations)
   
        
    end

    %%%%%%%%%%
    %Make genetic mutations, mutational supply is set to nugen/1 unit time
    %Distances from optimum are calculated individually
    %%%%%%%%%%
    
    %pause the clock
    c=clock;
    pause(c(6)*runnumber/100)
    rng('shuffle')
    
    %Vector of genetic mutations with each value being the magnitude of
    %their effect
    mut_effect=0+(mlimitgen*radius).*rand(1,nugen); 
    
    %pause the clock
    c=clock;
    pause(c(6)*runnumber/100)
    rng('shuffle')
    
    %Randomly sample 10 angles
    phi5=randsample(phi_values, nugen, true, prob_density5); 
    

    %%%Assign mutations to individuals%%%
    
    index_mut=1:N;
    
    %pause the clock
    c=clock;
    pause(c(6)*runnumber/100)
    rng('shuffle')
    
    %Randomly sample 10 individuals to genetically mutate without replacement (multiple hits are not allowed)
    index_mutated=randsample(index_mut,nugen,false); 
    
    %retrieve mutants from population 
    mutated_orgs=population.co(index_mutated,:);
   

    %%%calculate distance traveled after genetic mutation
    dist_travelled5 = [mutated_orgs(:,2)' - sqrt(mutated_orgs(:,2)'.^2 + ...
    mut_effect.^2 + 2.*mutated_orgs(:,2)'.*mut_effect.*sin(phi5))];
       

    if rep == 1 %always equals 1 in this model version
        

        %Calculate distance from the optimum phenotype
        dist5 = sqrt(mutated_orgs(:,2)'.^2 + mut_effect.^2 ...
        + 2.*sin(phi5).*mutated_orgs(:,2)'.*mut_effect);
        
        
        %%%%%%%%%%%%%%%%%%%%
        %Distance from the optimum phenotype based on oligo conditions;
        %%%%%%%%%%%%%%%%%%%%
        
        if oligo == 1 %n = 5 dimensions of the hypersphere
            
            %store distance from optimum for mutated individuals
            mutated_orgs(:,2) = dist5'; 
            
            %rename distance_travelled5 vector to dist_travelled_gen
            dist_travelled_gen = dist_travelled5;
            
        elseif oligo == 2 %n = 10
            
            %n = 10
            mutated_orgs(:,2) = dist10'; 
            
            %distance traveled after mutation
            dist_travelled_gen = dist_travelled10;
        
        elseif oligo == 3 %n = 15
            
            %n = 15
            mutated_orgs(:,2) = dist15'; 
            
            %distance traveled after mutation
            dist_travelled_gen = dist_travelled15;
            
        elseif oligo == 4 %n = 20
            
            %n = 20
            mutated_orgs(:,2) = dist20'; 
            
            %distance traveled after mutation
            dist_travelled_gen = dist_travelled20;
            
        elseif oligo == 5 %n = 25
            
            %n = 25
            mutated_orgs(:,2) = dist25'; 
            
            %distance traveled after mutation
            dist_travelled_gen = dist_travelled25;   
            
        elseif oligo == 6 %n = 30
            
            %n = 30
            mutated_orgs(:,2) = dist30'; 
            
            %distance traveled after mutation
            dist_travelled_gen = dist_travelled30; 
            
        elseif oligo == 7 %n = 35
            
            %n = 35
            mutated_orgs(:,2) = dist35'; 
            
            %distance traveled after mutation
            dist_travelled_gen = dist_travelled35;
            
        elseif oligo == 8 %n = 40
            
            %n = 40
            mutated_orgs(:,2) = dist40'; 
            
            %distance traveled after mutation
            dist_travelled_gen = dist_travelled40;
            
        elseif oligo == 9 %n = 5 + 20
  
            %n = 5 + 20 
            cat_5_20 = horzcat(dist5', dist20');
            dist_5_20 = mean(cat_5_20, 2)';
         
            mutated_orgs(:,2) = dist_5_20';
            
            %distance traveled after mutation
            dist_travelled_gen = dist_travelled_5_and_20;
            
        elseif oligo == 10 %n = 5 + 40
            
            %n = 5 + 20 
            cat_5_40 = horzcat(dist5', dist40');
            dist_5_40 = mean(cat_5_40, 2)';
         
            mutated_orgs(:,2) = dist_5_40';
            
            %distance traveled after mutation
            dist_travelled_gen = dist_travelled_5_and_40;
            
        end
        
        
    end
    
    
    %Save mean distance traveled toward optimum by genetic mutation
    results.co(i,25) = mean(mutated_orgs(:,2));
    
    %Save std of distance traveled toward optimum by genetic mutation
    results.co(i,26) = std(mutated_orgs(:,2));
    
    %Calculates fitness for the mutants after genetic mutation
    mutated_orgs(:,3) = exp((-mutated_orgs(:,2).^2)/2);
    
    %%%Cumulative distance travelled by genetic mutations
    mutated_orgs(:,4) = mutated_orgs(:,4) + dist_travelled_gen';
        
    %Increase the number of genetic mutations by one
    mutated_orgs(:,5) = mutated_orgs(:,5) + 1;
    
    %save distance travelled after genetic mutations
    mutated_orgs(:,6) = dist_travelled_gen';
    
    %NOT USED IN THIS VERSION OF MODEL CALCULATIONS - IGNORE
    %Using law of sines, calculate alternative angle. 
    fi5 = asin((sin(phi5 + (pi/2)).*(mut_effect))./dist5); 
    
    %Insert fi5
    mutated_orgs(:,10) = mutated_orgs(:,10) + fi5';

    %Storing the mutated individuals
    population.co(index_mutated,:)=mutated_orgs;

    %Storing distribution of genetic mutations
    for g=1:length(dist_travelled_gen')
        
        population_disgen{index_mutated(g), i} = strcat(num2str(dist_travelled_gen(g)));
    
    end
    
   
    %If the model run includes epigenetic effects, run the following
    tf = strcmp(epigenetics,'TRUE');

    if tf == 1
        
        %Make epigenetic mutations, mutational supply is set to Nu.epi / 1 unit of time
        
        %pause the clock
        c=clock;
        pause(c(6)*runnumber/100)
        rng('shuffle')
        
        %Mutational effects of epigenetics
        mut_effect_epi=0+(mlimitepi*radius).*rand(1, nuepi)';
        
        %Random angle for epigenetic mutations%%%%%%%%%%%%%%%%%%%%%
        
        %pause the clock
        c=clock;
        pause(c(6)*runnumber/100)
        rng('shuffle')
        
        %Randomly draw random angles for epigenetic mutations
        phi_epi5=randsample(phi_values, nuepi, true, prob_density5);


        %pause the clock
        c=clock;
        pause(c(6)*runnumber/100)
        rng('shuffle')
        
        %retrieve epimutants from population w/out replacement
        index_mutated_epi = randsample(index_mut, nuepi, false);
        mutated_orgs_epi = population.co(index_mutated_epi, :);
        
        
        %%%distance traveled after epimutations
        dist_travelled_epi5 = [mutated_orgs_epi(:,2)' - sqrt(mutated_orgs_epi(:,2)'.^2 + ...
        mut_effect_epi.^2' + 2.*mutated_orgs_epi(:,2)'.*mut_effect_epi'.*sin(phi_epi5))];

        if rep == 1 %always set to 1 for this model version

            %Distance from the optimum phenotype epimutations
            dist_epi5 = sqrt(mutated_orgs_epi(:,2)'.^2 + mut_effect_epi.^2' ...
            + 2.*sin(phi_epi5).*mutated_orgs_epi(:,2)'.*mut_effect_epi');

        end


         %Distance from the optimum phenotype based on oligo conditions;
         
        if oligo == 1 %5 dimensions of hypersphere
            
            %Store distance from the optimum phenotype 
            mutated_orgs_epi(:,2) = dist_epi5'; 

            %rename dist_travelled_epi5 to dist_travelled_e
            dist_travelled_e = dist_travelled_epi5;

        elseif oligo == 2 %n = 10
            
          
            %Distance from the optimum phenotype;
            mutated_orgs_epi(:,2) = dist_epi10';

            %distance traveled after epimutation
            dist_travelled_e = dist_travelled_epi10;
            
            
        elseif oligo == 3 %n = 15

            
            %Distance from the optimum phenotype 
            mutated_orgs_epi(:,2) = dist_epi15';

            %%distance traveled after epimutation
            dist_travelled_e = dist_travelled_epi15;
            
        elseif oligo == 4 %n = 20

            
            %Distance from the optimum phenotype 
            mutated_orgs_epi(:,2) = dist_epi20';

            %%distance traveled after epimutation
            dist_travelled_e = dist_travelled_epi20;
        
        elseif oligo == 5 %n = 25

            
            %Distance from the optimum phenotype 
            mutated_orgs_epi(:,2) = dist_epi25';

            %%distance traveled after epimutation
            dist_travelled_e = dist_travelled_epi25;
            
        elseif oligo == 6 %n = 30

            
            %Distance from the optimum phenotype 
            mutated_orgs_epi(:,2) = dist_epi30';

            %%distance traveled after epimutation
            dist_travelled_e = dist_travelled_epi30;
        
        elseif oligo == 7 %n = 35

            
            %Distance from the optimum phenotype 
            mutated_orgs_epi(:,2) = dist_epi35';

            %%distance traveled after epimutation
            dist_travelled_e = dist_travelled_epi35;
            
        elseif oligo == 8 %n = 40

            
            %Distance from the optimum phenotype 
            mutated_orgs_epi(:,2) = dist_epi40';

            %%distance traveled after epimutation
            dist_travelled_e = dist_travelled_epi40;

        elseif oligo == 9 %n = 5 + 20
            
         
            %n = 5 + 20
            cat_epi_5_20 = horzcat(dist_epi5', dist_epi20');
            dist_epi_5_20 = mean(cat_epi_5_20, 2)';
            
            %Distance from the optimum phenotype;
            mutated_orgs_epi(:,2) = dist_epi_5_20';

            %%Distance travelled by all - current distance after epigenetic mutations
            dist_travelled_epi_5_20 = horzcat(dist_travelled_epi5', dist_travelled_epi20');
            dist_travelled_epi_5_and_20 = mean(dist_travelled_epi_5_20, 2)';
            dist_travelled_e = dist_travelled_epi_5_and_20;
            
        elseif oligo == 10 %n = 5 + 40
            
         
            %n = 5 + 40
            cat_epi_5_40 = horzcat(dist_epi5', dist_epi40');
            dist_epi_5_40 = mean(cat_epi_5_40, 2)';
            
            %Distance from the optimum phenotype;
            mutated_orgs_epi(:,2) = dist_epi_5_40';

            %%Distance travelled by all - current distance after epigenetic mutations
            dist_travelled_epi_5_40 = horzcat(dist_travelled_epi5', dist_travelled_epi40');
            dist_travelled_epi_5_and_40 = mean(dist_travelled_epi_5_40, 2)';
            dist_travelled_e = dist_travelled_epi_5_and_40;

        end
        
        %Save mean distance travelled from optimum by epigenetic mutation
        results.co(i,27) = mean(mutated_orgs_epi(:,2));
    
        %Save std of distance travelled from optimum by epigenetic mutation
        results.co(i,28) = std(mutated_orgs_epi(:,2));
        
        %Calculate fitness for mutants
        mutated_orgs_epi(:,3) = exp((-mutated_orgs_epi(:,2).^2)/2);
        
        %Cumulative distance travelled by epigenetic mutations
        mutated_orgs_epi(:,7) = mutated_orgs_epi(:,7) + dist_travelled_e';
        
        %Increase the number of epigenetic mutations by one
        mutated_orgs_epi(:,8) = mutated_orgs_epi(:,8) + 1;
        
        %Store distance travelled
        mutated_orgs_epi(:,9) = dist_travelled_e';
        
        %CALCULATE ALTERNATIVE ANGLE - IGNORE THIS METRIC - NOT USED IN
        %THIS MODEL VERSION
        %Using law of sines
        fi_epi5 = asin((sin(phi_epi5 + (pi/2)).*(mut_effect_epi'))./dist_epi5);
        
        %Insert fi5 for epi
        mutated_orgs_epi(:,18) = mutated_orgs_epi(:,18) + fi_epi5';

        %Storing mutated individuals
        population.co(index_mutated_epi,:) = mutated_orgs_epi;
        
        %Distribution of epigenetic mutations
        %population_disepi(index_mutated_epi,i) = dist_travelled_e';
        
        for e=1:length(dist_travelled_e')
            
            population_disepi{index_mutated_epi(e), i} = dist_travelled_e(e);
        
        end
    
    end
    
    
    %As of now, multiple hits of either genetic and epigenetic mutations are not possible. 
    %However, a single individual can have two hits if one is genetic and one epigenetic mutation. 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Population regulation and reproduction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Soft selection is acting -> population size stays the same
    %Sampling individuals to the next time point with replacement
    
    ofsp_index = 1:N;
    
    if select == 0 %if true, population is in environment 2
        
        mode_type = strcmp(mode,'random'); %If true, sample without weighting
        %Sampling the next generation, no weighting
        if mode_type == 1
        
        %pause the clock
        c=clock;
        pause(c(6)*runnumber/100)
        rng('shuffle')
        
        sampled_index = randsample(ofsp_index, N, true);
        else
        %Sampling the next generation, negative weighting (i.e. the
        %reciprocal number of genetic mutations)
        reciprocal_weights = 1./population.co(:,5);
        
        %pause the clock
        c=clock;
        pause(c(6)*runnumber/100)
        rng('shuffle')
        
        sampled_index = randsample(ofsp_index, N, true, reciprocal_weights);
        
        end
        
    elseif select == 1 %if true, population is in environment 1
        
        %Sampling the next generation, weighting with fitness
        c=clock;
        pause(c(6)*runnumber/100)
        rng('shuffle')
        sampled_index = randsample(ofsp_index, N, true, population.co(:,3));
        
    end
    
    %Population after reproduction and regulation
    population.co = population.co(sampled_index,:);
    population_disgen = population_disgen(sampled_index,:);
    population_disepi = population_disepi(sampled_index,:);
    
    %Removed individuals from population
    pop_removed_index = setdiff(ofsp_index, sampled_index);
    population_removed = population.co(pop_removed_index,:);
    
    population_removed_disgen = population_disgen(pop_removed_index,:);
    population_removed_disepi = population_disepi(pop_removed_index,:);
    
    
    %%%%%%Store results removed population
   
    %Mean fitness of the population
    results_removed(i,2)= mean(population_removed(:,3)); 
 
    %Store the mean number of genetic mutations present in an individual
    results_removed(i,3)= mean(population_removed(:,5)); 
    
    %Store the standard deviation of mean genetic mutations present in an individual
    results_removed(i,23)= std(population_removed(:,5));
    
    %Store the mean number of epigenetic mutations present in an individual
    results_removed(i,4)= mean(population_removed(:,8));
    
    %Store the standard deviation of mean epigenetic mutations present in an individual
    results_removed(i,24)= std(population_removed(:,8));
    
    %Store the cumulative distance travelled towards the optimum by genetic mutations
    results_removed(i,5)= mean(population_removed(:,4));
    
    %Store the cumulative distance travelled towards the optimum by epigenetic mutations
    results_removed(i,6)= mean(population_removed(:,7));
    
    
    %%%%%%Store results sampled (non-removed) population

    %Mean fitness of the population
    results.co(i,2)= mean(population.co(:,3)); 

    %Store the mean number of genetic mutations present in an individual
    results.co(i,3)= mean(population.co(:,5)); 
    
    %Store the standard deviation of mean genetic mutations present in an individual
    results.co(i,23)= std(population.co(:,5));
    
    %Store the mean number of epigenetic mutations present in an individual
    results.co(i,4)= mean(population.co(:,8));
    
    %Store the standard deviation of mean epigenetic mutations present in an individual
    results.co(i,24)= std(population.co(:,8));
    
    %Store the cumulative distance travelled towards the optimum by genetic mutations
    results.co(i,5)= mean(population.co(:,4));
    
    %Store the cumulative distance travelled towards the optimum by epigenetic mutations
    results.co(i,6)= mean(population.co(:,7));
    
    %Store mean genetic mutation fi5
    results.co(i,7) = mean(population.co(:,10));

    %Store mean epigenetic mutation fi5
    results.co(i,15) = mean(population.co(:,18));
    
    
    if (i==check_count) %print results every 100 generations to check
        
        filename = ['results_check_' num2str(t_split) '_' mode '.txt'];
        %filename = strjoin(filename, {''});
        dlmwrite(filename, results.co, '-append', 'delimiter', '\t');
        check_count = check_count + 100;
    end
          
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Loss of epigenetic effects
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Next epigenetic changes are lost at a given rate 
    %(This needs to be slower than the forward mutation rate of epigenetic effects)
    %Epigenetic changes revert independently (always at a probability of: epi.backmutation in each generation)
    
    %Number of epigenetic mutations present
    no_epi_mut = sum(population.co(:,8));
    
    if (no_epi_mut>0)
        
        %Check which individuals have epigenetic mutation
        %[epi_p_index,col,v] = find(population.orgs(:,8));
        epi_p_index = find(population.co(:,8) > 0);
        epi_present = population.co(epi_p_index,:);
        epi_present_gen = population_disgen(epi_p_index,:);
        epi_present_epi = population_disepi(epi_p_index,:);
        
        %Generate backmutations probabilistically
        
        %epi_backmuts = binornd(1, epiback, sum(epi_present(:,8)),1); 
        if sum(epi_present(:,8)) > 0
            epi_backmuts = binornd(1, epiback, sum(epi_present(:,8)),1); 
        else
            continue
        end
        
        %Check for each backmutation whether reversion occurs or not
        %%%Mapping from individual mutations to individuals
        
        %Construct a helping vector with as many factors as there are inds with the equal 
        %length of mutation vector
        
        help_vec = repelem(1:size(epi_present(:,8)),epi_present(:,8));
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Sum by each factor
        %%%%%%%%%%%%%%%%%%%%Same function as "tapply" in R
        
        list_of_reversions = splitapply(@sum,epi_backmuts,help_vec');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Logical vector if backmutation occurs in a given individual (TRUE) or not
        
        epi_back_mut_index = find(list_of_reversions > 0);
        
        %Take only those individuals that had backmutations
        epi_back_mut_ind = epi_present(epi_back_mut_index,:);
        epi_back_mut_ind_gen = epi_present_gen(epi_back_mut_index,:);
        epi_back_mut_ind_epi = epi_present_epi(epi_back_mut_index,:);
        
        %How many backmutations occurred in each individual that had backmutations
        number_of_reversions = list_of_reversions(epi_back_mut_index);
        
        %if any backmutations occurred
        if any(list_of_reversions) > 0
            
            %pause clock
            c=clock;
            pause(c(6)*runnumber/100)
            rng('shuffle')
            
            
            for b=1:length(number_of_reversions)
                
                
                lost_mutation=zeros(1,number_of_reversions(b));
                
                current_org = cell2mat(epi_back_mut_ind_epi(b,:));
                current_org = current_org(current_org ~= 0); 
                %randomly sample epi effects from the distribution
                if length(current_org) >= number_of_reversions(b)
                    lost_mutation(1,:)=randsample(current_org, ...
                    number_of_reversions(b),false);
                end
               
                %Now individuals need to lose mutation
                if ~isempty(lost_mutation)
                    
                    %loop through all epimutations
                    for count=1:length(lost_mutation) 
                        %ensure current org contains epimutations to be
                        %lost
                        logic = any(current_org == lost_mutation(1,count)); 
                        
                        if logic == 1
                            %gather epimutation values to be lost and
                            %convert to strings in cell array r
                            r = cellfun(@num2str, epi_back_mut_ind_epi(b,:), 'UniformOutput', false);
                            
                            %find indicies of lost mutations in cell array
                            %r
                            idx2 = find(ismember(r, num2str(lost_mutation(1,count))));
                            
                            %if only 1 mutation is found,
                            %set to zero
                            if length(idx2) == 1 

                                epi_back_mut_ind_epi{b,idx2} = 0;
                            
                            %if more than 1 epimutation, loop through the
                            %number of epimutations to lose
                            elseif length(idx2) > 1

                                for l=1:length(idx2)

                                    epi_back_mut_ind_epi{b,idx2(l)} = 0;

                                end

                            end

                        end
                    end
                end

                %epi_back_mut_ind_epi(b,:)
                %lost_mutation(:,1)
              
               
               %Subtract the effect of the mutations from distance travelled by epigenetics
               epi_back_mut_ind(b,7)= epi_back_mut_ind(b,7)-sum(lost_mutation(1,:));
               
               %Add the effect of the mutation from distance to optimum
               epi_back_mut_ind(b,2)= epi_back_mut_ind(b,2)+sum(lost_mutation(1,:));
               
               %Subtract the number of epigenetic mutations
               epi_back_mut_ind(b,8)=epi_back_mut_ind(b,8)-number_of_reversions(b);
               
               %Recalculate fitness for backmutated individuals
               epi_back_mut_ind(b,3) = exp((-epi_back_mut_ind(b,2).^2)/2);
               

            end
                
            
          
            %Replace backmutated individuals in the population.co dataframe
            epi_present(epi_back_mut_index,:)=epi_back_mut_ind;
            population.co(epi_p_index,:)=epi_present;
            
            %Replace the epigenetic mutational distributions of the 
            %backmutated individuals in the population_disepi dataframe
            epi_present_epi(epi_back_mut_index,:)=epi_back_mut_ind_epi;
            population_disepi(epi_p_index,:)=epi_present_epi;
            
            
                
            
        end
        
        
        
    end
    
end

