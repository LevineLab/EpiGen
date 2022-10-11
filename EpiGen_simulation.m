function [results, population_disgen, population_disepi, results_removed, population_removed_disgen, population_removed_disepi, population_numgen]=EpiGen_simulation(N, time, radius, nugen, nuepi, epiback, mlimitgen, mlimitepi, epigenetics, patch, oligo, Driver, mode, runnumber, tfact)



clear
N=1000;
time=200;

radius=1;
nugen=10;
nuepi=90;
epiback=0.0205;
mlimitgen=2;
mlimitepi=0.3;
oligo = 1;
epigenetics = 'TRUE';
patch = 'TRUE';
t_split = 10;

Driver=ones(time,1);
for q=t_split(1):t_split(1)*2:time
    Driver(q:q+t_split(1)-1)=0; %generate time intervals to switch back and forth from
end
        
replicates=4; %always leave on for now

%FITNESS OF POPULATION 
%DIMENSIONS
population.co=zeros(N,25); %%create populations
population.co(:,1) = 1:N; %Index for individuals
population.co(:,2)=1; %initial distance from the optimum is radius
population.co(:,3)=exp((-population.co(:,2).^2)/2); %Using a Gaussian fitness function
%Column 4 is for distance traveled by genetic mutation
%Column 5 is for number of genetic mutations
%Column 6 is last time point of genetic mutations
%Column 7 is for distance traveled by epigenetic mutations
%Column 8 is for number of epigenetic mutations
%Column 9 is last time point of epigenetic mutations
%Column 10 is the net fi5 for genetic mutations
%Column 11 is the net fi10 for genetic mutations
%Column 12 is the net fi15 for genetic mutations
%Column 13 is the net fi20 for genetic mutations
%Column 14 is the net fi25 for genetic mutations
%Column 15 is the net fi30 for genetic mutations
%Column 16 is the net fi35 for genetic mutations
%Column 17 is the net fi40 for genetic mutations

%Column 18 is the net fi5 for epigenetic mutations
%Column 19 is the net fi10 for epigenetic mutations
%Column 20 is the net fi15 for epigenetic mutations
%Column 21 is the net fi20 for epigenetic mutations
%Column 22 is the net fi25 for epigenetic mutations
%Column 23 is the net fi30 for epigenetic mutations
%Column 24 is the net fi35 for epigenetic mutations
%Column 25 is the net fi40 for epigenetic mutations

%FITNESS OF REMOVED POPULATION 
%DIMENSIONS
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
population_disgen=cell(N, time);  %distribution of distance of genetic mutations of kept population
population_disepi=cell(N, time);  %distribution of epigenetic mutations of kept population
population_numgen=zeros(N, time);  %distribution of number genetic mutations of kept population
%population_numepi=zeros(N, time);  %distribution of number epigenetic mutations of kept population
population_removed_disgen = cell(N, time); %distribution of genetic mutations of removed population
population_removed_disepi = cell(N, time); %distribution of epigenetic mutations of removed population

%%
%Initialize geometric model
%dimensions of spherical phenotypic space 
%Z-values have been solved manually, dimensions supported are 3, 25, 30, 50
%If some other value of n is needed, need to solve the integral for Z and add the value to Z.mat
zmat=[5 0.75; 10 1.1641; 15 1.4663; 20 1.7162; 25 1.9342; 30 2.1299; 35 2.3092; 40 2.4755];

%Store Z-value for 5 dimension
Z5=zmat(1,2);  

%Store Z-value for 10 dimension
%Z10=zmat(2,2);  

%Store Z-value for 15 dimension
%Z15=zmat(3,2); 

%Store Z-value for 20 dimension
%Z20=zmat(4,2);  

%Store Z-value for 25 dimension
%Z25=zmat(5,2);  

%Store Z-value for 30 dimension
%Z30=zmat(6,2);

%Store Z-value for 35 dimension
%Z35=zmat(7,2);  

%Store Z-value for 40 dimension
%Z40=zmat(8,2);

%The random angle can have values from -pi/2 to pi/2
phi_values=linspace(-pi/2,pi/2,50000); 

%Calculate the n = 5 probability density function for the random angle
prob_density5=Z5*cos(phi_values).^(zmat(1,1) - 2);

%Calculate the n = 10 probability density function for the random angle
%prob_density10=Z10*cos(phi_values).^(zmat(2,1) - 2);

%Calculate the 15 probability density function for the random angle
%prob_density15=Z15*cos(phi_values).^(zmat(3,1) - 2);

%Calculate the n = 20 probability density function for the random angle
%prob_density20=Z20*cos(phi_values).^(zmat(4,1) - 2);

%Calculate the n = 25 probability density function for the random angle
%prob_density25=Z25*cos(phi_values).^(zmat(5,1) - 2);

%Calculate the 30 probability density function for the random angle
%prob_density30=Z30*cos(phi_values).^(zmat(6,1) - 2);

%Calculate the n = 35 probability density function for the random angle
%prob_density35=Z35*cos(phi_values).^(zmat(7,1) - 2);

%Calculate the 40 probability density function for the random angle
%prob_density40=Z40*cos(phi_values).^(zmat(8,1) - 2);

%Initialize Results of kept population
results.co=zeros(time,28);
results.co(:,1)=1:time;
%colname={'time','meanfitness','nogenmut','noepimut','distgen','distepi','net_fi_carbon_gen','net_fi_p_gen','net_fi_i_gen','net_fi_carbon_epi','net_fi_p_epi','net_fi_i_epi'};
%colindex={1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
%import bioma.data.*
%(http://www.mathworks.com/help/bioinfo/ug/representing-expression-data-values-in-datamatrix-objects.html)
%results.c=DataMatrix(result,1:time,colname);
%results.p=DataMatrix(result,1:time,colname);
%results.i=DataMatrix(result,1:time,colname);
%results.co=DataMatrix(result,1:time,colname);

%Initialize Results of removed population
results_removed=zeros(time,28);
results_removed(:,1)=1:time;


select = 1; %Switch to randomly resample population either not weighted (0) or weighted (1) by fitness 

%%%%%%%%%%%%%%%%%%%%%
%%Looping over time
%%%%%%%%%%%%%%%%%%%%%
for i=1:length(Driver)
    

    
    patch_tf = strcmp(patch,'FALSE'); %If true, only sample with weighting

    if patch_tf == 1
        
        select = 1; %only sample with weighting

    elseif (patch_tf == 0) && (Driver(i) == 1)
        
        
        select = 1; %resample population with weight by fitness
        %time_count = time_count + 1;
        
        
    elseif (patch_tf == 0) && (Driver(i) == 0)
    
        
        select = 0; %resample population with other weighting parameters
   
        
    end

    %%%%%%%%%%
    %Make genetic mutations, mutational supply is set to nugen/1 unit time
    %Distances from optimum are calculated individually
    %%%%%%%%%%
    
    %Magnitude of the vector
    c=clock;
    pause(c(6)*runnumber/100)
    rng('shuffle')
    mut_effect=0+(mlimitgen*radius).*rand(1,nugen); 
    
    %Angle of 5 vector
    c=clock;
    pause(c(6)*runnumber/100)
    rng('shuffle')
    phi5=randsample(phi_values, nugen, true, prob_density5); 
    
    %Angle of 10 vector
    %phi10=randsample(phi_values, nugen, true, prob_density10); 
    
    %Angle of 15 vector
    %phi15=randsample(phi_values, nugen, true, prob_density15); 
    
    %Angle of 20 vector
    %phi20=randsample(phi_values, nugen, true, prob_density20); 
    
    %Angle of 25 vector
    %phi25=randsample(phi_values, nugen, true, prob_density25); 
    
    %Angle of 30 vector
    %phi30=randsample(phi_values, nugen, true, prob_density30); 
    
    %Angle of 35 vector
    %phi35=randsample(phi_values, nugen, true, prob_density35); 
    
    %Angle of 40 vector
    %phi40=randsample(phi_values, nugen, true, prob_density40);
    
    %%%Assign mutations to individuals%%%
    
    index_mut=1:N;
    
    %Sampling mutated individuals without replacement (multiple hits are not allowed)
    c=clock;
    pause(c(6)*runnumber/100)
    rng('shuffle')
    index_mutated=randsample(index_mut,nugen,false); 
    
    %retrieve mutants from population 
    mutated_orgs=population.co(index_mutated,:);
    %mutated_orgs=population.c(index_mutated,:);
    %mutated_orgs=population.cp(index_mutated,:);
    %mutated_orgs=population.ci(index_mutated,:);
    %mut_orgs_gen=population_disgen(index_mutated,:);
    
    
    %%%distance traveled after n = 5 genetic mutation
    dist_travelled5 = [mutated_orgs(:,2)' - sqrt(mutated_orgs(:,2)'.^2 + ...
    mut_effect.^2 + 2.*mutated_orgs(:,2)'.*mut_effect.*sin(phi5))];
       
    %%%distance traveled after n = 10 genetic mutation
    %dist_travelled10 = [mutated_orgs(:,2)' - sqrt(mutated_orgs(:,2)'.^2 + ...
    %mut_effect.^2 + 2.*mutated_orgs(:,2)'.*mut_effect.*sin(phi10))];

    %%%distance traveled after n = 15 genetic mutation
    %dist_travelled15 = [mutated_orgs(:,2)' - sqrt(mutated_orgs(:,2)'.^2 + ...
    %mut_effect.^2 + 2.*mutated_orgs(:,2)'.*mut_effect.*sin(phi15))];

    %%%distance traveled after n = 20 genetic mutation
    %dist_travelled20 = [mutated_orgs(:,2)' - sqrt(mutated_orgs(:,2)'.^2 + ...
    %mut_effect.^2 + 2.*mutated_orgs(:,2)'.*mut_effect.*sin(phi20))];
       
    %%%distance traveled after n = 25 genetic mutation
    %dist_travelled25 = [mutated_orgs(:,2)' - sqrt(mutated_orgs(:,2)'.^2 + ...
    %mut_effect.^2 + 2.*mutated_orgs(:,2)'.*mut_effect.*sin(phi25))];

    %%%distance traveled after n = 30 genetic mutation
    %dist_travelled30 = [mutated_orgs(:,2)' - sqrt(mutated_orgs(:,2)'.^2 + ...
    %mut_effect.^2 + 2.*mutated_orgs(:,2)'.*mut_effect.*sin(phi30))];

    %%%distance traveled after n = 35 genetic mutation
    %dist_travelled35 = [mutated_orgs(:,2)' - sqrt(mutated_orgs(:,2)'.^2 + ...
    %mut_effect.^2 + 2.*mutated_orgs(:,2)'.*mut_effect.*sin(phi35))];

    %%%distance traveled after n = 40 genetic mutation
    %dist_travelled40 = [mutated_orgs(:,2)' - sqrt(mutated_orgs(:,2)'.^2 + ...
    %mut_effect.^2 + 2.*mutated_orgs(:,2)'.*mut_effect.*sin(phi40))];

    if replicates == 4 %carbon + p/i co-limitation
        

        %Distance from the optimum phenotype for n = 5 mutations
        dist5 = sqrt(mutated_orgs(:,2)'.^2 + mut_effect.^2 ...
        + 2.*sin(phi5).*mutated_orgs(:,2)'.*mut_effect);
        
        %Distance from the optimum phenotype for n = 10 mutations
        %dist10 = sqrt(mutated_orgs(:,2)'.^2 + mut_effect.^2 ...
        %+ 2.*sin(phi10).*mutated_orgs(:,2)'.*mut_effect);
    
        %Distance from the optimum phenotype for n = 15 mutations
        %dist15 = sqrt(mutated_orgs(:,2)'.^2 + mut_effect.^2 ...
        %+ 2.*sin(phi15).*mutated_orgs(:,2)'.*mut_effect);
    
        %Distance from the optimum phenotype for n = 20 mutations
        %dist20 = sqrt(mutated_orgs(:,2)'.^2 + mut_effect.^2 ...
        %+ 2.*sin(phi20).*mutated_orgs(:,2)'.*mut_effect);
        
        %Distance from the optimum phenotype for n = 25 mutations
        %dist25 = sqrt(mutated_orgs(:,2)'.^2 + mut_effect.^2 ...
        %+ 2.*sin(phi25).*mutated_orgs(:,2)'.*mut_effect);
    
        %Distance from the optimum phenotype for n = 30 mutations
        %dist30 = sqrt(mutated_orgs(:,2)'.^2 + mut_effect.^2 ...
        %+ 2.*sin(phi30).*mutated_orgs(:,2)'.*mut_effect);
        
        %Distance from the optimum phenotype for n = 35 mutations
        %dist35 = sqrt(mutated_orgs(:,2)'.^2 + mut_effect.^2 ...
        %+ 2.*sin(phi35).*mutated_orgs(:,2)'.*mut_effect);
    
        %Distance from the optimum phenotype for n = 40 mutations
        %dist40 = sqrt(mutated_orgs(:,2)'.^2 + mut_effect.^2 ...
        %+ 2.*sin(phi40).*mutated_orgs(:,2)'.*mut_effect);
   
        %%Distance travelled by n = 5 + 20 
        %dist_travelled_5_20 = horzcat(dist_travelled5', dist_travelled20');
        %dist_travelled_5_and_20 = mean(dist_travelled_5_20, 2)';
        
        %%Distance travelled by n = 5 + 40 
        %dist_travelled_5_40 = horzcat(dist_travelled5', dist_travelled40');
        %dist_travelled_5_and_40 = mean(dist_travelled_5_40, 2)';
        
        %%%%%%%%%%%%%%%%%%%%
        %Distance from the optimum phenotype based on oligo conditions;
        %%%%%%%%%%%%%%%%%%%%
        
        if oligo == 1 %n = 5
            
            %n = 5
            mutated_orgs(:,2) = dist5'; 
            
            %distance traveled after mutation
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
    
    %Last time point of distance travelled after genetic mutations
    mutated_orgs(:,6) = dist_travelled_gen';
    
    %Using law of sines
    %(https://www.mathsisfun.com/algebra/trig-sine-law.html) to
    %calculate fi n = 5
    %[c_tmp] = find(phi5<0);
    fi5 = asin((sin(phi5 + (pi/2)).*(mut_effect))./dist5); 
    %fi5(c_tmp) = -fi5(c_tmp);


    %Using law of sines
    %(https://www.mathsisfun.com/algebra/trig-sine-law.html) to
    %calculate fi n = 10
    %[p_tmp] = find(phi10<0);
    %fi10 = asin((sin(phi10 + (pi/2)).*(mut_effect))./dist10);
    %fi10(p_tmp) = -fi10(p_tmp);

    %Using law of sines
    %(https://www.mathsisfun.com/algebra/trig-sine-law.html) to
    %calculate fi for n = 15
    %[i_tmp] = find(phi15<0);
    %fi15 = asin((sin(phi15 + (pi/2)).*(mut_effect))./dist15);  
    %fi15(i_tmp) = -fi15(i_tmp);
    
    %Using law of sines
    %(https://www.mathsisfun.com/algebra/trig-sine-law.html) to
    %calculate fi for n = 20
    %[i_tmp] = find(phi15<0);
    %fi20 = asin((sin(phi20 + (pi/2)).*(mut_effect))./dist20);  
    %fi15(i_tmp) = -fi15(i_tmp);
    
    %Using law of sines
    %(https://www.mathsisfun.com/algebra/trig-sine-law.html) to
    %calculate fi for n = 25
    %[i_tmp] = find(phi15<0);
    %fi25 = asin((sin(phi25 + (pi/2)).*(mut_effect))./dist25);  
    %fi15(i_tmp) = -fi15(i_tmp);
    
    %Using law of sines
    %(https://www.mathsisfun.com/algebra/trig-sine-law.html) to
    %calculate fi for n = 30
    %[i_tmp] = find(phi15<0);
    %fi30 = asin((sin(phi30 + (pi/2)).*(mut_effect))./dist30);  
    %fi15(i_tmp) = -fi15(i_tmp);
    
    %Using law of sines
    %(https://www.mathsisfun.com/algebra/trig-sine-law.html) to
    %calculate fi for n = 35
    %[i_tmp] = find(phi15<0);
    %fi35 = asin((sin(phi35 + (pi/2)).*(mut_effect))./dist35);  
    %fi15(i_tmp) = -fi15(i_tmp);
    
    %Using law of sines
    %(https://www.mathsisfun.com/algebra/trig-sine-law.html) to
    %calculate fi for n = 40
    %[i_tmp] = find(phi15<0);
    %fi40 = asin((sin(phi40 + (pi/2)).*(mut_effect))./dist40);  
    %fi15(i_tmp) = -fi15(i_tmp);
    
    %Insert fi5
    mutated_orgs(:,10) = mutated_orgs(:,10) + fi5';
    
    %Insert fi10
    %mutated_orgs(:,11) = mutated_orgs(:,11) + fi10';
    
    %Insert fi15
    %mutated_orgs(:,12) = mutated_orgs(:,12) + fi15';
    
    %Insert fi20
    %mutated_orgs(:,13) = mutated_orgs(:,13) + fi20';
    
    %Insert fi25
    %mutated_orgs(:,14) = mutated_orgs(:,14) + fi25';
    
    %Insert fi30
    %mutated_orgs(:,15) = mutated_orgs(:,15) + fi30';
    
    %Insert fi35
    %mutated_orgs(:,16) = mutated_orgs(:,16) + fi35';
    
    %Insert fi40
    %mutated_orgs(:,17) = mutated_orgs(:,17) + fi40';
    
    
    %Storing the mutated individuals
    population.co(index_mutated,:)=mutated_orgs;
    %population(index_mutated,:) = mutated_index;
    
    %Distribution of genetic mutations
    %population_disgen(index_mutated) = dist_travelled_gen';
    
    for g=1:length(dist_travelled_gen')
        
        population_disgen{index_mutated(g), i} = strcat(num2str(dist_travelled_gen(g)));
    
    end
    
    population_numgen(:,i) = population.co(:,5); %save number of genetic mutations per generation
   
    %If the simulation includes epigenetic effects, run the following
    tf = strcmp(epigenetics,'TRUE');

    if tf == 1
        
        %Make epigenetic mutations, mutational supply is set to Nu.epi / 1 unit of time
        
        %Mutational effects of epigenetics
        c=clock;
        pause(c(6)*runnumber/100)
        rng('shuffle')
        mut_effect_epi=0+(mlimitepi*radius).*rand(1, nuepi)';
        
        %Random angle for epigenetic mutations%%%%%%%%%%%%%%%%%%%%%
        
        %Angle of n = 5 epi_vector
        c=clock;
        pause(c(6)*runnumber/100)
        rng('shuffle')
        phi_epi5=randsample(phi_values, nuepi, true, prob_density5);

        %Angle of n = 10 epi_vector
        %phi_epi10=randsample(phi_values, nuepi, true, prob_density10);

        %Angle of n = 15 epi_vector
        %phi_epi15=randsample(phi_values, nuepi, true, prob_density15);
        
        %Angle of n = 20 epi_vector
        %phi_epi20=randsample(phi_values, nuepi, true, prob_density20);

        %Angle of n = 25 epi_vector
        %phi_epi25=randsample(phi_values, nuepi, true, prob_density25);

        %Angle of n = 30 epi_vector
        %phi_epi30=randsample(phi_values, nuepi, true, prob_density30);
        
        %Angle of n = 35 epi_vector
        %phi_epi35=randsample(phi_values, nuepi, true, prob_density35);

        %Angle of n = 40 epi_vector
        %phi_epi40=randsample(phi_values, nuepi, true, prob_density40);
        
        %retrieve epimutants from population w/out replacement
        c=clock;
        pause(c(6)*runnumber/100)
        rng('shuffle')
        index_mutated_epi = randsample(index_mut, nuepi, false);
        mutated_orgs_epi = population.co(index_mutated_epi, :);
        %mut_orgs_epi=population_disepi(index_mutated_epi,:);
        
        %%%distance traveled after n = 5 epimutation
        dist_travelled_epi5 = [mutated_orgs_epi(:,2)' - sqrt(mutated_orgs_epi(:,2)'.^2 + ...
        mut_effect_epi.^2' + 2.*mutated_orgs_epi(:,2)'.*mut_effect_epi'.*sin(phi_epi5))];

        %%%distance traveled after n = 10 epimutation
        %dist_travelled_epi10 = [mutated_orgs_epi(:,2)' - sqrt(mutated_orgs_epi(:,2)'.^2 + ...
        %mut_effect_epi.^2' + 2.*mutated_orgs_epi(:,2)'.*mut_effect_epi'.*sin(phi_epi10))];

        %%%distance traveled after n = 15 epimutation
        %dist_travelled_epi15 = [mutated_orgs_epi(:,2)' - sqrt(mutated_orgs_epi(:,2)'.^2 + ...
        %mut_effect_epi.^2' + 2.*mutated_orgs_epi(:,2)'.*mut_effect_epi'.*sin(phi_epi15))];
        
        %%%distance traveled after n = 20 epimutation
        %dist_travelled_epi20 = [mutated_orgs_epi(:,2)' - sqrt(mutated_orgs_epi(:,2)'.^2 + ...
        %mut_effect_epi.^2' + 2.*mutated_orgs_epi(:,2)'.*mut_effect_epi'.*sin(phi_epi20))];

        %%%distance traveled after n = 25 epimutation
        %dist_travelled_epi25 = [mutated_orgs_epi(:,2)' - sqrt(mutated_orgs_epi(:,2)'.^2 + ...
        %mut_effect_epi.^2' + 2.*mutated_orgs_epi(:,2)'.*mut_effect_epi'.*sin(phi_epi25))];

        %%%distance traveled after n = 30 epimutation
        %dist_travelled_epi30 = [mutated_orgs_epi(:,2)' - sqrt(mutated_orgs_epi(:,2)'.^2 + ...
        %mut_effect_epi.^2' + 2.*mutated_orgs_epi(:,2)'.*mut_effect_epi'.*sin(phi_epi30))];
        
        %%%distance traveled after n = 35 epimutation
        %dist_travelled_epi35 = [mutated_orgs_epi(:,2)' - sqrt(mutated_orgs_epi(:,2)'.^2 + ...
        %mut_effect_epi.^2' + 2.*mutated_orgs_epi(:,2)'.*mut_effect_epi'.*sin(phi_epi35))];

        %%%distance traveled after n = 40 epimutation
        %dist_travelled_epi40 = [mutated_orgs_epi(:,2)' - sqrt(mutated_orgs_epi(:,2)'.^2 + ...
        %mut_effect_epi.^2' + 2.*mutated_orgs_epi(:,2)'.*mut_effect_epi'.*sin(phi_epi40))];
        
        if replicates == 4 %carbon + p/i co-limitation epimutation

            %Distance from the optimum phenotype for n = 5 epimutations
            dist_epi5 = sqrt(mutated_orgs_epi(:,2)'.^2 + mut_effect_epi.^2' ...
            + 2.*sin(phi_epi5).*mutated_orgs_epi(:,2)'.*mut_effect_epi');

            %Distance from the optimum phenotype for n = 10 epimutations
            %dist_epi10 = sqrt(mutated_orgs_epi(:,2)'.^2 + mut_effect_epi.^2' ...
            %+ 2.*sin(phi_epi10).*mutated_orgs_epi(:,2)'.*mut_effect_epi');

            %Distance from the optimum phenotype for n = 15 epimutations
            %dist_epi15 = sqrt(mutated_orgs_epi(:,2)'.^2 + mut_effect_epi.^2' ...
            %+ 2.*sin(phi_epi15).*mutated_orgs_epi(:,2)'.*mut_effect_epi');
            
            %Distance from the optimum phenotype for n = 20 epimutations
            %dist_epi20 = sqrt(mutated_orgs_epi(:,2)'.^2 + mut_effect_epi.^2' ...
            %+ 2.*sin(phi_epi20).*mutated_orgs_epi(:,2)'.*mut_effect_epi');

            %Distance from the optimum phenotype for n = 25 epimutations
            %dist_epi25 = sqrt(mutated_orgs_epi(:,2)'.^2 + mut_effect_epi.^2' ...
            %+ 2.*sin(phi_epi25).*mutated_orgs_epi(:,2)'.*mut_effect_epi');

            %Distance from the optimum phenotype for n = 30 epimutations
            %dist_epi30 = sqrt(mutated_orgs_epi(:,2)'.^2 + mut_effect_epi.^2' ...
            %+ 2.*sin(phi_epi30).*mutated_orgs_epi(:,2)'.*mut_effect_epi');
        
            %Distance from the optimum phenotype for n = 35 epimutations
            %dist_epi35 = sqrt(mutated_orgs_epi(:,2)'.^2 + mut_effect_epi.^2' ...
            %+ 2.*sin(phi_epi35).*mutated_orgs_epi(:,2)'.*mut_effect_epi');

            %Distance from the optimum phenotype for n = 40 epimutations
            %dist_epi40 = sqrt(mutated_orgs_epi(:,2)'.^2 + mut_effect_epi.^2' ...
            %+ 2.*sin(phi_epi40).*mutated_orgs_epi(:,2)'.*mut_effect_epi');

        end


         %Distance from the optimum phenotype based on oligo conditions;
         %These only affect fitness in certain nutrient-limited regimes
         
        if oligo == 1 %n = 5
            
            %Distance from the optimum phenotype 
            mutated_orgs_epi(:,2) = dist_epi5'; 

            %distance traveled after epimutation
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
        
        %Save mean distance travelled toward optimum by epigenetic mutation
        results.co(i,27) = mean(mutated_orgs_epi(:,2));
    
        %Save std of distance travelled toward optimum by epigenetic mutation
        results.co(i,28) = std(mutated_orgs_epi(:,2));
        
        %Calculate fitness for mutants
        mutated_orgs_epi(:,3) = exp((-mutated_orgs_epi(:,2).^2)/2);
        
        %Cumulative distance travelled by epigenetic mutations
        mutated_orgs_epi(:,7) = mutated_orgs_epi(:,7) + dist_travelled_e';
        
        %Increase the number of epigenetic mutations by one
        mutated_orgs_epi(:,8) = mutated_orgs_epi(:,8) + 1;
        
        %placeholder -- put most recent epimutations in this column
        mutated_orgs_epi(:,9) = dist_travelled_e';
        
        
        %Using law of sines
        %(https://www.mathsisfun.com/algebra/trig-sine-law.html) to
        %calculate fi n = 5
        %[c_epi_tmp] = find(phi_epi5<0);
        fi_epi5 = asin((sin(phi_epi5 + (pi/2)).*(mut_effect_epi'))./dist_epi5);
        %fi_epi5(c_epi_tmp) = -fi_epi5(c_epi_tmp);

        %Using law of sines
        %(https://www.mathsisfun.com/algebra/trig-sine-law.html) to
        %calculate fi n = 10
        %[p_epi_tmp] = find(phi_epi10<0);
        %fi_epi10 = asin((sin(phi_epi10 + (pi/2)).*(mut_effect_epi'))./dist_epi10);
        %fi_epi10(p_epi_tmp) = -fi_epi10(p_epi_tmp);

        %Using law of sines
        %(https://www.mathsisfun.com/algebra/trig-sine-law.html) to
        %calculate fi n = 15
        %[i_epi_tmp] = find(phi_epi15<0);
        %fi_epi15 = asin((sin(phi_epi15 + (pi/2)).*(mut_effect_epi'))./dist_epi15);
        %fi_epi15(i_epi_tmp) = -fi_epi15(i_epi_tmp);
        
        %Using law of sines
        %(https://www.mathsisfun.com/algebra/trig-sine-law.html) to
        %calculate fi n = 20
        %[c_epi_tmp] = find(phi_epi5<0);
        %fi_epi20 = asin((sin(phi_epi20 + (pi/2)).*(mut_effect_epi'))./dist_epi20);
        %fi_epi5(c_epi_tmp) = -fi_epi5(c_epi_tmp);

        %Using law of sines
        %(https://www.mathsisfun.com/algebra/trig-sine-law.html) to
        %calculate fi n = 25
        %[p_epi_tmp] = find(phi_epi10<0);
        %fi_epi25 = asin((sin(phi_epi25 + (pi/2)).*(mut_effect_epi'))./dist_epi25);
        %fi_epi10(p_epi_tmp) = -fi_epi10(p_epi_tmp);

        %Using law of sines
        %(https://www.mathsisfun.com/algebra/trig-sine-law.html) to
        %calculate fi n = 30
        %[i_epi_tmp] = find(phi_epi15<0);
        %fi_epi30 = asin((sin(phi_epi30 + (pi/2)).*(mut_effect_epi'))./dist_epi30);
        %fi_epi15(i_epi_tmp) = -fi_epi15(i_epi_tmp);
        
        %Using law of sines
        %(https://www.mathsisfun.com/algebra/trig-sine-law.html) to
        %calculate fi n = 35
        %[p_epi_tmp] = find(phi_epi10<0);
        %fi_epi35 = asin((sin(phi_epi35 + (pi/2)).*(mut_effect_epi'))./dist_epi35);
        %fi_epi10(p_epi_tmp) = -fi_epi10(p_epi_tmp);

        %Using law of sines
        %(https://www.mathsisfun.com/algebra/trig-sine-law.html) to
        %calculate fi n = 40
        %[i_epi_tmp] = find(phi_epi15<0);
        %fi_epi40 = asin((sin(phi_epi40 + (pi/2)).*(mut_effect_epi'))./dist_epi40);
        %fi_epi15(i_epi_tmp) = -fi_epi15(i_epi_tmp);
        
        %Insert fi5 for epi
        mutated_orgs_epi(:,18) = mutated_orgs_epi(:,18) + fi_epi5';
        
        %Insert fi10 for epi
        %mutated_orgs_epi(:,19) = mutated_orgs_epi(:,19) + fi_epi10';

        %Insert fi15 for epi
        %mutated_orgs_epi(:,20) = mutated_orgs_epi(:,20) + fi_epi15';
        
        %Insert fi20 for epi
        %mutated_orgs_epi(:,21) = mutated_orgs_epi(:,21) + fi_epi20';
        
        %Insert fi25 for epi
        %mutated_orgs_epi(:,22) = mutated_orgs_epi(:,22) + fi_epi25';

        %Insert fi30 for epi
        %mutated_orgs_epi(:,23) = mutated_orgs_epi(:,23) + fi_epi30';
        
        %Insert fi35 for epi
        %mutated_orgs_epi(:,24) = mutated_orgs_epi(:,24) + fi_epi35';

        %Insert fi40 for epi
        %mutated_orgs_epi(:,25) = mutated_orgs_epi(:,25) + fi_epi40';
        
        %Storing mutated individuals
        population.co(index_mutated_epi,:) = mutated_orgs_epi;
        
        %Distribution of epigenetic mutations
        %population_disepi(index_mutated_epi,i) = dist_travelled_e';
        
        for e=1:length(dist_travelled_e')
            
            population_disepi{index_mutated_epi(e), i} = dist_travelled_e(e);
        
        end
    
    end
    
    
    %As of now, multiple hits of either genetic and epigenetic mutations are not possible. 
    %However, a single individual can have two hits if other is genetic and other epigenetic mutation. 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Population regulation and reproduction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Soft selection is acting -> population size stays the same
    %Sampling individuals to the next time point with replacement
    
    ofsp_index = 1:N;
    
    if select == 0
        
        mode_type = strcmp(mode,'random'); %If true, only sample with weighting

        if mode_type == 1
            
            %Sampling the next generation, no weighting
            c=clock;
            pause(c(6)*runnumber/100)
            rng('shuffle')
            sampled_index = randsample(ofsp_index, N, true);
        
        else
            
            %Sampling the next generation, negative weighting
            tmp = 1./(population.co(:,5)+0.1);
            tmp2=(tmp-mean(tmp))./std(tmp);
            reciprocal_weights=(tmp2-min(tmp2))./max(tmp2-min(tmp2))*tfact;

            %w=tmp/sum(tmp);
            %max(w);

            c=clock;
            pause(c(6)*runnumber/100)
            rng('shuffle')

            sampled_index = randsample(ofsp_index, N, true, reciprocal_weights);
        
        end
        
    elseif select == 1
        
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
    population_numgen = population_numgen(sampled_index,:);
    
    
    %Removed individuals from population
    pop_removed_index = setdiff(ofsp_index, sampled_index);
    population_removed = population.co(pop_removed_index,:);
    
    population_removed_disgen = population_disgen(pop_removed_index,:);
    population_removed_disepi = population_disepi(pop_removed_index,:);
    
    
    %%%%%%Store results removed population
   
    %Mean fitness of the population
    results_removed(i,2)= mean(population_removed(:,3)); 
    %%Checking the fixation time of the first clone
    %if(any(population[,1] != population[1,1]) == FALSE) {stop(cat("Time to fixation was", i, "\n"))}
    
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
    
    
    %%%%%%Store results sampled population
    
    
    %Mean fitness of the population
    results.co(i,2)= mean(population.co(:,3)); 
    %%Checking the fixation time of the first clone
    %if(any(population[,1] != population[1,1]) == FALSE) {stop(cat("Time to fixation was", i, "\n"))}
    
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
    
    %Store mean genetic mutation fi10
    %results.co(i,8) = mean(population.co(:,11));
    
    %Store mean genetic mutation fi15
    %results.co(i,9) = mean(population.co(:,12));
    
    %Store mean genetic mutation fi20
    %results.co(i,10) = mean(population.co(:,13));
    
    %Store mean genetic mutation fi25
    %results.co(i,11) = mean(population.co(:,14));
    
    %Store mean genetic mutation fi30
    %results.co(i,12) = mean(population.co(:,15));
    
    %Store mean genetic mutation fi35
    %results.co(i,13) = mean(population.co(:,16));
    
    %Store mean genetic mutation fi40
    %results.co(i,14) = mean(population.co(:,17));

    %Store mean epigenetic mutation fi5
    results.co(i,15) = mean(population.co(:,18));
    
    %Store mean epigenetic mutation fi10
    %results.co(i,16) = mean(population.co(:,19));
    
    %Store mean epigenetic mutation fi15
    %results.co(i,17) = mean(population.co(:,20));
    
    %Store mean epigenetic mutation fi20
    %results.co(i,18) = mean(population.co(:,21));
    
    %Store mean epigenetic mutation fi25
    %results.co(i,19) = mean(population.co(:,22));
    
    %Store mean epigenetic mutation fi30
    %results.co(i,20) = mean(population.co(:,23));
    
    %Store mean epigenetic mutation fi35
    %results.co(i,21) = mean(population.co(:,24));
    
    %Store mean epigenetic mutation fi40
    %results.co(i,22) = mean(population.co(:,25));
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Loss of epigenetic effects
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Next epigenetic changes are lost at a given rate 
    %(This needs to be slower than the forward mutation rate of epigenetic effects)
    %In this new version of the algorithm, epigenetic changes revert independently
    %(always at a probability of: epi.backmutation in each generation)
    
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
        
        %%%%%%%%%%%%%%%%%same function as "rep" in R
        %x=0;
		%num_old=0;
		%help.vec=vector()
        %for j=1:length(epi_present(:,8))
		%	x=x+1;
		%	num=epi_present(j,8);
		%	if (x==1) 
				
		%		help_vec(j:num)=j;
		%		num_old=num_old+num+1;
		%	else 
		%		new_num=num_old+num;
		%		help_vec(num_old:new_num)=j;
		%		num_old=num_old+num;
        %   end
        %end
				
        %help_vec=help_vec(1:end-1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Sum by each factor
        %%%%%%%%%%%%%%%%%%%%Same function as "tapply" in R
        
        list_of_reversions = splitapply(@sum,epi_backmuts,help_vec');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Logical vector if backmutation occurs in a given individual (TRUE) or not
        %[epi_back_mut_index,col,v] = find(list_of_reversions(:,2));
        %epi_back_mut_index = list_of_reversions(:,2) > 0;
        
        epi_back_mut_index = find(list_of_reversions > 0);
        
        %Take only those individuals that had backmutations
        epi_back_mut_ind = epi_present(epi_back_mut_index,:);
        epi_back_mut_ind_gen = epi_present_gen(epi_back_mut_index,:);
        epi_back_mut_ind_epi = epi_present_epi(epi_back_mut_index,:);
        
        %How many backmutations occurred in each individual that had backmutations
        number_of_reversions = list_of_reversions(epi_back_mut_index);
        
        %Check if any backmutations occurred
        if any(list_of_reversions) > 0
            
        %epi_temp_ind=epi_back_mut_ind_epi(:,1);
        %epi_back_mut_ind_epi(:,1)=[];
        
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
                %C=[];
                %[c,epi_c,lost_c]=setxor(epi_back_mut_ind_epi(b,:),lost_mutation(:,:), 'stable');
                %[C,ia] = setdiff(current_org,lost_mutation(:,:));
                %[C,ia] = setdiff(epi_back_mut_ind_epi,lost_mutation_array);
                if ~isempty(lost_mutation)
                    
                    for count=1:length(lost_mutation)

                        logic = any(current_org == lost_mutation(1,count));

                        if logic == 1

                            r = cellfun(@num2str, epi_back_mut_ind_epi(b,:), 'UniformOutput', false);
                            idx2 = find(ismember(r, num2str(lost_mutation(1,count))));
                            if length(idx2) == 1 

                                epi_back_mut_ind_epi{b,idx2} = 0;

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
               
               %Put unlost" mutations into new matrix 
               %new_mat(b,:)=C;
               
               %Take the remaining epigenetic mutations, make a vector that
               %is as long as the difference between the number of
               %remaining mutations and the original time vector.
               
               %filler=0;
               %if filler<time
               %    filler=time-length(epi_c);
               %    filler=zeros(1,filler);
               %    C=[C filler];
               %end
               
               %Insert the new distribution of epigenetic mutations into
               %the individual that experienced the backmutation
               %epi_back_mut_ind_epi(ia)=epi_back_mut_ind_epi(ia) - C;
              
               
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

