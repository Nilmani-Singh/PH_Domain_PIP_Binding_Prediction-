%%
% this is the algorithm sript for RFC matrix generation and SRFC value
% calculation
% The fasta files contains aligned amino acid seqeunces, in order, 
% all PIP binding sequences
% non-binding sequences
% and sequences for prediction of binding 

tic 
clear;

%% Define variables 
% Read the .txt sequence file containing all aligned sequences 

Aligned_Seq = fastaread('Aligned_PH_domain_Seq_242.txt'); %read sequence file
amino_acids = 'ACDEFGHIKLMNPQRSTVWY-'; % string of all amino acids
aa_count = 20;  % number of amino acids

PIP_binders = 35;   %the number of proteins that bind PIPs 
PIP_Nbinders = 32;   %%the number of proteins that do not bind PIPs 

Len = length(Aligned_Seq(1).Sequence); % length of the aligned sequence

P_binding_mat  = zeros(aa_count,Len); %Prob matrix for binding PH domains
P_Nbinding_mat = zeros(aa_count,Len); %Prob matrix for non-binding PH domains
temp = [];

%%
% Create the probability matrix for the PIP binding proteins
clear Pos_string;
Pos_string = cell(1,Len); % for each position create a string that contains 
                     % all AA for that position 
for i = 1 : Len   % the length of the aligned sequences
    temp_1 = [];
    for j = 1 : PIP_binders  %number of PH domains that bind PIPs
        temp_1 = [temp_1 Aligned_Seq(j).Sequence(i)]; % create the 1x35 string
    end
    Pos_string{i} = temp_1;
%     if i == 17  %% manual result validation code
%         display(temp_1); 
%         length(temp_1);
%     end 
end

for j = 1 : Len
    for i = 1 : aa_count % Number of amino acids excluding gap position 21
        count = length(strfind(Pos_string{j},amino_acids(i)));
        P_binding_mat(i,j) = count/PIP_binders; % Divide by number of binding proteins
    end
end


%%
% Create the probability matrix for the non-PIP binding proteins
clear Pos_string;
Pos_string = cell(1,Len); % for each position create a string that contains 
                     % all AA for that position 
for i = 1 : Len   % the length of the aligned sequences
    temp_1 = [];
    for j = 1 : PIP_Nbinders  %number of PH domains that don't bind PIPs
        temp_1 = [temp_1 Aligned_Seq(j+PIP_binders).Sequence(i)]; % create the 1x32 string
        
        if i ==1 
            temp = [temp (j+PIP_binders)]; % to verify the correct seqeunces were read
        end 
    end
    Pos_string{i} = temp_1;
 
end

for j = 1 : Len
    for i = 1 : aa_count % Number of amino acids excluding - position at 21
        count = length(strfind(Pos_string{j},amino_acids(i)));
        P_Nbinding_mat(i,j) = count/PIP_Nbinders; % Divide by number of non-binding proteins
    end
end


%% 
% create the RFC matrix 

clear RFC;
RFC = [];
cond_1 = (1/PIP_binders) + (1/PIP_Nbinders);

for i = 1:aa_count             % amino acid positions
    for j = 1:Len              % alignement positions
      
        cond_2 = (abs(P_binding_mat(i,j))+ abs(P_Nbinding_mat(i,j)));
        
      if cond_2 > cond_1
      
              RFC(i,j) =  log(P_binding_mat(i,j)/P_Nbinding_mat(i,j));
           
           if isinf(RFC(i,j)) == 1 | isnan(RFC(i,j)) == 1
                  RFC(i,j) = 0; 
           end


      else 
           RFC(i,j) = 0;
              
      end
      
      
    end 
    
    
end

Filename = sprintf('RFC_matrix_%s.xlsx', datestr(now,'mm-dd-yyyy_HH-MM'));

xlswrite(Filename,RFC);
%%
% calculate the SRFC score
clear Name_list;
clear Predict;
Name_list = {};    % cell to store the SRFC scores

Len2 = length(Aligned_Seq);  % number of PH domains

for k = 1 : Len2                 
        Name_list{k,1} = Aligned_Seq(k).Header;
       % Name_list{k,3} = Aligned_Seq(k).Sequence;
        Seq = Aligned_Seq(k).Sequence;
        Sum_Predict = 0;  % for each protein, this will be re-initialized 
                       
        
    for i = 1:Len

       pos = strfind(amino_acids, Seq(i)); % 1-20 : amino acid code 
       
       if pos < 21
           Sum_Predict =  Sum_Predict + (RFC(pos,i));
       end 
     

    end 
    
    
   % Create a cell array conatining the name of the protein
   % followed by its RFC score
    Name_list{k,2} = Sum_Predict; 
   % Name_list{k,3} = Aligned_Seq(k).Sequence;
    
    % can write the cell array to excel file using xlswrite
    
end
%% Write to excel file

Filename_2 = sprintf('PIP_Binding_Prediction_%s.xlsx', datestr(now,'mm-dd-yy-HH-MM'));

num = PIP_binders+PIP_Nbinders;

xlswrite(Filename_2,Name_list(num+1:Len2,:),1);
xlswrite(Filename_2,Name_list(1:PIP_binders,:),2);
xlswrite(Filename_2,Name_list(PIP_binders+1:num,:),3);


toc
