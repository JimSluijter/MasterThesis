%% GRAPH GRAMMAR 
% By Jim Sluijter

% This script reproduces the graph grammar algorithm to search the topology
% design space as created by Zimmerman.
% The appropriate functions to execute rules r1 and r2, as well as their
% identifications M1 and M2, are added to the source code of the origami.m
% class in the matlabPTU source code.
% Further functions to check and filter isomorphic, symmetric, and
% semantically valid graphs are added as separate functions in the source
% code.
% An overview of the graph grammar is provided in the supplementary
% materials

%% SETUP
clear
close all

%function handles
r1 = {@r11, @r12, @r13, @r14, @r1r15, @r16};

%%USER DEFINED INPUTS
N_max = 3;              %max. number of internal vertices
filterOn = true;        %(used for diagnositcs/comparison with original results)

%%CREATE INITIAL GRAPH G_1
G_1 = origami('DoubleInput'); 
G_1.gen = 1;
G_1.g   = 1;

%%INITIATE G_SET
gen = 1;                %counter for READING generation fo G_SET

G_SET{gen,1} = G_1;     %put initial graph in set of graphs G_SET
clear G_1               %to keep workspace clean

%% GRAPH GRAMMAR

stopConditions = false;
while ~stopConditions

    %Count number of graphs in generation gen:
    N_g = nnz(~cellfun(@isempty,(G_SET(gen,:))));
    clear N_int_v M2_v
    g_next = 1; %index for input of G_nex in next generation of G_set
    
    %%% RULE 1: EXPANSION %%%
    for g = 1:N_g           %for every graph in this generation
        G = G_SET{gen,g};
        clear M1
        
        %%MATCHING
        G.matching1();      %call matching method
        M1 = G.M1;              %extract M1
            M1_V{g} = M1;       %storing M1 for stopCondition
            
        V_exp = find(M1);       %vector with all expandable nodes in graph G
        %number of interior nodes
        N_int = G.numberOfVertices('interior');
            N_int_v(g) = N_int; %Storing N_int for stopCondition
            
        %%EXPANSION LOOP
        if N_int < N_max
            for a = V_exp       %Perform separately for each expandable node
                clear G_temp
                G_temp = copy(G);   %copy needed for handle classes
                
                %Execute appropriate rule:
                G_temp.rule1(a)
                G_temp.gen = gen+1; %add appropriate labels
                G_temp.g = g_next;
                
                %Execture filters
                isIso = false;     
                isSym = false;
                if filterOn
                    isIso = filterIsomorphic(G_temp,G_SET);
                    if ~isIso
                        isSym = filterSymmetric(G_temp,G_SET);
                    end
                end
                % if G_temp is novel, add to G_SET
                if ~isIso && ~isSym
                     G_SET{gen+1,g_next} = G_temp;
                     g_next = g_next+1;
                end
            end
        end 
    end       

    %%% RULE 2: Merging %%%
    
    for g = 1:N_g
        G = G_SET{gen,g};
        if G.numberOfVertices()>3 %to prevent error at G1
         
          %%MATCHING
          G.matching2();
          M2 = G.M2;
          %store
          M2_v{g} = M2;
         
          %%apply rule 2
          for i = 1:size(M2,1)      %separately for every mergable node pair
              clear G_temp
              G_temp =  copy(G);
              c = M2(i,:);          %node pair to be merged
              G_temp.rule2(c);      %apply rule
              G_temp.gen = gen+1;
              G_temp.g = g_next;
              if filterOn
                  isIso = filterIsomorphic(G_temp,G_SET);
                  if ~isIso
                      isSym = filterSymmetric(G_temp,G_SET);
                  end
              else
                  isIso = false;
                  isSym = false;
              end
              if ~isIso && ~isSym
                  G_SET{gen+1,g_next} = G_temp;
                  g_next = g_next+1;
              end
          end
        end
    end
        
                
    %%% End of generation loop %%%
    %initiate next generation
    gen = gen+1;
    
    %Stop conditions
    %(N_max reached OR no M1 matches) AND no M2 matches    
    stopConditions = (nnz(N_int_v<N_max) == 0) && (nnz(~cellfun(@isempty,M2_v)) == 0);
end

%% SAVE
if filterOn == 1
    save('Results/sets/G_SET.mat','G_SET')
elseif filterOn == 0
    G_SET_UNF = G_SET;
    save('Results/sets/G_SET(unfiltered)','G_SET_UNF')
end