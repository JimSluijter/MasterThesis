%% Origami.m
% This code contains the main origami class, which is used to perform the
% graph grammar.
% This code is heavily based on the matlabPTU code by Andreas Walker, as
% provided on: https://doi.org/10.5905/ethz-1007-465
%%

classdef origami < matlab.mixin.Copyable
    properties (Access='public')
        vertices(1,:) vertex;   %Vector with all vertices in the origami
        edges(1,:) edge;        %Vector with all edges in the origami
        faces(1,:) face;        %Vector with all faces in the origami
        source1(1,1) {mustBeNonnegative, mustBeInteger, mustBeFinite} = 0;
                                %contains the index of the source vertex.
        source2(1,1) {mustBeNonnegative, mustBeInteger, mustBeFinite} = 0;
                                %(!) Added by Jim, contains index of second
                                %source vertex
        fixedFaceIndex(1,1) {mustBeNonnegative, mustBeInteger, mustBeFinite} = 0;
                                %Contains the index of the fixed face.
        status char {mustBeMember(status,{'editGraph','editFaces','doKinematics'})} = 'editGraph';
                                %Status of the origami. 
        symmetricFold(1,1) {mustBeNumericOrLogical} = false;
        %Always false in my code
                                
        %%%%%Handles for creating the plots of the deformed origamis%%%%%
        animationXLims(1,2) {mustBeNumeric, mustBeReal, mustBeFinite} = [-1, 5];
        animationYLims(1,2) {mustBeNumeric, mustBeReal, mustBeFinite} = [-3, 3];
        animationZLims(1,2) {mustBeNumeric, mustBeReal, mustBeFinite} = [-3, 3];
                                %limits of the x axis when plotting the
                                %deformed origami and/or creating the
                                %animation
        animationView(1,2) {mustBeNumeric, mustBeReal, mustBeFinite} = [-37.5, 30];
                                %View angles for the plotting of the
                                %deformed origami and/or creating the
                                %animation.
        yRotateAngle(1,1) {mustBeNumeric, mustBeReal, mustBeFinite} = 0;
                                %Angle by which deformed plots get rotate
                                %around the y-axis
        %%%%%Properties added by me%%%%%%
        label string = "";      %Label of rule operations performed   
        M1(1,:);                %Matching of all expandable nodes
        M2(:,2);                %Matching of all mergable nodes
        gen(1,1) uint8; %generation
        g(1,1) uint8;   %graph index in generation
        
        SemValid(1,1) {mustBeNumericOrLogical} = true; %If true: graph is semantically valid
        
        %%%%%For optimization%%%%%%
        grippingCandidates(1,:) {mustBeNumeric, mustBeNonnegative}; 
            %list of gripping candidates such that the graph remains
        vertexToGripWith(1,1) {mustBeNumeric, mustBeNonnegative} = 0;
            %index of the vertex that should reach the target point.
            
        %Variables
        R1(1,1)  {mustBePositive} = 1;
        R2(1,1)  {mustBePositive} = 0.8;
        a31(1,1) {mustBePositive, mustBeLessThan(a31,180)} = 80;
        a32(1,1) {mustBePositive, mustBeLessThan(a32,180)} = 40;
        a33(1,1) {mustBePositive, mustBeLessThan(a33,180)} = 50;
        a41(1,1) {mustBePositive, mustBeLessThan(a41,180)} = 120;
        a42(1,1) {mustBePositive, mustBeLessThan(a42,180)} = 90;
    end
    properties (Access='protected')
        isOptimizing(1,1) {mustBeNumericOrLogical} = false;
                                %If this property is "true", an
                                %optimization is in progress and all state
                                %checks get a free pass      
    end
%%%%PUBLIC METHODS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access='public')
        %%%%%Constructor%%%%%
        %for my project, I only use 'SingleCrease'
        function obj = origami(type,args)
            if strcmp(type,'SingleCrease') 
                len=1;          %default length is 1
                if nargin>=2    %When length is given
                    len=args;   %length = args (input lenght)
                end
                obj.setSourceVertex([0;0]);     %                          (see setSourceVertex)
                obj.addVertex([len;0],1,false); %sets second vertex:       (see addVertex)
                                                % at coords [len;0]
                                                % originVertex = 1
                                                % mode is false by default
                %set T
                obj.vertices(1).T = false;
                obj.vertices(2).T = true;
                obj.grippingCandidates = [2];
                
            %%%%%%%%% 2 inputs %%%%%%%%%
            elseif strcmp(type,'DoubleInput')
                
            %Create two source vertices
            obj.vertices=vertex([-1/sqrt(2);1/sqrt(2)], 1, false,'source');
            obj.source1=1; 
            obj.vertices= [obj.vertices vertex([-1/sqrt(2);-1/sqrt(2)], 2, false, 'source')];
            obj.source2=2;
            %Create third vertex and add extra crease
             obj.addVertex([0;0], 1,false); 
             obj.addCrease(2,3);
            
            %%Set expansion labels
            obj.vertices(1).T = false;
            obj.vertices(2).T = false;
%             obj.grippingCandidates = [3]; disabled for now 
            else
                error('specified type does not match any of the available options');
            end
        end
        
%%%%%%%%ADDED BY ME%%%%
        function [] = matching1(obj)
            
            obj.checkState('editGraph')
            
            numNodes = obj.numberOfVertices();
            %set all M1 to zero
            M1 = zeros(1,numNodes);
            for a = 1:numNodes        %for all nodes
                %Check if vertex deleted
                if obj.vertices(a).deleted
                    continue
                end

                %Check if node is expandable
                if obj.vertices(a).T 
                    if isempty(obj.vertices(a).outgoingNeighbours)
                        %r11 identified
                        obj.M1(a) = 1;
                    elseif length(obj.vertices(a).outgoingNeighbours) == 1
                        %r12 or r13
                        s = obj.vertices(a).outgoingEdges; %successors (only 1)
                        p1 = obj.vertices(a).incomingEdges(1); %predecessors
                        pm = obj.vertices(a).incomingEdges(end);
                        
                        FLs     = obj.edges(s).faceLeft;
                        FRs     = obj.edges(s).faceRight;
                        FLp1    = obj.edges(p1).faceLeft;
                        FRpm    = obj.edges(pm).faceRight;
                        if FLs == FLp1
                            obj.M1(a) = 2;%r12
                        elseif FRs == FRpm
                            obj.M1(a) = 3;%r13
                        else
                            error("error: neither r12 or r13")
                        end
                    elseif length(obj.vertices(a).outgoingNeighbours) == 2
                        %r14 r15 r16
                        s1 = obj.vertices(a).outgoingEdges(1);
                        s2 = obj.vertices(a).outgoingEdges(2);                        
                        p1 = obj.vertices(a).incomingEdges(1);
                        pm = obj.vertices(a).incomingEdges(end);
                        
                        FLs2 = obj.edges(s2).faceLeft;
                        FLp1 = obj.edges(p1).faceLeft;
                        FRs2 = obj.edges(s2).faceRight;
                        FLs1 = obj.edges(s1).faceLeft;
                        FRs1 = obj.edges(s1).faceRight;
                        FRpm = obj.edges(pm).faceRight;
                        
                        if FLs2 == FLp1 && FRs2 == FLs1
                            obj.M1(a) = 4; %r14
                        elseif FRa1 == FRpm && FLs2 == FLp1
                            obj.M1(a) = 5; %r15
                        elseif FRs1 == FRpm && FLs1 == FRs2
                            obj.M1(a) = 6; %r15
                        else
                            error("error: neither r14,r15,r16")
                        end
                    elseif length(obj.vertices(a).outgoingNeighbours) == 2
                        obj.vertices(a).T = false;
                    end
                end %if T=true
            end%for all a
        end %matching1
        function [] = rule1(obj,vertexIndex)
            obj.checkState('editGraph')
            %functions
            r1 = {@r11, @r12, @r13, @r14, @r15, @r16};
            %call appropriate rule1
            r1{obj.M1(vertexIndex)}(obj,vertexIndex)        
            %update T of vertex
            obj.vertices(vertexIndex).T = false;
            %update label of graph
            obj.label = obj.label+num2str(vertexIndex)+",";
            %clear matching
            obj.M1 = [];
            obj.M2 = [];
            
            %Update gripping candidates
            obj.updateGrippingCandidates(vertexIndex,'rule1');
            
        end
        function [] = r11(obj,vertexIndex)
            lengths = 1; %(!) make sure these lengths and angles change in generations in all rules
           
            ycoord = abs(obj.vertices(vertexIndex).position0(2));
            angle = double(70-4*ycoord);
            %%Deciding the direction of outgoing edges
            %bisection of the "inputs sector"
            firstIncomingEdge = obj.vertices(vertexIndex).incomingEdges(1);
            firstVec = obj.awayVector0(vertexIndex,firstIncomingEdge);
            lastIncomingEdge=obj.vertices(vertexIndex).incomingEdges(end);
            lastVec=obj.awayVector0(vertexIndex,lastIncomingEdge);
            incomingVec = unitVec(firstVec+lastVec);
            %vectors for directions of outgoing edges
            outgoingVec=-lengths*incomingVec;
            outgoingVecLeft=rotz(angle)*outgoingVec;
            outgoingVecRight=rotz(-angle)*outgoingVec;
            %position of expandable vertex
            pos = obj.vertices(vertexIndex).position0;
            
            %%Add actual vertices in counter-clockwise order
            obj.addVertex(pos+outgoingVecRight,vertexIndex,false);
            obj.addVertex(pos+outgoingVec,vertexIndex,false)  
            obj.addVertex(pos+outgoingVecLeft,vertexIndex,false);

        end
        
        function [] = matching2(obj)
            obj.checkState('editGraph')
            
            obj.M2 = []; %clear previous M2
            
            %1). Both vertices should be expandable
             %set of all expandable vertices:
             V_exp = find(cat(2,obj.vertices.T));
            %2). Border same facets
            obj.createFaces();
            k = 1; %counter for input index in cond2
            cond2 = [];         %prevents error in case no vertex pair meets condition
             for a1 = V_exp
             for a2 = V_exp
             if a1 ~= a2
                 %Edges
                 p11 = obj.vertices(a1).incomingEdges(1); %predecessors
                 p1m = obj.vertices(a1).incomingEdges(end);
                 p21 = obj.vertices(a2).incomingEdges(1);
                 p2m = obj.vertices(a2).incomingEdges(end);
                 %Extract facet labels
                 FLp11a1 = obj.edges(p11).faceLeft;
                 FRp1ma1 = obj.edges(p1m).faceRight;
                 FLp21a2 = obj.edges(p21).faceLeft;
                 FRp2ma2 = obj.edges(p2m).faceRight;
                 %Compare facets
                 %if: The rightmost facet of vertex 1 is the same as the
                  %leftmost facet of vertex 2; OR vice versa, cond2 is
                  %satisfied
                 if FRp1ma1 == FLp21a2 || FRp2ma2 == FLp11a1
                     cond2(k,:) = [a1 a2];
                     k = k+1;
                 end
             end
             end
             end%cond2 loop
             obj.goToEditGraph();
             %Check if cond2 is met by any pair, if not, return empty M2
             if isempty(cond2)
                 obj.M2 = [];
                 %"Cond2 not met"
                 return
             end
             
            %3). Cannot have the same predecessors
            ko=0; %counter for extraction index out cond2
            ki=1; %counter for input index cond3
            
            cond3 = []; %creates empty matrix, prevents error in case no
                        %vertex pair meets conditions
            for a1 = cond2(:,1)'
                ko = ko+1;
                a2 = cond2(ko,2);
                clear p1 p2
                p1 = obj.vertices(a1).incomingNeighbours; %vector with all predecessors of a1
                p2 = obj.vertices(a2).incomingNeighbours;
                if isempty(intersect(p1,p2))    %if NO the same predecessor
                    cond3(ki,:) = [a1 a2];
                    ki = ki+1;
                end
            end
            %Check if cond3 is met by any pair, if not, return empty M2
            if isempty(cond3)
                obj.M2 = [];
                %"Cond3 not met"
                return
            end
            
            %4.) Cannot both be in G_1
            ko=0; %counter for extraction index out cond3
            ki=1; %counter for input index cond4
            cond4 = zeros(1,2); %prevents error when cond4 doesnt exist yet
            %write condition 4...
            for a1 = cond3(:,1)'
                ko = ko+1;
                a2 = cond3(ko,2);
                %if cond4 is met: filter duplicates:
                if isempty(find(all(cond4 == [a2,a1],2))); %if pair is not in cond4 yet
                    cond4(ki,:) = [a1 a2];
                    ki = ki+1;
                end
            end
            if cond4(1) == 0 && cond4(2) == 0 %if zeros matrix not changed
                obj.M2 = [];
                %"Cond4 not met"
                return
            end
            
            obj.M2 = cond4;
            
                
        end%function matching2
        function [] = rule2(obj,c)
            obj.checkState('editGraph')
            c1 = c(1);
            c2 = c(2);


            %1. create edge from all p2 to c1
            p2 = obj.vertices(c2).incomingNeighbours;   %predecessors of c2
            for i = p2
                obj.addCrease(i,c1)
            end
            %2. delete edges from p2 to c2
            while size(obj.vertices(c2).incomingEdges)>0
                inEdge = obj.vertices(c2).incomingEdges(1);
                
                obj.removeCreaseByIndex(inEdge)
            end
            %Assert that no connections around the vertex remain
            if not(size(obj.vertices(c2).incomingEdges,2) == 0 && ...
                    size(obj.vertices(c2).outgoingEdges,2) == 0)
                error('Standard procedures were not able to completely disconnect the vertex from its surroundings');
            end
            %3. delete c2
            obj.vertices(c2).deleted = true;
            obj.vertices(c2).T = false;
            
            %update label of graph
            obj.label = obj.label+"("+num2str(c1)+","+num2str(c2)+"),";
            %clear matching
            obj.M1 = [];
            obj.M2 = [];
            
            %update gripping candidates
            obj.updateGrippingCandidates(c, 'rule2');
            
        end
        
        function [] = updateGrippingCandidates(obj,lastVertices, ruleType)
            % See file for changes in algorithm
            % ruleType can either be "rule1" or "rule2" 
            % (!) rule1 and rule2 are disabled for now so only the aux
            % vertices are gripping candidates!! (!)
            
            %%RULE1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if strcmp(ruleType,'rule1')
%                 %check size of input lastVertices
%                 if ~all(size(lastVertices) == [1, 1])
%                     error("Error: for rule1, size of input lastVertices should be (1 1). (currenly is ("+num2str(size(lastVertices))+"))")
%                 end
%                 
%                %if last vertex was also a gripping candidate--> new vertex
%                %also
%                 if any(obj.grippingCandidates == lastVertices)
%                    %delete all gripping candidates
%                    obj.grippingCandidates = [];
%                    %add new gripping candidates
%                    obj.grippingCandidates = obj.vertices(lastVertices).outgoingNeighbours;
%                 else %if last vertex was no gripping candidate --> new 
%                      %vertices are not either, and neither are any other vertices 
%                    obj.grippingCandidates = [];
%                 end
%             
            %%RULE2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif strcmp(ruleType, 'rule2')
%                 %check size of input lastVertices (two inputs for merger)
%                 if ~all(size(lastVertices) == [1 2])
%                     error("Error: for rule2, size of input lastVertices should be (1 2). (currently is ("+num2str(size(lastVertices))+"))")
%                 end
%                 
%                 %check if any of lastVertices is a grippingVertex
%                 if any(obj.grippingCandidates == lastVertices(:),'all')
%                    %scenario 1: only the non-dominant node will be taken
%                    %out
%                    obj.grippingCandidates = [obj.grippingCandidates lastVertices(1)];
%                    obj.grippingCandidates = obj.grippingCandidates(find(obj.grippingCandidates~=lastVertices(2))); 
%                    obj.grippingCandidates = sort(obj.grippingCandidates);
%                     %deletes the second number of lastVertices from gripCan
%                 else
%                    %check if for the new merged node, lastVertices(1) all
%                    %internal nodes are present in the predecessors
%                    allPreds = obj.allPredecessors(lastVertices(1));
%                    
%                    %all internal nodes
%                    intNodes = [];
%                    for i = cat(2,obj.vertices.index)
%                         if strcmp(obj.vertices(i).type, 'interior')
%                             intNodes = [intNodes i];
%                         end
%                    end
%                    
%                    %See if all internal nodes are present in allPreds
%                    if all(any(allPreds == intNodes(:),2))
%                        %scenario 3
%                        obj.grippingCandidates = [obj.grippingCandidates lastVertices(1)];
%                        obj.grippingCandidates = unique(obj.grippingCandidates);
%      
%                    end
%                 end
            %%AuxVertices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif strcmp(ruleType,'aux')
                auxVertices = lastVertices;
                
                %for all new auxVertices
                for thisVertex = auxVertices
                    %find all predecessors
                    allPreds = obj.allPredecessors(thisVertex);
                    %find all internal nodes
                   intNodes = [];
                   for i = cat(2,obj.vertices.index)
                        if strcmp(obj.vertices(i).type, 'interior')
                            intNodes = [intNodes i];
                        end
                   end
                   %check if all internal nodes are present in allPreds
                   if all(any(allPreds == intNodes(:),2))
                       obj.grippingCandidates = [obj.grippingCandidates thisVertex];
                       obj.grippingCandidates = unique(obj.grippingCandidates);
                   end 
                end
            else
                error('specified rule type does not match any of the available options')
            end

            
        end
        function [] = findGrippingCandidatesOLD(obj)

                %find all aux vertices
                N_aux = obj.numberOfVertices('aux');
                auxVertices = zeros(1,N_aux);
                k=1;
                for j =1:obj.numberOfVertices()
                    if strcmp(obj.vertices(j).type,'aux')
                        auxVertices(k) = j;
                        k = k+1;
                    end
                end
                
                %clear grippingCandidates property
                obj.grippingCandidates = [];
                %for all new auxVertices
                for thisVertex = auxVertices
                    %find all predecessors
                    allPreds = obj.allPredecessors(thisVertex);
                    %find all internal nodes
                   intNodes = [];
                   for i = cat(2,obj.vertices.index)
                        if strcmp(obj.vertices(i).type, 'interior')
                            intNodes = [intNodes i];
                        end
                   end
                   %check if all internal nodes are present in allPreds
                   if all(any(allPreds == intNodes(:),2))
                       obj.grippingCandidates = [obj.grippingCandidates thisVertex];
                       obj.grippingCandidates = unique(obj.grippingCandidates);
                   end 
                end
        end
        
        function [] = findGrippingCandidates(obj)
            %clear grippingCandidates property
            obj.grippingCandidates = [];
            
            %find all internal nodes (used in step 2)
            intNodes = [];
            for i = cat(2,obj.vertices.index)
                if strcmp(obj.vertices(i).type, 'interior')
                    intNodes = [intNodes i];
                end
            end
            
            %adjacent facets to first internal node
            fixedIntNode = intNodes(1);
            adjFacets = obj.vertices(fixedIntNode).adjacentFaces;
            
            %For every vertex:
            numVers = obj.numberOfVertices;
            for i = 1:numVers
                %%STEP 1: Check if deleted --> clearly not a suitable EE
                if obj.vertices(i).deleted
                    continue
                end
                %%STEP 2: Check if they have all internal nodes in their
                %predecessors (semantic validity)
                
                %find all predecessors
                allPreds = obj.allPredecessors(i);
                
                %check if all internal nodes are present in allPreds, if
                %not --> continue
                if ~(all(any(allPreds == intNodes(:),2)))
                    continue
                end 
                %%STEP 3: Check if the facet they lie on borders the fixed
                %internal node. (These can only move in spherical surface
                %i.e. cannot form a straight horizontal line)
                if any(intersect(obj.vertices(i).adjacentFaces, adjFacets))
                    continue
                end
                
                %%If vertex i has survived the three steps, it can be
                %%considered a grippingCandidate
                obj.grippingCandidates = [obj.grippingCandidates i];    
            end
            %sort in the end
            obj.grippingCandidates = unique(obj.grippingCandidates);
        end
        
        function allPreds = allPredecessors(obj, vertexIndex)
            allPreds = [vertexIndex]; %total set of predecessors
            genPreds = vertexIndex; %predecessors of this gen
            
            while ~isempty(genPreds)
                newNodes = genPreds;
                genPreds = [];
                for newNode = newNodes
                    nodePreds = obj.vertices(newNode).incomingNeighbours;
                        %predecessors of current newNode
                    genPreds = [genPreds nodePreds];
                    allPreds = [allPreds nodePreds];
                end
            end
            allPreds = unique(allPreds);
        end
        
        %%Scripts for variable inputs of Parallel.mat
        
        function [] = updateDimensions(obj,varargin)
            %(!) Updates node positions based on defined variables.
            %Everything mus be updated at once in order to prevent
            %non-synchronized node positions
            %To increase speed, only the necessary node updates are
            %performed

            %create flag container
            flag = false(1,7);
            %% Update all given variables
            if numel(varargin) >0
                for ii = 1:2:numel(varargin)
                    varName = varargin{ii};
                    if ii+1 > numel(varargin)
                        error("no value assigned to input variable" +varName)
                    end
                    %R1
                    if strcmp(varName,'R1')
                        %"R1 updated"
                        flag(1) = true;
                        obj.R1 = varargin{ii+1};
                    elseif strcmp(varName,'R2')
                        %"R2 updated"
                        flag(2) = true;
                        obj.R2 = varargin{ii+1};
                    elseif strcmp(varName,'a31')
                        %"a31 updated"
                        flag(3) = true;
                        obj.a31 = varargin{ii+1};
                    elseif strcmp(varName,'a32')
                        %"a32 updated"
                        flag(4) = true;
                        obj.a32 = varargin{ii+1};
                    elseif strcmp(varName,'a33')
                        %"a33 updated"
                        flag(5) = true;
                        obj.a33 = varargin{ii+1};
                    elseif strcmp(varName,'a41')
                        %"a41 updated"
                        flag(6) = true;
                        obj.a41 = varargin{ii+1};
                    elseif strcmp(varName,'a42')
                        %"a42 updated"
                        flag(7) = true;
                        obj.a42 = varargin{ii+1};
                    else 
                        error('invalid input')
                    end
                end
            end
            
            %% Check if new variable values are allowed:
            a34 = 270 - (obj.a31 + obj.a32 + obj.a33);
            % Sector angles in node 3 should not exceed 260 deg
            % Sector angles in node 5 should be greated than 90 deg
%             if a34<=10 || a34>= 160    %160 deg otherwise node 4 will be very high up
%                 error("sector angles in 3 not allowed: "+num2str([obj.a31 obj.a32 obj.a33]))   
%             end
            % Sector angles in node 4 should not exceed 350 deg
            % Sector angles in node 4 should be greated than 180 deg
            a44 = 180 - a34;
            a43 = 360 - (a44 + obj.a41 + obj.a42);
            if a43 <= 20 || a43 >= 160          
                error("sector angles in 4 not allowed: "+num2str([a44 obj.a41 obj.a42]))
            end
            %Update the affected vertices (seperated for speed)
            if flag(3) %a31
                obj.vertex4
                obj.vertex5
                obj.vertex6
                obj.vertex8
                obj.vertex9 
                obj.resetAdjacentAuxiliaries([4,5,6,8,9]);
            elseif flag(4) %a32
                obj.vertex4
                obj.vertex5
                obj.vertex8
                obj.vertex9
                obj.resetAdjacentAuxiliaries([4,5,8,9]);
            elseif flag(5) || flag(1) %a33 or R1
                obj.vertex4
                obj.vertex8
                obj.vertex9
                obj.resetAdjacentAuxiliaries([4,8,9]);
            elseif flag(6) %a41
                obj.vertex8
                obj.vertex9
                obj.resetAdjacentAuxiliaries([8,9]);
            elseif flag(7) || flag(2) %a42 or R2
                obj.vertex8
                obj.resetAdjacentAuxiliaries(8);
            end    
            
            
            
            
            
        end
        function [] = vertex4(obj)
        x4 = obj.R1;
        y4 = tand(180-(obj.a31+obj.a32+obj.a33))*obj.R1;

        %prevent conflicts: make sure v4 stays above v7
        if y4 <obj.vertices(7).position0(2)+0.1
            obj.vertices(7).position0(2) = y4-1;
        end

        %move vertex
        obj.vertices(4).position0  = [x4;y4;0];
        obj.vertices(4).position   = [x4;y4;0];
%         obj.resetAdjacentAuxiliaries(4);
        end
        function [] = vertex5(obj)
        y5 = sind(180-(obj.a31+obj.a32));
        x5 = cosd(180-(obj.a31+obj.a32));
        %old:
%         y5 = 1;
%         x5 = y5/tand(180-(obj.a31+obj.a32));

%         %prevent conflicts: make sure v5 stays left of v9
%         while x5 > obj.vertices(9).position0(1)-0.1
%             %decrease crease length
%             y5 = 0.8*y5;
%             x5 = y5/tand(180-(obj.a31+obj.a32));
%         end

        %move vertex
        obj.vertices(5).position0 = [x5;y5;0];
        obj.vertices(5).position  = [x5;y5;0];
%         obj.resetAdjacentAuxiliaries(5);
        end
        function [] = vertex6(obj)
        x6 = cosd(180-obj.a31);
        y6 = sind(180-obj.a31);
        %old:
%         y6 = 1;
%         x6 = -y6/tand(obj.a31); 

        %no conflicts to be prevented
        %move vertex
        obj.vertices(6).position0 = [x6;y6;0];
        obj.vertices(6).position  = [x6;y6;0];
%         obj.resetAdjacentAuxiliaries(6);
        end
        function [] = vertex8(obj)
        x8 = obj.R1+obj.R2;

        y4   = tand(180-(obj.a31+obj.a32+obj.a33))*obj.R1;
        dy48 = tand( 180 - (obj.a41+obj.a42-(180-(obj.a31+obj.a32+obj.a33))))*obj.R2; 
        y8 = y4+dy48;
        %no conflicts to be prevented
        %move vertex
        obj.vertices(8).position0 = [x8;y8;0];
        obj.vertices(8).position  = [x8;y8;0];
%         obj.resetAdjacentAuxiliaries(8);
        end
        function [] = vertex9(obj)
        a34 = 270-(obj.a31+obj.a32+obj.a33);
        a44 = 180-a34;
        x9 = obj.vertices(4).position0(1) - sind(a44 + obj.a41);
        y9 = obj.vertices(4).position0(2) - cosd(a44 + obj.a41);
        %old
%         y9 = 1;
%         x4 = obj.R1;
%         y4 = tand(180-(obj.a31+obj.a32+obj.a33))*obj.R1;
%         dy49 = y9-y4;
% 
%         dx49 = -dy49/tand(obj.a41-(180-(obj.a31+obj.a32+obj.a33)));
%         x9 = x4+dx49;

%         %prevent conflict: v9 should stay right of v5
%         while x9 < obj.vertices(5).position0(1)+0.1
%             %decrease crease length
%             y9 = 0.8*y9;
%             x9 = x4+0.8*dx49;
%         end
        %move vertex
        obj.vertices(9).position0 = [x9;y9;0];
        obj.vertices(9).position  = [x9;y9;0];
        %obj.resetAdjacentAuxiliaries(9);
        end
        
        %Plotting
        function [] = drawRangeOfMotionGIF(obj,angleRange1,angleRange2,filename,varargin)
            %% Created by Jim (!)
            % *creates a drawing of all points the 2DOF mechanism can reach
            %  with the given range of input angles
            % *Basically runs an animation of many combinations of angles in
            %  order to show what the range of motion is 
            % *Very similar to animateOrigami
            % *angleRange can be a 1x1, 1x2 or >1x2 row vector. See below
            % for different interpretations
            % *varargin can be:
            %  # 'nSteps': number of steps in a single one way stroke
            %  # phaseShift: number 
            % (!) Implement function where phase shift can be added
            % (!) Add check for even number of varargin
            
            %Requires status doKinematics
            obj.checkState('doKinematics');
            
%             %Set number of steps per angle if not specified in input
%             %! this means the total number of steps will be (2*nStep-2)^2
%             %! 
%             if nargin<5
%                 nSteps = 16;
%             end

%             %Get varargin names
%             varNames = varargin(1:2:end)
% 
%             %get nSteps
%             if ~isempty(findstr(varNames,'nSteps'))
%                     nSteps = findstr(varNames,'nSteps');
%             else 
            
            %default variables for varargin:
            nSteps = 16;
            phaseShift = 90;

            %%Get varargins
            if numel(varargin) >0
                for ii = 1:2:numel(varargin)
                    varName = varargin{ii};
                    
                    %check if there is also a value assigned
                    if ii+1 > numel(varargin)
                        error('no value assigned to input variable')
                    end
                    %Number of steps
                    if strcmp(varName, 'nSteps')
                        "nSteps detected :)"
                        nSteps = varargin{ii+1};
                    %Phase shift
                    elseif strcmp(varName, 'nSteps')
                        "phaseShift detected :)"
                        phaseShift = varargin{ii+1};
                    else
                        error('variable name not recognized')
                    end
                end
            end
            
            
            
            %give warning if rounding error occurs when shifting phase
            if round(nSteps*phaseShift/360) ~= nSteps*phaseShift/360
                warning('Phase shift may be off due to innacurracy in rounding. HINT: make sure nSteps*phaseShift/360 is an integer')
            end
            phaseShift = round(nSteps*phaseShift/360);
            
            %%Calculate input angles for each step
            % angleRange 1x1:
            %  Input is taken as max: the path is taken from 0 to
            %  angleRange. Currently, a phase shift of 90 deg is created
            %  between both angles (angle2 lags with 90 deg)
            % angleRange 1x2:
            %  Input is taken as range. Same phaseShift
            % angelRange > 1x2:
            %  This is interpreted as custom mode. The input is now
            %  considered to represent the full 
            
            %check if size of both inputs are the same
            if size(angleRange1) ~= size(angleRange2)
                error('Make sure angleRange1 and anlgeRange2 are the same size')
            end
            
           
            if size(angleRange1,2) == 1 
                %input1
                halfCycle = linspace(1e-3*sign(angleRange1), angleRange1,nSteps);
                dihAngles1 = [halfCycle fliplr(halfCycle(2:end-1))];
                %input2
                halfCycle = linspace(1e-3*sign(angleRange2), angleRange2,nSteps);
                dihAngles2 = [halfCycle fliplr(halfCycle(2:end-1))];
                %create phaseShift:
                dihAngles2 = [dihAngles2(phaseShift+1:end) dihAngles2(1:phaseShift)];
                
            elseif size(angleRange1,2) == 2
                %input1
                halfCycle = linspace(angleRange1(1), angleRange1(2),nSteps);
                dihAngles1 = [halfCycle fliplr(halfCycle(2:end-1))];
                %input2
                halfCycle = linspace(angleRange2(1), angleRange2(2),nSteps);
                dihAngles2 = [halfCycle fliplr(halfCycle(2:end-1))];
                %create phaseShift:
                dihAngles2 = [dihAngles2(phaseShift+1:end) dihAngles2(1:phaseShift)];
                
            elseif size(angleRange1,2) >2
                dihAngles1 = angleRange1;
                dihAngles2 = angleRange2;
            else
                error('angleRange1 should be a 1x1 or 1x2 vector')
            end
            
            %make sure no trivial fold occurs
                index0 = find(dihAngles1 == 0);
                dihAngles1(index0) = 1e-3;
                index0 = find(dihAngles2 == 0);
                dihAngles2(index0) = 1e-3;                

            pauseTime = 0.1;
           
            EEcoords = [];
            for j1 = 1:nSteps
            for j2 = 1:nSteps

                %%Create the deformed origami
                obj.setDihedralAngles([dihAngles1(j1), dihAngles2(j2)]);
                obj.determineOrigami();
                %get coordinates of end effector & add to set of all coordinates
                EEcoords = [EEcoords [obj.vertices(obj.vertexToGripWith).position;dihAngles1(j1);dihAngles2(j2)] ];
                    %%EEcoords is a matrix with:
                    % every column is an end effector position
                    % first three rows are its X,Y,Z coords
                    % 4th and 5th rows are the inputs angles for which this
                    % position occurs
                
                
                %plot origami in the current axis
                if j1 == 1 && j2 == 1
                    axAni = obj.plotOrigamiDeformed();
                else
                    axAni = obj.plotOrigamiDeformed(axAni);
                end
                hold(axAni,'on')
                %%plot all coordinates of end effector
                EEgrid = plot3(EEcoords(1,:),EEcoords(2,:),EEcoords(3,:),'o','Color','b','MarkerSize',4,...
    'MarkerFaceColor','#D9FFFF');
                %%add data tips
                dtRows = [dataTipTextRow("\rho_1",EEcoords(4,:))
                          dataTipTextRow("\rho_2",EEcoords(5,:))];
                EEgrid.DataTipTemplate.DataTipRows(end+1:end+2) = dtRows;

                %plot z-surface
                fmesh(1,'FaceAlpha',0.1)
                s.EdgeColor = 'none';
                obj.resetToUndetermined();
%                 hold(ax,'off')
                
%             
                %get movie frame
                frame = getframe(gcf);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);
                
                %Write to the GIF File
                if j1 == 1 &j2 == 1
                    testcounter=1
                    imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
                else
                    testcounter = testcounter+1
                    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime', pauseTime);
                end
            end
            end
        end
        function [] = drawRangeOfMotion(obj,angleRange1,angleRange2,filename,varargin)
            %% Created by Jim (!)
            % *creates a drawing of all points the 2DOF mechanism can reach
            %  with the given range of input angles
            % *Basically runs an animation of many combinations of angles in
            %  order to show what the range of motion is 
            % *Very similar to animateOrigami
            % *angleRange can be a 1x1, 1x2 or >1x2 row vector. See below
            % for different interpretations
            % *varargin can be:
            %  # 'nSteps': number of steps in a single one way stroke
            %  # phaseShift: number 
            % (!) Implement function where phase shift can be added
            % (!) Add check for even number of varargin
            
            %Requires status doKinematics
            obj.checkState('doKinematics');
            
%             %Set number of steps per angle if not specified in input
%             %! this means the total number of steps will be (2*nStep-2)^2
%             %! 
%             if nargin<5
%                 nSteps = 16;
%             end

%             %Get varargin names
%             varNames = varargin(1:2:end)
% 
%             %get nSteps
%             if ~isempty(findstr(varNames,'nSteps'))
%                     nSteps = findstr(varNames,'nSteps');
%             else 
            
            %default variables for varargin:
            nSteps = 16;
            phaseShift = 90;

            %%Get varargins
            if numel(varargin) >0
                for ii = 1:2:numel(varargin)
                    varName = varargin{ii}
                    
                    %check if there is also a value assigned
                    if ii+1 > numel(varargin)
                        error('no value assigned to input variable')
                    end
                    %Number of steps
                    if strcmp(varName, 'nSteps')
                        "nSteps detected :)"
                        nSteps = varargin{ii+1};
                    %Phase shift
                    elseif strcmp(varName, 'nSteps')
                        "phaseShift detected :)"
                        phaseShift = varargin{ii+1};
                    else
                        error('variable name not recognized')
                    end
                end
            end
            
            
            
            %give warning if rounding error occurs when shifting phase
            if round(nSteps*phaseShift/360) ~= nSteps*phaseShift/360
                warning('Phase shift may be off due to innacurracy in rounding. HINT: make sure nSteps*phaseShift/360 is an integer')
            end
            phaseShift = round(nSteps*phaseShift/360);
            
            %%Calculate input angles for each step
            % angleRange 1x1:
            %  Input is taken as max: the path is taken from 0 to
            %  angleRange. Currently, a phase shift of 90 deg is created
            %  between both angles (angle2 lags with 90 deg)
            % angleRange 1x2:
            %  Input is taken as range. Same phaseShift
            % angelRange > 1x2:
            %  This is interpreted as custom mode. The input is now
            %  considered to represent the full 
            
            %check if size of both inputs are the same
            if size(angleRange1) ~= size(angleRange2)
                error('Make sure angleRange1 and anlgeRange2 are the same size')
            end
            
           
            if size(angleRange1,2) == 1 
                %input1
                halfCycle = linspace(1e-3*sign(angleRange1), angleRange1,nSteps);
                dihAngles1 = [halfCycle fliplr(halfCycle(2:end-1))];
                %input2
                halfCycle = linspace(1e-3*sign(angleRange2), angleRange2,nSteps);
                dihAngles2 = [halfCycle fliplr(halfCycle(2:end-1))];
                %create phaseShift:
                dihAngles2 = [dihAngles2(phaseShift+1:end) dihAngles2(1:phaseShift)];
                
            elseif size(angleRange1,2) == 2
                %input1
                halfCycle = linspace(angleRange1(1), angleRange1(2),nSteps);
                dihAngles1 = [halfCycle fliplr(halfCycle(2:end-1))];
                %input2
                halfCycle = linspace(angleRange2(1), angleRange2(2),nSteps);
                dihAngles2 = [halfCycle fliplr(halfCycle(2:end-1))];
                %create phaseShift:
                dihAngles2 = [dihAngles2(phaseShift+1:end) dihAngles2(1:phaseShift)];
                
            elseif size(angleRange1,2) >2
                dihAngles1 = angleRange1;
                dihAngles2 = angleRange2;
            else
                error('angleRange1 should be a 1x1 or 1x2 vector')
            end
            
            %make sure no trivial fold occurs
                index0 = find(dihAngles1 == 0);
                dihAngles1(index0) = 1e-3;
                index0 = find(dihAngles2 == 0);
                dihAngles2(index0) = 1e-3;                

            pauseTime = 0.1;
           
            %% Calculate all coordinates
            
            %for all grippingCandidates

                
                
            %create empty containers
            for k = 1:length(obj.grippingCandidates)
                EEcoords{k} = [];
            end
                        
            for j1 = 1:nSteps
            for j2 = 1:nSteps

                %%Create the deformed origami
                obj.setDihedralAngles([dihAngles1(j1), dihAngles2(j2)]);
                obj.determineOrigami();
                for k = 1:length(obj.grippingCandidates)
                    EE = obj.grippingCandidates(k);
                    %get coordinates of end effector & add to set of all coordinates
                    EEcoords{k} = [EEcoords{k} [obj.vertices(EE).position;dihAngles1(j1);dihAngles2(j2)] ]
                        %%EEcoords is a matrix with:
                        % every column is an end effector position
                        % first three rows are its X,Y,Z coords
                        % 4th and 5th rows are the inputs angles for which this
                        % position occurs
                end
            obj.resetToUndetermined();
            end
            end

            
            
            % plot deformed origami for visual purposes
            obj.setDihedralAngles([15, 45]);
            obj.determineOrigami();
            axAni = obj.plotOrigamiDeformed();

            hold(axAni,'on')
            
            %%plot all coordinates of end effector
            for k = 1:length(obj.grippingCandidates)
                colors = rand(1,3);
                EEgrid{k} = plot3(EEcoords{k}(1,:),EEcoords{k}(2,:),EEcoords{k}(3,:),'o','Color','k','MarkerSize',4,...
                'MarkerFaceColor',rand(1,3));
                %%add data tips
                dtRows = [dataTipTextRow("\rho_1",EEcoords{k}(4,:))
                          dataTipTextRow("\rho_2",EEcoords{k}(5,:))];
                EEgrid{k}.DataTipTemplate.DataTipRows(end+1:end+2) = dtRows;
            end
            %plot z-surface
            fmesh(1,'FaceAlpha',0.1)
            s.EdgeColor = 'none';

            
        end
        function [] = drawEEpathWithOrigami(obj,angleRange1,angleRange2,varargin)
            %(!)JIM
%             default variables for varargin:
            nSteps = 128;
            phaseShift = 90;

            %%Get varargins
            if numel(varargin) >0
                for ii = 1:2:numel(varargin)
                    varName = varargin{ii}
                    
                    %check if there is also a value assigned
                    if ii+1 > numel(varargin)
                        error('no value assigned to input variable')
                    end
                    %Number of steps
                    if strcmp(varName, 'nSteps')
                        "nSteps detected :)"
                        nSteps = varargin{ii+1};
                    %Phase shift
                    elseif strcmp(varName, 'nSteps')
                        "phaseShift detected :)"
                        phaseShift = varargin{ii+1};
                    else
                        error('variable name not recognized')
                    end
                end
            end
            
            
            inputAngles = generateInputProfile(angleRange1,nSteps);
            dihAngles1 = inputAngles(1,:);
            dihAngles2 = inputAngles(2,:);

            
            %% Calculate all EE coordinates
            
            %for all grippingCandidates
            
            %create empty containers
            for k = 1%:length(obj.grippingCandidates)
                EEcoords{k} = [];
            end
            
            for j = 1:size(dihAngles1,2)
                obj.setDihedralAngles([dihAngles1(j), dihAngles2(j)]);
                obj.determineOrigami();
                
                for k = 1%:length(obj.grippingCandidates)
                    EE = 8;%obj.grippingCandidates(k);
                    EEcoords{k} = [EEcoords{k} [obj.vertices(EE).position;dihAngles1(j);dihAngles2(j)] ];
                    %%EEcoords{k} is a matrix with:
                    % every column is an end effector position
                    % first three rows are its X,Y,Z coords
                    % 4th and 5th rows are the inputs angles for which this
                    % position occurs
                end
               obj.resetToUndetermined();
            end
            
            %% Set view 
            allvertexcoords = cat(2,obj.vertices(:).position0);
            obj.animationXLims = [min(allvertexcoords(1,:))-0.3 max(allvertexcoords(1,:))+0.3];
            obj.animationYLims = [min(allvertexcoords(2,:))-0.3 max(allvertexcoords(2,:))+0.3]; 
            obj.animationZLims = [min([EEcoords{k}(3,:) 0])-0.3 max([EEcoords{k}(3,:) 0])+0.3];
            
                
            %% Plot deformed origami for visual purposes
            obj.setDihedralAngles([dihAngles1(ceil(end/2)), dihAngles2(ceil(end/2))]);
                obj.determineOrigami();
            ax = obj.plotOrigamiDeformed();
            hold(ax,'on')
            
            %% plot coordinates
            for k = 1:length(obj.grippingCandidates)
                colors = rand(1,3);
                EEpath = plot3(EEcoords{k}(1,:),EEcoords{k}(2,:),EEcoords{k}(3,:),'-o','Color','k','MarkerSize',4,...
                'MarkerFaceColor',colors);
                %%add data tips
                dtRows = [dataTipTextRow("\rho_1",EEcoords{k}(4,:))
                          dataTipTextRow("\rho_2",EEcoords{k}(5,:))];
                EEpath.DataTipTemplate.DataTipRows(end+1:end+2) = dtRows;
            end
        end
        function [] = drawEEpath(obj,angleRange1,angleRange2,varargin)
            %(!)JIM
            %default variables for varargin:
            nSteps = 32;
            phaseShift = 90;

            %%Get varargins
            if numel(varargin) >0
                for ii = 1:2:numel(varargin)
                    varName = varargin{ii}
                    
                    %check if there is also a value assigned
                    if ii+1 > numel(varargin)
                        error('no value assigned to input variable')
                    end
                    %Number of steps
                    if strcmp(varName, 'nSteps')
                        "nSteps detected :)"
                        nSteps = varargin{ii+1};
                    %Phase shift
                    elseif strcmp(varName, 'nSteps')
                        "phaseShift detected :)"
                        phaseShift = varargin{ii+1};
                    else
                        error('variable name not recognized')
                    end
                end
            end
            
            %give warning if rounding error occurs when shifting phase
            if round(nSteps*phaseShift/360) ~= nSteps*phaseShift/360
                warning('Phase shift may be off due to innacurracy in rounding. HINT: make sure nSteps*phaseShift/360 is an integer')
            end
            phaseShift = round(nSteps*phaseShift/360);
            
            %%Calculate input angles for each step
            % angleRange 1x1:
            %  Input is taken as max: the path is taken from 0 to
            %  angleRange. Currently, a phase shift of 90 deg is created
            %  between both angles (angle2 lags with 90 deg)
            % angleRange 1x2:
            %  Input is taken as range. Same phaseShift
            % angelRange > 1x2:
            %  This is interpreted as custom mode. The input is now
            %  considered to represent the full 
            
            %check if size of both inputs are the same
            if size(angleRange1) ~= size(angleRange2)
                error('Make sure angleRange1 and anlgeRange2 are the same size')
            end
            
           
            if size(angleRange1,2) == 1 
                %input1
                halfCycle = linspace(1e-3*sign(angleRange1), angleRange1,nSteps);
                dihAngles1 = [halfCycle fliplr(halfCycle(2:end-1))];
                %input2
                halfCycle = linspace(1e-3*sign(angleRange2), angleRange2,nSteps);
                dihAngles2 = [halfCycle fliplr(halfCycle(2:end-1))];
                %create phaseShift:
                dihAngles2 = [dihAngles2(phaseShift+1:end) dihAngles2(1:phaseShift)];
                
            elseif size(angleRange1,2) == 2
                %input1
                halfCycle = linspace(angleRange1(1), angleRange1(2),nSteps);
                dihAngles1 = [halfCycle fliplr(halfCycle(2:end-1))];
                %input2
                halfCycle = linspace(angleRange2(1), angleRange2(2),nSteps);
                dihAngles2 = [halfCycle fliplr(halfCycle(2:end-1))];
                %create phaseShift:
                dihAngles2 = [dihAngles2(phaseShift+1:end) dihAngles2(1:phaseShift)];
                
            elseif size(angleRange1,2) >2
                dihAngles1 = angleRange1;
                dihAngles2 = angleRange2;
            else
                error('angleRange1 should be a 1x1 or 1x2 vector')
            end
            
            %% Calculate all EE coordinates
            EEcoords = [];
            for j = 1:size(dihAngles1,2)
                obj.setDihedralAngles([dihAngles1(j), dihAngles2(j)]);
                obj.determineOrigami();
                EEcoords = [EEcoords [obj.vertices(obj.vertexToGripWith).position;dihAngles1(j);dihAngles2(j)] ];
                    %%EEcoords is a matrix with:
                    % every column is an end effector position
                    % first three rows are its X,Y,Z coords
                    % 4th and 5th rows are the inputs angles for which this
                    % position occurs
               obj.resetToUndetermined();
            end
            
            %% Set view 
            allvertexcoords = cat(2,obj.vertices(:).position0);
            obj.animationXLims = [min(allvertexcoords(1,:))-0.3 max(allvertexcoords(1,:))+0.3];
            obj.animationYLims = [min(allvertexcoords(2,:))-0.3 max(allvertexcoords(2,:))+0.3]; 
            obj.animationZLims = [min([EEcoords(3,:) 0])-0.3 max([EEcoords(3,:) 0])+0.3];
            
            
            %% plot coordinates
            EEpath = plot3(EEcoords(1,:),EEcoords(2,:),EEcoords(3,:),'-o','Color','b','MarkerSize',4,...
            'MarkerFaceColor','#D9FFFF');
            %%add data tips
            dtRows = [dataTipTextRow("\rho_1",EEcoords(4,:))
                      dataTipTextRow("\rho_2",EEcoords(5,:))];
            EEpath.DataTipTemplate.DataTipRows(end+1:end+2) = dtRows;
         
        end        
        function [] = animateEEpath(obj,angleRange1,angleRange2,filename,xlimit,ylimit,zlimit,varargin)
            %% Created by Jim (!)
            % *creates a drawing of all points the 2DOF mechanism can reach
            %  with the given range of input angles
            % *Basically runs an animation of many combinations of angles in
            %  order to show what the range of motion is 
            % *Very similar to animateOrigami
            % *angleRange can be a 1x1, 1x2 or >1x2 row vector. See below
            % for different interpretations
            % *varargin can be:
            %  # 'nSteps': number of steps in a single one way stroke
            %  # phaseShift: number 
            % (!) Implement function where phase shift can be added
            % (!) Add check for even number of varargin
            
            %Requires status doKinematics
            obj.checkState('doKinematics');
            
%             %Set number of steps per angle if not specified in input
%             %! this means the total number of steps will be (2*nStep-2)^2
%             %! 
%             if nargin<5
%                 nSteps = 16;
%             end

%             %Get varargin names
%             varNames = varargin(1:2:end)
% 
%             %get nSteps
%             if ~isempty(findstr(varNames,'nSteps'))
%                     nSteps = findstr(varNames,'nSteps');
%             else 
            
            %default variables for varargin:
            nSteps = 64;
            phaseShift = 90;

            %%Get varargins
            if numel(varargin) >0
                for ii = 1:2:numel(varargin)
                    varName = varargin{ii}
                    
                    %check if there is also a value assigned
                    if ii+1 > numel(varargin)
                        error('no value assigned to input variable')
                    end
                    %Number of steps
                    if strcmp(varName, 'nSteps')
                        "nSteps detected :)"
                        nSteps = varargin{ii+1};
                    %Phase shift
                    elseif strcmp(varName, 'phaseShift')
                        "phaseShift detected :)"
                        phaseShift = varargin{ii+1};
                    else
                        error('variable name not recognized')
                    end
                end
            end
            
            
            
            %give warning if rounding error occurs when shifting phase
            if round(nSteps*phaseShift/360) ~= nSteps*phaseShift/360
                warning('Phase shift may be off due to innacurracy in rounding. HINT: make sure nSteps*phaseShift/360 is an integer')
            end
            phaseShift = round(nSteps*phaseShift/360);
            
            %%Calculate input angles for each step
            % angleRange 1x1:
            %  Input is taken as max: the path is taken from 0 to
            %  angleRange. Currently, a phase shift of 90 deg is created
            %  between both angles (angle2 lags with 90 deg)
            % angleRange 1x2:
            %  Input is taken as range. Same phaseShift
            % angelRange > 1x2:
            %  This is interpreted as custom mode. The input is now
            %  considered to represent the full 
            
            %check if size of both inputs are the same
            if size(angleRange1) ~= size(angleRange2)
                error('Make sure angleRange1 and anlgeRange2 are the same size')
            end
            
            inputAngles = generateInputProfile(angleRange1,nSteps);
            dihAngles1 = inputAngles(1,:);
            dihAngles2 = inputAngles(2,:);

            pauseTime = 0.1;
           
            EEcoords = [];

            for j = 1:size(dihAngles1,2)

                %%Create the deformed origami
                obj.setDihedralAngles([dihAngles1(j), dihAngles2(j)]);
                obj.determineOrigami();
                %get coordinates of end effector & add to set of all coordinates
                EEcoords = [EEcoords [obj.vertices(obj.vertexToGripWith).position;dihAngles1(j);dihAngles2(j)] ];
                    %%EEcoords is a matrix with:
                    % every column is an end effector position
                    % first three rows are its X,Y,Z coords
                    % 4th and 5th rows are the inputs angles for which this
                    % position occurs
                
                
                %plot origami in the current axis
                if j == 1 
                    axAni = obj.plotOrigamiDeformed();
                    f = gcf;
                    f.Position = [0 0 1920 1080];
                else
                    axAni = obj.plotOrigamiDeformed(axAni);
                end
                hold(axAni,'on')
                view(obj.animationView);
                xlim(xlimit)
                ylim(ylimit)
                zlim(zlimit)
                
                
                %%plot all coordinates of end effector
                EEgrid = plot3(EEcoords(1,:),EEcoords(2,:),EEcoords(3,:),'-o','Color','b','MarkerSize',4,...
                'MarkerFaceColor','#D9FFFF');
                %%add data tips
                dtRows = [dataTipTextRow("\rho_1",EEcoords(4,:))
                          dataTipTextRow("\rho_2",EEcoords(5,:))];
                EEgrid.DataTipTemplate.DataTipRows(end+1:end+2) = dtRows;

                %plot z-surface
                s.EdgeColor = 'none';
                obj.resetToUndetermined();
%                 hold(ax,'off')
                
%             
                %get movie frame
                frame = getframe(gcf);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);
                
                %Write to the GIF File
                if j == 1 &j == 1
                    
                    imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
                else
                    
                    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime', pauseTime);
                end
            end    
        end     
                
        function [] = plotOrigamiUndeformed3D(obj, ax, plotEdges)
            %can be called at any state but editGraph
            if strcmp(obj.status,'editGraph')
                obj.checkState('editFaces');
            end
            
            if size(obj.faces,2)==0
                error('No faces defined');
            end
            
            edgeCol=[0 0 0];
            
            if nargin==1
                ax2=figure;
                ax2=axes('Parent',ax2);
                hold(ax2,'on');
            elseif nargin >= 2
                ax2=ax;
                hold(ax2,'on');
                if nargin==3
                    if not(plotEdges)
                        edgeCol='none';
                    end
                end
            end
            
            %%%%%%%%%
            %iterate through all vertices and extract all positions
            %(position0)
            
            positions=zeros(obj.numberOfVertices(),3);
            for j=1:obj.numberOfVertices()
                positions(j,:)=obj.vertices(j).position0';
            end
            %generate random colors between 0.8 and 1 for the faces
            rng('default');
            col=rand([obj.numberOfFaces(),1])+0.8;
            
            
            
            
        end
        
        
            
%%%%%%%%TOPOLOGY SECTION%%%%
        function [] = addVertex(obj, pos, origin, mode)
            %ADDVERTEX(OBJ, POS, ORIGIN, MODE) adds vertex to origami
            %object
            %   INPUT:
            %   pos: 2-by-1 vector with position in undeformed
            %   configuration.
            %   origin: index of "parent" vertex to the new vertex
            %   mode (optional): kinematic mode of new vertex. true for up,
            %   false for down. false is default.
            
            %Requires state editGraph
            obj.checkState('editGraph');                                   %see checkState
            %check if there is already a vertex at this position
            obj.checkIsPositionOccupied(pos);            
            
            %If mode is not given, set mode to default: false
            if nargin<3
                mode=false;
            end
            
            %%Adding the new vertex to the origami object.
            index=obj.numberOfVertices()+1;
                                %automatically add next 
                                                    %index available
            zAxis=obj.vertices(origin).unfoldedZ;   %set zAxis same as 
                                                    %zAxis of origin vertex
            obj.vertices=[obj.vertices vertex(pos,...
                index, mode, 'boundary',zAxis)];    %add vertex to object
                                                    %with constructor of
                                                    %vertex.m class
                                                    %(!) regionOfOrigin is
                                                    %omitted

            %%Adding the new edge
            edgeIndex=obj.numberOfEdges()+1;        
            obj.edges=[obj.edges edge(origin,index,edgeIndex)];
            
            %Updating some parameters...
            obj.addOutgoingCrease(origin,edgeIndex);                       %see addOutgoingCrease
            obj.addIncomingCrease(index,edgeIndex);                        %see addIncomingCrease
            
        end %function addVertex      
        function [] = removeVertex(obj, vertexIndex)
            %Requires state editGraph
            obj.checkState('editGraph');
            %Check if the vertex has outputs
            if size(obj.vertices(vertexIndex).outgoingEdges,2)~=0
                error('Trying to delete a vertex with one or several outputs. This is not allowed.');
            end
            %Go through incoming edges and delete them
            while size(obj.vertices(vertexIndex).incomingEdges,2)>0
                edgeIndex=obj.vertices(vertexIndex).incomingEdges(1); 
                %delete the crease
                obj.removeCreaseByIndex(edgeIndex);
            end
            %Assert that no connections around the vertex remain
            if not(size(obj.vertices(vertexIndex).incomingEdges,2) == 0 && ...
                    size(obj.vertices(vertexIndex).outgoingEdges,2) == 0)
                error('Standard procedures were not able to completely disconnect the vertex from its surroundings');
            end
            %delete the vertex
            obj.vertices(vertexIndex).deleted=true;
            obj.vertices(vertexIndex).type = 'deleted'; %ADDED
            %update T
            obj.vertices(vertexIndex).T = false;
            
        end %function removeVertex
        
        function [] = addCrease(obj, origin, endIndex)
            %ADDCREASE Adds a crease to an origami object
            %   Use this function to connect two existing vertices
            
            obj.checkState('editGraph');

            %Check if the two vertices are already connected
            if ismember(endIndex,obj.vertices(origin).outgoingNeighbours)
                error('These two vertices are already connected! Check vertices nr. %d and %d.',origin,endIndex);
            end
            
            %First, add another crease to the origami object
            index=obj.numberOfEdges()+1;
            obj.edges=[obj.edges edge(origin,endIndex,index)];
            
            %Now link the new edge to the existing vertices
            obj.addOutgoingCrease(origin,index);
            obj.addIncomingCrease(endIndex,index);
        end
        
        function [] = removeCreaseByIndex(obj,edgeIndex)
            %REMOVECREASEBYINDEX Removes a crease given its index. Returns
            %an error if the crease can't be found or is already deleted.
            
            %Requires state editGraph if the crease is internal. If the
            %crease is on the boundary, editFaces and editGraph are both
            %permitted. It is required that the deletion of a boundary in
            %the state editGraph works in order for the function
            %goToEditGraph to work.
            
            if strcmp(obj.edges(edgeIndex).type,'internal')
                obj.checkState('editGraph');
            elseif strcmp(obj.edges(edgeIndex).type,'boundary')
                if strcmp(obj.status,'editFaces')
                    obj.checkState('editFaces');
                elseif strcmp(obj.status,'editGraph')
                    obj.checkState('editGraph');
                else
                    obj.checkState('editFaces'); % throws error
                end
            end
            
            %Check if the calling origami even contains an edge with the
            %given index
            if obj.numberOfEdges<edgeIndex
                error('Edge Index not found among origami object');
            end
            
            %Check whether the edge was already marked as deleted (can't
            %delete objects twice)
            if obj.edges(edgeIndex).deleted
                error('Edge was already deleted!');
            end
            
            %Now mark the edge as deleted
            obj.edges(edgeIndex).deleted=true;
            
            %And adapt the connected vertices
            obj.removeOutgoingCrease(obj.edges(edgeIndex).sourceVertex,...
                edgeIndex);
            obj.removeIncomingCrease(obj.edges(edgeIndex).endVertex,...
                edgeIndex);
            
        end
        
        %%Counting functions
        function numVertices = numberOfVertices(obj, opts)
            %NUMBEROFVERTICES returns the number of vertices in the origami
            %object. If the second argument is 'omitDeleted', only the
            %number of non-deleted vertices is returned.
            
            if nargin==1
                numVertices=size(obj.vertices,2);
            elseif strcmp(opts,'omitDeleted')
                counter=0;
                for j=1:obj.numberOfVertices()
                    if not(obj.vertices(j).deleted)
                        counter=counter+1;
                    end
                end
                numVertices=counter;
            elseif strcmp(opts,'interior') %(!) added internal vertices counter
                counter=0;
                for j=1:obj.numberOfVertices()
                    if strcmp(obj.vertices(j).type,'interior')
                        counter=counter+1;
                    end
                end
                numVertices = counter;
            elseif strcmp(opts,'source') %(!) added internal vertices counter
                counter=0;
                for j=1:obj.numberOfVertices()
                    if strcmp(obj.vertices(j).type,'source')
                        counter=counter+1;
                    end
                end
                numVertices = counter;
            elseif strcmp(opts,'aux') %(!) added internal vertices counter
                counter=0;
                for j=1:obj.numberOfVertices()
                    if strcmp(obj.vertices(j).type,'aux')
                        counter=counter+1;
                    end
                end
                numVertices = counter;
                
            end
        end
        function numEdges = numberOfEdges(obj,opts)
            %NUMBEROFEDGES returns the number of edges/creases in the
            %origami object.
            
            if nargin==1
                numEdges=size(obj.edges,2);
            elseif strcmp(opts,'omitDeleted')
                %In that case, only count those edges that aren't marked as
                %deleted.
                counter=0;
                for j=1:size(obj.edges,2)
                    if not(obj.edges(j).deleted)
                        counter=counter+1;
                    end
                end
                numEdges=counter;
            end
        end
        function numFaces = numberOfFaces(obj)
            %NUMBEROFFACES returns the number of faces in the origami
            %object.
            numFaces=size(obj.faces,2);
        end
    
%%%%%%%%FACE CREATION%%%%
        %(unreviewed)
        function [] = createFaces(obj)
            %CREATEFACES(OBJ) creates all faces in the origmi (and links them to
            %the vertices and edges, of course)
            %Also creates default-faces for all subOrigamis.
            
            %It is possible that this function is called when the faces
            %have already been created. In that case, just go to the part
            %where the subOrigamis are called.
            if not(strcmp(obj.status,'editFaces'))
                
                %Can be executed from state editGraph. Changes origami object
                %state to editFaces.
                obj.checkState('editGraph');
                
                %Loop over vertices
                nOV=obj.numberOfVertices();
                for j=1:nOV
                    j
                    
                    %Check if neither deleted nor auxiliary - in that case, skip
                    if obj.vertices(j).deleted || strcmp(obj.vertices(j).type,'aux') %|| strcmp(obj.vertices(j).type,'source')
                        continue;
                    end
                    
                    %Loop over incoming edges
                    for k=1:size(obj.vertices(j).incomingEdges,2)
                        
                        edgeIndex=obj.vertices(j).incomingEdges(k);
                        
                        %Double-check to ensure they aren't deleted - in that case,
                        %throw error
                        if obj.edges(edgeIndex).deleted
                            error('Deleted edge was found among connected edges. This is not allowed. Did you remove the deleted edges correctly?');
                        end
                        "current |node: "+num2str(j)+"|edge:"+num2str(edgeIndex)
                        %if there is no face to this edge's right, define a face there
                        
                        if obj.edges(edgeIndex).faceRight==0
                            "No Face right of edge "+num2str(k)+" going to node "+num2str(j)
                            obj.defineFaceFromVertexAndEdge(j,edgeIndex);
                        else
                            "edge: "+num2str(edgeIndex)+"|faceRight: "+num2str(obj.edges(edgeIndex).faceRight)
                        end
                        
                    end %loop over incoming edges
                    
                    %Loop over outgoing edges
                    for k=1:size(obj.vertices(j).outgoingEdges,2)
                        
                        edgeIndex=obj.vertices(j).outgoingEdges(k);
                        
                        %if it's already marked as boundary edge: skip it
                        if strcmp(obj.edges(edgeIndex).type,'boundary')
                            continue;
                        end
                        "current |node: "+num2str(j)+"|edge:"+num2str(edgeIndex)
                        %Double-check to ensure they aren't deleted - in that case,
                        %throw error
                        if obj.edges(edgeIndex).deleted
                            error('Deleted edge was found among connected edges. This is not allowed. Did you remove the deleted edges correctly?');
                        end
                        
                        %if there is no face to this edge's left, define a face there
                        if obj.edges(edgeIndex).faceLeft==0
                            "No FaceLeft found yet"
                            obj.defineFaceFromVertexAndEdge(j,edgeIndex);
                            
                        else
                            "edge: "+num2str(edgeIndex)+"|faceLeft: "+num2str(obj.edges(edgeIndex).faceLeft)
                        end
                        
                    end %loop over outgoing edges
                    
                end % loop over vertices
                
                %find gripping candidates (added by Jim (!))
                obj.findGrippingCandidates();
%                 
%                 
%                 %JIM (!):
%                 %Add an extra boundary crease between the two source nodes
%                 numSource = obj.numberOfVertices('source');
%                 if numSource==2
%                     obj.addAuxiliaryVertex(obj.source1, obj.source2)
%                 end
                
                
                %Changing status to editFaces
                obj.status='editFaces';
                

                
            end %if-condition "only create faces if they haven't been created yet"
          
        end
        function [] = goToEditGraph(obj)
            %GOTOEDITGRAPH(OBJ) Resets the origami (and all associated
            %containers) to the state editGraph.
            
            %Requires any state but editGraph
            if strcmp(obj.status,'editGraph')
                obj.checkState('editFaces'); %throws error
            end
            
            %If origami object is currently in the state doKinematics,
            %first reset to editFaces
            if strcmp(obj.status,'doKinematics')
                obj.goToEditFaces();
            end
            
            %update status. Needs to happen here to allow for the
            %operations that come later on (removing vertices etc.)
            obj.status='editGraph';
            
            %Empty the faces container
            obj.faces=face.empty(1,0);
            
            %get all vertices to delete
            toDelete=obj.verticesOnBoundary();
            
            while not(all([obj.vertices(toDelete).deleted]))
                for j=toDelete
                    
                    
                    if obj.vertices(j).deleted
                        %We don't care about deleted vertices
                        continue;
                    end
                    
                    
                    %check if it can be deleted
                    if isempty(obj.vertices(j).outgoingEdges)
                        %Delete the vertex
                        obj.removeVertex(j);
                        %Delete the adjacent faces vector
                        obj.vertices(j).adjacentFaces(:)=0;
                        %NOTE BY JIM: we only delete the adjacentFaces of
                        %the boundary creases. This becomes a problem later
                        %with copying
                    end
                    
                end
            end
            
            for j=1:size(obj.edges,2)
                if obj.edges(j).deleted
                    %We don't care about deleted edges
                    continue;
                end
                
                %Delete the face to the left and to the right of the crease
                obj.edges(j).faceLeft=0;
                obj.edges(j).faceRight=0;
                
            end
            
            %clear gripping candidates
            obj.grippingCandidates = [];
        end
        function [] = goToEditFaces(obj)
            %GOTOEDITFACES(OBJ) Resets the containers of the origami object s.t.
            %the status editFaces is restored.
            
            %Requires status doKinematics
            obj.checkState('doKinematics');
            
            %Clear all sector angle assignments to the vertices
            for j=1:size(obj.vertices,2)
                
                if obj.vertices(j).deleted
                    %We don't care about deleted vertices
                    continue;
                end
                
                %set all angles to 0
                obj.vertices(j).sectorAngles(:)=0;
                
                %reset all positions to position0
                obj.vertices(j).position=obj.vertices(j).position0;
                
                %reset all determinacies and determinabilities
                obj.vertices(j).determinable=false;
                obj.vertices(j).determined=false;
                
            end
            
            %Clear all length assignments to the edges
            for j=1:size(obj.edges,2)
                
                if obj.edges(j).deleted
                    %We don't care about deleted edges
                    continue;
                end
                
                obj.edges(j).length=0;
                
                %Reset the direction vector
                obj.edges(j).directionVector.vec=[0;0;0];
                
                %Reset the dihedral angle unless the source is the source
                %vertex of the origami
                if not(obj.edges(j).sourceVertex==obj.source1 || obj.edges(j).sourceVertex==obj.source2)
                    obj.edges(j).dihedralAngle.ang=0;
                end
                
            end
            
            %Now the faces
            for j=1:size(obj.faces,2)
                %Set the faces to 'not determined'
                obj.faces(j).determined=false;
               
                %Reset the normal vectors to the faces
                obj.faces(j).normalVector.vec(:)=0;
                
            end
            
            %change the state of the origami object
            obj.status='editFaces';
            
            
        end
        
        function [] = fixFace(obj, faceIndex)
            %FIXFACE(OBJ, FACEINDEX) fixes the face with the specified
            %index and unfixes all other faces.
            
            %requires state editFaces
            obj.checkState('editFaces');
            
            %Sanity check on the specified faceIndex
            if faceIndex<1 || faceIndex>obj.numberOfFaces()
                error('Invalid Face Index');
            end
            
            %Mark the specified face as fixed
            obj.faces(faceIndex).fixed=true;
            
            %Mark the currently fixed face as not fixed
            if faceIndex~=obj.fixedFaceIndex && obj.fixedFaceIndex~=0
                %The fixed face is being changed
                obj.faces(obj.fixedFaceIndex).fixed=false;
            end
            
            %Update the index of the fixed face of the current origami
            %object
            obj.fixedFaceIndex=faceIndex;
            
        end 
        
%%%%%%%%EVALUATION%%%%
        function [] = setDihedralAnglesOld(obj, angle)
            %% OLD VERSION ONLY WORKS FOR A SINGLE SOURCE VERTEX
            %SETDIHEDRALANGLES(OBJ, VERTEXINDEX, ANGLE) Sets the dihedral 
            %angle or, if several internal edges emanate from the source
            %vertex, several dihedral angles.
            %All angles in degrees.
            %Can be executed in any state
            
            %Check if enough angles are given
            if size(angle,2)~=obj.vertices(obj.source1).degree
                error('Number of specified dihedral angles does not match with degree of source vertex');
            end
            
            counter=1;
            %j is a vector spanning all outgoing edges of the source vertex
            for j=1:size(obj.vertices(obj.source1).outgoingEdges,2)
                edgeIndex=obj.vertices(obj.source1).outgoingEdges(j); %appoint edge
                if strcmp(obj.edges(edgeIndex).type,'internal')                 %if the edge is internal
                    obj.edges(edgeIndex).dihedralAngle.ang=angle(counter);      %set angle of edge
                end
            end
            
            obj.resetToUndetermined();  
        end %setDihedralAngles
        function [] = setDihedralAngles(obj,angle)
            %% NEW VERSION BY JIM (!)
            % sets the two dihedral angles of double input mechanism
            
            %check if amount of angles given is same as amount of crease
            %emanating from source vertices
            if size(angle,2) ~= obj.vertices(obj.source1).degree + obj.vertices(obj.source2).degree
            error('Number of dihedral angles given does not match amount of input creases')
            end
            
            counter=1;
            %j1 is a vector spanning all outgoing edges of the 1st source
            %index
            for j1=1:size(obj.vertices(obj.source1).outgoingEdges,2)
                edgeIndex=obj.vertices(obj.source1).outgoingEdges(j1); %appoint edge
                if strcmp(obj.edges(edgeIndex).type,'internal')                 %if the edge is internal
                    obj.edges(edgeIndex).dihedralAngle.ang=angle(counter);      %set angle of edge
                    counter = counter+1;
                end
            end
            
            %j2 is a vector spanning all outgoing edges of the 1st source
            %index
            for j2=1:size(obj.vertices(obj.source2).outgoingEdges,2)
                edgeIndex=obj.vertices(obj.source2).outgoingEdges(j2); %appoint edge
                if strcmp(obj.edges(edgeIndex).type,'internal')                 %if the edge is internal
                    obj.edges(edgeIndex).dihedralAngle.ang=angle(counter);      %set angle of edge
                    counter = counter+1;
                end
            end
            
        end
        
        %%SET MODES%%
        function [] = setModeTo(obj, vertexIndex, inputMode)
            %SETMODETO Sets the mode of the indicated vertex to the
            %specified mode
            
            %Can be used in any state
            
            obj.vertices(vertexIndex).mode=inputMode;
            
            if not(obj.isOptimizing)
                if strcmp(obj.status,'doKinematics')
                    obj.resetToUndetermined();
                end
            end
            
        end
        function [] = switchMode(obj,vertexIndex)
            %SWITCHMODE Switches Mode at indicated vertex (between
            %true/false)
            
            %Can be used in any state
            
            obj.vertices(vertexIndex).mode=not(obj.vertices(vertexIndex).mode);
            
            if not(obj.isOptimizing)
                if strcmp(obj.status,'doKinematics')
                    obj.resetToUndetermined();
                end
            end
        end
        
     %%%Determinging origami
        function [] = calculateSectorLengths(obj)
         
     %Requires state editFaces
     obj.checkState('editFaces');
     
     %Calculate sector angles
     obj.calculateSectorAngles();
     
     %Calculate edge lengths
     obj.calculateEdgeLengths();
     
     %Before updating the state, make sure that a face has been fixed
     if obj.fixedFaceIndex==0 
        error('A valid face must be fixed before proceeding to the doKinematics status');
     end
     
     %Set the status to doKinematics
     obj.status='doKinematics';
     end %calculateSectorLengths
     
        function [] = determineOrigami(obj)
        %Requires status doKinematics
        obj.checkState('doKinematics');
            
        %Introducing a counter to avoid infinite loops
        counter=0;
        
        %While loop, checking whether the computation is complete
            while not(obj.isOrigamiDetermined())
                
                oldcounter=counter;
                
                %Looping over vertices and determining whatever can be
                %determined...
                for j=1:size(obj.vertices,2)
                    
                    %Skip deleted vertices
                    if obj.vertices(j).deleted
                        continue;
                    end
                    "current vertex = "+num2str(j);
                    
                    %determine everything that can be determined and hasn't
                    %been determined yet
                    if obj.isVertexDeterminable(j) && ...
                            not(obj.vertices(j).determined)
                        obj.determineVertex(j);
                        counter=counter+1;
                    end
                    
                end
                
                if oldcounter==counter
                    error('Evaluation of origami kinematics failed. Are you sure there are no loops in the crease graph?');
                end
                
            end
            
            %Check if everything is determined
            if not(obj.isOrigamiCompletelyDetermined())
                error('Failed to determine the entire origami. Probably, not all subOrigamis were determined.');
            end
    
     end
     
             %%%%%Function for resetting the whole origami to undetermined%%%%%
        
        function [] = resetToUndetermined(obj)
            %RESETTOUNDETERMINED(OBJ) resets all edges, faces and vertices
            %to undetermined and undetermined s.t. the determination can
            %start anew.
            
            [obj.vertices(:).determined]=deal(false);
            
            [obj.edges(:).determined]=deal(false);
            
            [obj.faces(:).determined]=deal(false);
            
        end
     
        
 %%%%%%%GEOMETRY%%%%
      %%Crease length%%
        function [] = expandCreaseToLength(obj,edgeIndex,length)
            %EXPANDCREASETOLENGTH(OBJ,EDGEINDEX,FACTOR) Moves the end 
            %vertex of the specified crease s.t. the specified crease is
            %elongated to the specified length.
            
            %Check if the specified edge is valid
            if edgeIndex>obj.numberOfEdges()
                error('The origami currently only has %d edges, but the specified edge index is %d. This is invalid.',obj.numberOfEdges(),edgeIndex);
            end
            if obj.edges(edgeIndex).deleted
                error('The specified edge (nr. %d) has been deleted and can not be elongated.',edgeIndex);
            end
            
            %Check if the length is valid
            if length<1e-6
                error('The specified new length of crease nr. %d of %.4e is invalid.',edgeIndex,length);
            end
            
            %Get the source vertex and the direction vector of the crease
            sourceVertex=obj.edges(edgeIndex).sourceVertex;
            dirVec=obj.awayVector0(sourceVertex,edgeIndex); %returns a unit vector
            
            %Get new position
            newPos=obj.vertices(sourceVertex).position0+length*dirVec;
            
            %Get end position (vertex to move)
            endVertex=obj.edges(edgeIndex).endVertex;
            
            %move vertex
            obj.moveVertexTo(endVertex,newPos);           
        end
        function [] = expandCreaseByFactor(obj,edgeIndex,factor)
            %EXPANDCREASEBYFACTOR(OBJ,EDGEINDEX,FACTOR) Moves the end 
            %vertex of the specified crease s.t. the specified crease is
            %elongated by the specified factor.
            %Mainly calls the method expandCreaseToLength.
            
            %Check if factor is valid
            
            if factor<0 || factor>1e6
                error('The specified extensionFactor of crease nr. %d of %.4e is invalid.',edgeIndex,factor);
            end
            
            %Get the current length of the edge (also checks whether index
            %is reasonable)
            oldLength=obj.calculateEdgeLength(edgeIndex);
            
            %stretch/contract the edge
            obj.expandCreaseToLength(edgeIndex,factor*oldLength);
            
        end
        function [] = expandCreaseTillXEquals(obj,edgeIndex,newX)
            %EXPANDCREASETILLXEQUALS(OBJ,EDGEINDEX,NEWX) Expands or
            %contracts the specified edge (moving the end point while
            %holding the source vertex) until the x-coordinate of the
            %endpoint is equal to the input argument.
            
            %Check if edge specification is valid
            if edgeIndex>obj.numberOfEdges()
                error('The origami currently only has %d edges, but the specified edge index is %d. This is invalid.',...
                    obj.numberOfEdges(),edgeIndex);
            end
            
            if obj.edges(edgeIndex).deleted
                error('The specified edge (nr. %d) has been deleted and can not be elongated.',...
                    edgeIndex);
            end
            
            %Get the source position
            sourcePos=obj.vertices(obj.edges(edgeIndex).sourceVertex).position0;
            
            %Get the direction vector
            dirVec=obj.sourceTargetVec0(edgeIndex); %returns a unit vector
            
            %Check if such an extension is possible
            deltaX=newX-sourcePos(1);
            if sign(deltaX)~=sign(dirVec(1))
                error('The edge %d can not be extended to reach the specified x-coordinate.',edgeIndex);
            end
            
            %Get the new length
            newLength=deltaX/dirVec(1);
            
            %Perform the stretch
            obj.expandCreaseToLength(edgeIndex,newLength);
            
        end
        function [] = expandCreaseTillYEquals(obj,edgeIndex,newY)
            %EXPANDCREASETILLYEQUALS(OBJ,EDGEINDEX,NEWX) Expands or
            %contracts the specified edge (moving the end point while
            %holding the source vertex) until the y-coordinate of the
            %endpoint is equal to the input argument.
            
            %Check if edge specification is valid
            if edgeIndex>obj.numberOfEdges()
                error('The origami currently only has %d edges, but the specified edge index is %d. This is invalid.',...
                    obj.numberOfEdges(),edgeIndex);
            end
            
            if obj.edges(edgeIndex).deleted
                error('The specified edge (nr. %d) has been deleted and can not be elongated.',...
                    edgeIndex);
            end
            
            %Get the source position
            sourcePos=obj.vertices(obj.edges(edgeIndex).sourceVertex).position0;
            
            %Get the direction vector
            dirVec=obj.sourceTargetVec0(edgeIndex); %returns a unit vector
            
            %Check if such an extension is possible
            deltaY=newY-sourcePos(2);
            if sign(deltaY)~=sign(dirVec(2))
                error('The edge %d can not be extended to reach the specified y-coordinate.',edgeIndex);
            end
            
            %Get the new length
            newLength=deltaY/dirVec(2);
            
            %Perform the stretch
            obj.expandCreaseToLength(edgeIndex,newLength);
        end
      %%Moving Vertices%%
        function [] = moveVertexBy(obj, vertexIndex, translationVec)
            %MOVEVERTEXBY(OBJ, VERTEXINDEX, TRANSLATIONVEC) Moves vertex by
            %specified translationVec. Briefly, this is a shortcut to
            %moveVertexTo.
            
            %If the vertex is only adjacent to boundary edges, editFaces is
            %sufficient. Otherwise, editGraph is required.
            
            %find out if there are adjacent edges that are not on the
            %boundary
            isBoundary=true;
            adjEdg=[obj.vertices(vertexIndex).incomingEdges,...
                obj.vertices(vertexIndex).outgoingEdges];
            for j=adjEdg
                if not(strcmp(obj.edges(j).type,'boundary'))
                    %this edge is not on the boundary
                    isBoundary=false;
                    break;
                end
            end
            
            if isBoundary
                obj.checkState('editFaces');
            else
                obj.checkState('editGraph');
            end
            
            %new position needs to be 2-by-1 vector. If it's a 3-by-1
            %vector, the z-component will be set to 0.
            if not(size(translationVec,1)<=3 && size(translationVec,1)>=2)...
                    || size(translationVec,2)~=1
                error('Translation vector must be specified as 2-by-1 or 3-by-1 (last entry ignored) vector');
            end
            
            newPos=obj.vertices(vertexIndex).position0 + [translationVec(1);...
                translationVec(2); 0];
            obj.moveVertexTo(vertexIndex,newPos);
            
        end
        function [] = moveVertexTo(obj, vertexIndex, newPos)
            %MOVEVERTEXTO(OBJ, VERTEXINDEX, NEWPOS) moves vertex
            %(vertexIndex) to the specified new position. The function
            %considers the possibility that a vertex is moved "to the other
            %side" of an edge and resorts entities accordingly.
            
            %If the vertex is only adjacent to boundary edges, editFaces is
            %sufficient. Otherwise, editGraph is required.
            
            %find out if there are adjacent edges that are not on the
            %boundary
            isBoundary=true;
            adjEdg=[obj.vertices(vertexIndex).incomingEdges,...
                obj.vertices(vertexIndex).outgoingEdges];
            for j=adjEdg
                if not(strcmp(obj.edges(j).type,'boundary'))
                    %this edge is not on the boundary
                    isBoundary=false;
                    break;
                end
            end
            
            if isBoundary
                obj.checkState('editFaces');
            else
                obj.checkState('editGraph');
            end

            %new position needs to be 2-by-1 vector. If it's a 3-by-1
            %vector, the z-component will be set to 0.
            if not(size(newPos,1)<=3 && size(newPos,1)>=2) || size(newPos,2)~=1
                error('New Position must be specified as 2-by-1 or 3-by-1 (last entry ignored) vector');
            end
            
            %Save the containers to temporary containers.
            sourcesIncoming=obj.vertices(vertexIndex).incomingNeighbours;
            incomingEdges=obj.vertices(vertexIndex).incomingEdges;
            endsOutgoing=obj.vertices(vertexIndex).outgoingNeighbours;
            outgoingEdges=obj.vertices(vertexIndex).outgoingEdges;
            
            %Remove all outgoing edges, one by one.
            
            while size(obj.vertices(vertexIndex).outgoingEdges,2)>0
                edgeIndexTemp=obj.vertices(vertexIndex).outgoingEdges(1);
                obj.removeCreaseByIndex(edgeIndexTemp);
                %Note: This also sets the edges to deleted. Reverse this
                %when re-adding them
            end
            
            %Remove all incoming edges, one by one
            
            while size(obj.vertices(vertexIndex).incomingEdges,2)>0
                edgeIndexTemp=obj.vertices(vertexIndex).incomingEdges(1);
                obj.removeCreaseByIndex(edgeIndexTemp);
                %Once more, this sets the edges to deleted...
            end
            
            %Reset position of vertex
            if length(newPos)==2
                newPos=[newPos;0];
            end
            obj.vertices(vertexIndex).position0=newPos;
            if not(obj.isOptimizing)
                obj.vertices(vertexIndex).position=newPos;
            end
            
            %Re-add first incoming edge - if existing (source vertex, e.g.)
            
            if size(incomingEdges,2)~=0
                sourceVertex=sourcesIncoming(1);
                %adding the edge manually to source and end
                obj.addOutgoingCrease(sourceVertex,incomingEdges(1));
                obj.addIncomingCrease(vertexIndex,incomingEdges(1));
                obj.edges(incomingEdges(1)).deleted=false;  %reversing the marking as deleted from before
            end
            
            %Re-add outgoing edges (if existing)
            
            for j=1:size(outgoingEdges,2)
                endIndex=endsOutgoing(j);
                %Adding the edges manually
                obj.addOutgoingCrease(vertexIndex,outgoingEdges(j));
                obj.addIncomingCrease(endIndex,outgoingEdges(j));
                obj.edges(outgoingEdges(j)).deleted=false;  %reversing the marking as deleted from before
            end
            
            %Re-add remaining incoming edges
            
            for j=2:size(incomingEdges,2)   %The first edge was already re-added
                sourceVertex=sourcesIncoming(j);
                %adding the edge manually to source and end
                obj.addOutgoingCrease(sourceVertex,incomingEdges(j));
                obj.addIncomingCrease(vertexIndex,incomingEdges(j));
                obj.edges(incomingEdges(j)).deleted=false;  %reversing the marking as deleted from before
            end
            
            
            %if we're optimizing, recalculate everything around this vertex
            if obj.isOptimizing
                affectedEdges=[obj.vertices(vertexIndex).incomingEdges,...
                    obj.vertices(vertexIndex).outgoingEdges];
                affectedVertices=[vertexIndex, ...
                    obj.vertices(vertexIndex).incomingNeighbours,...
                    obj.vertices(vertexIndex).outgoingNeighbours];
                
                for j=affectedEdges
                    obj.calculateEdgeLength(j);
                end
                
                for j=affectedVertices
                    obj.calculateSectorAngle(j);
                end
            end
            
         
        end
        function [] = resetAdjacentAuxiliaries(obj,verticesToReset)
            %RESETADJACENTAUXILIARIES(OBJ,VERTICESTORESET) Goes through all
            %vertices specified in verticesToReset and resets all adjacent
            %auxiliary vertices. After the reset, these auxiliary vertices
            %are in the middle between their two neighbours.
            
            %State should be editFaces - but it's expected that this
            %function is usually called with isOptimizing=true
            obj.checkState('editFaces');
            
            %Iterate through vertices
            for j=1:size(verticesToReset,2)
                currentIndex=verticesToReset(j);
                
                %Check whether the entry is valid
                if obj.vertices(currentIndex).deleted
                    warning('Vertex nr. %d was specified but is marked as deleted. Skipping this one.',currentIndex);
                    continue;
                end
                
                if currentIndex>obj.numberOfVertices()
                    error('The specified vertices for resetting of aux-vertices include nr. %d, but the origami object only has %d vertices. Aborting.',currentIndex, obj.numberOfVertices());
                end
                
                %Go through outgoing neighbours and check for aux-vertices
                outNeighbours=obj.vertices(currentIndex).outgoingNeighbours;
                %check the first one
                if strcmp(obj.vertices(outNeighbours(1)).type,'aux')
                    obj.resetAuxVertex(outNeighbours(1));
                end
                %...and the last one
                if strcmp(obj.vertices(outNeighbours(end)).type,'aux')
                    obj.resetAuxVertex(outNeighbours(end));
                end
                
            end
            
            %second loop to recalculate the edge lengths and sector angles
            for j=1:size(verticesToReset,2)
                currentIndex=verticesToReset(j);
                
                %Skipping tests this time: already passed in first loop
                
                %Go through outgoing neighbours and check for aux-vertices
                outNeighbours=obj.vertices(currentIndex).outgoingNeighbours;
                %check the first one
                if strcmp(obj.vertices(outNeighbours(1)).type,'aux')
                    obj.resetAuxVertexRecalculate(outNeighbours(1));
                end
                %...and the last one
                if strcmp(obj.vertices(outNeighbours(end)).type,'aux')
                    obj.resetAuxVertexRecalculate(outNeighbours(end));
                end
                
            end
            
        end
        
%%%%%%%%PLOT%%       
        function [ax2] = plotOrigamiGraph(obj,ax)
            %PLOTORIGAMIGRAPH Plots the base graph of the origami object.
            %Deleted crease lines and vertices are not shown.
            %   INPUTS:
            %   ax (optional): Axis handle in which the graph is plotted.
            %   If no axis is provided, the graph is plotted to a new
            %   figure.
            %   OUTPUTS:
            %   ax2: Handle to the figure the graph was plotted to
            %Can be called at any status
            
            if nargin==1
                ax2=figure;
                ax2=axes('Parent',ax2);
                hold(ax2,'on');
            else
                ax2=ax;
            end
            
            %Set the axis to invisible to plot "in the background"
            ax2.Visible='off';
            
            %This parameter controls how far away from the vertices and
            %edges their labels are 
            distLabel=0.03;
            
            %Plotting the vertices
            for j=1:obj.numberOfVertices()
                if not(obj.vertices(j).deleted)
                    format=obj.vertexTypeFormatFromIndex(j);
                    plot(ax2,obj.vertices(j).position0(1),...
                        obj.vertices(j).position0(2),format,'MarkerSize',30);
                    
                    %Add a label
                    text(obj.vertices(j).position0(1)+distLabel,...
                        obj.vertices(j).position0(2)-distLabel,num2str(j));
                    hold on;
                    if strcmp(obj.vertices(j).type,'interior') ...  %also plot the mode for interior vertices
                         || (strcmp(obj.vertices(j).type,'source') && isa(obj,'subOrigami'))
                        format=obj.vertexModeFormatFromIndex(j);
                        plot(ax2,obj.vertices(j).position0(1),...
                            obj.vertices(j).position0(2),format,'MarkerSize',30);
                    end
                end
            end
            
            %Plotting the creases
            for j=1:obj.numberOfEdges()
                if not(obj.edges(j).deleted)
                    format=obj.creaseFormatFromIndex(j);
                    origin=obj.vertices(obj.edges(j).sourceVertex).position0;
                    target=obj.vertices(obj.edges(j).endVertex).position0;
                    difference=target-origin;
                    q=quiver(origin(1),origin(2),difference(1),difference(2),format,'AutoScale','off');
                    %q.MaxHeadSize=3/norm(difference,2);
                    %Add a label
                    labelPos=origin+0.5*difference+[distLabel;-distLabel;0];
                    text(labelPos(1),labelPos(2),num2str(j),'Color','red');
                end
            end
            
            grid minor;
            daspect([1 1 1]);
            
            %Make the axis visible again
            ax2.Visible='on';
            
            %BY JIM: add title
            title(obj.label);
        end %plotOrigamiGraph
        function [ax2] = plotOrigamiUndeformed(obj, ax, plotEdges)
            %PLOTORIGAMIUNDEFORMED(OBJ, AX, PLOTEDGES) plots the origami in 
            %undeformed configuration. Uses the face definitions and 
            %initial positions of the vertices.
            %Use the second argument to specify an axis object in which to
            %plot
            %The third argument can be used to pass a flag specifying
            %whether or not to plot the edges.
            
            %can be called at any state but editGraph
            if strcmp(obj.status,'editGraph')
                obj.checkState('editFaces');
            end
            
            if size(obj.faces,2)==0
                error('No faces defined');
            end
            
            edgeCol=[0 0 0];
            
            if nargin==1
                ax2=figure;
                ax2=axes('Parent',ax2);
                hold(ax2,'on');
            elseif nargin >= 2
                ax2=ax;
                hold(ax2,'on');
                if nargin==3
                    if not(plotEdges)
                        edgeCol='none';
                    end
                end
            end
            
            %iterate through all vertices and extract all positions
            %(position0)
            
            positions=zeros(obj.numberOfVertices(),3);
            for j=1:obj.numberOfVertices()
                positions(j,:)=obj.vertices(j).position0';
            end
            
            %generate random colors between 0.8 and 1 for the faces
            rng('default');
            col=rand([obj.numberOfFaces(),1])+0.8;
            
            %iterate through all faces
            for j=1:obj.numberOfFaces()
                currentFace=obj.faces(j);
                
                %skip the face if it has been deleted
                if currentFace.deleted
                    continue;
                end
                
                %patch the face
                if currentFace.fixed
                    col(j)=0;
                end
                
                %obtain the face center from its edge vertices
                labelPos=zeros(3,1);
                nCorners=size(currentFace.connectedVertices,2);
                for k=1:nCorners
                    labelPos=labelPos+(1/nCorners)*obj.vertices(...
                        currentFace.connectedVertices(k)).position0;
                end
                
                patch('Faces',currentFace.connectedVertices, ...
                    'Vertices', positions,'FaceAlpha',0.5,...
                    'FaceVertexCData',col(j),'FaceColor','flat',...
                    'EdgeColor',edgeCol);
                
                %add a label with the id of the face
                text(labelPos(1), labelPos(2), num2str(j));
                
                hold on;
                
            end %loop over faces
            
            camlight;
            daspect([1,1,1]);
            
        end %function plotOrigamiUndeformed

        %Functions for plotting the deformed origami
        function [ax2] = plotOrigamiDeformed(obj, ax, plotEdges)
            %PLOTORIGAMIUNDEFORMED(OBJ, AX, PLOTEDGES) plots the origami in 
            %undeformed configuration. Uses the face definitions and 
            %initial positions of the vertices.
            %Use the second argument to specify an axis object in which to
            %plot
            %The third argument can be used to pass a flag specifying
            %whether or not to plot the edges.
            
            %Requires status doKinematics
            obj.checkState('doKinematics');
            
            if size(obj.faces,2)==0
                error('No faces defined');
            end
            
            %(!) makeshift to prevent error in plotting 2.1445 for i=2834 | sqrt(3) for 22586
            idx = 2834;
            if idx == 2834 
                obj.vertices(1).position = [-1;2.1445;0];
            elseif idx == 22586
                obj.vertices(1).position = [-1;sqrt(3);0];
            end
            
            edgeCol=[0 0 0];
            
            
            if nargin==1
                ax2=figure;
                ax2=axes('Parent',ax2);
                hold(ax2,'on');
            elseif nargin >= 2
                ax2=ax;
                
                cla; %clear anything that was on that axis
                
                hold(ax2,'on');
                if nargin==3
                    if not(plotEdges)
                        edgeCol='none';
                    end
                end
            end
            
            %iterate through all vertices and extract all positions
            %(position0)
            
            positions=zeros(obj.numberOfVertices(),3);
            for j=1:obj.numberOfVertices()
                positions(j,:)=obj.vertices(j).position';
            end
            
            %generate random colors between 0.8 and 1 for the faces
            rng('default');
            col=rand([obj.numberOfFaces(),1])+0.8;
            
            %iterate through all faces
            for j=1:obj.numberOfFaces()
                currentFace=obj.faces(j);
                
                %skip the face if it has been deleted
                if currentFace.deleted
                    continue;
                end
                
                %patch the face
                if currentFace.fixed
                    col(j)=0;
                end
                
                %obtain the face center from its edge vertices
                labelPos=zeros(3,1);
                nCorners=size(currentFace.connectedVertices,2);
                for k=1:nCorners
                    labelPos=labelPos+(1/nCorners)*obj.vertices(...
                        currentFace.connectedVertices(k)).position0;
                end
                
                patch('Faces',currentFace.connectedVertices, ...
                    'Vertices', positions,'FaceAlpha',0.5,...
                    'FaceVertexCData',col(j),'FaceColor','flat',...
                    'EdgeColor',edgeCol);
                
                %add a label with the id of the face
%                 text(labelPos(1), labelPos(2), num2str(j));
                
                hold on;
                
            end %loop over faces
            
            %set aspect ratio
            daspect([1 1 1]);
            
            %some more plot options
            camlight;
            %adjust limits manually
            zlim(obj.animationZLims);
            xlim(obj.animationXLims);
            ylim(obj.animationYLims);
            %view position
            view(obj.animationView);
            
        end %function plotOrigamiDeformed
        

            
        
        function [] = animateOrigami(obj,endAngle,filename,nSteps)
            %ANIMATEORIGAMI(OBJ, ENDANGLE, FILENAME) Creates an animation
            %of an origami as it folds from the initial configuration to
            %the state when the input dihedral angle is equal to endAngle.
            %Saves the animation in filename (as a gif file). Uses nSteps
            %frames (default: 20)
            %If endAngle is a vector of size 2, the origami is animated
            %between the two stated by these two dihedral angles.
            
            %Requires status doKinematics
            obj.checkState('doKinematics');
            
            %Setting the number of steps, if not specified by input.
            if nargin<4
                nSteps=20;
            end
            
            %Calculate Input angles for each step
            
            %specify the dihedral angles. Do not use 0 as starting
            %angle, because 0 is not recognized as valid starting
            %angle by the "isReady"-functions
            if size(endAngle,2)<2
                startAngle=1e-3*sign(endAngle);
            else
                startAngle=endAngle(1);
            end
            
            dihedralAngles=linspace(startAngle,endAngle(end),nSteps)';
            
            %define pauseTime (time between frames)
            pauseTime=0.1;
            
            for j=1:nSteps
                
                %Create the deformed origami
                obj.setDihedralAngles(dihedralAngles(j,:));
                obj.determineOrigami();
                
                %plot the origami in the current axis - if available
                if j==1
                    ax=obj.plotOrigamiDeformed();
                else
                    ax=obj.plotOrigamiDeformed(ax);
                end
                
                %get a movie frame
                frame = getframe(gcf);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);
                
                % Write to the GIF File
                if j == 1
                    imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
                else
                    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime', pauseTime);
                end
            end
            
            
        end
        function [] = animateOrigamiWithZoomIn(obj,endAngle,filename,nSteps,...
                endXLims, endYLims, endZLims, endViewAngle)
            %ANIMATEORIGAMIWITHZOOMIN(OBJ,ENDANGLE,FILENAME,NSTEPS,ENDXLIMS,ENDYLIMS,ENDZLIMS,ENDVIEWANGLE)
            %Is similar to animateOrigami, but in the and also adds a
            %zoom-in during which the frame zooms in to, in the end, have
            %the limits (and view angle) specified in the input argument.
            
            %execute animateOrigami, as usual
            obj.animateOrigami(endAngle,filename,nSteps);
            
            %number of steps during zoom-in
            zoomInSteps=nSteps;
            
            pauseTime=0.1;
            
            %get the x-limits
            xLimits=[linspace(obj.animationXLims(1),endXLims(1),zoomInSteps);...
                linspace(obj.animationXLims(2),endXLims(2),zoomInSteps)];
            %y-limits
            yLimits=[linspace(obj.animationYLims(1),endYLims(1),zoomInSteps);...
                linspace(obj.animationYLims(2),endYLims(2),zoomInSteps)];
            %z-limits
            zLimits=[linspace(obj.animationZLims(1),endZLims(1),zoomInSteps);...
                linspace(obj.animationZLims(2),endZLims(2),zoomInSteps)];
            %view-angles
            viewAngles=[linspace(obj.animationView(1),endViewAngle(1),zoomInSteps);...
                linspace(obj.animationView(2),endViewAngle(2),zoomInSteps)];
            
            %loop over steps
            for j=1:zoomInSteps
                
                %reset the limits and the view angle
                xlim(xLimits(:,j)');
                ylim(yLimits(:,j)');
                zlim(zLimits(:,j)');
                view(viewAngles(1,j),viewAngles(2,j));
                
                %get frame and attach to the gif
                
                %get a movie frame
                frame = getframe(gcf);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);
                
                % Write to the GIF File
                imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime', pauseTime);
                
            end
            
        end
        %Exporting the folded origami to an stl file
        function [] = writeToSTL(obj,filename,writeEdges)
            %WRITETOSTL(OBJ,FILENAME) Creates an stl file with the
            %specified name that contains the origami in the folded
            %configuration. This function requires that the origami is
            %determined (throws error if it isn't).
            %This function requires the geom3d toolbox. This toolbox can be
            %installed from here: https://www.mathworks.com/matlabcentral/fileexchange/24484-geom3d
            %if the second input argument is true, two additional files are
            %created: One containing all the coordinates of the points and
            %one with the edge connectivity. Both of these files derive
            %their names from the name of the stl-file and are in
            %csv-format.
            
            %requires state doKinematics
            obj.checkState('doKinematics');
            
            %origami must be determined
            if not(obj.isOrigamiDetermined())
                error('Origami must be determined in order to create an stl-file.');
            end
            
            %filename must be a char array
            if not(ischar(filename))
                error('Filename must be specified as char array.');
            end
            
            %gather connectivity
            [points,edgeConn,faceConn,~]=obj.gatherTriangulation([],[],{},0);
            
            %connectivity of triangulated faces
            tri=triangulateFaces(faceConn);
            
            %ensure that both sides will be visible in the stl
            tri=[tri;fliplr(tri)];
            
            %triangulation
            tr=triangulation(tri,points(1,:)',points(2,:)',points(3,:)');
            
            %write triangulation to stl
            stlwrite(tr,filename);
            
            if nargin==3 && writeEdges
                %get an array with all vertex positions

                basename=filename(1:end-4);
                pointName=strcat(basename,'points.csv');
                edgeName=strcat(basename,'edges.csv');
                
                %write to csv
                csvwrite(pointName,points');
                csvwrite(edgeName,edgeConn');
            end
            
        end
        %Function for exporting the crease pattern as pd
        function [] = exportCreasePatternToPDF(obj,filename,paperOrientation,paperType)
            %EXPORTCREASEPATTERNTOPDF(OBJ,FILENAME,FORMAT) Exports the
            %crease pattern of the origami object to a pdf. This pdf can be
            %printed and used to fold the origami.
            %The second input (filename) is a string with the name of the
            %file.
            %The last two inputs are to specify the paper format. They are
            %directly passed to the figure. Cf. the MATLAB manual for
            %reference: https://www.mathworks.com/help/matlab/ref/matlab.ui.figure-properties.html#buiwuyk-1-PaperOrientation
            
            %Start new figure and axes
            fig=figure();
            ax=axes();
            
            %if specified, forward paper orientation
            if nargin>=3
                fig.PaperOrientation=paperOrientation;
            end
            
            %if specified, forward paper type/size
            if nargin>=4
                fig.PaperType=paperType;
            end
            
            %plot the crease pattern
            obj.plotCreasePattern(ax);
            
            %set aspect ratio and turn off the axes
            daspect([1 1 1]);
            axis off;
            
            %scale the figure to the paper
            
            %set the units of the paper
            fig.PaperUnits='points';
            %get the size of the paper and the figure position
            size=fig.PaperSize;
            pos=fig.Position;
            %expand figure to paper size
            fig.Position=[pos(1:2),size];
            
            %scale axis to fill whole figure
            ax.Position=[0,0,1,1];
            
            %save to file
            saveas(fig,filename);
            
            
        end
      
    end %public methods section
%%%%PROTECTED METHODS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = 'protected')
        
     %%%My auxiliary methods

        
%%%%%%%%COPYING ORIGAMI OBJECT%%
        %Needed in handle class
        function cp = copyElement(obj)
            %CP=COPYELEMENT(OBJ) creates a deep copy of the origami object.
            %This overrides the default copy function of the mixin.copyable
            %class.
            
            cp=origami('SingleCrease'); %something. Will be destroyed and rebuilt later; this is just for creating the constructor.
            
            props=properties(obj);
            for j=1:length(props)
                thisProperty=props{j};
                if isa(obj.(thisProperty),'handle')
                    %handle->use copy routine
                    cp.(thisProperty)=copy(obj.(thisProperty));
                else
                    cp.(thisProperty)=obj.(thisProperty);
                end
            end
            %reconnect everything
            cp.reconnectOrigami();
            
end
        function [] = reconnectOrigami(obj)
            %RECONNECTORIGAMI(OBJ) Redoes all "pointer-like" connections:
            %Reconnects the edges with the faces and vertices (normal
            %vectors, dihedral angles) and the faces with the vertices
            %(normal vectors).
            
            
            %reconnect the vectors and the dihedral angles
            for j=1:obj.numberOfVertices()
                
                if obj.vertices(j).deleted
                    continue;
                end
                
                %connect to the normal vectors
                for k=1:length(obj.vertices(j).adjacentFaces)
                    %EDIT: added if statement to prevent error when
                    %returning to editGraph from editFaces. adjacentFaces
                    %of internal creases are not deleted, but faces are.
                    %Therefore an error occurs when indexing
                    %obj.faces(connectedFaces)
                    
                    if ~isempty(obj.faces) %only run if faces exist
                        connectedFace=obj.vertices(j).adjacentFaces(k);
                        obj.vertices(j).normalVectors(k)=obj.faces(connectedFace).normalVector;
                    end
                end
                
                %connect to the edges
                %first the incoming ones
                for k=1:length(obj.vertices(j).incomingEdges)
                    thisEdge=obj.vertices(j).incomingEdges(k);
                    indexHere=obj.edges(thisEdge).indexAtEnd;
                    obj.vertices(j).creaseVectors(indexHere)=obj.edges(thisEdge).directionVector;
                    obj.vertices(j).dihedralAngles(indexHere)=obj.edges(thisEdge).dihedralAngle;
                end
                
                %now the outgoing ones
                for k=1:length(obj.vertices(j).outgoingEdges)
                    thisEdge=obj.vertices(j).outgoingEdges(k);
                    indexHere=obj.edges(thisEdge).indexAtSource;
                    obj.vertices(j).creaseVectors(indexHere)=obj.edges(thisEdge).directionVector;
                    obj.vertices(j).dihedralAngles(indexHere)=obj.edges(thisEdge).dihedralAngle;
                end
                
            end
            
            %reconnect the edges to their adjacent faces (i.e. normals)
            for j=1:obj.numberOfEdges()
                if obj.edges(j).deleted
                    continue;
                end
                
                if obj.edges(j).faceLeft~=0
                    obj.edges(j).normalVecLeft=obj.faces(obj.edges(j).faceLeft).normalVector;
                end
                
                if obj.edges(j).faceRight~=0
                    obj.edges(j).normalVecRight=obj.faces(obj.edges(j).faceRight).normalVector;
                end
                
            end
            
            
                end
        
%%%%%%%%Functions for bookkeeping around vertices%%
        %set first index of origami to 'source'
        function [] = setSourceVertex(obj,pos0)
            %SETSOURCEVERTEX Sets the first vertex in the origami and sets
            %its type to 'source'.
            
            if obj.numberOfVertices~=0
                error('Source vertex must be the first vertex to be set');
            end
            
            %by default the region with the source vertex is the first one
            obj.vertices=vertex(pos0,1,false,'source',[0;0;1]);
            obj.source1=1;    
        end %setSourceVertex
    
        %Bookkeeping of parameters with extension/removal of edges
        function [] = addOutgoingCrease(obj,vertexIndex,edgeIndex)
            %ADDOUTGOINGCREASE Adjusts the involved vertex and edge
            %objects s.t. the new outgoing edge is registered.
            
            %First check whether we're dealing with an auxiliary node
            if strcmp(obj.vertices(vertexIndex).type,'aux')
                error('Extension of auxiliary vertices is not allowed!');
            end
            
            %Check if the new crease has a near-0 angle to an existing edge
            obj.checkIsNewCreaseAllowed(vertexIndex, edgeIndex);
            
            %Get number of edges that are already registered as "outgoing"
            numberOfOutgoingEdges=size(obj.vertices(vertexIndex).outgoingEdges,2);
            numberOfIncomingEdges=size(obj.vertices(vertexIndex).incomingEdges,2);
            deg=size([obj.vertices(vertexIndex).incomingEdges, ...
                obj.vertices(vertexIndex).outgoingEdges],2);
            corrDeg=deg;
            
            %Get edge to the left and to the right of the new edge
            [left, right, jl, jr]=obj.findLeftAndRight(vertexIndex,edgeIndex);
            jl=jl-numberOfIncomingEdges;
            if jr==1    %the new edge is the right-most edge
                jr=numberOfOutgoingEdges+1;
            elseif jr==0    %nothing at this vertex so far
                jr=1;
            else
                jr=jr-numberOfIncomingEdges;
            end
            
            %Performing some checks
            
            %only the last incoming edge is allowed to be the next edge to
            %the left of the new edge
            forbiddenLefts=obj.vertices(vertexIndex).incomingEdges(1:end-1);
            if ismember(left,forbiddenLefts) && size(obj.vertices(vertexIndex).outgoingEdges,2)==0
                obj.resortEdges(vertexIndex,edgeIndex);
                return; %the "resortEdges" function is doing all the work. Quit this function to avoid double-doing the work
            elseif ismember(left,forbiddenLefts) && size(obj.vertices(vertexIndex).outgoingEdges,2)>0
                error('New outgoing edge appears to be pointing towards the incoming ones. This is not allowed.');
            end
            
            %only the first incoming edge is allowed to be the next edge to
            %the right of the new edge
            forbiddenRights=obj.vertices(vertexIndex).incomingEdges(2:end);
            if ismember(right,forbiddenRights) && size(obj.vertices(vertexIndex).outgoingEdges,2)==0
                obj.resortEdges(vertexIndex, edgeIndex);
                return; %the "resortEdges" function is doing all the work. Quit this function to avoid double-doing the work
            elseif ismember(right,forbiddenRights) && size(obj.vertices(vertexIndex).outgoingEdges,2)>0
                error('New outgoing edge appears to be pointing towards the incoming ones. This is not allowed.');
            end
            
            %Check if outgoing crease is added to a vertex with no
            %inputs and whether the new crease is not "within existing
            %ones"
            
            criterion1=numberOfOutgoingEdges>0; %not adding the first one - then, things are trivial
            criterion2=deg==numberOfOutgoingEdges; %no incoming edges
            criterion3=(jr==1) || (jl==numberOfOutgoingEdges);
            
            if all([criterion1, criterion2, criterion3])
                %The findLeftAndRight function, in this case, can't
                %properly determine whether to add the crease to the left
                %or to the right of the existing crease(s). We therefore
                %use this function to adapt accordingly.
                
                %check if we're dealing with a source vertex or a local 
                %(region-wise) source vertex - the only cases when this is 
                %allowed
                if not(strcmp(obj.vertices(vertexIndex).type,'source')) && ...
                        not(obj.isLocalSource(vertexIndex))
                    error('Only source vertices may have outputs without having inputs');
                end
                
                [jl, jr] = obj.decideLeftRight(vertexIndex,edgeIndex);
                if jl~=0
                    corrDeg=deg+1;  %for the construction of arrayleft. 
                    %If this is not incremented, arrayleft will stay empty 
                    %but should contain everything.
                end
            end
            
            %Create array with indices of outgoing edges coming before the
            %new edge
            arrayleft=1:(mod(jl,corrDeg));
            
            %Create array with indices of outgoing edges coming after the
            %new edge
            arrayright=jr:numberOfOutgoingEdges;
            
            %update indexAtSource for outgoing edges
            for j=1:size(arrayright,2)
                edgeIndextemp=obj.vertices(vertexIndex).outgoingEdges(arrayright(j));
                obj.edges(edgeIndextemp).indexAtSource=...
                    obj.edges(edgeIndextemp).indexAtSource+1;
            end
            
            %update indexAtSource of "this" edge
            obj.edges(edgeIndex).indexAtSource=jl+1+numberOfIncomingEdges;
            
            %IMPORTANT: Update the edges before updating the vertices. The
            %edge update is written s.t. it refers to the old/unchanged
            %edge list at the vertex.
            
            %re-write list of outgoing edges at vertex
            obj.vertices(vertexIndex).outgoingEdges=[...
                obj.vertices(vertexIndex).outgoingEdges(arrayleft),...
                edgeIndex,...
                obj.vertices(vertexIndex).outgoingEdges(arrayright)];
            
            %re-write list of "child vertices"
            obj.vertices(vertexIndex).outgoingNeighbours=[...
                obj.vertices(vertexIndex).outgoingNeighbours(arrayleft),...
                obj.edges(edgeIndex).endVertex,...
                obj.vertices(vertexIndex).outgoingNeighbours(arrayright)];
            
            %update degree of vertex
            obj.updateDegree(vertexIndex);
            
            %rearrange crease vectors and dihedral angles
            obj.vertices(vertexIndex).dihedralAngles=[...
                obj.vertices(vertexIndex).dihedralAngles(1:numberOfIncomingEdges),...
                obj.vertices(vertexIndex).dihedralAngles(arrayleft+numberOfIncomingEdges),...
                obj.edges(edgeIndex).dihedralAngle,...
                obj.vertices(vertexIndex).dihedralAngles(arrayright+numberOfIncomingEdges)];
            
            obj.vertices(vertexIndex).creaseVectors=[...
                obj.vertices(vertexIndex).creaseVectors(1:numberOfIncomingEdges),...
                obj.vertices(vertexIndex).creaseVectors(arrayleft+numberOfIncomingEdges),...
                obj.edges(edgeIndex).directionVector,...
                obj.vertices(vertexIndex).creaseVectors(arrayright+numberOfIncomingEdges)];
            
            %update type of vertex
            if strcmp(obj.vertices(vertexIndex).type,'boundary') && ...
                    not(strcmp(obj.edges(edgeIndex).type,'boundary')) && ...
                    size(obj.vertices(vertexIndex).incomingNeighbours,2)>=1
                obj.vertices(vertexIndex).type='interior';
            end
            
        end     %function addOutgoingCrease
        function [] = removeOutgoingCrease(obj,vertexIndex,edgeIndex)
            %REMOVEOUTGOINGCREASE Removes crease with index edgeIndex
            %(last/third/second input argument) from the container lists of
            %vertex with index vertexIndex (second-last input argument).
            %Does not change edge indices. Does not delete the crease from
            %the list or change the attributes of the crease itself. Set
            %Edge to deleted (deleted=true) before using this function.
            
            %Check if edge was set to deleted already and throw error
            %otherwise.
            if not(obj.edges(edgeIndex).deleted)
                error('Set Edge to deleted before using this function');
            end
            
            %Find out where among the outgoing edges the "edge to delete"
            %is.
            index=find(obj.vertices(vertexIndex).outgoingEdges==edgeIndex);
            
            %Check whether result is meaningful
            if size(index,2)~=1 %Gives error if the "index" wasn't found at all or several times. In the latter case, something went terribly wrong.
                error('Could not identify the location of the specified outgoing edge on the specified vertex');
            end
            
            %Get the number of (currently) outgoing and incoming edges
            numberOfOutgoingCreases=length(obj.vertices(vertexIndex).outgoingEdges);
            numberOfIncomingCreases=length(obj.vertices(vertexIndex).incomingEdges);
            
            %Array with local indices of outgoing edges before the one to
            %be deleted
            arrayleft=1:(index-1);
            
            %Array with local indices of outgoing edges after the one to be
            %deleted
            arrayright=(index+1):numberOfOutgoingCreases;
            
            %Update list of "child vertices"
            obj.vertices(vertexIndex).outgoingNeighbours=[...
                obj.vertices(vertexIndex).outgoingNeighbours(arrayleft),...
                obj.vertices(vertexIndex).outgoingNeighbours(arrayright)];
            
            %Update indexAtSource for the edges registered in arrayright
            for j=1:size(arrayright,2)
                tempEdgeIndex=obj.vertices(vertexIndex).outgoingEdges(arrayright(j));
                obj.edges(tempEdgeIndex).indexAtSource=...
                    obj.edges(tempEdgeIndex).indexAtSource-1;
            end
            
            %Now re-write the vector of outgoing edges. From here on, the
            %indexation used above is no longer valid
            obj.vertices(vertexIndex).outgoingEdges=[...
                obj.vertices(vertexIndex).outgoingEdges(arrayleft),...
                obj.vertices(vertexIndex).outgoingEdges(arrayright)];
            
            %Rearrange the dihedral angles and the crease vectors
            obj.vertices(vertexIndex).dihedralAngles=[...
                obj.vertices(vertexIndex).dihedralAngles(1:numberOfIncomingCreases),...
                obj.vertices(vertexIndex).dihedralAngles(numberOfIncomingCreases+arrayleft),...
                obj.vertices(vertexIndex).dihedralAngles(numberOfIncomingCreases+arrayright)];
            
            obj.vertices(vertexIndex).creaseVectors=[...
                obj.vertices(vertexIndex).creaseVectors(1:numberOfIncomingCreases),...
                obj.vertices(vertexIndex).creaseVectors(arrayleft+numberOfIncomingCreases),...
                obj.vertices(vertexIndex).creaseVectors(arrayright+numberOfIncomingCreases)];
            
            %Calculating new number of outgoing Edges
            numberOfOutgoingCreasesNew=size(obj.vertices(vertexIndex).outgoingEdges,2);
            
            %Verifying that the deletion was successful
            if numberOfOutgoingCreasesNew~=numberOfOutgoingCreases-1
                error('Deletion of outgoing edge not successful');
            end
            
            %Update degree of vertex
            obj.updateDegree(vertexIndex);
            
            %Possibly the type needs to be updated too
            if strcmp(obj.vertices(vertexIndex).type,'interior') && ...
                    numberOfOutgoingCreasesNew==0
                obj.vertices(vertexIndex).type='boundary';
            end
            
        end
        
        function [] = addIncomingCrease(obj,vertexIndex,edgeIndex)
            %ADDINCOMINGCREASE Adapts the involved vertex, edge and face
            %objects to consider the new incoming vertex.
            
            %Check if the new crease has a near-0 angle to an existing edge
            obj.checkIsNewCreaseAllowed(vertexIndex, edgeIndex);
            
            %Find out how the situation at the vertex looks like
            numberOfIncomingEdges=length(obj.vertices(vertexIndex).incomingEdges);
            numberOfOutgoingEdges=length(obj.vertices(vertexIndex).outgoingEdges);
            deg=size([obj.vertices(vertexIndex).incomingEdges, ...
                obj.vertices(vertexIndex).outgoingEdges],2);
            corrDeg=deg;
            
            %Find the edge to the left and to the right of the new incoming
            %edge
            
            [left, right, jl, jr] = obj.findLeftAndRight(vertexIndex,edgeIndex);
            jr=max([1 jr]); %Correction for "new" vertices where 0 is a possible (but undesired) output.
            
            %Check if incoming creases are added to a vertex with no
            %outputs and whether the new creases are not "within existing
            %ones"
            criterion1=numberOfIncomingEdges>0;
            criterion2=deg==numberOfIncomingEdges; %no outgoing edges
            criterion3=(jr==1) || (jl==numberOfIncomingEdges);
            
            if all([criterion1, criterion2, criterion3])
                %The findLeftAndRight function, in this case, can't
                %properly determine whether to add the crease to the left
                %or to the right of the existing crease(s). We therefore
                %use this function to adapt accordingly.
                [jl, jr] = obj.decideLeftRight(vertexIndex,edgeIndex);
                if jl~=0
                    corrDeg=deg+1;  %for the construction of arrayleft. 
                    %If this is not incremented, arrayleft will stay empty 
                    %but should contain everything.
                end
            end
            
            %Quick check to ensure that the new incoming vector is within
            %the allowable sector.
            forbiddenLefts=obj.vertices(vertexIndex).outgoingEdges(1:end-1);
            if ismember(left,forbiddenLefts)
                error('New incoming vector appears to be coming in between outgoing edges. This is not allowed.');
            end
            
            forbiddenRights=obj.vertices(vertexIndex).outgoingEdges(2:end);
            if ismember(right,forbiddenRights)
                error('New incoming vector appears to be coming in between outgoing edges. This is not allowed.');
            end
            
            %Get array with indices of incoming edges before the new one
            arrayleft=1:(mod(jl,corrDeg));  %stays empty if "left" is the last outgoing edge.
            
            
            %Get array with indices of incoming edges after the new one
            arrayright=jr:numberOfIncomingEdges;    %stays empty if "right" is the first outgoing edge.
            
            %Update indexAtEnd for all incoming edges after the new one (and the new one)
            for j=1:size(arrayright,2)
                edgeIndextemp=obj.vertices(vertexIndex).incomingEdges(arrayright(j));
                obj.edges(edgeIndextemp).indexAtEnd=...
                    obj.edges(edgeIndextemp).indexAtEnd+1;
            end
            
            %Assign the correct value to the new edge
            obj.edges(edgeIndex).indexAtEnd=jr;
            
            %Since an incoming crease has been added, also the
            %indexAtSource entries of the outgoing edges have to be
            %incremented...
            for j=1:size(obj.vertices(vertexIndex).outgoingEdges,2)
                edgeIndextemp=obj.vertices(vertexIndex).outgoingEdges(j);
                obj.edges(edgeIndextemp).indexAtSource=...
                    obj.edges(edgeIndextemp).indexAtSource+1;
            end
            
            %IMPORTANT: Again, update the edge indices before updating the
            %vertex indices.
            
            %Re-write vector of incoming edges to the vertex
            obj.vertices(vertexIndex).incomingEdges=[...
                obj.vertices(vertexIndex).incomingEdges(arrayleft),...
                edgeIndex,...
                obj.vertices(vertexIndex).incomingEdges(arrayright)];
            
            %Update the list of incoming Neighbours to the vertex
            obj.vertices(vertexIndex).incomingNeighbours=[...
                obj.vertices(vertexIndex).incomingNeighbours(arrayleft),...
                obj.edges(edgeIndex).sourceVertex,...
                obj.vertices(vertexIndex).incomingNeighbours(arrayright)];
            
            %update the list of dihedral angles and crease vectors
            obj.vertices(vertexIndex).dihedralAngles=[...
                obj.vertices(vertexIndex).dihedralAngles(arrayleft),...
                obj.edges(edgeIndex).dihedralAngle,...
                obj.vertices(vertexIndex).dihedralAngles(arrayright),...
                obj.vertices(vertexIndex).dihedralAngles(numberOfIncomingEdges+1:end)];
            
            %and the crease vectors...
            obj.vertices(vertexIndex).creaseVectors=[...
                obj.vertices(vertexIndex).creaseVectors(arrayleft),...
                obj.edges(edgeIndex).directionVector,...
                obj.vertices(vertexIndex).creaseVectors(arrayright),...
                obj.vertices(vertexIndex).creaseVectors(numberOfIncomingEdges+1:end)];
            
            %updating the target vertex degree
            obj.updateDegree(vertexIndex);
        end %function addIncomingEdge   
        function [] = removeIncomingCrease(obj, vertexIndex, edgeIndex)
            %REMOVEINCOMINGCREASE Removes incoming crease with index
            %edgeIndex (last/third/second input argument) from lists
            %associated to vertex with index vertexIndex (second-last input
            %argument). Does not delete the edge from the edge container in
            %obj. Does not change any global indices. Set edge to "deleted"
            %(deleted=true) before using this function.
            
            %Check whether edge was set to deleted already and throw error
            %otherwise
            if not(obj.edges(edgeIndex).deleted)
                error('Set edge to be deleted before using this function.');
            end
            
            %Check if the crease is removed to stay deleted or just to move
            %a vertex
            
            %The first item is this function itself; the second one would
            %be removeCreaseByIndex
            moving=false;
            callingstack=dbstack();
            if length(callingstack)>=3
                callingFcn=callingstack(3).name;
                if strcmp(callingFcn,'origami.moveVertexTo')
                    moving=true;
                end
            end
            
            %check if an auxiliary boundary vertex is being added
            splitBoundary=false;
            if length(callingstack)>=2
                callingFcn=callingstack(2).name;
                if strcmp(callingFcn,'origami.splitBoundaryEdge')
                    splitBoundary=true;
                end
            end
            
            %Find out where in the vector of incoming creases the edge to
            %be deleted is located.
            index=find(obj.vertices(vertexIndex).incomingEdges==edgeIndex);
            
            %Check whether result is meaningful and throw error otherwise
            if size(index,2)~=1     %Throws an error if the edge wasn't 
                %found at all or was found several times. In the latter
                %case, something went terribly wrong.
                error('Could not locate specified edge among incoming edges.');
            end
            
            %Get number of incoming creases
            numberOfIncomingCreases=size(obj.vertices(vertexIndex).incomingEdges,2);
            
            %Array with local indices of incoming creases "before" the
            %crease to be deleted
            arrayleft=1:(index-1);
            
            %Array with local indices of incoming creases "after" the
            %crease to be deleted
            arrayright=(index+1):numberOfIncomingCreases;
            
            %Check whether removing the edge results in a "disconnected"
            %graph. I.e. it is not allowed to have vertices with one or
            %several outputs but no inputs.
            if size([arrayleft, arrayright],2)==0 &&...
                    size(obj.vertices(vertexIndex).outgoingEdges,2)~=0 &&...
                    not(moving) && not(splitBoundary)
                error('Removing the last incoming edge to a vertex with one or several outputs. This is not allowed.');
            end
            
            %Update list of "parent vertices"
            obj.vertices(vertexIndex).incomingNeighbours=[...
                obj.vertices(vertexIndex).incomingNeighbours(arrayleft),...
                obj.vertices(vertexIndex).incomingNeighbours(arrayright)];
            
            %Update the "indexAtEnd" entry of the edges coming in after the
            %edge to be deleted
            for j=1:size(arrayright,2)
                tempEdgeIndex=obj.vertices(vertexIndex).incomingEdges(arrayright(j));
                obj.edges(tempEdgeIndex).indexAtEnd=...
                    obj.edges(tempEdgeIndex).indexAtEnd-1;
            end
            
            %Update the "indexAtSource" entry of the outgoing edges
            for j=1:size(obj.vertices(vertexIndex).outgoingEdges,2)
                tempEdgeIndex=obj.vertices(vertexIndex).outgoingEdges(j);
                obj.edges(tempEdgeIndex).indexAtSource=...
                    obj.edges(tempEdgeIndex).indexAtSource-1;
            end
            
            %Now re-write vector of incoming edges. The above-used
            %indexation is not valid any more from here on.
            obj.vertices(vertexIndex).incomingEdges=[...
                obj.vertices(vertexIndex).incomingEdges(arrayleft),...
                obj.vertices(vertexIndex).incomingEdges(arrayright)];
            
            %update the dihedral angles and the crease vectors
            obj.vertices(vertexIndex).dihedralAngles=[...
                obj.vertices(vertexIndex).dihedralAngles(arrayleft),...
                obj.vertices(vertexIndex).dihedralAngles(arrayright),...
                obj.vertices(vertexIndex).dihedralAngles(numberOfIncomingCreases+1:end)];
            
            obj.vertices(vertexIndex).creaseVectors=[...
                obj.vertices(vertexIndex).creaseVectors(arrayleft),...
                obj.vertices(vertexIndex).creaseVectors(arrayright),...
                obj.vertices(vertexIndex).creaseVectors(numberOfIncomingCreases+1:end)];
            
            %Calculate new number of incoming creases
            numberOfIncomingCreasesNew=size(obj.vertices(vertexIndex).incomingEdges,2);
            
            %Check if deletion was successful
            if numberOfIncomingCreasesNew~=numberOfIncomingCreases-1
                error('Deletion of incoming edge not successful');
            end
            
            %Update the vertex degree
            obj.updateDegree(vertexIndex);
            
        end

        %Create vectors in edge directions (folded)
        function edgeVec = awayVector(obj, vertexIndex, edgeIndex)
            %AWAYVECTOR(OBJ, VERTEXINDEX, EDGEINDEX) Returns a vector
            %pointing along edge with index edgeIndex, away from vertex
            %with index vertexIndex. Returns an error if edgeIndex is not
            %incident to vertexIndex.
            %Depending on situation at vertex, either targetSourceVec or
            %sourceTargetVec is called.
            
            isSource=false;
            isEnd=false;
            
            if obj.edges(edgeIndex).sourceVertex==vertexIndex
                isSource=true;
            end
            if obj.edges(edgeIndex).endVertex==vertexIndex
                isEnd=true;
            end
            
            if not(isSource) && not(isEnd)
                error('Edge does not seem to be connected to Vertex');
            end
            
            if isSource && isEnd
                error('Edge seems to connect vertex with itself. This is not expected.');
            end
            
            if isSource
                edgeVec=obj.sourceTargetVec(edgeIndex);
            end
            if isEnd
                edgeVec=obj.targetSourceVec(edgeIndex);
            end
            
        end
        function edgeVec = sourceTargetVec(obj, edgeIndex,opts)
            %SOURCETARGETVEC Gives a unit vector pointing in the direction
            %of the edge, from source to end
            
            sourceIndex=obj.edges(edgeIndex).sourceVertex;
            targetIndex=obj.edges(edgeIndex).endVertex;
            
            edgeVec=obj.vertices(targetIndex).position - ...
                obj.vertices(sourceIndex).position;
            
            edgeVec=unitVec(edgeVec);
            
            if nargin==3
                if strcmp(opts,'2D')
                    edgeVec=edgeVec(1:2);   %The used demanded a 2D-vector
                end
            end
        end     
        function edgeVec = targetSourceVec(obj, edgeIndex,opts)
            %TARGETSOURCEVEC Gives a unit vector pointing in the direction
            %of the edge, from end to source
            
            if nargin==3
                edgeVec=-sourceTargetVec(obj,edgeIndex,opts);
            else
                edgeVec=-sourceTargetVec(obj,edgeIndex);
            end
        end
        %Create vectors in edge directions (folded)
        function edgeVec = targetSourceVec0(obj, edgeIndex,opts)
            %TARGETSOURCEVEC0 Gives a unit vector pointing in the direction
            %of the edge, from end to source
            
            if nargin==3
                edgeVec=-sourceTargetVec0(obj,edgeIndex,opts);
            else
                edgeVec=-sourceTargetVec0(obj,edgeIndex);
            end
        end
        function edgeVec = sourceTargetVec0(obj, edgeIndex,opts)
            %SOURCETARGETVEC Gives a unit vector pointing in the direction
            %of the edge, from source to end, in the unfolded configuration
            
            sourceIndex=obj.edges(edgeIndex).sourceVertex;
            targetIndex=obj.edges(edgeIndex).endVertex;
            
            edgeVec=obj.vertices(targetIndex).position0 - ...
                obj.vertices(sourceIndex).position0;
            
            edgeVec=unitVec(edgeVec);
            
            if nargin==3
                if strcmp(opts,'2D')
                    edgeVec=unitVec(edgeVec(1:2));   %The used demanded a 2D-vector
                end
            end
        end
        function edgeVec = awayVector0(obj, vertexIndex, edgeIndex)
            %AWAYVECTOR0(OBJ, VERTEXINDEX, EDGEINDEX) Returns a vector
            %pointing along edge with index edgeIndex, away from vertex
            %with index vertexIndex. Returns an error if edgeIndex is not
            %incident to vertexIndex.
            %Depending on situation at vertex, either targetSourceVec or
            %sourceTargetVec is called.
            %unlike awayVector, this function operates in the unfolded
            %configuration
            
            isSource=false;
            isEnd=false;
            
            if obj.edges(edgeIndex).sourceVertex==vertexIndex
                isSource=true;
            end
            if obj.edges(edgeIndex).endVertex==vertexIndex
                isEnd=true;
            end
            
            if not(isSource) && not(isEnd)
                error('Edge does not seem to be connected to Vertex');
            end
            
            if isSource && isEnd
                error('Edge seems to connect vertex with itself. This is not expected.');
            end
            
            if isSource
                edgeVec=obj.sourceTargetVec0(edgeIndex);
            end
            if isEnd
                edgeVec=obj.targetSourceVec0(edgeIndex);
            end
            
        end
        %Find edges left and right of current edge
        function [left, right, internalIndexLeft, internalIndexRight] = ...
                findLeftAndRight(obj,vertexIndex,edgeIndex)
            %FINDLEFTANDRIGHT Given an existing vertex, returns the edge
            %index of the next edge in clockwise direction (left) and in
            %counterclockwise direction (right).
            
            %Check if the edge was already registered at the central vertex
            registeredInputs=obj.vertices(vertexIndex).incomingEdges;
            registeredOutputs=obj.vertices(vertexIndex).outgoingEdges;
            alreadyregistered=ismember(edgeIndex,[registeredInputs, ...
                registeredOutputs]);
            
            deg=obj.vertices(vertexIndex).degree;
            if alreadyregistered
                deg=deg-1;
            end
            
            %normal vector for the computation of sector angles
            zAxis=obj.vertices(vertexIndex).unfoldedZ;
            
            %countedDeg specifies the degree if boundary edges are
            %considered as well
            countedDeg=size([obj.vertices(vertexIndex).incomingEdges, ...
                obj.vertices(vertexIndex).outgoingEdges],2);
            if alreadyregistered
                countedDeg=countedDeg-1;
            end
            
            %Preparing empty containers
            
            creaseIndices=zeros(1,countedDeg);
            creaseVecs=zeros(3,countedDeg);
            
            counter=0;
            controlCounter=0;
            
            %Gathering vectors and indices from incoming edges.
            for j=1:size(registeredInputs,2)    %iterating through already registerd input edges
                
                if registeredInputs(j)~=edgeIndex   %Making sure not to register the "input" edge twice
                    
                    counter=counter+1;  %increment counter
                    if not(strcmp(obj.edges(registeredInputs(j)).type,'boundary'))
                        controlCounter=controlCounter+1;
                    end
                    
                    %add this edge to the "collection"
                    creaseIndices(counter)=registeredInputs(j);
                    creaseVecs(:,counter)=obj.targetSourceVec0(...
                        registeredInputs(j)); %incoming vector -> get vector pointing away from vertex
                    
                end
                
            end
            
            %Gathering vectors and indices from outgoing edges
            for j=1:size(registeredOutputs,2)   %Iterating through already registered outgoing edges
                
                if registeredOutputs(j)~=edgeIndex %
                    
                    if not(strcmp(obj.edges(registeredOutputs(j)).type,'boundary'))
                        controlCounter=controlCounter+1;
                    end
                    counter=counter+1;  %increment counter
                    creaseIndices(counter)=registeredOutputs(j);
                    creaseVecs(:,counter)=obj.sourceTargetVec0(...
                        registeredOutputs(j)); %outgoing edge -> vector pointing away from vertex is sourceTarget.
                    
                end
                
            end
            
            %At this point, deg and counter should have the same value
            if controlCounter~=deg || counter~=countedDeg
                error('number of adjacent edges not equal to degree of vertex');
            end
            
            if counter==0
                %Nothing has been registered so far. Assigning zeros to all
                %output vars to indicate this.
                left=0;
                right=0;
                internalIndexLeft=0;
                internalIndexRight=0;
            end
            
            %Calculating thisVector, a vector pointing away from the
            %current vector in the direction of the crease in question
            thisVector=obj.awayVector0(vertexIndex,edgeIndex);
            
            for j=1:counter     %iterating through all adjacent edges
                
                potLeft=j;
                potRight=mod(j,countedDeg)+1;
                
                leftVec=creaseVecs(:,potLeft);
                rightVec=creaseVecs(:,potRight);
                
                %Check if this sector is "the right one"
                angle1=counterClockWiseAngle(leftVec,thisVector,zAxis);
                angle2=counterClockWiseAngle(leftVec,rightVec,zAxis);
                
                if angle1<angle2
                    left=creaseIndices(potLeft);
                    internalIndexLeft=potLeft;
                    right=creaseIndices(potRight);
                    internalIndexRight=potRight;
                    break;
                end
    
            end %for-loop
            
        end %function findleftandright
        function [internalIndexLeft, internalIndexRight] = ...
                decideLeftRight(obj,vertexIndex,edgeIndex)
            %DECIDELEFTRIGHT Is used to handle the special case of a vertex
            %with several inputs but no output or several outputs and no
            %input
            
            %z-axis for the calculation of sector angles
            zAxis=obj.vertices(vertexIndex).unfoldedZ;
            
            if size(obj.vertices(vertexIndex).outgoingEdges,2)==0
                %here, we treat the case of ordering incoming edges
                
                
                %First fetch all angles relevant for this analysis
                thisVector=obj.targetSourceVec0(edgeIndex);
                
                %Get the existing incoming vectors
                firstInputVec=obj.targetSourceVec0(...
                    obj.vertices(vertexIndex).incomingEdges(1));
                lastInputVec=obj.targetSourceVec0(...
                    obj.vertices(vertexIndex).incomingEdges(end));
                
                %Calculate the angle between new edge and first incoming edge
                ang1=counterClockWiseAngle(thisVector,firstInputVec,zAxis);
                
                %...and the angle between the new edge and the last incoming
                %edge
                ang2=counterClockWiseAngle(lastInputVec,thisVector,zAxis);
                
                if ang1<ang2
                    %new edge is closer to the first incoming vector ->
                    %it's the new first one
                    internalIndexLeft=0;
                    internalIndexRight=1;
                else
                    %new edge is closer to the last incoming edge -> it's
                    %the new last one
                    internalIndexLeft=size(obj.vertices(vertexIndex).incomingEdges,2);
                    internalIndexRight=internalIndexLeft+1;
                end
                
            elseif size(obj.vertices(vertexIndex).incomingEdges,2)==0
                %Here, we're dealing with ordering edges around a source
                %vertex
                
                %Treat the case of a source vertex where the last outgoing
                %boundary edge is added.
                if strcmp(obj.edges(edgeIndex).type,'boundary') && ...
                        (strcmp(obj.edges(obj.vertices(vertexIndex).outgoingEdges(1)).type,'boundary') || ...
                        strcmp(obj.edges(obj.vertices(vertexIndex).outgoingEdges(end)).type,'boundary'))
                    
                    %A boundary edge is added to a vertex that
                    %already has a boundary edge registered.
                    %This can apply to source vertices or local sources.
                    
                    if strcmp(obj.edges(obj.vertices(vertexIndex).outgoingEdges(1)).type,'boundary')
                        %first outgoing crease is a boundary crease. We are
                        %adding another outgoing crease -> must come at the
                        %end
                        internalIndexLeft=size(obj.vertices(vertexIndex).outgoingEdges,2);
                        internalIndexRight=internalIndexLeft+1;
                        return;
                    elseif strcmp(obj.edges(obj.vertices(vertexIndex).outgoingEdges(end)).type,'boundary')
                        %the last registered incoming one is already a
                        %boundary edge, so this one will be the first one
                        internalIndexLeft=0;
                        internalIndexRight=1;
                        return;
                    end
                        
                end
                
                %First fetch all angles relevant for this analysis
                thisVector=obj.sourceTargetVec0(edgeIndex);
                
                %Get the existing outgoing vectors
                firstOutputVec=obj.sourceTargetVec0(...
                    obj.vertices(vertexIndex).outgoingEdges(1));
                lastOutputVec=obj.sourceTargetVec(...
                    obj.vertices(vertexIndex).outgoingEdges(end));
                
                %Calculate the angle between new edge and first outgoing edge
                ang1=counterClockWiseAngle(thisVector,firstOutputVec,zAxis);
                
                %...and the angle between the new edge and the last
                %outgoing edge
                ang2=counterClockWiseAngle(lastOutputVec,thisVector,zAxis);
                
                if ang1<ang2
                    internalIndexLeft=0;
                    internalIndexRight=1;
                else
                    internalIndexLeft=size(obj.vertices(vertexIndex).outgoingEdges,2);
                    internalIndexRight=internalIndexLeft+1;
                end
                
            else
                error('Vertex has incoming and outgoing edges. Only use this function if it has only incoming ones or only outgoing ones');
            end
            
        end
        %Update degree of vertex
        function [newDeg] = updateDegree(obj, vertexIndex)
            %UPDATEDEGREE(OBJ, VERTEXINDEX) Updates the degree of
            %vertexIndex by looking at its incident edges.
            
            %Create a temporary vector with all incident edges
            tempVec=[obj.vertices(vertexIndex).incomingEdges, ...
                obj.vertices(vertexIndex).outgoingEdges];
            
            deg=0;
            
            %Now get the degree by counting the incident interior edges
            for j=1:size(tempVec,2)
                if not(strcmp(obj.edges(tempVec(j)).type,'boundary'))
                    deg=deg+1;
                end
            end
            
            obj.vertices(vertexIndex).degree=deg;
            
            newDeg=deg;
            
        end
        
        %%%%%Auxiliary Functions for Creation of Faces%%%%%
        %(unreviewed)
        function verdict = isOnBoundary(obj,vertexIndex)
            %VERDICT=ISONBOUNDARY(OBJ,VERTEXINDEX) Returns true if all of
            %the edges incident to the vertex are boundary edges.
            
            %verify the input
            obj.verifyVertexIndex(vertexIndex);
            
            %get the adjacent edges
            adjEdg=[obj.vertices(vertexIndex).incomingEdges,...
                obj.vertices(vertexIndex).outgoingEdges];
            
            if length(adjEdg)>2
                %if the edge has more than two incident edges, it can't be
                %on the boundary 
                verdict=false;
                return;
            end
            
            %here it is verified that there are two adjacent edges
            for j=adjEdg
                if not(strcmp(obj.edges(j).type,'boundary'))
                    verdict=false;
                    return;
                end
            end
            
            %if control reaches this point, the vertex must be on the
            %boundary
            verdict=true;
            
        end
        function boundaryVertices = verticesOnBoundary(obj)
            %BOUNDARYVERTICES = VERTICESONBOUNDARY(OBJ) returns all
            %vertices for which isOnBoundary is true. 
            
            boundaryVertices=false(1,obj.numberOfVertices());
            for v=obj.vertices
                %remember all vertices that are on the boundary
                if v.deleted
                    continue;
                end
                
                if obj.isOnBoundary(v.index)
                    boundaryVertices(v.index)=true;
                end
            end
            
            %only return those that turned out to be on the boundary
            allVertices=1:obj.numberOfVertices();
            boundaryVertices=allVertices(boundaryVertices);
            
        end
        function [] = defineFaceFromVertexAndEdge(obj, startVertex, startEdge)
            %DEFINEFACEFROMVERTEXANDEDGE(OBJ, STARTVERTEX, STARTEDGE)
            %starts at startvertex and adds the face located next to
            %startEdge in counterclockwise direction (rotating around
            %startVertex).
            %If the face turns out to be a boundary face, an auxiliary
            %vertex is added to close the face.
            
            %Get the counterclockwise marching and check whether it ends at
            %a boundary node
            [verticesCCW, edgesCCW, isLeftsCCW, isRightsCCW, faceIndicesCCW] = ...
                obj.cCWMarching(startVertex, startEdge);
            isBoundaryCCW=isBoundary(edgesCCW);
            
            %Get the clockwise marching and check whether it ends at a
            %boundary node
            [verticesCW, edgesCW, isLeftsCW, isRightsCW, faceIndicesCW] = ...
               obj.cWMarching(startVertex, startEdge);
           isBoundaryCW=isBoundary(edgesCW);
            
            %If one marching ends at a boundary node and the other one
            %doesn't, throw an error
            
            %NOTE FROM JIM: isBoundary doesnt identify if the last node in
            %particular has the type 'boundary'. It identifies if the last
            %entry to edgesCW or edgesCCW is zero. 
            %Therefore I think it can
            %also identify source vertices as boundary vertices????
            
            if not((isBoundaryCCW&&isBoundaryCW) || ...
                    (not(isBoundaryCCW)&&not(isBoundaryCW)))
                error('Contradictory Results from marching in counter- and -clockwise direction');
            end
            
            %If both end at boundary nodes, add an auxiliary node and
            %assemble the containers accordingly.
            if isBoundaryCCW && isBoundaryCW
                
               
                
                endVertexCCW=verticesCCW(end);
                endVertexCW=verticesCW(end);
                
                %Construct the auxiliary vertex
                [auxVertexIndex, edgeIndex1, edgeIndex2] = ...
                    obj.addAuxiliaryVertex(endVertexCCW, endVertexCW);
                
                %consider the auxiliary vertex when assembling the list of
                %connected vertices
                faceVertices=[verticesCW, auxVertexIndex, fliplr(verticesCCW)];
                
                %Consider that the edge-sets end on 0 (omit) and that two
                %boundary edges were added
                faceEdges=[edgesCW(1:end-1), edgeIndex2, edgeIndex1, ...
                    fliplr(edgesCCW(1:end-1))];
                
                %The clockwise marching output contains the output for the
                %start vertex (omit accordingly for the counterclockwise
                %result).
                faceLefts=[isLeftsCW, false, true, fliplr(isLeftsCCW(2:end))];
                faceRights=[isRightsCW, true, false, fliplr(isRightsCCW(2:end))];
                
                %The indices at the vertices work similar to the vertices
                %themselves. By definition, the face (only one present) at
                %the auxiliary vertices always gets index 1.
                
                faceIndices=[faceIndicesCW, 1, fliplr(faceIndicesCCW)];
                
            elseif not(isBoundaryCW) && not(isBoundaryCCW)
                
                
                
            %Else (if neither ends at a boundary node), check that the
            %marchings reveal identical containers. Use result from
            %clockwise marching to construct containers
            
                if not(all(verticesCCW==fliplr(verticesCW)))
                    error('The two marchings both resulted in a closed path, but a different one');
                end
                
                %The containers can simply be taken from the
                %clockwise marching - result.
                faceVertices=verticesCW;
                faceEdges=edgesCW(1:end-1);
                faceLefts=isLeftsCW;
                faceRights=isRightsCW;
                faceIndices=faceIndicesCW;
            
            else
                %One of the marchings goes into the open and the other one
                %goes "around", i.e. closes itself
                error('Contradictory Results from marching in counter- and -clockwise direction');
            end
            
            
            %Use face-constructor to construct a new face
            faceIndex=obj.numberOfFaces() + 1;
            obj.faces=[obj.faces face(faceVertices, ...
                faceEdges, faceIndices, faceIndex)];
            
            %Loop over involved vertices and register the new face
            
            for j=1:size(faceVertices,2)
                vertexIndex=faceVertices(j);
                obj.vertices(vertexIndex).adjacentFaces(faceIndices(j))=faceIndex;
                obj.vertices(vertexIndex).normalVectors(faceIndices(j))=...
                    obj.faces(faceIndex).normalVector;
            end
            
            %Loop over all involved edges and register the new face
            
            for j=1:size(faceEdges,2)
                edgeIndex=faceEdges(j);
                
                if not(faceLefts(j)) && not(faceRights(j))
                    error('Both faceLeft and faceRight are true. This is invalid.');
                end
                
                if faceLefts(j)
                    obj.edges(edgeIndex).faceLeft=faceIndex;
                    obj.edges(edgeIndex).normalVecLeft=obj.faces(faceIndex).normalVector;
                elseif faceRights(j)
                    obj.edges(edgeIndex).faceRight=faceIndex;
                    obj.edges(edgeIndex).normalVecRight=obj.faces(faceIndex).normalVector;
                else
                    error('Neither isLeft nor isRight is true. This is invalid.');
                end
            end
            
        end
        function [vertexIndex, edgeIndex1, edgeIndex2] = ... 
                addAuxiliaryVertex(obj, vertex1, vertex2)
            %ADDAUXILIARYVERTEX(OBJ, VERTEX1, VERTEX2) Adds a new auxiliary
            %vertex halfway between vertices vertex1 and vertex2 and
            %connects the three vertices. The purpose of such vertices is
            %the definition of the boundary facet. These vertices, or, more
            %precisely, their positions, can be used to define the
            %geometries of the boundary facets.
            %[vertexIndex, edgeIndex1, edgeIndex2] =
            %addAuxiliaryVertex(obj, vertex1, vertex2) returns the index of
            %the auxiliary vertex as well as the index of the edge between
            %vertex1 and the auxiliary vertex (edgeIndex1) and the index of
            %the edge between vertex2 and the auxiliary Vertex.
            %Situation Sketch:
            %
            % |
            % O vertex1
            %  \
            %   \ edgeIndex1
            %    \
            %     _|
            %       O_vertexIndex
            %       |\
            %         \
            %          \ edgeIndex2
            %           \ 
            %            O vertex2 -----------
            
            %Run initial checks: are the two vertices not identical and are
            %they not at the same position?
            
            if vertex1==vertex2
                error('It looks like the clockwise and counterclockwise marching procedures terminated at the same vertex. Check the known issues section of the manual.');
            end
            
            if norm(obj.vertices(vertex1).position0-...
                    obj.vertices(vertex2).position0,2)<1e-9
                error('The two specified vertices for the construction of an auxiliary vertex are too close to each other. Try merging them instead.');
            end
            
                
            %Check whether both vertices are already marked as boundary
            %vertices or source vertex.
            
            if not(strcmp(obj.vertices(vertex1).type,'boundary') ...
                    || strcmp(obj.vertices(vertex1).type,'source'))
                error('Only boundary vertices and the source vertex can be used as basis to construct an auxiliary vertex.');
            end
            
            if not(strcmp(obj.vertices(vertex2).type,'boundary') ...
                    || strcmp(obj.vertices(vertex2).type,'source'))
                error('Only boundary vertices and the source vertex can be used as basis to construct an auxiliary vertex.');
            end
            
            %Consider the following: Maybe the two are already connected.
            %In that case, we just return the usual data and create
            %nothing.
            out1=obj.vertices(vertex1).outgoingNeighbours;
            out2=obj.vertices(vertex2).outgoingNeighbours;
                        
            %%search for common members         
           
            [isCommon,loc]=ismember(out1,out2);
            if sum(isCommon)==1 && ...
                    not(strcmp(obj.vertices(vertex1).type,'source') &&...
                    strcmp(obj.vertices(vertex2).type,'source'))
                %ADDED BY JIM (!) above:
                % This creates a bug where the two source vertices are not
                % connected. Therefore an excemption is made in the if
                % statement
               
                %a common entry was found. It is...
                vertexIndex=out2(loc(isCommon));
                edgeIndex2=obj.vertices(vertex2).outgoingEdges(loc(isCommon));
                localIndexAt1=find(obj.vertices(vertex1).outgoingNeighbours==vertexIndex);
                edgeIndex1=obj.vertices(vertex1).outgoingEdges(localIndexAt1);
                return;
            elseif sum(isCommon)>1
                error('More than one connecting vertex found between vertices %d and %d.',...
                    vertex1,vertex2);
            end
                        
            %Reading out the vertex types at the two outgoing vertices
            
            type1=obj.vertices(vertex1).type;
            type2=obj.vertices(vertex2).type;
            
            %Calculate position and index of new vertex
            %CHANGED BY JIM (!)
            %old:
            newPos=0.5*(obj.vertices(vertex1).position0 +...
                obj.vertices(vertex2).position0);
            %new:
%              %Check if face added left or right
%             
%             
%              %define vector from first to last vertex
%              vect = obj.vertices(vertex2).position0 -...
%                  obj.vertices(vertex1).position0;
%              %rot. matrix 90 deg (on left of edge) and -90 deg (on right of
%              %edge), and scale by 1/2
%              Rl = [0 -1/2; 1/2 0];
%              Rr = [0 1/2; -1/2 0];
%              %coords of midpoint
%              m = 0.5*(obj.vertices(vertex1).position0 +...
%              obj.vertices(vertex2).position0);
%             %if left or right: (Added by jim)
%             if     obj.edges(edgeIndex).faceLeft == 0
%                 newPos = m + Rl*vect;
%             elseif obj.edges(edgeIndex).faceRight == 0
%                 newPos = m + Rr*vect;
%             else
%                error("ah oh something went wrong")
%             end
                        
            
            vertexIndex=obj.numberOfVertices()+1;
            
            zAxis=obj.vertices(vertex1).unfoldedZ;
            
            %Create new vertex
            obj.vertices=[obj.vertices vertex(newPos, ...
                vertexIndex, false, 'aux',zAxis)];
            
            %Create new edges from vertex 1 and 2 to new vertex.
            edgeIndex1=obj.numberOfEdges()+1;
            obj.edges=[obj.edges, edge(vertex1, vertexIndex, edgeIndex1)];
            
            edgeIndex2=obj.numberOfEdges()+1;
            obj.edges=[obj.edges, edge(vertex2, vertexIndex, edgeIndex2)];
            
            %edges aren't truly connected yet...
            
            %Assign correct type to new edges
            
            obj.edges(edgeIndex1).type='boundary';
            obj.edges(edgeIndex2).type='boundary';
            
            %Add new edges as outgoing to vertex1 and vertex2
            
            obj.addOutgoingCrease(vertex1,edgeIndex1);
            obj.addOutgoingCrease(vertex2,edgeIndex2);
            
            %Add new edges as incoming to new vertex
            
            obj.addIncomingCrease(vertexIndex,edgeIndex1);
            obj.addIncomingCrease(vertexIndex,edgeIndex2);
            
            %Make sure the type of vertex1 and vertex2 remains unchanged
            
            obj.vertices(vertex1).type=type1;
            obj.vertices(vertex2).type=type2;
            
        end
        function [vertices, edges, isLefts, isRights, faceIndices] = ...
                cCWMarching(obj, startVertex, startEdge)
            %CCWMARCHING(OBJ, STARTVERTEX, STARTEDGE) Starts at startvertex
            % and marches, starting by going in the direction of startEdge,
            % in counterclockwise direction around the next face (that is,
            % the face seen on the left of startEdge when looking in the
            % direction of startEdge from startVertex.
            %When using counterclockwise marching, the starting vertex is
            %not included in the output. This is different for clockwise
            %marching.
            %OUTPUTS:
            %   vertices: indices of vertices that are encountered during
            %   the marching.
            %   edges: indices of edges that are followed during the
            %   marching.
            %   isLefts: Boolean values specifying whether the "face of
            %   interest" is seen to the left of the "current" edge - when
            %   looking in the edge direction.
            %   isRights: Boolean values specifying whether the "face of
            %   interest" is seen to the right of the "current" edge - when
            %   looking in the edge direction.
            %   faceIndices: Indices specifying the local index of the
            %   "face of interest" on the encountered vertices.
            
            %Before starting anything: Make sure we're not dealing with
            %deleted objects
            if obj.edges(startEdge).deleted
                error('Edge is deleted. This function is not intended to be used with deleted edges.');
            end
            
            if obj.vertices(startVertex).deleted
                error('Vertex is deleted. This function is not intended to be used with deleted edges.');
            end
            
            %prepare empty containers
            vertices=[];
            edges=[];
            isLefts=[];
            isRights=[];
            faceIndices=[];
           
            %initialize iterators
            currentVertex=startVertex;
            currentEdge=startEdge;
            
            while true
                
                %getting the next vertex
                currentVertex=obj.getNextVertex(currentEdge, currentVertex);
                
                %march ahead
                [currentEdge, isLeft, isRight, localFaceIndex] = ...
                    obj.getNextEdgeCCW(currentVertex,currentEdge);
                
                %update containers
                vertices=[vertices, currentVertex];
                edges=[edges, currentEdge];
                isLefts=[isLefts, isLeft];
                isRights=[isRights, isRight];
                faceIndices=[faceIndices, localFaceIndex];
                
                %Checking whether to abort
                if currentEdge==0 || currentVertex==startVertex
                    break;
                end
                
            end
            
        end
        function [vertices, edges, isLefts, isRights, faceIndices] = ...
                cWMarching(obj, startVertex, startEdge)
            %CWMARCHING(OBJ, STARTVERTEX, STARTEDGE) Starts at startvertex
            % and marches, starting by going in the direction of startEdge,
            % in clockwise direction around the next face (that is,
            % the face seen on the right of startEdge when looking at
            % startVertex from startEdge.
            %When using clockwise marching, the starting vertex is
            %included in the output. This is different for counterclockwise
            %marching.
            %OUTPUTS:
            %   vertices: indices of vertices that are encountered during
            %   the marching.
            %   edges: indices of edges that are followed during the
            %   marching.
            %   isLefts: Boolean values specifying whether the "face of
            %   interest" is seen to the left of the "current" edge - when
            %   looking in the edge direction.
            %   isRights: Boolean values specifying whether the "face of
            %   interest" is seen to the right of the "current" edge - when
            %   looking in the edge direction.
            %   faceIndices: Indices specifying the local index of the
            %   "face of interest" on the encountered vertices.
            
            %Before starting anything: Make sure we're not dealing with
            %deleted objects
            if obj.edges(startEdge).deleted
                error('Edge is deleted. This function is not intended to be used with deleted edges.');
            end
            
            if obj.vertices(startVertex).deleted
                error('Vertex is deleted. This function is not intended to be used with deleted edges.');
            end
            
            %Check if the starting edge is connected to the starting vertex
            if not(obj.edges(startEdge).sourceVertex==startVertex || ...
                    obj.edges(startEdge).endVertex==startVertex)
                error('Specified start edge and start vertex are not connected');
            end
            
            %prepare empty containers
            vertices=[];
            edges=[];
            isLefts=[];
            isRights=[];
            faceIndices=[];
          
            %initialize iterators
            currentVertex=startVertex;
            currentEdge=startEdge;
            
            while true
                
                %Marching step
                
                [nextEdge, isLeft, isRight, localFaceIndex] = ...
                    obj.getNextEdgeCW(currentVertex, currentEdge);
                
                %Write to containers
                edges=[edges, currentEdge];
                vertices=[vertices, currentVertex];
                isLefts=[isLefts, isLeft];
                isRights=[isRights, isRight];
                faceIndices=[faceIndices, localFaceIndex];
                
                %Check termination criteria
                if nextEdge==0 || nextEdge==startEdge
                    edges=[edges, nextEdge]; %add the zero s.t. it becomes obvious from the output why the marching ended.
                    break;
                end
                
                %Update iterators
                currentEdge=nextEdge;
                currentVertex=obj.getNextVertex(currentEdge,currentVertex);
                
            end
            
        end
        function [nextEdge, isLeft, isRight, localFaceIndex] = ... 
                getNextEdgeCCW(obj, vertexIndex, edgeIndex)
            %GETNEXTEDGECCW(OBJ,VERTEXINDEX,EDGEINDEX) Gives the next edge
            %(i.e. its index) when walking around a face in
            %counterclockwise direction. edgeIndex is the current edge
            %index and vertexIndex is the next vertex on the face.
            %[nextEdge] = getNextEdgeCCW(...) returns only the index of the
            %next Edge.
            %[nextEdge, isLeft, isRight] = getNextEdgeCCW(...) Also returns
            %the bools isLeft and isRight. If, when looking in the
            %direction of edgeIndex, the face that is "circumnavigated", is
            %seen on the left, isLeft is true. Otherwise (when the face is
            %seen on the right), isRight is true.
            %[nextEdge, isLeft, isRight, localFaceIndex] =
            %getNextEdgeCCW(...) also returns the local Index of the face
            %at the vertex. The first index belongs to the sector between
            %first incoming and last outgoing edge and is then enumerated
            %in counterclockwise direction.
            
            %Before starting anything: Make sure we're not dealing with
            %deleted objects
            if obj.edges(edgeIndex).deleted
                error('Edge is deleted. This function is not intended to be used with deleted edges.');
            end
            
            obj.verifyVertexIndex(vertexIndex);
            
            %Check whether the edge is the source or the end of the edge
            isSource=false;
            isEnd=false;
            if obj.edges(edgeIndex).sourceVertex==vertexIndex
                isSource=true;
            end
            if obj.edges(edgeIndex).endVertex==vertexIndex
                isEnd=true;
            end
            
            %Verifying that the vertex is actually connected to the edge
            if not(isSource) && not(isEnd)
                error('Vertex not connected to specified edge');
            end
            
            %Just for safety: ensure that self-connection does not occur.
            %This would correspond to some weird edge defect.
            if isSource && isEnd
                error('Vertex seems to be connected to itself. This is not allowed.');
            end
            
            %If the vertex is the end of the edge, the circumnavigated face
            %is seen on the left
            isLeft=isEnd;
            isRight=not(isLeft);
            
            %See whether we arrived at a boundary node
            onBoundary=false;
            if strcmp(obj.vertices(vertexIndex).type,'boundary') || ...
                    strcmp(obj.vertices(vertexIndex).type,'source')
                if obj.vertices(vertexIndex).degree == 1 || ...
                        obj.isAuxiliaryBoundary(vertexIndex)
                    %before aborting, we check if there is an auxiliary
                    %boundary vertex ahead somewhere
                    if strcmp(obj.vertices(vertexIndex).type,'source')
                        %source vertex: the boundary vertex would be the
                        %first outgoing edge
                        if not(isempty(obj.vertices(vertexIndex).outgoingEdges))
                            %hopefully, this condition is always true (this
                            %is just a double-check to really avoid
                            %accessing nonexistent entries
                            potentialEdge=obj.vertices(vertexIndex).outgoingEdges(1);
                            potentialVertex=obj.vertices(vertexIndex).outgoingNeighbours(1);
                            if strcmp(obj.edges(potentialEdge).type,'boundary') && ...
                                    strcmp(obj.vertices(potentialVertex).type,'boundary')
                                nextEdge=potentialEdge;
                                localFaceIndex=1;
                                return;
                            end
                        end
                    elseif strcmp(obj.vertices(vertexIndex).type,'boundary')
                        %boundary vertex: potential auxiliary boundary
                        %vertex would be the last outgoing one
                        if not(isempty(obj.vertices(vertexIndex).outgoingEdges))
                            potentialEdge=obj.vertices(vertexIndex).outgoingEdges(end);
                            potentialVertex=obj.vertices(vertexIndex).outgoingNeighbours(end);
                            if strcmp(obj.edges(potentialEdge).type,'boundary') && ...
                                    strcmp(obj.vertices(potentialVertex).type,'boundary')
                                nextEdge=potentialEdge;
                                localFaceIndex=1;
                                return;
                            end
                            
                        end
                    end
                    
                    %Set the output to 0
                    nextEdge=0;
                    %the local index is 1, in this case.
                    localFaceIndex=1;
                    %End evaluation here
                    return;
                else
                    %Degree is 2 or higher. We should still check
                    %convexity(!)
                    %remember for later that we have a special case here.
                    onBoundary=true;
                end
            end
            
            %Calculate the index of the next edge
            vertexDeg=obj.vertices(vertexIndex).degree;
            
            %Ensure that there is no degree 1 internal vertex
            if vertexDeg<2 && not(onBoundary)
                error('Degree 1 vertices not allowed to have any type other than boundary or source. Check vertex nr. %d.',...
                    vertexIndex);
            end
            
            %Execute the actual computation: find the next edge
            
            %Create a vector that contains all edges incident to the
            %vertex.
            localEdgeVec=[obj.vertices(vertexIndex).incomingEdges,...
                obj.vertices(vertexIndex).outgoingEdges];
            if isSource
                %this is an outgoing edge
                localIndex=obj.edges(edgeIndex).indexAtSource;
                nextIndex=localIndex-1; %the next index in counterclockwise
                %direction is just the next one in clockwise direction on 
                %the vertex.
                
                %If this results in an invalid local index, this is a sign
                %that there is no incoming edge to this vertex. This is
                %only allowed (only for the source vertex, in some cases).
                if nextIndex<=0 && not(strcmp(obj.vertices(vertexIndex).type,'source'))
                    error('Error obtaining next edge. Make sure there are no non-source vertices without inputs.');
                elseif nextIndex==0 && strcmp(obj.vertices(vertexIndex).type,'source')
                    %In this case, there is a multi-output source vertex.
                    %This case is untested.
                    %In this case, it is allowed to set the next index
                    %manually.
                    nextIndex=vertexDeg;
                end
            elseif isEnd
                %this is an incoming edge
                localIndex=obj.edges(edgeIndex).indexAtEnd;
                %Caution: subtraction of 1 may result in an index of 0.
                %This needs to be corrected: 0 is "translated" to the last
                %index.
                nextIndex=mod(localIndex+vertexDeg-2,vertexDeg)+1;
            end
            
            nextEdge=localEdgeVec(nextIndex);
            
            
            %Face Index is equal to local Index.
            localFaceIndex=localIndex;
            
            %Before terminating, check the convexity of the face
            currentVec=obj.awayVector0(vertexIndex,edgeIndex);
            nextVec=obj.awayVector0(vertexIndex,nextEdge);
            zAxis=obj.vertices(vertexIndex).unfoldedZ;
            ang=counterClockWiseAngle(nextVec,currentVec,zAxis);
            
            if ang>=180  %Face is not strictly convex.
                if not(onBoundary)
                    %internal vertex
                    %->sector angles like this are not allowed.
                    %error/abort.
                    error('Sector angles of 180° or more are not allowed. Check Vertex nr. %d.',vertexIndex);
                end
            end
            
            if (onBoundary && ang>=180) ||...
                    (onBoundary && (nextIndex>localIndex))
                    %angle found on boundary
                    %such a sector angle demands abortion of the marching
                    %scheme.
                    
                    %before returning the signal to abort, we check if
                    %there is an auxiliary boundary vertex ahead
                    if strcmp(obj.vertices(vertexIndex).type,'boundary')
                        %potential vertex is last outgoing one
                        if not(isempty(obj.vertices(vertexIndex).outgoingNeighbours))
                            potentialVertex=obj.vertices(vertexIndex)...
                                .outgoingNeighbours(end);
                            if obj.isAuxiliaryBoundary(potentialVertex)
                                nextEdge=obj.vertices(vertexIndex).outgoingEdges(end);
                                localFaceIndex=1;
                                return;
                            end
                        end
                        
                    elseif strcmp(obj.vertices(vertexIndex).type,'source')
                        %potential vertex is first outgoing one
                        
                        %there are always outgoing vertices to the source
                        %vertex
                        potentialVertex=obj.vertices(vertexIndex).outgoingNeighbours(1);
                        if obj.isAuxiliaryBoundary(potentialVertex)
                            nextEdge=obj.vertices(vertexIndex).outgoingEdges(1);
                            localFaceIndex=1;
                            return;
                        end
                        
                    end
                    
                    nextEdge=0;
                    localFaceIndex=1; %whenever this case occurs, a new 
                    %("boundary") face will be created. The face will have
                    %local index 1.
            end 
            
        end
        function [nextEdge, isLeft, isRight, localFaceIndex] = ... 
                getNextEdgeCW(obj, vertexIndex, edgeIndex)
            %GETNEXTEDGECW(OBJ,VERTEXINDEX,EDGEINDEX) Gives the next edge
            %(i.e. its index) when walking around a face in
            %clockwise direction. vertexIndex is the current edge
            %index and vertexIndex is the next edge on the face.
            %[nextEdge] = getNextEdgeCW(...) returns only the index of the
            %next Edge.
            %[nextEdge, isLeft, isRight] = getNextEdgeCW(...) Also returns
            %the bools isLeft and isRight. If, when looking in the
            %direction of edgeIndex, the face that is "circumnavigated", is
            %seen on the left, isLeft is true. Otherwise (when the face is
            %seen on the right), isRight is true.
            %[nextEdge, isLeft, isRight, localFaceIndex] =
            %getNextEdgeCCW(...) also returns the local Index of the face
            %at the vertex. The first index belongs to the sector between
            %first incoming and last outgoing edge and is then enumerated
            %in counterclockwise direction.
            
            %Check whether the edge is the source or the end of the edge
            isSource=false;
            isEnd=false;
            if obj.edges(edgeIndex).sourceVertex==vertexIndex
                isSource=true;
            end
            if obj.edges(edgeIndex).endVertex==vertexIndex
                isEnd=true;
            end
            
            %Verifying that the vertex is actually connected to the edge
            if not(isSource) && not(isEnd)
                error('Vertex not connected to specified edge');
            end
            
            %Just for safety: ensure that self-connection does not occur.
            %This would correspond to some weird edge defect.
            if isSource && isEnd
                error('Vertex seems to be connected to itself. This is not allowed.');
            end
            
            %If the vertex is the end of the edge, the circumnavigated face
            %is seen on the left
            isLeft=isSource;
            isRight=not(isLeft);
            
            %See whether we arrived at a boundary node
            onBoundary=false;
            if strcmp(obj.vertices(vertexIndex).type,'boundary') || ...
                    strcmp(obj.vertices(vertexIndex).type,'source')
                if obj.vertices(vertexIndex).degree == 1 || ...
                        obj.isAuxiliaryBoundary(vertexIndex)
                    %before aborting, we check if there is an auxiliary
                    %boundary vertex ahead somewhere
                    if strcmp(obj.vertices(vertexIndex).type,'source')
                        %source vertex: the boundary vertex would be the
                        %last outgoing edge
                        if not(isempty(obj.vertices(vertexIndex).outgoingEdges))
                            %hopefully, this condition is always true (this
                            %is just a double-check to really avoid
                            %accessing nonexistent entries
                            potentialEdge=obj.vertices(vertexIndex).outgoingEdges(end);
                            potentialVertex=obj.vertices(vertexIndex).outgoingNeighbours(end);
                            if strcmp(obj.edges(potentialEdge).type,'boundary') && ...
                                    strcmp(obj.vertices(potentialVertex).type,'boundary')
                                nextEdge=potentialEdge;
                                localFaceIndex=2;
                                return;
                            end
                        end
                    elseif strcmp(obj.vertices(vertexIndex).type,'boundary')
                        %boundary vertex: potential auxiliary boundary
                        %vertex would be the first outgoing one
                        if not(isempty(obj.vertices(vertexIndex).outgoingEdges))
                            potentialEdge=obj.vertices(vertexIndex).outgoingEdges(1);
                            potentialVertex=obj.vertices(vertexIndex).outgoingNeighbours(1);
                            if strcmp(obj.edges(potentialEdge).type,'boundary') && ...
                                    strcmp(obj.vertices(potentialVertex).type,'boundary')
                                nextEdge=potentialEdge;
                                localFaceIndex=2;
                                return;
                            end
                            
                        end
                    end
                    
                    %Set main output to 0 (signal end of marching scheme)
                    nextEdge=0;
                    %Set local Face Index
                    localFaceIndex=2; %2 in this case: index 1 goes to the "opposite" face.
                    %exception: auxiliary boundary edges
                    if obj.isAuxiliaryBoundary(vertexIndex)
                        localFaceIndex=1;
                    end
                    %End evaluation here
                    return;
                else
                    %Degree is 2 or higher. We should still check
                    %convexity(!)
                    %remember for later that we have a special case here.
                    onBoundary=true;
                end
            end
            
            %Calculate the vertex degree
            vertexDeg=obj.vertices(vertexIndex).degree;
            
            %Ensure that there is no degree 1 internal vertex
            if vertexDeg<2 && not(onBoundary)
                error('Degree 1 vertices not allowed to have any type other than boundary, source or aux. Check vertex nr. %d.',...
                    vertexIndex);
            end
            
            %Execute the actual computation: find the next edge
            
            %Create a vector that contains all edges incident to the
            %vertex.
            localEdgeVec=[obj.vertices(vertexIndex).incomingEdges,...
                obj.vertices(vertexIndex).outgoingEdges];
            if isSource
                %this is an outgoing edge
                localIndex=obj.edges(edgeIndex).indexAtSource;
                %Caution: Addition of 1 may result in an index of deg+1.
                %This needs to be corrected to 1
                nextIndex=mod(localIndex,vertexDeg)+1;
            elseif isEnd
                %this is an incoming edge
                localIndex=obj.edges(edgeIndex).indexAtEnd;
                nextIndex=mod(localIndex,vertexDeg)+1; %This has to be
                %corrected too: A vertex can have no outputs; in that case
                %simply adding one would lead to senseless results.
                %most usually, this result is changed again in one of the
                %very last lines of this function.
            end
            
            nextEdge=localEdgeVec(nextIndex);
            %set the local face index
            localFaceIndex=nextIndex;
            
            %Before terminating, check the convexity of the face
            currentVec=obj.awayVector0(vertexIndex,edgeIndex);
            nextVec=obj.awayVector0(vertexIndex,nextEdge);
            zAxis=obj.vertices(vertexIndex).unfoldedZ;
            ang=counterClockWiseAngle(currentVec,nextVec,zAxis);
            
            if ang>=180  %Face is not strictly convex.
                if not(onBoundary)
                    %internal vertex
                    %->sector angles like this are not allowed.
                    %error/abort.
                    error('Sector angles of 180° or more are not allowed. Check vertex nr. %d.',vertexIndex);
                end
            end
            
            if (onBoundary && ang>=180) || (onBoundary && (localIndex>nextIndex))
                    %angle found on boundary
                    %such a sector angle demands abortion of the marching
                    %scheme.
                    
                    
                    %before returning the signal to abort, we check if
                    %there is an auxiliary boundary vertex ahead
                    if strcmp(obj.vertices(vertexIndex).type,'boundary')
                        %potential vertex is first outgoing one
                        if not(isempty(obj.vertices(vertexIndex).outgoingNeighbours))
                            potentialVertex=obj.vertices(vertexIndex)...
                                .outgoingNeighbours(1);
                            if obj.isAuxiliaryBoundary(potentialVertex)
                                nextEdge=obj.vertices(vertexIndex).outgoingEdges(1);
                                localFaceIndex=vertexDeg+1;
                                return;
                            end
                        end
                        
                    elseif strcmp(obj.vertices(vertexIndex).type,'source')
                        %potential vertex is last outgoing one
                        
                        %there are always outgoing vertices to the source
                        %vertex
                        potentialVertex=obj.vertices(vertexIndex).outgoingNeighbours(end);
                        if obj.isAuxiliaryBoundary(potentialVertex)
                            nextEdge=obj.vertices(vertexIndex).outgoingEdges(end);
                            localFaceIndex=vertexDeg+1;
                            return;
                        end
                        
                    end
                    
                    nextEdge=0;
                    %Set the local face index
                    localFaceIndex=vertexDeg+1;
            end
            
        end
        function [verdict] = isAuxiliaryBoundary(obj, vertexIndex)
            %VERDICT=ISAUXILIARYBOUNDARY(OBJ,VERTEXINDEX) checks whether
            %the specified vertex index is an auxiliary boundary vertex,
            %i.e. a special boundary vertex that helps defining the surface
            
            %make sure we are dealing with a valid vertex index
            obj.verifyVertexIndex(vertexIndex);
            
            %type must be boundary
            if not(strcmp(obj.vertices(vertexIndex).type,'boundary'))
                verdict=false;
                return;
            end
            
            %exactly one input
            if not(length(obj.vertices(vertexIndex).incomingEdges)==1)
                verdict=false;
                return;
            end
            
            %the one input must be a boundary edge
            if not(strcmp(obj.edges(obj.vertices(vertexIndex)...
                    .incomingEdges(1)).type,'boundary'))
                verdict=false;
                return;
            end
            
            %exactly one output
            if not(length(obj.vertices(vertexIndex).outgoingEdges)==1)
                verdict=false;
                return;
            end
            
            %the one output must be a boundary edge
            if not(strcmp(obj.edges(obj.vertices(vertexIndex)...
                    .outgoingEdges(1)).type,'boundary'))
                verdict=false;
                return;
            end
            
            %all tests passed
            verdict=true;
            
        end
        function [nextVertex] = getNextVertex(obj, edgeIndex, prevVertex)
            %GETNEXTVERTEX(EDGEINDEX, PREVVERTEX) Returns the next vertex
            %when walking along edge with edgeIndex. In other words,
            %returns the vertex at the opposite end of prevVertex on
            %edgeIndex.
            
            if obj.edges(edgeIndex).sourceVertex==prevVertex
                nextVertex=obj.edges(edgeIndex).endVertex;
            elseif obj.edges(edgeIndex).endVertex==prevVertex
                nextVertex=obj.edges(edgeIndex).sourceVertex;
            else
                error('Edge seems to be unconnected to the previous vertex');
            end
            
        end
        function [] = setFaceNormal(obj, faceIndex, normalVector)
            %SETFACENORMAL Sets the face normal on a specified face and
            %also sets the corresponding quantity on all connected vertices
            %and edges.
            
            %check whether normal vector is of length 1
            if abs(norm(normalVector)-1)>1e-6
                error('The provided normal vector must be of unit length');
            end
            
            %set the normal vector on the face
            obj.faces(faceIndex).normalVector.vec=normalVector;
            
            %Set the face to "determined"
            obj.faces(faceIndex).determined=true;
            
        end 
                %%%%%Input validation functions%%%% (unreviewed)
        function [] = verifyVertexIndex(obj,vertexIndex)
            %VERIFYVERTEXINDEX(OBJ,VERTEXINDEX) Checks whether the
            %specified vertex index is deleted, nonpositive or too large
            %and throws an error if any of this is the case.
            
            %can be used for multilevel input
            if length(vertexIndex)>1
                %dealing with a multilevel input
                if vertexIndex(1)<1 || vertexIndex(1)>length(obj.subOrigamis)
                    error('The specified subOrigami index of %d is invalid. The current level has %d subOrigamis.',...
                        vertexIndex(1),length(obj.subOrigamis));
                end
                %recursive call: do the same on the next-lower level
                obj.subOrigamis{vertexIndex(1)}.verifyVertexIndex(vertexIndex(2:end));
                return;
            end
            
            if vertexIndex<1
                error('The specified vertex index (nr. %d) is not strictly positive and therefore invalid.',...
                    vertexIndex);
            end
            
            if vertexIndex>obj.numberOfVertices()
                error('The specified vertex index (nr. %d) is larger than the number of vertices in the origami (%d) and is therefore invalid.',...
                    vertexIndex,obj.numberOfVertices());
            end
            
            if obj.vertices(vertexIndex).deleted
                error('The specified vertex index (nr. %d) is deleted and therefore invalid.',...
                    vertexIndex);
            end
            
        end
        
%%%%Auxiliary Functions for the preparation and evaluation of the origami %%%%%
        
        function [] = calculateSectorAngles(obj)
            %CALCULATESECTORANGLES(OBJ) calculates all sector angles at all
            %vertices except for auxiliary vertices
            
            %iterate over all vertices
            for j=1:obj.numberOfVertices()
                
                %calculate sector angle of current vertex
                obj.calculateSectorAngle(j);
                
            end %loop over vertices
            
        end
        function [] = calculateSectorAngle(obj,vertexIndex)
            %CALCULATESECTORANGLE(OBJ,VERTEXINDEX) Calculates the sector
            %angles around the specified vertex and saves them to the
            %respective containers.
            currentV=obj.vertices(vertexIndex);
            
            %make sure they're not deleted or auxiliary -> skip
            if currentV.deleted || strcmp(currentV.type,'aux')
                return;
            end
            
            %shift parameter to deal with source vector
            shift=0;
            
            %if the vertex is of type source or boundary,
            if strcmp(currentV.type,'source') || strcmp(currentV.type,'boundary')
                %initialize container to degree+1
                obj.vertices(vertexIndex).sectorAngles=zeros(1,currentV.degree+1);
                
                %set localDeg to deg+2
                localDeg=currentV.degree+2;
                
                %and ensure that the last outgoing edge is omitted when going
                %through the outgoing edges
                endIterator=size(currentV.outgoingEdges,2)-1;
                
                %When dealing with the source vertex, consider that the
                %first face is "invalid": according to common
                %indexation convention, the first face is located
                %between the last and first outgoing vertex (which is
                %invalid). This is fixed by adding a shifting parameter
                %to "relocate" the first face to what is actually the
                %second.
                %Shift is also set to 1 if the vertex is a boundary vertex
                %with two adjacent boundary edges. In that case, the
                %convention is to store the sector angle that is required
                %to transform the input crease vector to the output crease
                %vector. This is the sector angle from the incoming to the
                %outgoing crease. Therefore, we skip the first face (which
                %would be from outgoing crease to incoming crease).
                if strcmp(currentV.type,'source') || ...
                        obj.isLocalSource(vertexIndex) || ...
                        obj.isOnBoundary(vertexIndex)
                    shift=1;
                end
                
            elseif strcmp(currentV.type,'interior')
                %if the vertex is internal, set the container size to the
                %degree of the vertex
                obj.vertices(vertexIndex).sectorAngles=zeros(1,currentV.degree);
                
                %set localDeg to deg
                localDeg=currentV.degree;
                
                %iterate through all adjacent edges
                endIterator=size(currentV.outgoingEdges,2);
                
            end
            
            incidentEdges=[currentV.incomingEdges, currentV.outgoingEdges];
            nIncoming=size(currentV.incomingEdges,2);
            
            %iterate through all incident edges
            %shift parameter considers the situation at the source
            %vertex.
            for k=(1+shift):(nIncoming+endIterator+shift)
                
                %find the next edge to the left of the crease
                leftIndex=mod((k+localDeg-2),localDeg)+1;
                leftIndex=incidentEdges(leftIndex);
                %find the index of the current crease
                rightIndex=incidentEdges(k);
                
                %find the angle in clockwise direction to the current
                %edge
                outVec1=obj.awayVector0(vertexIndex,leftIndex);
                outVec2=obj.awayVector0(vertexIndex,rightIndex);
                
                %see whether either edge is on the boundary
                isLeftBoundary=strcmp(obj.edges(leftIndex).type,'boundary');
                isRightBoundary=strcmp(obj.edges(rightIndex).type,'boundary');
                
                zAxis=obj.vertices(vertexIndex).unfoldedZ;
                
                angle=counterClockWiseAngle(outVec1,outVec2,zAxis);
                if angle<1e-9
                    error('Very small angle detected at vertex %d at sector nr. %d',vertexIndex,k);
                elseif (angle>=180 || imag(angle)~=0) && ...
                        not(isLeftBoundary) && not(isRightBoundary) && ...
                        not(strcmp(obj.vertices(vertexIndex).type,'boundary'))
                    error('Sector angles of 180° or more are not allowed. Check Vertex %d, sector nr. %d', vertexIndex, k-shift);
                end
                
                obj.vertices(vertexIndex).sectorAngles(k-shift)=angle;
            end
            
        end
        function [] = calculateEdgeLengths(obj)
            %CALCULATEEDGELENGTHS(OBJ) Calculates the lengths of all edges
            %in the origami.
            
            %Loop over all edges
            for j=1:obj.numberOfEdges()
                
                %calculate length
                obj.calculateEdgeLength(j);
                
            end %loop over all edges
            
        end
        function length = calculateEdgeLength(obj,edgeIndex)
            %CALCULATEEDGELENGTH(OBJ,EDGEINDEX) Calculates the length of a
            %specified edge and saves the length to the respective
            %container. (and also returns that value)
            
            if edgeIndex>obj.numberOfEdges()
                error('The specified edge has index %d but the origami currently only has %d vertices. This is invalid.',edgeIndex,obj.numberOfEdges);
            end
            
            currentE=obj.edges(edgeIndex);
            
            %Check whether edge is deleted -> skip
            if currentE.deleted
                return;
            end
            
            %Get position of source vertex
            sourcePos=obj.vertices(currentE.sourceVertex).position0;
            
            %Get position of end vertex
            endPos=obj.vertices(currentE.endVertex).position0;
            
            %Calculate length
            obj.edges(edgeIndex).length=norm(endPos-sourcePos,2);
            
            %Check if length is valid
            if obj.edges(edgeIndex).length < 1e-9
                error('Very short edge detected');
            end
            
            length=obj.edges(edgeIndex).length;
            
        end
        
        function [allDet] = isOrigamiDetermined(obj)
            %ISORIGAMIDETERMINED(OBJ) Checks whether all vertices within an
            %origami are determined
            
            %First getting all vertices that are not deleted, then making a
            %vector of whether they are determined or not. Then, 'all' can
            %be used to check if all of the entries are true.
            allDet=all([obj.vertices(not([obj.vertices(:).deleted])).determined]);
            
        end
        function [allDet,allSubsDet] = isOrigamiCompletelyDetermined(obj)
            %[ALLDET,ALLSUBSDET]=ISORIGAMICOMPLETELYDETERMINED(OBJ) Specifies whether
            %the origami is completely determined, i.e. the origami and all
            %of its subOrigamis (and their subOrigamis etc.) are fully
            %determined. The second output specifies whether all
            %subOrigamis have been determined.
            
            %collecting all criteria in a vector
            allDet=false(1,1);
            allDet(1)=obj.isOrigamiDetermined();       
        end 
        function [] = determineSubOrigamis0(obj)
            %DETERMINESUBORIGAMIS0(OBJ) Goes through all subOrigamis
            %and passes the interface information (vertex positions and
            %crease vectors) to the subOrigamis and calls determineOrigami
            %with the subOrigamis. As interface, the positions and vectors
            %from the undeformed origami are passed. The positions obtained
            %this way are required when using a subOrigami vertex as
            %reference for a mirroring or translation operation.
            
            %requires state doKinematics or editFaces
            if strcmp(obj.status,'editGraph')
                obj.checkState('editFaces');
            end
            
            %when called from a subOrigami, the subOrigami must be
            %determined
            if isa(obj,'subOrigami')
                if not(obj.isOrigamiDetermined())
                    error('SubOrigami must be completely determined before determining lower-level subOrigamis.');
                end
            end
            
            %pass information to all subOrigamis
            for j=1:length(obj.subOrigamis)
                %caution: this function is called recursively, i.e. also
                %from the subOrigamis.
                
                if isa(obj,'subOrigami')
                    obj.passInformationToSubOrigami(j);
                else
                    obj.passInformationToSubOrigami0(j);
                end
                
                %calculate sector angles and lengths
                obj.subOrigamis{j}.calculateSectorLengths(); %this call is recursive
                
                %determine the subOrigami - and everything on lower levels
                %too.
                obj.subOrigamis{j}.determineOrigami(); %this call is recursive
                
            end
            
        end            
        function [] = determineSubOrigamis(obj)
            %DETERMINESUBORIGAMIS(OBJ) Goes through all subOrigamis
            %and passes the interface information (vertex positions and
            %crease vectors) to the subOrigamis and calls determineOrigami
            %with the subOrigamis.
            
            %requires state doKinematics
            obj.checkState('doKinematics');
            
            %try to determine all subOrigamis
            for j=1:length(obj.subOrigamis)
                
                %skip subOrigamis that are already determined
                if obj.subOrigamis{j}.isOrigamiCompletelyDetermined()
                    continue;
                end
                
                %find out whether the subOrigami can be determined
                sourceVert=obj.subOrigamis{j}.sourceParent;
                sourceDet=obj.vertices(sourceVert).determined;
                
                firstParentEdge=obj.subOrigamis{j}.firstEdgeParent;
                firstEdgDet=obj.edges(firstParentEdge).determined;
                
                secondParentEdge=obj.subOrigamis{j}.secondEdgeParent;
                secondEdgDet=obj.edges(secondParentEdge).determined;
                
                %if all information is ready, we can march on
                if sourceDet && firstEdgDet && secondEdgDet
                    obj.passInformationToSubOrigami(j);
                    obj.subOrigamis{j}.determineOrigami();
                end
            end
            
        end
        function [] = passInformationToSubOrigami0(obj,index)
            %PASSINFORMATIONTOSUBORIGAMI0(OBJ,INDEX) Passes the information
            %to the subOrigami with the specified index. After this, the
            %source vertex of the subOrigami should be determinable. This
            %function passes down information obtained from the undeformed
            %configuration. (For passing infomation from the deformed
            %configuration, use passInformationToSubOrigami().
            
            %make sure the index is meaningful
            if index<1 || index>length(obj.subOrigamis)
                error('The specified subOrigami index is invalid. It is either smaller than 1 or larger than the number of subOrigamis.');
            end
            
            %set the position of the parent vertex
            parentVertex=obj.subOrigamis{index}.sourceParent;
            obj.subOrigamis{index}.sourcePos=obj.vertices(parentVertex).position0;
            
            %get the first crease vector and write it to the subOrigami
            creaseIndex=obj.subOrigamis{index}.firstEdgeParent;
            creaseVec=obj.awayVector0(parentVertex,creaseIndex);
            obj.subOrigamis{index}.firstCreaseVec=creaseVec;
            
            %same story with second crease
            creaseIndex=obj.subOrigamis{index}.secondEdgeParent;
            creaseVec=obj.awayVector0(parentVertex,creaseIndex);
            obj.subOrigamis{index}.secondCreaseVec=creaseVec;
            
            %and write to the subOrigami that this information was passed
            %down
            obj.subOrigamis{index}.inputInformationObtained=true;
            
            %resetting the subOrigami to undetermined (most likely, this
            %was new information -> does not agree with potentially
            %previously saved information
            obj.subOrigamis{index}.resetToUndetermined();
            
        end
        function [] = passInformationToSubOrigami(obj,index)
            %PASSINFORMATIONTOSUBORIGAMI(OBJ,INDEX) Passes the information
            %to the subOrigami with the specified index. After this, the
            %source vertex of the subOrigami should be determinable.
            
            %make sure the index is meaningful
            if index<1 || index>length(obj.subOrigamis)
                error('The specified subOrigami index is invalid. It is either smaller than 1 or larger than the number of subOrigamis.');
            end
            
            %set the position of the parent vertex
            parentVertex=obj.subOrigamis{index}.sourceParent;
            obj.subOrigamis{index}.sourcePos=obj.vertices(parentVertex).position;
            
            %get the first crease vector and write it to the subOrigami
            creaseIndex=obj.subOrigamis{index}.firstEdgeParent;
            creaseVec=obj.awayVector(parentVertex,creaseIndex);
            obj.subOrigamis{index}.firstCreaseVec=creaseVec;
            
            %same story with second crease
            creaseIndex=obj.subOrigamis{index}.secondEdgeParent;
            creaseVec=obj.awayVector(parentVertex,creaseIndex);
            obj.subOrigamis{index}.secondCreaseVec=creaseVec;
            
            %and write to the subOrigami that this information was passed
            %down
            obj.subOrigamis{index}.inputInformationObtained=true;
            
            %resetting the subOrigami to undetermined (most likely, this
            %was new information -> does not agree with potentially
            %previously saved information
            obj.subOrigamis{index}.resetToUndetermined();
            
        end
        %
        function [ready] = isVertexDeterminable(obj, vertexIndex)
            %ISVERTEXDETERMINABLE(OBJ, VERTEXINDEX) Determines whether a
            %vertex can be determined/evaluated. Briefly, calls the
            %function appropriate to the vertex type to determine this.
            
            
            if strcmp(obj.vertices(vertexIndex).type,'source')
                ready=obj.isSourceVertexDeterminable(vertexIndex);
            elseif strcmp(obj.vertices(vertexIndex).type,'aux')
                ready=obj.isAuxVertexDeterminable(vertexIndex);
            elseif strcmp(obj.vertices(vertexIndex).type,'boundary')
                ready=obj.isBoundaryVertexDeterminable(vertexIndex);
            elseif strcmp(obj.vertices(vertexIndex).type,'interior')
                ready=obj.isInternalVertexDeterminable(vertexIndex);
            else
                error('Unknown vertex type');
            end
            
        end        
        function [ready] = isSourceVertexDeterminable(obj, vertexIndex)
            %ISSOURCEVERTEXDETERMINABLE(OBJ,VERTEXINDEX) Checks whether
            %source vertex is ready to be evaluated/determinable
            
            if not(strcmp(obj.vertices(vertexIndex).type,'source'))
                error('Wrong determinability function');
            end
            
            %Going through all criteria, if one of them is not given,
            %return false.
            ready=true;
            
            %criterion1: At least one outgoing edge must be internal. 
            %criterion2: Among the internal outgoing edges, at least one
            %must have a non-trivial dihedral angle.
            criterion1=false;
            criterion2=false;
            for j=1:size(obj.vertices(vertexIndex).outgoingEdges,2)
                edgeIndex=obj.vertices(vertexIndex).outgoingEdges(j);
                if strcmp(obj.edges(edgeIndex).type,'internal')
                    criterion1=true;
                    if obj.edges(edgeIndex).dihedralAngle.ang~=0
                        criterion2=true;
                    end
                end
            end
            
            %Check if these criteria are met
            if not(criterion1 && criterion2)
                error('Can not determine the source vertex. Either no outgoing internal edge was found or all creases have a trivial dihedral angle. HINT: check if the input angle is nonzero');
            end
            
            %Exactly one adjacent face must be 
            %(!) is this even needed?
            numFixed=0;
            for j=1:size(obj.vertices(vertexIndex).adjacentFaces,2)
                faceIndex=obj.vertices(vertexIndex).adjacentFaces(j);
                if obj.faces(faceIndex).fixed
                    numFixed=numFixed+1;
                end
            end
            
            %Check if criterion is met
            if numFixed~=1 && not(obj.symmetricFold)
                error('Can not determine the source vertex. Did you fix exactly one face?');
            end
            
            
            if ready
                obj.vertices(vertexIndex).determinable=true;
            end
            
        end
        function [ready] = isAuxVertexDeterminable(obj, vertexIndex)
            %ISAUXVERTEXDETERMINABLE(OBJ, VERTEXINDEX) Checks whether an
            %auxiliary vertex is ready to be evaluated/determined
            
            %Check if the vertex actually is an auxiliary vertex
            if not(strcmp(obj.vertices(vertexIndex).type,'aux'))
                error('Wrong Determinability Function');
            end
            
            ready=true;
            
            %Check if the vertex has two incoming edges and no outgoing
            %ones
            criterion1=size(obj.vertices(vertexIndex).incomingEdges,2)==2;
            criterion2=size(obj.vertices(vertexIndex).outgoingEdges,2)==0;
            if not(criterion1 && criterion2)
                error('Auxiliary vertex with index %d does not have exactly two inputs and no outputs as required.',vertexIndex);
            end
            
            %Check whether one of the input edges is determined
            edge1=obj.vertices(vertexIndex).incomingEdges(1);
            edge2=obj.vertices(vertexIndex).incomingEdges(2);
            criterion1=obj.edges(edge1).determined;
            criterion2=obj.edges(edge2).determined;
            if not(criterion1 || criterion2)
                ready=false;
                return;
            end
            
            %Set vertex to determinable
            if ready
                obj.vertices(vertexIndex).determinable=true;
            end
            
        end    
        function [ready] = isBoundaryVertexDeterminable(obj, vertexIndex)
            %ISBOUNDARYVERTEXDETERMINABLE(OBJ, VERTEXINDEX) Determines
            %whether boundary vertex is ready to be evaluated/determined
            
            if not(strcmp(obj.vertices(vertexIndex).type,'boundary'))
                error('Wrong determinability function');
            end
            
            ready=true;
            
            %Vertex should have at least one input
            criterion=size(obj.vertices(vertexIndex).incomingEdges,2)>=1;
            if not(criterion)
                error('Boundary Vertex with index %d appears to not have at least one input.',vertexIndex);
            end
            
            %All inputs should be determined
            numInputs=size(obj.vertices(vertexIndex).incomingEdges,2);
            for j=1:numInputs
                edgeIndex=obj.vertices(vertexIndex).incomingEdges(j);
                if not(obj.edges(edgeIndex).determined)
                    ready=false;
                    return;
                end
            end
            
            %Vertex should have exactly no more than two outputs
            criterion=size(obj.vertices(vertexIndex).outgoingEdges,2)<=2;
            if not(criterion)
                error('Boundary vertex with index %d appears to have more than two outputs that are boundary edges.', vertexIndex);
            end
            
            %All adjacent faces should be determined
            adjFaces=obj.vertices(vertexIndex).adjacentFaces;
            numFaces=size(adjFaces,2);
            for j=1:numFaces
                currentFace=adjFaces(j);
                if not(obj.faces(currentFace).determined)
                    ready=false;
                    return;
                end
            end
            
            %Both outputs should be boundary edges
            edge1=obj.vertices(vertexIndex).outgoingEdges(1);
            criterion1=strcmp(obj.edges(edge1).type,'boundary');
            edge2=obj.vertices(vertexIndex).outgoingEdges(end);
            criterion2=strcmp(obj.edges(edge2).type,'boundary');
            
            %if the vertex defines the shape of a boundary face, the two
            %adjacent edges should both be boundary edges
            if obj.isOnBoundary(vertexIndex)
                adjEdg=[obj.vertices(vertexIndex).incomingEdges,...
                    obj.vertices(vertexIndex).outgoingEdges];
                for j=adjEdg
                    if not(strcmp(obj.edges(j).type,'boundary'))
                        criterion2=false;
                    end
                end
            end
            
            if not(criterion1 && criterion2)
                error('Boundary vertex with index %d appears to not have exactly two outputs that are boundary edges.', vertexIndex);
            end
            
            %Set vertex to determinable
            if ready
                obj.vertices(vertexIndex).determinable=true;
            end
            
        end
        function [ready] = isInternalVertexDeterminable(obj, vertexIndex)
            %ISINTERNALVERTEXDETERMINABLE(OBJ, VERTEXINDEX) Checks whether
            %an internal vertex is ready to be evaluated/determined
            
            %Check vertex type
            if not(strcmp(obj.vertices(vertexIndex).type,'interior'))
                error('Wrong determinability function');
            end
            
            ready=true;
            
            %Vertex should have at least one input
            criterion=size(obj.vertices(vertexIndex).incomingEdges,2)>=1;
            if not(criterion)
                error('Internal vertex with index %d appears to not have at least one input.',vertexIndex);
            end
            
            %All inputs should be determined
            numInputs=size(obj.vertices(vertexIndex).incomingEdges,2);
            for j=1:numInputs
                edgeIndex=obj.vertices(vertexIndex).incomingEdges(j);
                if not(obj.edges(edgeIndex).determined)
                    ready=false;
                    return;
                end
            end
            
            %The faces between incoming edges should be determined, as well
            %as the one "before" the first edge
            adjFaces=obj.vertices(vertexIndex).adjacentFaces(1:numInputs);
            numFaces=size(adjFaces,2);
            for j=1:numFaces
                faceIndex=adjFaces(j);
                if not(obj.faces(faceIndex).determined)
                    ready=false;
                    return;
                end
            end
            
            %Vertex should have exactly three outputs
            criterion=size(obj.vertices(vertexIndex).outgoingEdges,2)==3;
            if not(criterion)
                error('Internal vertex with index %d does not have exactly three outputs.',vertexIndex);
            end
            
            %All incident edges should be of type 'internal'
            incidentEdges=[obj.vertices(vertexIndex).incomingEdges, ...
                obj.vertices(vertexIndex).outgoingEdges];
            numIncident=size(incidentEdges,2);
            for j=1:numIncident
                edgeIndex=incidentEdges(j);
                if not(strcmp(obj.edges(edgeIndex).type,'internal'))
                    error('Internal vertex with index %d has incident edge with index %d which is not an internal edge. Only internal edges are allowed as incident edges to internal vertices.',vertexIndex,edgeIndex);
                end
            end
            
            %Set vertex to determinable
            obj.vertices(vertexIndex).determinable=true;
            
        end
        %
        function [] = determineVertex(obj, vertexIndex)
            %DETERMINEVERTEX(OBJ, VERTEXINDEX) Determines/evaluates vertex.
            %Briefly, calls the determination function according to the
            %type of the vertex.
            
            if strcmp(obj.vertices(vertexIndex).type,'source')
                obj.determineSourceVertex(vertexIndex);
            elseif strcmp(obj.vertices(vertexIndex).type,'aux')
                obj.determineAuxVertex(vertexIndex);
            elseif strcmp(obj.vertices(vertexIndex).type, 'boundary')
                obj.determineBoundaryVertex(vertexIndex);
            elseif strcmp(obj.vertices(vertexIndex).type, 'interior')
                obj.determineInternalVertex(vertexIndex);
            else
                error('Unknown vertex type');
            end
            
        end
        function [] = determineSourceVertex(obj, vertexIndex)
            %DETERMINESOURCEVERTEX(OBJ, VERTEXINDEX) Determines a source
            %vertex and calculates all dihedral angles, direction vectors
            %of adjacent creases and normal vectors of adjacent faces.
            
            %Initial check
            if not(obj.isSourceVertexDeterminable(vertexIndex))
                error('Can not determine non-determinable vertex.');
            end
            
            %Find out where the fixed face is
            fFaceIndex=obj.fixedFaceIndex;
            
            faceIndex=find(obj.vertices(vertexIndex).adjacentFaces==fFaceIndex);
            
            
            if size(faceIndex,2)~=1 && not(obj.symmetricFold)
                error('Could not find the fixed face near the source vertex');
            end
            
            %Set the properties of the related face and edges (and the
            %vertex)
                zAxis=[0;0;1];
                
            %zAxis is the normal vector to the fixed face
            obj.setFaceNormal(fFaceIndex,zAxis);
            
            %Set the edge "before" the fixed edge
            rightEdgeLocal=faceIndex;
            rightEdge=obj.vertices(vertexIndex).outgoingEdges(rightEdgeLocal);
            
            %same story with other edge
            leftEdgeLocal=faceIndex+1;
            leftEdge=obj.vertices(vertexIndex).outgoingEdges(leftEdgeLocal);
            
            %define the dihedral angle to send to the edge - setup
            %function
            if strcmp(obj.edges(rightEdge).type,'boundary')
                dihAngRight=0;
            else 
                dihAngRight=obj.vertices(vertexIndex).dihedralAngles(faceIndex).ang;
            end
            
            if strcmp(obj.edges(leftEdge).type,'boundary')
                dihAngLeft=0;
            else 
                dihAngLeft=obj.vertices(vertexIndex).dihedralAngles(faceIndex+1).ang;
            end
            
                %no symmetric folding -> just use the vector between the
                %connected vertices in undeformed configuration
                xAxisRight=obj.awayVector0(vertexIndex,rightEdge);
                
                leftEdgeLocal=faceIndex+1;
                leftEdge=obj.vertices(vertexIndex).outgoingEdges(leftEdgeLocal);
                
                xAxisLeft=obj.awayVector0(vertexIndex,leftEdge);
            
            %setup the "right" edge
            obj.propagateAlongEdge(rightEdge,xAxisRight,dihAngRight);
            
            %...and the "left" edge
            obj.propagateAlongEdge(leftEdge,xAxisLeft,dihAngLeft);
            
            %Create an array with faces "to be determined" in
            %counterclockwise direction
            arrayleft=(faceIndex+1):obj.vertices(vertexIndex).degree+1;
            
            %same in clockwise direction
            arrayright=faceIndex-1:-1:1;
            
            %save a copy of the z-axis for later
            zAxisOld=zAxis;
            
            %Go counterclockwise and set all encountered faces accordingly
            xAxis=xAxisLeft;
            for j=1:size(arrayleft,2)
                %arrayleft(j) is the "sector" index
                %arrayleft(j)-1 is the index of the dihedral angle "before" the
                %sector face
                %arrayleft(j) is the index of the dihedral angle "after"
                %the sector face
                %arrayleft(j) is the index of the edge "before" the sector
                %face in outgoingEdges
                %arrayleft(j)+1 is the index of the edge "after" the sector
                %face
                
                %Get the dihedral and sector angle in this sector
                dihedral=obj.vertices(vertexIndex).dihedralAngles(arrayleft(j)).ang;
                sector=obj.vertices(vertexIndex).sectorAngles(arrayleft(j));
                
                %Get the normal vector of the next face
                zAxis=rotateAboutArbitrary(dihedral,xAxis)*zAxis;
                obj.setFaceNormal(obj.vertices(vertexIndex).adjacentFaces(...
                    arrayleft(j)),zAxis);
                
                %Get the crease vector of the next crease
                xAxis=rotateAboutArbitrary(sector,zAxis)*xAxis;
                
                %Get the dihedral angle of the next crease
                if j==size(arrayleft,2)
                    nextDihAng=0;
                else
                    nextDihAng=obj.vertices(vertexIndex).dihedralAngles(arrayleft(j)+1).ang;
                end
                
                nextEdgeIndexLocal=arrayleft(j)+1;
                nextEdge=obj.vertices(vertexIndex).outgoingEdges(nextEdgeIndexLocal);
                obj.propagateAlongEdge(nextEdge,xAxis,nextDihAng);
            end
            
            %Go clockwise and set all encountered faces accordingly.
            zAxis=zAxisOld;
            xAxis=xAxisRight;
            
            for j=1:size(arrayright,2)
                %Get the dihedral and sector angle in this sector
                dihedral=obj.vertices(vertexIndex).dihedralAngles(arrayright(j)+1).ang;
                sector=obj.vertices(vertexIndex).sectorAngles(arrayright(j)+1);
                
                %Get the normal vector of the next face
                zAxis=rotateAboutArbitrary(-dihedral,xAxis)*zAxis;
                obj.setFaceNormal(obj.vertices(vertexIndex).adjacentFaces(...
                    arrayright(j)),zAxis);
                
                %Get the crease vector of the next crease
                xAxis=rotateAboutArbitrary(-sector,zAxis)*xAxis;
                
                %Get the dihedral angle of the next crease
                if j==size(arrayright,2)
                    nextDihAng=0;
                else
                    nextDihAng=obj.vertices(vertexIndex).dihedralAngles(arrayright(j+1)).ang;
                end
                
                nextEdgeIndexLocal=arrayright(j); %the local face index is equal to the local edge index of the edge "before" it
                nextEdge=obj.vertices(vertexIndex).outgoingEdges(nextEdgeIndexLocal);
                obj.propagateAlongEdge(nextEdge,xAxis,nextDihAng);
            end
            
            %set vertex to determined
            obj.vertices(vertexIndex).determined=true;
            
        end %function determinesourcevector        
        function [] = determineAuxVertex(obj, vertexIndex)
            %DETERMINEAUXVERTEX(OBJ, VERTEXINDEX) Determines/evaluates an
            %auxiliary vertex.
            
            %Initial check (for double-checking)
            if not(obj.isAuxVertexDeterminable(vertexIndex))
                error('Cannot determine non-determinable auxiliary vertex');
            end
            
            %Nothing to evaluate about an auxiliary vertex...
            
            %Set vertex to determined
            obj.vertices(vertexIndex).determined=true;
            
        end
        function [] = determineBoundaryVertex(obj, vertexIndex)
            %DETERMINEBOUNDARYVERTEX(OBJ, VERTEXINDEX) Determines/evaluates
            %a boundary vertex. Evaluates the two adjacent boundary
            %creases.
            
            %Initial check
            if not(obj.isBoundaryVertexDeterminable(vertexIndex))
                error('Cannot determine non-determinable vertex');
            end
            
            %sector angles
            sectors=obj.vertices(vertexIndex).sectorAngles;
            
            %Deal with the first outgoing edge
            
            %Get direction vector of last incoming crease
            creaseIndex=obj.vertices(vertexIndex).incomingEdges(end);
            xAxis=obj.awayVector(vertexIndex,creaseIndex);
            
            %Get normal vector to last adjacent face
            zAxis=obj.vertices(vertexIndex).normalVectors(end).vec;
            
            %Get direction vector of outgoing crease from last sector angle
            sectorAngle=sectors(end);
            xAxis=rotateAboutArbitrary(sectorAngle,zAxis)*xAxis;
            
            %Propagate along first outgoing edge
            creaseIndex=obj.vertices(vertexIndex).outgoingEdges(1);
            obj.propagateAlongEdge(creaseIndex,xAxis,0);
            
            %if the vertex is on the boundary (i.e. an auxiliary-boundary 
            %vertex), there is nothing else to do
            if obj.isOnBoundary(vertexIndex)
                obj.vertices(vertexIndex).determined=true;
                return;
            end
            
            %Same story with second outgoing edge. This time, obtain
            %information from first incoming crease (mind rotation
            %direction!)
            creaseIndex=obj.vertices(vertexIndex).incomingEdges(1);
            xAxis=obj.awayVector(vertexIndex,creaseIndex);
            
            zAxis=obj.vertices(vertexIndex).normalVectors(1).vec;
            
            sectorAngle=sectors(1);
            xAxis=rotateAboutArbitrary(-sectorAngle,zAxis)*xAxis;
            
            creaseIndex=obj.vertices(vertexIndex).outgoingEdges(end);
            obj.propagateAlongEdge(creaseIndex,xAxis,0);
            
            %Set vertex to determined
            obj.vertices(vertexIndex).determined = true;
            
        end
        function [] = determineInternalVertex(obj, vertexIndex)
            %DETERMINEINTERNALVERTEX(OBJ, VERTEXINDEX) Determines an
            %internal vertex by calculating, for all adjacent creases, the
            %direction vector and dihedral angle and the normal vectors to
            %the adjacent faces.
            %Throws an error if not rigidly foldable
            
            %Initial Check
            if not(obj.isInternalVertexDeterminable(vertexIndex))
                error('Cannot determine non-determinable vertex');
            end
            
            %some containers needed later on
            deg=obj.vertices(vertexIndex).degree;
            incidentEdges=[obj.vertices(vertexIndex).incomingEdges, ...
                obj.vertices(vertexIndex).outgoingEdges];
            
            %Call PTU to get all dihedral angles
            numInputs=size(obj.vertices(vertexIndex).incomingEdges,2);
            
            sectors=obj.vertices(vertexIndex).sectorAngles;
            
            sectors1=sectors(1:numInputs+1);
            sectors2=sectors(numInputs+2);
            sectors3=sectors(numInputs+3);
            dihedrals1=[obj.vertices(vertexIndex).dihedralAngles(1:numInputs).ang];
            try
            [d1, d2, d3]=PTU(obj.vertices(vertexIndex).mode,sectors1,...
                dihedrals1,sectors2,[],sectors3,[]);
            catch MExc
                MExc
                error('Can not fold Vertex nr. %d',vertexIndex);
            end
            dihedrals=[dihedrals1, d3, d1, d2];
            
            %Get initial x-axis
            lastInput=obj.vertices(vertexIndex).incomingEdges(end);
            xAxis=obj.awayVector(vertexIndex,lastInput);
            
            %Get initial z-axis
            zAxis=obj.vertices(vertexIndex).normalVectors(numInputs).vec;
            
            %iterate over all outgoing edges
            for j=numInputs:deg-1
                %j is equivalent to local index of the crease (the one at the
                %beginning of this sector)
                %j+1 is equivalent to the local index of the sector
                
                %find next z-Axis
                zAxis=rotateAboutArbitrary(dihedrals(j),xAxis)*zAxis;
                
                %set normal to current sector
                faceIndex=obj.vertices(vertexIndex).adjacentFaces(j+1);
                obj.setFaceNormal(faceIndex,zAxis);
                
                %find next x-Axis
                xAxis=rotateAboutArbitrary(sectors(j+1),zAxis)*xAxis;
                
                %propagate along next edge
                edgeIndex=incidentEdges(j+1);
                obj.propagateAlongEdge(edgeIndex,xAxis,dihedrals(j+1));
                
            end
            
            %set vertex to determined
            obj.vertices(vertexIndex).determined=true;
            
        end
        function [] = propagateAlongEdge(obj, edgeIndex, creaseVec, dihedralAngle)
            %PROPAGATEALONGEDGE(OBJ, EDGEINDEX, CREASEVEC, DIHEDRALANGLE)
            %Uses the information given to the edge (direction vector and
            %dihedral angle) to set all associated properties on the edge
            %itself and the two connected vertices.
            
            if strcmp(obj.edges(edgeIndex).type,'boundary')
                %this is a boundary edge
                if nargin<3
                    error('Not enough input arguments');
                end
                
                obj.edges(edgeIndex).directionVector.vec=creaseVec;
                obj.propagateAlongBoundaryEdge(edgeIndex);
                
            elseif strcmp(obj.edges(edgeIndex).type,'internal')
                %this is an internal edge
                if nargin<4
                    error('Not enough input arguments');
                end
                
                obj.edges(edgeIndex).directionVector.vec=creaseVec;
                obj.edges(edgeIndex).dihedralAngle.ang=dihedralAngle;
                obj.propagateAlongInternalEdge(edgeIndex);
                
            else
                error('Unknown edge type');
            end
            
        end
        function [] = propagateAlongInternalEdge(obj, edgeIndex)
            %PROPAGATEALONGINTERNALEDGE(OBJ, EDGEINDEX) Uses the information given
            %to an edge (i.e. direction vector, length and dihedral angle)
            %to pass information along it. Determines location of target
            %vertex and the input dihedral angle corresponding to this
            %edge.
            
            currentEdge=obj.edges(edgeIndex);
            
            if strcmp(currentEdge.type,'boundary')
                error('Wrong function for propagation along edge');
            end
            
            %Check if direction vector is available and is a unit vector
            if abs(norm(currentEdge.directionVector.vec,2)-1)>1e-6
                error('Direction vector must be specified and be a unit vector.');
            end
            
            %calculate end position
            originPosition=obj.vertices(currentEdge.sourceVertex).position;
            endPosition=originPosition+currentEdge.length*currentEdge.directionVector.vec;
            
            %write position to end-vertex
            obj.vertices(currentEdge.endVertex).position=endPosition;
            
            %Set the current edge to "determined"
            obj.edges(edgeIndex).determined=true;
            
        end
        function [] = propagateAlongBoundaryEdge(obj, edgeIndex)
            %PROPAGATEALONGBOUNDARYEDGE(OBJ, EDGEINDEX) Uses the information given
            %to an edge (i.e. direction vector and length) to pass 
            %information along it. Determines location of target vertex
            %corresponding to this edge.
            
            currentEdge=obj.edges(edgeIndex);
            
            if strcmp(currentEdge.type,'internal')
                error('Wrong function for propagation along edge');
            end
            
            %Check if direction vector is available and is a unit vector
            if abs(norm(currentEdge.directionVector.vec,2)-1)>1e-6
                error('Direction vector must be specified and be a unit vector.');
            end
            
            %calculate end position
            originPosition=obj.vertices(currentEdge.sourceVertex).position;
            endPosition=originPosition+currentEdge.length*currentEdge.directionVector.vec;
            
            %write position to end-vertex
            obj.vertices(currentEdge.endVertex).position=endPosition;
            
            %Set the current edge to "determined"
            obj.edges(edgeIndex).determined=true;
            
        end
    

        
%%%%%%%%%Checks%%
        function [] = checkState(obj,statestr)
            %CHECKSTATE(OBJ, STATESTR) Checks whether the origami is in the
            %state defined by statestr. Returns error if this is not the
            %case.
            
            %free pass if isOptimizing is active
            if obj.isOptimizing
                return;
            end
            
            if strcmp(obj.status,statestr)
                %The origami is in the correct state
                return;
            elseif strcmp(obj.status,'editGraph')
                
                %Is-state: editGraph
                if strcmp(statestr,'editFaces')
                    error('The origami object has status editGraph but the requested action requires status editFaces. Did you call the function createFaces?');
                elseif strcmp(statestr,'doKinematics')
                    error('The origami object has status editGraph but the requested action requires status doKinematics. Use the function createFaces to get to the status editFaces and from there use calculateSectorLengths to obtain status doKinematics.');
                end
                    
            elseif strcmp(obj.status,'editFaces')
                
                %Is-state: editFaces
                if strcmp(statestr,'editGraph')
                    error('The origami object has status editFaces but the requested action requires status editGraph. Use the function goToEditGraph to return to editGraph.');
                elseif strcmp(statestr,'doKinematics')
                    error('The origami object has status editFaces but the requested action requires status doKinematics. Make sure to fix one face (function fixFace) and use the function calculateSectorLengths to proceed to status doKinematics.');
                end
                
            elseif strcmp(obj.status, 'doKinematics')
                
                %Is-state: doKinematics
                if strcmp(statestr, 'editGraph')
                    error('The origami object has status doKinematics but the requested action requires status editGraph. Use the function goToEditGraph to reset the origami to the status editGraph.');
                elseif strcmp(statestr, 'editFaces')
                    error('The origami object has status doKinematics but the requested action requires status editFaces. Use the function goToEditFaces to reset the origami to the status editFaces.');
                end
                
            end
            
        end %function checkState
        %check to avoid duplicate vertices
        function [] = checkIsPositionOccupied(obj, pos)
            %CHECKISPOSITIONOCCUPIED(OBJ, POS, REGION) verifies that if a 
            %vertex was added at pos, this vertex would not be "on top of" 
            %or very close to an existing vertex. Searches within the
            %specified region (vertices from different regions may
            %overlap).
            
            %Loop over all existing vertices
            
            for j=1:obj.numberOfVertices()
                
                %skip deleted vertices
                if obj.vertices(j).deleted
                    continue;
                end
                
                %calculate the length between them
                dis=obj.vertices(j).position0(1:2)-pos(1:2);
                
                if norm(dis,2)<1e-9
                    error('matlabPTU:otherVertexTooClose',...
                        'Vertex Nr. %d is too close to new vertex',j);
                end
                
            end
            
        end
        %check to avoid duplicate edges
        function [] = checkIsNewCreaseAllowed(obj, vertexIndex, edgeIndex)
            %CHECKISOUTGOINGCREASEALLOWED(OBJ, VERTEXINDEX, EDGEINDEX)
            %checks whether the new crease would lie at a very low (near-0)
            %angle to an existing crease
            
            %Get existing creases around the central vertex
            incident=[obj.vertices(vertexIndex).incomingEdges, ...
                obj.vertices(vertexIndex).outgoingEdges];
            
            %Get this edge vector, pointing away from the vertex
            thisVec=obj.awayVector0(vertexIndex, edgeIndex);
            thisZAxis=obj.vertices(vertexIndex).unfoldedZ;
            
            %loop over them
            for j=1:size(incident,2)
                thisEdgeIndex=incident(j);
                
                %ignore deleted creases
                if obj.edges(thisEdgeIndex).deleted
                    continue;
                end
                
                %calculate angle inbetween
                thisEdgeVec=obj.awayVector0(vertexIndex,thisEdgeIndex);
                
                ccwAngle=counterClockWiseAngle(thisVec,...
                    thisEdgeVec,thisZAxis);
            
                if abs(mod(real(ccwAngle),360))<1e-6
                    error('The new crease would form a very small angle with crease Nr. %d and is therefore not allowed',thisEdgeIndex);
                end
                
            end
            
        end
        
        function verdict = isLocalSource(obj,vertexIndex)
            %VERDICT=ISLOCALSOURCE(OBJ,VERTEXINDEX) Checks whether the
            %vertex is a local source.
            
            %verify input
            obj.verifyVertexIndex(vertexIndex);
            
            %If the vertex has type source, then it's all clear-> true
            if strcmp(obj.vertices(vertexIndex).type,'source')
                verdict=true;
            else 
                verdict=false;
            end
        end
        
        
%%%%%%%%%Auxiliary functions for plotting%%
        function format = vertexTypeFormatFromIndex(obj,vertexIndex)
            %VERTEXFORMATFROMINDEX Auxiliary function to plotOrigamiGraph.
            %Returns the format entry to plot the vertex at vertexIndex.
            
            type=obj.vertices(vertexIndex).type;
            
            if strcmp(type,'interior')
                format='.b';
            elseif strcmp(type,'boundary')
                format='.r';
            elseif strcmp(type,'source')
                format='.k';
            elseif strcmp(type,'aux')
                format='.g';
            else
                error('invalid vertex type!');
            end
        end
        function format = vertexModeFormatFromIndex(obj,vertexIndex)
            %VERTEXFORMATFROMINDEX Auxiliary function to plotOrigamiGraph.
            %Returns the format entry to plot the vertex at vertexIndex.
            
            if obj.vertices(vertexIndex).mode
                format='ro';
            else
                format='rx';
            end
            
        end
        function format = creaseFormatFromIndex(obj,edgeIndex)
            %CREASEFORMATFROMINDEX Auxiliary function to plotOrigamiGraph.
            %Returns format of a solid line for internal creases and dashed
            %lines for boundary creases.
            
            type=obj.edges(edgeIndex).type;
            if strcmp(type,'internal')
                format='r-';
            elseif strcmp(type,'boundary')
                format='r--';
            else
                error('Invalid Crease Type!');
            end
            
        end
    %%%%%Auxiliary functions for printing the crease pattern%%%%%
        function [] = plotCreasePattern(obj,ax)
            %PLOTCREASEPATTERN(OBJ,AX) Plots the crease pattern of the origami
            %s.t. it can be exported to pdf and printed out.
            %Plots into the specified axis.
            %Mountain creases are shown as strong lines and valley creases
            %as dashed lines.
            
            %origami object must be in state doKinematics and fully
            %determined.
            obj.checkState('doKinematics');
            
            if not(obj.isOrigamiDetermined())
                error('Determine the origami completely before plotting the crease pattern.');
            end
            
            %hold on to axis
            hold(ax,'on');
            
            %loop over all edges
            for e=obj.edges
                
                if e.deleted
                    %skip deleted edges
                    continue;
                end
                
                startPos=obj.vertices(e.sourceVertex).position0(1:2);
                endPos=obj.vertices(e.endVertex).position0(1:2);
                format=obj.creaseFormatForPattern(e.index);
                
                Xcoords=[startPos(1),endPos(1)];
                Ycoords=[startPos(2),endPos(2)];
                
                %plot the line
                line(Xcoords,Ycoords,'LineStyle',format,'Color','k');
                hold on;
                
            end
            
        end
        function [] = plotCreasePatternDouble(obj,ax,dist)
            %PLOTCREASEPATTERNDOUBLE(OBJ,AX,DIST) is similar to
            %plotCreasePattern, but plots each line twice, with a distance
            %of dist where the line that plotCreasePattern would plot is in
            %their middle. All lines are plotted as strong lines
            
            %origami object must be in state doKinematics and fully
            %determined.
            obj.checkState('doKinematics');
            
            if not(obj.isOrigamiDetermined())
                error('Determine the origami completely before plotting the crease pattern.');
            end
            
            %hold on to axis
            hold(ax,'on');
            
            %loop over all edges
            for e=obj.edges
                
                if e.deleted
                    %skip deleted edges
                    continue;
                end
                
                dirVec=obj.awayVector0(e.sourceVertex,e.index);
                leftwards=rotz(90)*[dirVec(1:2);0];
                startPos1=[obj.vertices(e.sourceVertex).position0(1:2);0]+0.5*dist*leftwards;
                endPos1=[obj.vertices(e.endVertex).position0(1:2);0]+0.5*dist*leftwards;
                
                Xcoords=[startPos1(1),endPos1(1)];
                Ycoords=[startPos1(2),endPos1(2)];
                
                %plot the line
                line(Xcoords,Ycoords,'LineStyle','-','Color','k');
                hold on;
                
                startPos2=[obj.vertices(e.sourceVertex).position0(1:2);0]-0.5*dist*leftwards;
                endPos2=[obj.vertices(e.endVertex).position0(1:2);0]-0.5*dist*leftwards;
                
                Xcoords=[startPos2(1),endPos2(1)];
                Ycoords=[startPos2(2),endPos2(2)];
                
                %plot the line
                line(Xcoords,Ycoords,'LineStyle','-','Color','k');
                
            end
            
        end
        function format = creaseFormatForPattern(obj,edgeIndex)
            %FORMAT=CREASEFORMATFORPATTERN(OBJ,EDGEINDEX) Returns the
            %correct format of the specified crease for plotting the crease
            %pattern. Legend:
            %0 fold angle: 'none' (not shown in crease pattern)
            %dihedral angle > 0: '--' (dashed, valley fold)
            %dihedral angle < 0: '-' (solid line, mountain fold)
            %boundary edge: '-.' (dash-dotted)
            
            if not(obj.edges(edgeIndex).determined)
                error('Undetermined edge encountered. Check edge nr. %d.',...
                    edgeIndex);
            end
            
            %return correct format
            if strcmp(obj.edges(edgeIndex).type,'boundary')
                format='-.';
                return;
            elseif strcmp(obj.edges(edgeIndex).type,'internal')
                
                if abs(obj.edges(edgeIndex).dihedralAngle.ang)<1e-6
                    format='none';
                    return;
                elseif obj.edges(edgeIndex).dihedralAngle.ang>0
                    format='--';
                    return;
                elseif obj.edges(edgeIndex).dihedralAngle.ang<0
                    format='-';
                    return;
                else
                    error('Failed to obtain format from edge index. Check edge nr. %d.',...
                        edgeIndex);
                end
            else
                error('Unknown edge type. Check edge nr. %d.',edgeIndex);
            end
            
        end
        
    %%%%%Auxiliary function to gather triangulation of all subOrigamis%%%%%
        function [points, connections, facets, lastVertexIndex] = ...
                gatherTriangulation(obj,points,connections,facets,currentLast)
            %[POINTS,CONNECTIONS,FACETS,LASTVERTEXINDEX]=
            %GATERHTRIANGULATION(POINTS,CONNECTIONS,FACETS,CURRENTLAST)
            %Returns the positions of all points, all connections (i.e.
            %edges), all facets and the number of involved vertices for the
            %complete origami, i.e. including the source vertex
            
            %get the indices of the undeleted vertices in this origami
            undeleted=[obj.vertices(not([obj.vertices(:).deleted])).index];
            
            %add the undeleted vertices to the positions vector
            points=[points,[obj.vertices(undeleted).position]];
            
            %get the connectivitiy of the faces
            connectivity={obj.faces(not([obj.faces(:).deleted]))...
                .connectedVertices};
            
            %update the indices
            for j=1:size(connectivity,2)
                [~,correctedIndices]=ismember(connectivity{j},undeleted);
                connectivity{j}=correctedIndices+currentLast;
            end
            
            %add to the facets container
            facets=cat(2,facets,connectivity);
            
            %get the connectivity of all edges
            edgeConn=[obj.edges(not([obj.edges(:).deleted])).sourceVertex;...
                obj.edges(not([obj.edges(:).deleted])).endVertex];
            
            [~,correctedEdgeConn]=ismember(edgeConn,undeleted);
            correctedEdgeConn=correctedEdgeConn+currentLast;
            
            connections=[connections,correctedEdgeConn];
            
            lastVertexIndex=currentLast+length(undeleted);
            
            
        end
        
    %%%%%Arbitrary utilities for auxiliary vertices%%%%%
        
        function [] = moveAuxAway(obj,auxIndex,vertexIndex)
            %MOVEAUXAWAY(OBJ,AUXINDEX,VERTEXINDEX) moves the vertex with
            %index auxIndex away from the vertex with index vertexIndex.
            %Potentially useful during optimization.
            
            %Intended to be used with status editFaces
            obj.checkState('editFaces');
            
            %find the other vertex that's a neighbour to the specified aux
            %vertex
            otherLocalIndex=find(obj.vertices(auxIndex).incomingNeighbours...
                ~=vertexIndex);
            if size(otherLocalIndex,2)~=1
                error('Identification of other adjacent vertex not successful');
            end
            
            otherVertex=obj.vertices(auxIndex).incomingNeighbours(otherLocalIndex);
            
            %Get the vector that connects the otherVertex with vertexIndex
            connVec0=obj.vertices(vertexIndex).position0-...
                obj.vertices(otherVertex).position0;
            connVec=unitVec(connVec0);
            
            %New position for the aux Vertex - a small distance away from
            %the other vertex
            if norm(connVec0,2)>0.02 %the two vertices are far enough apart
                newPos=obj.vertices(otherVertex).position0+0.01*connVec;
            else %the two vertices are very close to each other -> move vertex to middle
                newPos=obj.vertices(otherVertex).position0+0.5*connVec0;
            end
            
            %Move the aux vertex to the new position
            obj.moveVertexTo(auxIndex,newPos);
            
            %Recalculate all affected sector angles and edge lengths
            obj.resetAuxVertexRecalculate(auxIndex);
            
        end

        function [] = resetAuxVertex(obj,auxIndex)
            %RESETAUXVERTEX(OBJ,AUXINDEX) Resets the auxiliary vertex
            %specified in the input to be exactly between the two adjacent
            %vertices.
            
            if not(strcmp(obj.vertices(auxIndex).type,'aux'))
                error('The vertex with index %d is not an auxiliary vertex. Can not reset this vertex.',auxIndex);
            end
            
            adjVertices=obj.vertices(auxIndex).incomingNeighbours;
            if not(size(adjVertices,2)==2)
                error('Vertex %d is marked as auxiliary vertex but does not have exactly two incoming vertices (found: %d).',auxIndex,size(adjVertices,2));
            end
            
            adjacent1Pos=obj.vertices(adjVertices(1)).position0(1:2);
            adjacent2Pos=obj.vertices(adjVertices(2)).position0(1:2);
            newPos=0.5*(adjacent1Pos+adjacent2Pos);
            
            %move the vertex
            obj.moveVertexTo(auxIndex,newPos);            
        end
        function [] = resetAuxVertexRecalculate(obj,auxIndex)
            %RESETAUXVERTEXRECALCULATE(OBJ,AUXINDEX) Recalculates all edge
            %lengths and sector angles that could have changed while moving
            %the auxiliary vertex
            %Call this function after moving auxiliary vertices with
            %resetAuxVertex
            
            adjVertices=obj.vertices(auxIndex).incomingNeighbours;
            
            obj.calculateSectorAngle(auxIndex);
            obj.calculateSectorAngle(adjVertices(1));
            obj.calculateSectorAngle(adjVertices(2));
            
            obj.calculateEdgeLength(obj.vertices(auxIndex).incomingEdges(1));
            obj.calculateEdgeLength(obj.vertices(auxIndex).incomingEdges(2));
            
        end
        
    end
end