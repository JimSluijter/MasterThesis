%% License information

%This file is part of the matlabPTU framework.

% © 2021 ETH Zurich, Andreas Walker, D-MAVT, EDAC

%Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

%The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

%%
classdef edge < matlab.mixin.Copyable
    %EDGE Represents an edge/crease in an origami
    %   Part of an origami object
    
    properties
        sourceVertex(1,1) {mustBeInteger, mustBeNonnegative, mustBeFinite};
        %Index of the source Vertex. Must be specified by constructor.
        
        endVertex(1,1) {mustBeInteger, mustBeNonnegative, mustBeFinite};
        %Index of the end vertex. Must be specified by constructor.
        
        regionNumber(1,1) {mustBeInteger, mustBeNonnegative, mustBeFinite};
        %Region to which this edge belongs. Must be specified by
        %constructor.
        
        directionVector(1,1) oriVec = oriVec();
        %Unit vector of vector pointing from source to target. Is assigned
        %when source vertex is determined/"evaluated"
        
        length(1,1) {mustBeReal, mustBeNonnegative} = 0;
        %Length of Crease. Assigned during construction phase.
        
        index(1,1) {mustBeInteger, mustBeNonnegative, mustBeFinite};
        %Index of this crease in the origami. Must be specified in 
        %constructor.
        
        dihedralAngle(1,1) oriAng;
        %Dihedral angle between the two adjacent faces in degrees. Positive
        %values correspond to valley folds, negative values to mountain
        %folds.
        
        faceLeft(1,1) {mustBeInteger, mustBeNonnegative, mustBeFinite} = 0;
        %Index of face to the left of the crease (as seen when standing at
        %the source and looking at the target). Value must be specified
        %during the construction phase.
        
        faceRight(1,1) {mustBeInteger, mustBeNonnegative, mustBeFinite} = 0;
        %Index of face to the right of the crease (as seen when standing at
        %the source and looking at the target). Value must be specified
        %during the construction phase.
        
        normalVecLeft(1,1) oriVec;
        %Normal Vector of the face to the left of the crease (same
        %convention as for the indices of the adjacent faces). Must have
        %length 1 and point "upwards".
        
        normalVecRight(1,1) oriVec;
        %Normal Vector of the face to the right of the crease (same
        %convention as for the indices of the adjacent faces). Must have
        %length 1 and point "upwards".
        
        determined(1,1) {mustBeNumericOrLogical} = false;
        %Specifies whether the crease/edge has been "evaluated", i.e. its
        %direction vector, dihedral angle and adjacent normal vectors have
        %been set.
        
        indexAtSource(1,1) {mustBeInteger, mustBeNonnegative, mustBeFinite} = 0;
        %Specifies the index among the outgoing edges at the source. Cf.
        %the vertex class for enumeration/indexation convention. 
        
        indexAtEnd(1,1) {mustBeInteger, mustBeNonnegative, mustBeFinite} = 0;
        %Specifies the index among the incoming edges at the end. Cf. the
        %vertex class for enumeration/indexation convention.
        
        type char {mustBeMember(type,{'internal','boundary'})} = 'internal';
        %Specifies whether the edge is internal or on the boundary. Default
        %is internal.
        
        deleted(1,1) {mustBeNumericOrLogical} = false;
        %Is set to true to mark that the edge has been deleted. After
        %deletion, the object remains in the containers (and is marked as
        %deleted). The deletion can be reverted or the edge can be
        %overwritten. This strategy helps to avoid indexation problems.
        
    end
    
    methods (Access='public')
        function obj = edge(sourceIndex, endIndex, edgeIndex)
            %EDGE(SOURCEINDEX,ENDINDEX,EDGEINDEX,REGION) Constructs edge object.
            %   Requires Index of source vertex, index of end vertex, edge
            %   index and region number as input.
            obj.sourceVertex = sourceIndex;
            obj.endVertex = endIndex;
            obj.index = edgeIndex;
            
            %explicitly call the constructor to ensure that a new instance
            %is created with each constructor call.
            obj.dihedralAngle=oriAng();
            obj.directionVector=oriVec();
            
        end
        
    end
    
    methods (Access='protected')
    end
    
end

