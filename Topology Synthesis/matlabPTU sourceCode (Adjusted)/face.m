%% License information

%This file is part of the matlabPTU framework.

% © 2021 ETH Zurich, Andreas Walker, D-MAVT, EDAC

%Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

%The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

%% 
classdef face < matlab.mixin.Copyable
    %FACE Represents a facet in an origami
    %   Part of an origami object
    
    properties
        connectedVertices(1,:) {mustBeInteger, mustBeNonnegative, mustBeFinite} = [];
        %Vector with indices of all connected Vertices. The indices
        %correspond to walking around the face in clockwise
        %direction.
        
        boundaryEdges(1,:) {mustBeInteger, mustBeNonnegative, mustBeFinite} = [];
        %Vector with indices of all edges at the boundary of this facet.
        %Order: counterclockwise.
        
        %regionNumber(1,1) {mustBeInteger,mustBeNonnegative,mustBeFinite};
  
        indicesAtVertices(1,:) {mustBeInteger, mustBeNonnegative, mustBeFinite} = [];
        %Vector where each entry corresponds to one entry in
        %connectedVertices. The entry in this vector states the index of
        %the facet in the adjacentFaces vector of the corresponding vertex.
        %Cf. the vertex class for the indexation convention of facets at
        %vertices.
        
        normalVector(1,1) oriVec;
        %Upwards-pointing normal unit vector to the plane. Refers to the
        %folded/deformed configuration. Set during evaluation/folding
        %phase.
        
        determined(1,1) {mustBeNumericOrLogical} = false;
        %Specifies whether all connected vertices have been determined.
        
        index(1,1) {mustBeInteger, mustBeNonnegative, mustBeFinite} = 0;
        %Index of this face object within the global origami.
        
        fixed(1,1) {mustBeNumericOrLogical} = false;
        %Specifies whether this face is fixed. Only one face in a whole
        %origami is allowed to be fixed.
        
        deleted(1,1) {mustBeNumericOrLogical} = false;
        %Indicates whether the face has been deleted. When a face gets
        %deleted, it is not actually deleted but rather marked as deleted
        %(this is to avoid confusion with indices).
        
    end
    
    methods (Access='public')
        function obj = face(vertexList, edgeList,...
                indicesAtVertices, index)
            %FACE(REGIONNUMBER, VERTEXLIST, EDGELIST, INDICESATVERTICES, 
            %   INDEX) 
            %   Constructs a face - object
            %   First input is mandatory (region number). Rest is not
            %   mandatory.
            %   vertexList (optional, first argument) contains list of the
            %   connected vertices (their indices, resp.)
            %   edgeList (optional, second argument) contains list of the
            %   boundary edges (their indices, resp.)
            
            
            if nargin>=2
                obj.connectedVertices=vertexList;
            end
            
            if nargin>=3
                %ensure the number of edges is identical to the number of
                %vertices
                if size(edgeList, 2) ~= size(vertexList,2)
                    error('Number of edges and number of vertices must be identical');
                end
                obj.boundaryEdges=edgeList;
            end
            
            if nargin>=4
                %ensure that number of indices is identical to the number
                %of faces
                if size(indicesAtVertices,2) ~= size(vertexList,2)
                    error('Number of indices at vertices must be identical to number of vertices (one index per vertex)');
                end
                obj.indicesAtVertices=indicesAtVertices;
            end
            
            if nargin>=5
                obj.index=index;
            end
            
            %explicitly call the constructor to ensure that a new oriVec is
            %created
            obj.normalVector=oriVec();
            
        end
        
    end
    
    methods (Access='protected')
        %%%%%Copying the face%%%%%
        
        function cp = copyElement(obj)
            %CP=COPYELEMENT(OBJ) creates a deep copy of the face object.
            %This overrides the default copy function of the mixin.copyable
            %class. This is necessary in order to create independent normal
            %vectors when copying an origami. Important: after copying, the
            %edges and vertices of the top-origami must be reconnected to
            %this face normal.
            
            cp=face([],[],[],1); %something. Will be destroyed and rebuilt later; this is just for calling the constructor.
            
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
            
        end
    end
end

