%%
classdef vertex < matlab.mixin.Copyable
    %VERTEX Describes a vertex in an origami
    %   is part of an origami object, holds all relevant information on a
    %   vertex.
    properties
        degree(1,1) {mustBeInteger, mustBeFinite, mustBeNonnegative} = 0;
        incomingEdges(1,:) {mustBeInteger, mustBePositive, mustBeFinite} = [];
        outgoingEdges(1,:) {mustBeInteger, mustBePositive, mustBeFinite} = [];
        adjacentFaces(1,:) {mustBeInteger, mustBeNonnegative, mustBeFinite} = [];
        position0(3,1) {mustBeReal, mustBeFinite};
        unfoldedZ(3,1) {mustBeReal, mustBeFinite} = [0;0;1];
        position(3,1) {mustBeReal, mustBeFinite};
        determinable(1,1) {mustBeNumericOrLogical} = false;
        determined(1,1) {mustBeNumericOrLogical} = false;
        normalVectors(1,:) oriVec = oriVec.empty();
        creaseVectors(1,:) oriVec = oriVec.empty();
        mode(1,1) {mustBeNumericOrLogical} = false;
        dihedralAngles(1,:) oriAng = oriAng.empty();
        sectorAngles(1,:) {mustBeReal, mustBeFinite} = [];
        incomingNeighbours(1,:) {mustBeInteger, mustBeFinite} = [];
        outgoingNeighbours(1,:) {mustBeInteger, mustBeFinite} = [];
        type char {mustBeMember(type,{'source','interior','boundary',...
            'aux','deleted'})} = 'boundary';
        %(!) added 'deleted' to types to avoid confusion
        index(1,1) {mustBeInteger, mustBeNonnegative, mustBeFinite};
        deleted(1,1) {mustBeNumericOrLogical} = false;
        
        %%%ADDED BY JIM%%%
        %Add a property with amount of DOF added?
        addedDOF(1,1)   =   0;
        %Add a property about the expandability? T
        T(1,1) {mustBeNumericOrLogical} =   true;
        
        
    end
    methods
        %%%Constructor%%%
        function obj = vertex(pos0, index, mode, type,zVector)
            %regionNumber is omitted
            %(!) add input of addedDOF (if so, change nargin counters in if
            %                           statement)
            
            %position vector will be 3 coords:
            if length(pos0)==2
                pos=[pos0;0];
            elseif length(pos0)==3
                pos=pos0;
            else
                error('The position of the vertex must be specified either as a 2-by-1 or as a 3-by-1 vertex.');
            end
            obj.position0=pos;
            obj.position=pos;   %initializing the deformed position to the 
                                %reference position
            obj.index=index;    %assigning an index
            
            %(!) add input addedDOF?
            %obj.addedDOF = addedDOF
            
            if nargin>=3    %mode specified by third argument
                obj.mode=mode;
            end
            if nargin>=4    %type specified by fourth argument
                obj.type=type;
            end
            if nargin>=5    %normal direction specified in fifth argument
                obj.unfoldedZ=zVector;
            end
        end
    end
end
            