%% License information

%This file is part of the matlabPTU framework.

% © 2021 ETH Zurich, Andreas Walker, D-MAVT, EDAC

%Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

%The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

%%
function [isBoundary] = isBoundary(setOfEdges)
%ISBOUNDARY(SETOFEDGES) determines wheter setOfEdges contains a boundary
%edge. setOfEdge is typically an output of (counter)clockwiseMarching.

%Is there a 0 before the end?
if ismember(0,setOfEdges(1:end-1))
    error('The provided set of edges contains a zero at a forbidden location');
end

isBoundary=setOfEdges(end)==0;
end

