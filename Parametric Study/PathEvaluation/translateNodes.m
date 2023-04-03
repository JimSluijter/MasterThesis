function newNode = translateNodes(oldNode,factor)
% by Jim Sluijter
% This function translates the node number of a specific node when the
% nSteps is changed. 
% factor = nStepsNew/nStepsOld
newNode = (oldNode-1)*factor+1;
end