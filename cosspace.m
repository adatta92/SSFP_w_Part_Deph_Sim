function x = cosspace(startPoint, endPoint, varargin)
%COSSPACE Cosine spaced vector.
%  COSSPACE(X1, X2) generates a row vector of 100 cosine-spaced
%  points between X1 and X2.
%
%  COSSPACE(X1, X2, N) generates N cosine-spaced points between X1 and X2.
%
%  This method of spacing concentrates samples at the ends while
%  producing fewer sample points in the center.
%
%             ** *  *   *    *    *   *  * * **
%             1...                            N
%
%  See also LOGSPACE, LINSPACE


%% See how many points we need.
if ~isempty(varargin)
    numPoints=varargin{1};
else
    numPoints = 100;
end
    
%% We already know the endpoints, so save a couple of calculations.
%We also benefit by preallocating this fixed-length variable.
x(numPoints)=endPoint;
x(1)=startPoint;

%% Find the mid point of the data set
%Do some checking here to save us issues downstream.
if endPoint <= startPoint
    disp('End point must be greater than the start point');
    x = [];
    return;
else
    midPoint = (endPoint-startPoint)/2;
end

%% Calculate the iteration increment
%Each point in the array corresponds to an angle.
angleInc = pi/(numPoints-1);

%% Brute Force way...
%For reference only.  The algorithm is easier to see this way,
%but slower because it doesn't take advantage of symmetry.
% curAngle=angleInc;
% for idx = 2:numPoints-1
%     x(idx)=startPoint+midPoint*(1-cos(curAngle));
%     curAngle=curAngle+angleInc;
% end

%% Alternative way (1/2 the "for loop" iterations)
%The spacing before and after the midpoint is the same.
%We can save some calculations by making use of this symmetry.
curAngle=angleInc;
for idx = 2:ceil(numPoints/2)
    x(idx)=startPoint+midPoint*(1-cos(curAngle));
    x(end-(idx-1))=x(end-(idx-2))-(x(idx)-x(idx-1));
    curAngle=curAngle+angleInc;
end

% by Cole Stephens 
% 14 Jul 2004 

% Copyright (c) 2016, The MathWorks, Inc.
% All rights reserved.

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%    * Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in
%      the documentation and/or other materials provided with the distribution.
%    * In all cases, the software is, and all modifications and derivatives
%      of the software shall be, licensed to you solely for use in conjunction
%      with MathWorks products and service offerings.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
