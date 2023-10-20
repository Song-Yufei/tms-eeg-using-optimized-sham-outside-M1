% EEG = rename_events(EEG, oldmarker, newmarker)
%
%       Renames all events of (a) given type(s) to the indicated new type.
%
% In:
%       EEG - EEGLAB dataset containing event information
%       oldmarker - string or cell of strings giving the type of old events
%                   to be renamed
%       newmarker - string giving the type of the new event type
%
% Out:
%       EEG - EEGLAB dataset with indicated events renamed
%
% Usage example:
%       >> EEG = rename_events(EEG, {'event1a', 'event1b'}, 'event1');
% 
%                       Laurens R. Krol
%                       Neuroadaptive Human-Computer Interaction
%                       Brandenburg University of Technology

% Copyright 2020 Laurens R. Krol

% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the
% following conditions:

% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.

function EEG = rename_events(EEG, oldmarker, newmarker)

idx = ismember({EEG.event.type}, oldmarker);
[EEG.event(idx).type] = deal(newmarker);

fprintf('renamed %d events to ''%s''\n', sum(idx), newmarker);

end