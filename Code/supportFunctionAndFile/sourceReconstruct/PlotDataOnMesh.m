function hp=PlotDataOnMesh(mesh,data,varargin)
% function hp=PlotDataOnMesh(mesh,data,varargin)
% Plots scalar field on a mesh.
% mesh: a struct with fields mesh.p for points and mesh.e for element
%       descriptions
% data: a Nx1 array, where N is the number of points in the mesh.
%
% varargin = options and settings for the plot
%   - a struct with optional fields:
%       - 'colormap', <colormap array> or <number of colors for jet colormap>
%       - 'caxis', <colormap axis, either [min, max] or [absmax]>
%       - 'figure', <figure number or handle>
%       - 'position', <figure position as in set(hf, 'position')>
%       - 'inflated', <use inflated mesh? 1 or 0 (mesh.pinf)>
%       - 'colorbar', <plot colorbar? 1 or 0>
%       - 'view',   <view angles in a format taken by view()-function>
%       - 'pointset', <[N x 3]-array of points or [N x 1]-array of indexes
%           to vertices to be plotted>
%       - 'pointstyle', <How the plotted points should look like, e.g., 'k.'>
%       - 'pointsize', <Size of the plotted points
%    - Or, if no struct given, parameter-value pairs as described above.
% If any struct is given as input, no other arguments are parsed.
%
% version 160923 (stand alone)
% (c) Matti Stenroos (matti.stenroos@aalto.fi)

%if arguments are passed as varargin from another function, they need to be
%converted...
if iscell(varargin) && length(varargin)==1 && iscell(varargin{1})
    varargin=varargin{1};
end
plotoptions=ParseOptions(varargin);

%make sure data is given the correct way
if size(data,1)==1;
    data=data';
end
%check that the data is of correct size...
if size(data,1)~=size(mesh.p,1),
    error('PlotDataOnMesh: The sizes of mesh and datavector do not match.');
end
%OK, let's start...

%check some settings
if isempty(plotoptions.figure),plotoptions.figure=gcf;
else figure(plotoptions.figure);end
if ~isempty(plotoptions.position),set(figure(plotoption.figure),'position',plotoption.position);end
if plotoptions.inflated,p=mesh.pinf;else p=mesh.p;end

%this allow good-quality isocolor surface visualization:
set(gcf,'renderer','zbuffer');if ~ishold,cla;end
%...make the plot
hp = patch('faces',mesh.e,'vertices',p,'facevertexcdata',data,'facecolor','interp','edgecolor',plotoptions.edgecolor);
%set colormap
if size(plotoptions.colormap,2)==3
    colormap(plotoptions.colormap)
else
    colormap(jet(plotoptions.colormap));
end
caxis(SetPlotScale(data,plotoptions.caxis));%set scale

if ~isempty(plotoptions.pointset)%deal with the pointset
    holdstate=ishold;
    if ~holdstate
        hold on;
    end
    if min(size(plotoptions.pointset))==1
        plotoptions.pointset=p(plotoptions.pointset,:);
    end
    hpo=plot3(plotoptions.pointset(:,1),plotoptions.pointset(:,2),plotoptions.pointset(:,3),plotoptions.pointstyle);
    set(hpo,'MarkerSize',plotoptions.pointsize);
    if ~holdstate
        hold off;
    end
end

if plotoptions.colorbar%add colorbar
    colorbar;
end
view(plotoptions.view);%set view

%and make everything look nice...
axis tight equal off;material dull;lighting gouraud;
if isempty(findobj(gca,'Type','light'));
    camlight
end

function plotoptions=ParseOptions(varargin)
if iscell(varargin) && length(varargin)==1 && iscell(varargin{1})
    varargin=varargin{1};
end
%defaults
plotdef.caxis=[];
plotdef.colorbar=1;
plotdef.colormap=10;
plotdef.edgecolor='none';
plotdef.figure=[];
plotdef.inflated=0;
plotdef.pointset=[];
plotdef.pointsize=10;
plotdef.pointstyle='k.';
plotdef.position=[];
plotdef.view=[-90 0];

%parse options
plotoptions=GetStruct(varargin);%try, if there is any struct...?
if ~isstruct(plotoptions)%if not, check the parameters...
    plotoptions=struct('kind','plotoptions');
    if IsParameter(varargin,'caxis'),plotoptions.caxis=GetValue(varargin,'caxis');end
    if IsParameter(varargin,'colorbar'),plotoptions.colorbar=GetValue(varargin,'colorbar');end
    if IsParameter(varargin,'colormap'),plotoptions.colormap=GetValue(varargin,'colormap');end
    if IsParameter(varargin,'edgecolor'),plotoptions.edgecolor=GetValue(varargin,'edgecolor');end
    if IsParameter(varargin,'figure'),plotoptions.figure=GetValue(varargin,'figure');end
    if IsParameter(varargin,'inflated'),plotoptions.inflated=GetValue(varargin,'inflated');end
    if IsParameter(varargin,'pointset'),plotoptions.pointset=GetValue(varargin,'pointset');end
    if IsParameter(varargin,'pointsize'),plotoptions.pointsize=GetValue(varargin,'pointsize');end
    if IsParameter(varargin,'pointstyle'),plotoptions.pointstyle=GetValue(varargin,'pointstyle');end
    if IsParameter(varargin,'position'),plotoptions.position=GetValue(varargin,'position');end
    if IsParameter(varargin,'view'),plotoptions.view=GetValue(varargin,'view');end
end
%fill missing option values with defaults.
if ~isfield(plotoptions,'caxis'),plotoptions.caxis=plotdef.caxis;end
if ~isfield(plotoptions,'colorbar'),plotoptions.colorbar=plotdef.colorbar;end
if ~isfield(plotoptions,'colormap'),plotoptions.colormap=plotdef.colormap;end
if ~isfield(plotoptions,'edgecolor'),plotoptions.edgecolor=plotdef.edgecolor;end
if ~isfield(plotoptions,'figure'),plotoptions.figure=plotdef.figure;end
if ~isfield(plotoptions,'inflated'),plotoptions.inflated=plotdef.inflated;end
if ~isfield(plotoptions,'pointset'),plotoptions.pointset=plotdef.pointset;end
if ~isfield(plotoptions,'pointsize'),plotoptions.pointsize=plotdef.pointsize;end
if ~isfield(plotoptions,'pointstyle'),plotoptions.pointstyle=plotdef.pointstyle;end
if ~isfield(plotoptions,'position'),plotoptions.position=plotdef.position;end
if ~isfield(plotoptions,'view'),plotoptions.view=plotdef.view;end

function cscale=SetPlotScale(data,cscale)
% function cscale=SetPlotScale(data,cscale)
% Compute scale for data to be plotted.
% cscale (optional): if empty or omitted, autoscale to the abs(max(plotdata))

if ~isempty(cscale) && length(cscale)==2,
    return;
end
if ~isempty(cscale),
    if all(data>=0)
        cscale=[0 1]*abs(cscale);
    elseif all(data<=0)
        cscale=[-1 0]*abs(cscale);
    else
        cscale=[-1 1]*abs(cscale);
    end
elseif all(data>=0)
    cscale=[0 1]*max(data);
elseif all(data<=0)
    cscale=[1 0]*min(data);
else
    cscale=[-1 1]*max(abs(data));
end

function [res,index]=GetStruct(paralist)
% function [res,index]=GetStruct(paralist)
% Checks, whether there is a struct in the input list and
% returns it and the index to the struct (or 0 and []). If there are many
% structs, this returns only the first one.
N=length(paralist);
res=0;index=[];
for I=1:N,
    if isstruct(paralist{I}),
        res=paralist{I};
        index=I;
        break
    end
end
function [val,IsItThere,Index]=GetValue(input,parameter)
% function [val,IsItThere,Index]=GetValue(input,parameter)
% Gets a value from a parameter/value list
% input: parameter/value list; X1,Y1,X2,Y2,...
% parameter: the parameter to find
% val: value of the parameter

paralist=input;
test=strcmpi(paralist,parameter);
if ~any(test)
    val=[];
    IsItThere=0;
    Index=[];
else
    ind=find(test);
    val=paralist{ind+1};
    IsItThere=1;
    Index=ind;
end

function res=IsParameter(input,parameter)
% function res=IsParameter(input,parameter)
% Checks, whether a value exists in a parameter/value list
% input: parameter/value list; X1,Y1,X2,Y2,...
% parameter: the parameter to find
% res: TRUE or FALSE

% paralist=input(1:2:end);
paralist=input;
test=strcmpi(paralist,parameter);
res=any(test);
    
