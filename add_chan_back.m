
function EEG = add_chan_back(EEG,chanLabel)

% [EEG.chanlocs.ref] = deal('Cz');
% EEG = pop_chanedit(EEG, 'append',63,'changefield',{64 'labels' 'Cz'},'lookup','eeglabroot/plugins/dipfit3.3/standard_BEM/elec/standard_1005.elc',...
%                    'eval','chans = pop_chancenter( chans, [],[]);','changefield',{64 'type' 'REF'});
% EEG = pop_reref( EEG, [],'refloc',struct('labels',{'Cz'},'type',{'REF'},'theta',{-94.2457},'radius',{0.023785},'X',{-0.54445},'Y',{7.3339},'Z',{98.2349},...
%                     'sph_theta',{94.2457},'sph_phi',{85.7187},'sph_radius',{98.5098},'urchan',{64},'ref',{''},'datachan',{0}));

% Read BEM template coordinates
locs = readeetraklocs('standard_1005.elc');
locs = convertlocs(locs,'cart2all');
locs = rmfield(locs, {'sph_theta' 'sph_theta_besa'}); % for the conversion below

nElecs = length(chanLabel);

% find elec to add back in template list
for iElec = 1:nElecs
    idx(iElec,:) = find(strcmpi({locs.labels},chanLabel{iElec}));

    if isempty(idx(iElec,:))
        warning('Could not find coorindates for this electrode'); return; 
    elseif length(idx) > 1
        error('Several electrodes with th esam elabel were found in the BEM coordinate template. There should be just one.')
    end

    idx


end

%% Subfunctions

% Read 3-D location files saved using the EETrak digitizing software.
% Author: Arnaud Delorme, CNL / Salk Institute, Nov 2003

function chanlocs = readeetraklocs( filename )
    
% read file
locs  = loadtxt( filename );

% get label names
indlabels = [];
indpos    = [];
for ind = 1:size(locs,1)
    if ischar(locs{ind,1})
        if strcmpi(locs{ind,1}, 'Labels')
            indlabels = ind;
        end
        if strcmpi(locs{ind,1}, 'Positions')
            indpos = ind;
        end
    end
end
if isempty(indpos) || isempty(indlabels)
    error('Could not find ''Labels'' or ''Position'' tag in electrode file');
end

% get positions
if strcmp(locs(indpos+1,2),':')
    positions = locs(indpos+1:indlabels-1,3:5);
else
    positions = locs(indpos+1:indlabels-1,1:3);
end
labels    = locs(indlabels+1:end,:);

for index = 1:length(labels)
    chanlocs(index).labels = labels{index};
    chanlocs(index).X      = positions{index,1};
    chanlocs(index).Y      = positions{index,2};
    chanlocs(index).Z      = positions{index,3};
end

chanlocs = convertlocs(chanlocs, 'cart2all');


%% subfunction: convertlocs
% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 22 Dec 2002, arno@salk.edu

function chans = convertlocs(chans, command, varargin)

if nargin < 1
   help convertlocs;
   return;
end

if ~isfield(chans, 'theta') && ~isfield(chans, 'X') && ~isfield(chans, 'radius') && ~isfield(chans, 'sph_theta_besa')
    return
end

if nargin < 2
   command = 'auto';
end
if nargin == 4 && strcmpi(varargin{2}, 'on')
    verbose = 1;
else
    verbose = 0; % off
end

% test if value exists for default
% --------------------------------
if strcmp(command, 'auto')
    if isfield(chans, 'X') && any(~cellfun(@isempty, { chans.X }))
        command = 'cart2all';
        if verbose
            disp('Make all coordinate frames uniform using Cartesian coords');
        end
    else
        if isfield(chans, 'sph_theta') && ~isempty(chans(1).sph_theta)
            command = 'sph2all';
            if verbose
                disp('Make all coordinate frames uniform using spherical coords');
            end
        else
            if isfield(chans, 'sph_theta_besa') && ~isempty(chans(1).sph_theta_besa)
                command = 'sphbesa2all';
                if verbose
                    disp('Make all coordinate frames uniform using BESA spherical coords');
                end
            else
                command = 'topo2all';
                if verbose
                    disp('Make all coordinate frames uniform using polar coords');
                end
            end
        end
    end
end

% convert
% -------         
switch command
 case 'topo2sph'
   theta  = {chans.theta};
   radius = {chans.radius};
   indices = find(~cellfun('isempty', theta));
   [sph_phi, sph_theta] = topo2sph( [ [ theta{indices} ]' [ radius{indices}]' ] );
   if verbose
       disp('Warning: electrodes forced to lie on a sphere for polar to 3-D conversion');
   end
   for index = 1:length(indices)
      chans(indices(index)).sph_theta  = sph_theta(index);
      chans(indices(index)).sph_phi    = sph_phi  (index);
   end
   if isfield(chans, 'sph_radius')
       meanrad = mean([ chans(indices).sph_radius ]);
       if isempty(meanrad)
           [chans(indices).sph_radius] = deal(85);
           meanrad = 85; 
       end
   else
        [chans(indices).sph_radius] = deal(85);
        meanrad = 85;
   end
   sph_radius(1:length(indices)) = {meanrad};
case 'topo2sphbesa'
   chans = convertlocs(chans, 'topo2sph', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2sphbesa', varargin{:}); % search for spherical coords
case 'topo2cart'
   chans = convertlocs(chans, 'topo2sph', varargin{:}); % search for spherical coords
   if verbose
       disp('Warning: spherical coordinates automatically updated');
   end
   chans = convertlocs(chans, 'sph2cart', varargin{:}); % search for spherical coords
case 'topo2all'
   chans = convertlocs(chans, 'topo2sph', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2sphbesa', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2cart', varargin{:}); % search for spherical coords
case 'sph2cart'
   sph_theta  = {chans.sph_theta};
   sph_phi    = {chans.sph_phi};
   indices = find(~cellfun('isempty', sph_theta));
   if ~isfield(chans, 'sph_radius')
        [chans(indices).sph_radius] = deal(85);
       sph_radius(1:length(indices)) = {85};
   else                              
       sph_radius = {chans.sph_radius};
   end
   inde = find(cellfun('isempty', sph_radius));
   if ~isempty(inde)
       meanrad = mean( [ sph_radius{:} ]);
       sph_radius(inde) = { meanrad };
   end
   [x, y, z] = sph2cart([ sph_theta{indices} ]'/180*pi, [ sph_phi{indices} ]'/180*pi, [ sph_radius{indices} ]');
   for index = 1:length(indices)
      chans(indices(index)).X = x(index);
      chans(indices(index)).Y = y(index);
      chans(indices(index)).Z = z(index);
   end
case 'sph2topo'
 if verbose
     % disp('Warning: all radii constrained to one for spherical to topo transformation');
 end
 sph_theta  = {chans.sph_theta};
 sph_phi    = {chans.sph_phi};
 indices = find(~cellfun('isempty', sph_theta));
 [chan_num,angle,radius] = sph2topo([ ones(length(indices),1)  [ sph_phi{indices} ]' [ sph_theta{indices} ]' ], 1, 2); % using method 2
 for index = 1:length(indices)
     chans(indices(index)).theta  = angle(index);
     chans(indices(index)).radius = radius(index);
     if ~isfield(chans, 'sph_radius') || isempty(chans(indices(index)).sph_radius)
         chans(indices(index)).sph_radius = 85;
     end
 end
case 'sph2sphbesa'
   % using polar coordinates
   sph_theta  = {chans.sph_theta};
   sph_phi    = {chans.sph_phi};
   indices = find(~cellfun('isempty', sph_theta));
   [chan_num,angle,radius] = sph2topo([ones(length(indices),1)  [ sph_phi{indices} ]' [ sph_theta{indices} ]' ], 1, 2);
   [sph_theta_besa, sph_phi_besa] = topo2sph([angle radius], 1, 1);
   for index = 1:length(indices)
      chans(indices(index)).sph_theta_besa  = sph_theta_besa(index);
      chans(indices(index)).sph_phi_besa    = sph_phi_besa(index);
   end
case 'sph2all'
   chans = convertlocs(chans, 'sph2topo', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2sphbesa', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2cart', varargin{:}); % search for spherical coords
case 'sphbesa2sph'
   % using polar coordinates
   sph_theta_besa  = {chans.sph_theta_besa};
   sph_phi_besa    = {chans.sph_phi_besa};
   indices = find(~cellfun('isempty', sph_theta_besa));
   [chan_num,angle,radius] = sph2topo([ones(length(indices),1)  [ sph_theta_besa{indices} ]' [ sph_phi_besa{indices} ]' ], 1, 1);
   %for index = 1:length(chans)
   %   chans(indices(index)).theta  = angle(index);
   %   chans(indices(index)).radius = radius(index);
   %   chans(indices(index)).labels = int2str(index);
   %end;   
   %figure; topoplot([],chans, 'style', 'blank', 'electrodes', 'labelpoint');
   
   [sph_phi, sph_theta] = topo2sph([angle radius], 2);
   for index = 1:length(indices)
      chans(indices(index)).sph_theta  = sph_theta(index);
      chans(indices(index)).sph_phi    = sph_phi  (index);      
   end
case 'sphbesa2topo'
   chans = convertlocs(chans, 'sphbesa2sph', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2topo', varargin{:}); % search for spherical coords
case 'sphbesa2cart'
   chans = convertlocs(chans, 'sphbesa2sph', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2cart', varargin{:}); % search for spherical coords   
case 'sphbesa2all'
   chans = convertlocs(chans, 'sphbesa2sph', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2all', varargin{:}); % search for spherical coords
case 'cart2topo'
   chans = convertlocs(chans, 'cart2sph', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2topo', varargin{:}); % search for spherical coords
case 'cart2sphbesa'
   chans = convertlocs(chans, 'cart2sph', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2sphbesa', varargin{:}); % search for spherical coords
case 'cart2sph'
    if verbose
        disp('WARNING: If XYZ center has not been optimized, optimize it using Edit > Channel Locations');
    end
    X  = {chans.X};
    Y  = {chans.Y};
    Z  = {chans.Z};
    indices = find(~cellfun('isempty', X));
    [th, phi, radius] = cart2sph( [ X{indices} ], [ Y{indices} ], [ Z{indices} ]);
	for index = 1:length(indices)
		 chans(indices(index)).sph_theta     = th(index)/pi*180;
		 chans(indices(index)).sph_phi       = phi(index)/pi*180;
		 chans(indices(index)).sph_radius    = radius(index);
	end
case 'cart2all'
   chans = convertlocs(chans, 'cart2sph', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2all', varargin{:}); % search for spherical coords
end
