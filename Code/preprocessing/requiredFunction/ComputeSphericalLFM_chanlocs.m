
function LFM_sphere=ComputeSphericalLFM_chanlocs(chanlocs, refChLabel)

%
% This function extract infromation from chanlocs struct from EEGlab data
% struct, and the label of reference channel
% Then LFM_sphere is computed by leadfield1 function
%
% .........................................................................
% 13 October 2023 : Johanna Metsomaa, Aalto university  
% .........................................................................


% THE SPHERICAL HEAD MODEL
r0 = 76*1e-3;
r3=88*1e-3;
rad=[81,85,88]*1e-3;
sig=[.33,.33/100,.33];

     x=cat(1,chanlocs(:).X);
     y=cat(1,chanlocs(:).Y);
     z=cat(1,chanlocs(:).Z);
     elec_coords=[x,y,z];

R=r3*elec_coords./repmat(sqrt(sum(elec_coords.^2,2)), [1,3]); % electrode positions,

DipN=5000;
P=randn(3,DipN);
P=r0*P./(ones(3,1)*sqrt(sum(P.^2))); 
Pns=[P(1,:);P(2,:);(P(3,:))]; % dipole positions on the upper hemisphere
                               % with radius r0,
%Moments of the unit dipoles                               
Q = Pns;                              
Qns = Q./repmat(sqrt(sum(Q.^2,1)),[3,1]);

%Computing the leadfield accoring to [1]
LFM_sphere = leadfield1(R',Pns,Qns,rad,sig,30);

if refChLabel
    chanlookupfile = fullfile(fileparts(which('eeglab')), 'plugins', 'dipfit', 'standard_BESA', 'standard-10-5-cap385.elp');
    chanlocs_temp.labels=refChLabel;
    chanlocs_temp = pop_chanedit(chanlocs_temp, 'lookup', chanlookupfile);
     x=cat(1,chanlocs_temp(:).X);
     y=cat(1,chanlocs_temp(:).Y);
     z=cat(1,chanlocs_temp(:).Z);
     elec_coords=[x,y,z];
     R=r3*elec_coords./repmat(sqrt(sum(elec_coords.^2,2)), [1,3]); % electrode positions,
    LFM_sphere_ref = leadfield1(R',Pns,Qns,rad,sig,30);
    LFM_sphere=LFM_sphere-LFM_sphere_ref;
    
end