function tissue = makeTissueList(nm)
%function tissueProps = makeTissueList(nm)
%   Returns the tissue optical properties at the wavelength nm:
%       tissueProps = [mua; mus; g]';
%   Uses 
%       SpectralLIB.mat
%
% Steven L. Jacques, March 3, 2013.

%% Load spectral library
load SpectralLIB
%   muadeoxy      701x1              5608  double              
%   muamel        701x1              5608  double              
%   muaoxy        701x1              5608  double              
%   muawater      701x1              5608  double              
%   musp          701x1              5608  double              
%   nmLIB         701x1              5608  double              
MU(:,1) = interp1(nmLIB,muaoxy,nm);
MU(:,2) = interp1(nmLIB,muadeoxy,nm);
MU(:,3) = interp1(nmLIB,muawater,nm);
MU(:,4) = interp1(nmLIB,muamel,nm);
LOADED = 1;

%% Create tissueList

j=1;
tissue(j).name = 'air';
tissue(j).mua = 0.1;
tissue(j).mus = 10.0;
tissue(j).g   = 1.0;

j=2;
tissue(j).name = 'water';
tissue(j).mua = 0.1;
tissue(j).mus = 10.0;
tissue(j).g   = 1.0;

j=3;
tissue(j).name = 'blood';
B       = 1.00;
S       = 0.75;
W       = 0.95;
M       = 0;
musp500 = 10;
fray    = 0.0;
bmie    = 1.0;
gg      = 0.90;
X = [B*S B*(1-S) W M]';
musp = musp500*(fray*(nm/500).^-4 + (1-fray)*(nm/500).^-bmie);
tissue(j).mua = MU*X;
tissue(j).mus = musp/(1-gg);
tissue(j).g   = gg;

j = 4;
tissue(j).name = 'dermis';
B = 0.002; 
S = 0.67;
W = 0.65;
M = 0;
musp500 = 42.4;
fray    = 0.62;
bmie    = 1.0;
X = [B*S B*(1-S) W M]';
musp = musp500*(fray*(nm/500).^-4 + (1-fray)*(nm/500).^-bmie);
gg = 0.90;
tissue(j).mua = MU*X;
tissue(j).mus = musp/(1-gg);
tissue(j).g   = gg;

j=5;
tissue(j).name = 'epidermis';
B = 0;
S = 0.75;
W = 0.75;
M = 0.03;
musp500 = 40;
fray    = 0.0;
bmie    = 1.0;
X = [B*S B*(1-S) W M]';
musp = musp500*(fray*(nm/500).^-4 + (1-fray)*(nm/500).^-bmie);
gg = 0.90;
tissue(j).mua = MU*X;
tissue(j).mus = musp/(1-gg);
tissue(j).g   = gg;

j=6;
tissue(j).name = 'skull';
B = 0.0005;
S = 0.75;
W = 0.35;
M = 0;
musp500 = 30;
fray    = 0.0;
bmie    = 1.0;
X = [B*S B*(1-S) W M]';
musp = musp500*(fray*(nm/500).^-4 + (1-fray)*(nm/500).^-bmie);
gg = 0.90;
tissue(j).mua = MU*X;
tissue(j).mus = musp/(1-gg);
tissue(j).g   = gg;

j=7;
tissue(j).name = 'gray matter';
B = 0.01;
S = 0.75;
W = 0.75;
M = 0;
musp500 = 20;
fray    = 0.2;
bmie    = 1.0;
X = [B*S B*(1-S) W M]';
musp = musp500*(fray*(nm/500).^-4 + (1-fray)*(nm/500).^-bmie);
gg = 0.90;
tissue(j).mua = MU*X;
tissue(j).mus = musp/(1-gg);
tissue(j).g   = gg;

j=8;
tissue(j).name = 'white matter';
B = 0.01;
S = 0.75;
W = 0.75;
M = 0;
musp500 = 20;
fray    = 0.2;
bmie    = 1.0;
X = [B*S B*(1-S) W M]';
musp = musp500*(fray*(nm/500).^-4 + (1-fray)*(nm/500).^-bmie);
gg = 0.90;
tissue(j).mua = MU*X;
tissue(j).mus = musp/(1-gg);
tissue(j).g   = gg;

Nt = length(tissue);

for j=1:Nt, tissue(j).nm = nm; end % wavelength

disp(sprintf('---- tissueList ------ \t%10s   \t%s  \t%s  \t%s','mua','mus','g','musp'))
for i=1:Nt
    musp = tissue(i).mus*(1-tissue(i).g);
    disp(sprintf('%d\t%15s\t%10.4f\t%0.1f\t%0.3f\t%0.1f',...
        i,tissue(i).name, tissue(i).mua,tissue(i).mus,tissue(i).g,musp))
end
disp(' ')

