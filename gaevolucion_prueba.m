function gaDat=gaevolucion_prueba(gaDat,nint)
% One generation -------
Chrom=gaDat.Chrom;
nind=size(Chrom,1);
ObjV=inf(nind,1);
%add edwin
parametros=gaDat.ObjfunPar;
for i=1:nind
    if isempty(gaDat.ObjfunPar)
        ObjV(i)=feval(gaDat.Objfun,Chrom(i,:));
        
    else
        parametros.filename=['mc_genet_'  num2str(i)];
        %filenm=parametros.filename;
         copyfile('*.m',['arquivos_his/generacion' num2str(nint) ]);
           cd(['arquivos_his/generacion' num2str(nint) ])  
          
        ObjV(i)=feval(gaDat.Objfun,Chrom(i,:),parametros);
        % movefile([filenm '.his'],['arquivos_his/generacion' num2str(nint) ])
%          movefile([filenm '_F.bin'],['arquivos_F/generacion' num2str(nint) ])
%          movefile([filenm '_T.bin'],['arquivos_T/generacion' num2str(nint) ])
%          movefile([filenm '_H.mci'],['arquivos_mci/generacion' num2str(nint) ])
        clear parametros.filename
        clear filenm;
       
    end
end
gaDat.ObjV=ObjV;

% Best individual of the generation -------------------------
[v,p]=min(gaDat.ObjV);
if v<=gaDat.fxmin
    gaDat.xmin=Chrom(p,:);
    gaDat.fxmin=v;
end
% Next generation
% RANKING -------------------------------------------------
FitnV = ranking(gaDat.ObjV,gaDat.rf);
% SELECTION -----------------------------------------------
% Stochastic Universal Sampling (SUS).
SelCh = select('sus',Chrom,FitnV,1);
% CROSSOVER ---------------------------------------------------
% Uniform crossover.
SelCh = lxov(SelCh,gaDat.Pc,gaDat.alfa);
% MUTATION ------------------------------------------------
Chrom = mutbga(SelCh,gaDat.FieldD,[gaDat.Pm 1]); % Codificación Real.
% Reinsert the best individual  ---------------------------
Chrom(round(gaDat.NIND/2),:) = gaDat.xmin;
gaDat.Chrom=Chrom;
% Optional additional task required by user
gaiteration(gaDat)

%% ---------------------------------------------------------
function FitV=ranking(ObjV,RFun)
% Ranking function
if nargin==1
    error('Ranking function needs two parameters');
end

if ~(length(ObjV)==length(RFun))
    error('RFun have to be of the same size than ObjV.');
end

[~,pos]=sort(ObjV);
FitV(pos)=flipud(RFun);
FitV=FitV';

%% ---------------------------------------------------------
function [SelCh]=select(SEL_F, Chrom, FitnV, GGAP)
% Selection Function
if (nargin==3) %  No overlap -------------------
    if strcmp(SEL_F,'rws')
        % Roulette wheel selection method
        indices=rws(FitnV,length(FitnV));
        SelCh=Chrom(indices,:);
    elseif strcmp(SEL_F,'sus')
        % Stochastic unversal sampling selection
        indices=sus(FitnV,length(FitnV));
        SelCh=Chrom(indices,:);
    else
        error('Incorrect selection method');
    end
elseif (nargin==4) % With overlap -----------------------------
	% Indexes of new individuals
    if strcmp(SEL_F,'rws')
        indices=rws(FitnV,round(length(FitnV)*GGAP));
    elseif strcmp(SEL_F,'sus')
        indices=sus2(FitnV,round(length(FitnV)*GGAP));
    else
        error('Incorrect selection method');
    end

    if (GGAP<1) % there is overlap
        % Members of the population to overlap
        oldpos=(1:length(FitnV))';
        for k=1:length(FitnV)
            pos=round(rand*length(FitnV)+0.5);
            % exchange indexes
            oldpos([pos k])=oldpos([k pos]);
        end
        oldpos=oldpos(1:round(length(FitnV)*GGAP));
        SelCh=Chrom;
        SelCh(oldpos,:)=Chrom(indices,:);
    else % more childs than parents
        SelCh=Chrom(indices,:);
    end
else
    error('Incorrect number of paramenters');
end

% Disorder the population.
[~,indi]=sort(rand(length(FitnV),1));
SelCh=SelCh(indi,:);

%% ------------------------------------------------------------------
function NewChrom =lxov(OldChrom, XOVR, alpha)
% Linear crossover
% Produce a~ new population by linear crossover and XOVR crossover probability
%   NewChroms =lxov(OldChrom, XOVR, alpha, FieldDR)
%
% Linear recombination.
% Parameters 'beta1' and 'beta2' are randomly obtained inside [-alpha, 1+alpha]
% interval
%   Child1 = beta1*Parent1+(1-beta1)*Parent2
%   Child2 = beta2*Parent1+(1-beta2)*Parent2

if nargin==1
    XOVR = 0.7;
    alpha = 0;
elseif nargin==2
    alpha = 0;
end

n = size(OldChrom,1);   % Number of individuals and chromosome length
npares = floor(n/2);    % Number of pairs
cruzar = rand(npares,1)<= XOVR;    % Pairs to crossover
NewChrom=OldChrom;

for i=1:npares
    pin = (i-1)*2+1;
    if ~(cruzar(i)==0)
        betas=rand(2,1)*(1+2*alpha)-(0.5+alpha);
        A=[betas(1) 1-betas(1); 1-betas(2) betas(2)];
        NewChrom(pin:pin+1,:)=A*OldChrom(pin:pin+1,:);
    end
end

% Coerce points outside search space
% aux = ones(n,1);
% auxf1=aux*FieldDR(1,:);
% auxf2=aux*FieldDR(2,:);
% NewChrom = (NewChrom>auxf2).*auxf2+(NewChrom<auxf1).*auxf1+(NewChrom<=auxf2 & NewChrom>=auxf1).*NewChrom;

%% -------------------------------------------------------------
function NewChrom=mutbga(OldChrom,FieldDR,MutOpt)
% Mutation function
% Real coded mutation. 
% Mutation is produced adding a low random value
% OldChrom: Initial population.
% FieldChrom: Upper and lower bounds.
% MutOpt: mutation options,
%         MutOpt(1)=mutation probability (0 to 1).
%         MutOpt(2)=compression of the mutation value (0 to 1).
%         default MutOpt(1)=1/Nvar y MutOpt(2)=1

if (nargin==3)
    pm=MutOpt(1);
    shr=MutOpt(2);
elseif (nargin==2)
    pm=1/size(FieldDR,2);
    shr=1;
else
    error('Incorrect number of parameters');
end

Nind=size(OldChrom,1);
m1=0.5-(1-pm)*0.5;
m2=0.5+(1-pm)*0.5;
aux=rand(size(OldChrom));
MutMx=(aux>m2)-(aux<m1);
range=[-1 1]*FieldDR*0.5*shr;
range=ones(Nind,1)*range;
index=find(MutMx);
m=20;
alpha=rand(m,length(index))<(1/m);
xx=2.^(0:-1:(1-m));
aux2=xx*alpha;
delta=zeros(size(MutMx));
delta(index)=aux2;
NewChrom=OldChrom+(MutMx.*range.*delta);

% Coerce points outside bounds
aux = ones(Nind,1);
auxf1=aux*FieldDR(1,:);
auxf2=aux*FieldDR(2,:);
NewChrom = (NewChrom>auxf2).*auxf2+(NewChrom<auxf1).*auxf1+(NewChrom<=auxf2 & NewChrom>=auxf1).*NewChrom;

%% ----------------------------------------------------------
function NewChrIx=sus2(FitnV, Nsel)
suma=sum(FitnV);     
% Position of the roulette pointers
j=0;
sumfit=0; 
paso=suma/Nsel; % distance between pointers
flecha=rand*paso; % offset of the first pointer
NewChrIx(Nsel,1)=0; 
for i=1:Nsel
    sumfit=sumfit+FitnV(i);
    while (sumfit>=flecha)
        j=j+1;
        NewChrIx(j)=i;
        flecha=flecha+paso;
    end
end