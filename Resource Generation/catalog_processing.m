% load protein information from uniprot
P = readtable('uniprot_hs_rev_nameseq.tab','FileType','text');
P.Properties.VariableNames{1} = 'uid';
P.Properties.VariableNames{2} = 'name';
P.Properties.VariableNames{3} = 'sequence';
P.name = strhead(P.name,'_');
P.length = cellfun(@length,P.sequence);

% load lists of condensate-forming proteins
llpsdbs = readtable('uidllpsdbs.txt');
psap = readtable('uidpsap.txt');
dp = readtable('uiddp.txt');

llps = unique(subvertmerge(subvertmerge(llpsdbs,dp(dp.dp>=.9,:)),psap(psap.psap>quantile(psap.psap,.9),:)));

P.llps = ismember(P.uid,llps.uid);

% load disorder scores from metapredict
load('211004_metapredict_scores.mat');
P = mergetabl(P,scores);

docutoff = 0.2;
P.dmask = cellfun(@(x)x>=docutoff,P.scores,'UniformOutput',false);
P.smask = cellfun(@(x)~x,P.dmask,'UniformOutput',false);

% define idr regions
R = subtab(P,{'uid','dmask','sequence'});
opfcn = @(x)struct2cell(regionprops(x,'PixelIdxList'))';

R = tabexpand(R,'dmask',opfcn,'regions');
R.start = cellfun(@min,R.regions);
R.stop = cellfun(@max,R.regions);

R.sequence = cellfun(@(x,y,z)x(y:z),R.sequence,num2cell(R.start),num2cell(R.stop),'uniformoutput',false);
R.length = cellfun(@length,R.sequence);
R.id = [1:height(R)]';

R = remtabvars(R,{'dmask','regions'},1);

% define non-idr regions
SR = subtab(P,{'uid','smask','sequence'});
opfcn = @(x)struct2cell(regionprops(x,'PixelIdxList'))';

SR = tabexpand(SR,'smask',opfcn,'regions');
SR.start = cellfun(@min,SR.regions);
SR.stop = cellfun(@max,SR.regions);

SR.sequence = cellfun(@(x,y,z)x(y:z),SR.sequence,num2cell(SR.start),num2cell(SR.stop),'uniformoutput',false);
SR.length = cellfun(@length,SR.sequence);
SR.id = [1:height(SR)]';

SR = remtabvars(SR,{'smask','regions'},1);

P.do = cellfun(@sum,P.dmask); % idr content
P.st = cellfun(@sum,P.smask);

% read in domain information from interpro and mid subset
D = readtable('uiddom.txt','Delimiter','\t');
mid = readtable('uidmid.txt','Delimiter','\t');
D = unique(D); mid = unique(mid);

mid = mergetabl(mid,countclass(mid.uid,'classname','uid','countsname','midval')); % mid valency
mid.di = strcat(mid.uid,':',mid.ipro);
mid = mergetabl(mid,countclass(mid.di,'classname','di','countsname','domval'));

D.di = strcat(D.uid,':',D.ipro); 
D = mergetabl(D,countclass(D.di,'classname','di','countsname','domval')); % all domain valency

% read in various lcs coordinates
patches = readtable('/Users/salmanbanani/Documents/Computational/Projects/patches/211011_patches_filtered.txt');
prd = readtable('uidprd.txt');
lark = readtable('uidlark.txt');
pscore = readtable('uidpscore.txt');

% filter out lcs < 5 amino acids
patches((patches.stop-patches.start+1)<5,:) = [];

patches.pi = strcat(patches.uid,':',patches.type);

patchmask = coords2mask(mergetabl(patches,subtab(P,{'uid','length'})),'rmask');
	patchmask = mergemasks(patchmask,'pi','rmask');
		patchmask.uid = strhead(patchmask.pi,':');
		patchmask.type = strtail(patchmask.pi,':');

prdmask = coords2mask(mergetabl(prd,subtab(P,{'uid','length'})),'rmask');
	prdmask = mergemasks(prdmask,'uid','rmask');
	prdmask.type = repmat({'prd'},[height(prdmask),1]);

larkmask = coords2mask(mergetabl(lark,subtab(P,{'uid','length'})),'rmask');
	larkmask = mergemasks(larkmask,'uid','rmask');
	larkmask.type = repmat({'lark'},[height(larkmask),1]);
	
pscoremask = coords2mask(mergetabl(pscore,subtab(P,{'uid','length'})),'rmask');
	pscoremask = mergemasks(pscoremask,'uid','rmask');
	pscoremask.type = repmat({'pscore'},[height(pscoremask),1]);
	
alllcsmasks = subvertmerge(subvertmerge(subvertmerge(patchmask,prdmask),larkmask),pscoremask);

lcsmask = mergemasks(alllcsmasks,'uid','rmask');
	t = table(P.uid(~ismember(P.uid,lcsmask.uid)),repmat({0},[sum(~ismember(P.uid,lcsmask.uid)),1]),'VariableNames',{'uid','rmask'});
	t.rmask = cellfun(@(x)logical(zeros([1,filtertable(P,'uid',x).length])),t.uid,'uniformoutput',false);
	lcsmask = [lcsmask;t]; clear t;
lcsmask.lcsval = cellfun(@sum,lcsmask.rmask);

midmask = subtab(mid,{'uid','start','stop'});
midmask = mergetabi(midmask,subtab(P,{'uid','length'}));
midmask = coords2mask(midmask,'rmask');
midmask = mergemasks(midmask,'uid','rmask');
	t = table(P.uid(~ismember(P.uid,midmask.uid)),repmat({0},[sum(~ismember(P.uid,midmask.uid)),1]),'VariableNames',{'uid','rmask'});
	t.rmask = cellfun(@(x)logical(zeros([1,filtertable(P,'uid',x).length])),t.uid,'uniformoutput',false);
	midmask = [midmask;t]; clear t;
	
lcs = subvertmerge(subvertmerge(subvertmerge(patches,prd),lark),pscore);
lcs.fxlen = lcs.stop-lcs.start+1;

% additional valency metrics
P = mergetabl(P,subtab(mid,{'uid','midval'},1));
P.midval(isnan(P.midval)) = 0;
P = mergetabl(P,subtab(lcsmask,{'uid','lcsval'}));
P = mergetabl(P,tvf2(mid,'uid','domval',@max));
P.domval_max(isnan(P.domval_max)) = 0;

% read in mutation information from vep
Tcv = readtable([inpath 'Tcv.txt'],'Delimiter','\t','MultipleDelimsAsOne',false); % clinvar
Tgn = readtable([inpath 'Tgn.txt'],'Delimiter','\t','MultipleDelimsAsOne',false); % genie
Tcb = readtable([inpath 'Tcb.txt'],'Delimiter','\t','MultipleDelimsAsOne',false); % cbioportal
	Tcb = Tcb(~Tcb.msk,:); %remove msk impact since data alrady in present in genie dataset
Thg = readtable([inpath 'Thg.txt'],'Delimiter','\t','MultipleDelimsAsOne',false); % hgmd

Tcv = tabvarlogical(Tcv,{'av','b','p','u','m','n','f','msk','psa','trunc'});
Tgn = tabvarlogical(Tgn,{'av','b','p','u','m','n','f','msk','psa','trunc'});
Tcb = tabvarlogical(Tcb,{'av','b','p','u','m','n','f','msk','psa','trunc'});
Thg = tabvarlogical(Thg,{'av','b','p','u','m','n','f','msk','psa','trunc'});

% mutation filtering for mutations with ambiguous or problematic properties

filt = @(T)T(~T.syn,:);
%^filter out synonymous variants
	Tcv = filt(Tcv);
	Tgn = filt(Tgn);
	Tcb = filt(Tcb);
	Thg = filt(Thg);

filt = @(T)T(cellstrfind(T.pmi,'p.'),:);
%^filter out cases where protein change not mapped
	Tcv = filt(Tcv);
	Tgn = filt(Tgn);
	Tcb = filt(Tcb);
	Thg = filt(Thg);

filt = @(T)T(~cellstrfind(T.pmi,{'Ter=','Ter?'}),:);
	%^filter to remove stop codon readthroughs and stop->stop cases since
	%these are not truncations per se
	Tcv = filt(Tcv);
	Tgn = filt(Tgn);
	Tcb = filt(Tcb);
	Thg = filt(Thg);

filt = @(T)T(~cellisempty(T.exon_rank),:);
% filt = @(T)T(~isempty(T.exon_rank),:);
	Tcv = filt(Tcv);
	Tgn = filt(Tgn);
	Tcb = filt(Tcb);
	Thg = filt(Thg);

% working tables
Tpsa = subvertmerge(subvertmerge(subvertmerge(Tcv(Tcv.psa,:),Tgn(Tgn.psa,:)),Tcb(Tcb.psa,:)),Thg(Thg.psa,:)); % protein sequence altering (missense and indels)

Ttrunc = subvertmerge(subvertmerge(subvertmerge(Tcv(Tcv.trunc,:),Tgn(Tgn.trunc,:)),Tcb(Tcb.trunc,:)),Thg(Thg.trunc,:)); % truncations (nonsense and frameshift)

Tpall = subvertmerge(subvertmerge(subvertmerge(Tcv(Tcv.p,:),Tgn(Tgn.p,:)),Tcb(Tcb.p,:)),Thg(Thg.p,:)); % all pathogenic mutations

% additional filter for mutations with ambiguous position
t = countclass(subtab(Tpsa,{'pmi','protein_position'},1).pmi,'classname','pmi');
t = t(t.counts>1,:);
Tpsa(ismember(Tpsa.pmi,t.pmi),:) = [];
clear t;

t = countclass(subtab(Ttrunc,{'pmi','protein_position'},1).pmi,'classname','pmi');
t = t(t.counts>1,:);
Ttrunc(ismember(Ttrunc.pmi,t.pmi),:) = [];
clear t;

t = countclass(subtab(Tpall,{'pmi','protein_position'},1).pmi,'classname','pmi');
t = t(t.counts>1,:);
Tpall(ismember(Tpall.pmi,t.pmi),:) = [];
clear t;

% additional filter for mutations with ambiguous/multiple leading consequences
t = subtab(Tpall,{'pmi','consequence'},1);
t = mergetabl(t,countclass(t.pmi,'classname','pmi'));
t = t(t.counts>1,:);
Tpall(ismember(Tpall.pmi,t.pmi),:) = [];
clear t;

% define disease-associated genes
dzg = subtab(Tpall,'uid',1);
P.dzg = ismember(P.uid,dzg.uid);

% define mendelian and cancer mutations
Tpall.mend = cellstrfind(Tpall.id,{'cv','hg'});
Tpall.canc = cellstrfind(Tpall.id,{'gn','cb'});

Tpsa.mend = cellstrfind(Tpsa.id,{'cv','hg'});
Tpsa.canc = cellstrfind(Tpsa.id,{'gn','cb'});

Ttrunc.mend = cellstrfind(Ttrunc.id,{'cv','hg'});
Ttrunc.canc = cellstrfind(Ttrunc.id,{'gn','cb'});

% remove version from ensemble transcript id
Tpall.enst = ensnodot(Tpall.enst);
Tpsa.enst = ensnodot(Tpsa.enst);
Ttrunc.enst = ensnodot(Ttrunc.enst);

% define mutations in condensate-forming proteins
Tpsa.llps = ismember(Tpsa.uid,llps.uid);
Ttrunc.llps = ismember(Ttrunc.uid,llps.uid);
Tpall.llps = ismember(Tpall.uid,llps.uid);

% read in transcript coordinates for nmd predictions
ens = readtable('ens.txt'); 
ens = ens(ismember(ens.enst,Ttrunc.enst),:);
ens = sortrows(ens,{'chr','enst_start','enst','exon_rank'},'ascend');

% identify last exon
ens = mergetabl(ens,...
		renametabcol(mergetabi(subtab(ens,{'enst','exon_rank'},1),...
			tvf(ens,'enst','exon_rank',@(x)max(x),'colname','exon_rank')),...
			'exon_rank','lastex'));

ens = sortrows(ens,{'chr','ensg_start','exon_rank'},'ascend');

% identify second-to-last exon
penultex = removevars(...
	mergetabi(subtab(ens,{'enst','exon_cdna_stop','exon_rank'},1),...
	tvf(ens,'enst','exon_rank',@(x)max(x)-1,'colname','exon_rank')),...
	'exon_rank');
penultex = renametabcol(penultex,'exon_cdna_stop','penultex_stop');
ens = mergetabl(ens,penultex);

% identify start codon position
enstaug = mergetabi(subtab(ens,{'enst','exon_cdna_coding_start'},1),...
	tvf(ens,'enst','exon_cdna_coding_start',@min,'colname','exon_cdna_coding_start'));
enstaug = renametabcol(enstaug,'exon_cdna_coding_start','enstaug_start');
ens = mergetabl(ens,enstaug);

nmdparams = subtab(ens,{'enst','exon_rank','ense_length','penultex_stop','enstaug_start','lastex'},1);

% nmd prediction table
nmd = Ttrunc;
nmd.enst = ensnodot(nmd.enst);

% stop codon position
nmd.stopcod = regexp(strtail(nmd.pmi,':'),'(?<!(fsTer|\d))\d++','match');
nmd.stopoffset = regexp(strtail(nmd.pmi,'p.'),'(?<=fsTer)\d++','match');
f = cellisempty(nmd.stopoffset);
nmd(f,:).stopoffset = repmat({{'0'}},[sum(f),1]);
clear f;

f = cellfun(@length,nmd.stopcod);
nmd.stopcod(f>1) = cellfun(@(x)join(x(:),','),nmd.stopcod(f>1),'UniformOutput',false);
nmd.stopcod = cellfun(@(x)x{:},nmd.stopcod,'UniformOutput',false);
clear f;

nmd.stopoffset = cellfun(@(x)x{:},nmd.stopoffset,'UniformOutput',false);

nmd.stopcod = cellstr2num(nmd.stopcod,'celloutput',true);
nmd.stopoffset = cellstr2num(nmd.stopoffset,'celloutput',true);
nmd.stopcod = cellfun(@mean,nmd.stopcod);
nmd.stopoffset = cellfun(@mean,nmd.stopoffset);

nmd.exon_rank = cellstr2num(nmd.exon_rank);
nmd = mergetabl(nmd,nmdparams);

nmd.stoppos = nmd.enstaug_start-1+(nmd.stopcod+nmd.stopoffset)*3-1;
toc;

% nmd prediction rules
nmd.nmd = (nmd.exon_rank~=nmd.lastex)&...
	nmd.stoppos<(nmd.penultex_stop-50)&...
	nmd.stoppos>(nmd.enstaug_start+200)&...
	nmd.ense_length<=400;

nmd = subtab(nmd,{'id','pmi',...
	'uid','ensg','enst','ensp','exon_rank','b','p','u',...
	'stopcod','stopoffset','ense_length','lastex','penultex_stop','enstaug_start','stoppos','nmd'},1);

Ttrunc = mergetabl(Ttrunc,subtab(nmd,{'enst','pmi','nmd'},1));
Tpall = mergetabl(Tpall,subtab(nmd,{'enst','pmi','nmd'},1));

% filter out mutations with ambiguous nmd predictions due to mutation
% mapping to multiple transcript isoforms
t = subtab(nmd,{'pmi','nmd'},1);
t = mergetabl(t,countclass(t.pmi,'classname','pmi'));
Tpall(ismember(Tpall.pmi,t(t.counts>1,:).pmi),:) = [];
clear t;

% add lcs-related metrics to truncation mutations
T = mergetabi(subtab(Ttrunc,{'uid','pmi','protein_position'},1),subtab(P,{'uid','smask'}));
T.protein_position = cellstr2num(T.protein_position,'celloutput',true);
T.protein_position = cellfun(@min,T.protein_position);
T.stlost = cellfun(@(x,y)sum(x(y:end)),T.smask,num2cell(T.protein_position));
Ttrunc = mergetabl(Ttrunc,subtab(T,{'pmi','stlost'},1));
clear T;

% loss of idr from truncation
T = mergetabi(subtab(Ttrunc,{'uid','pmi','protein_position'},1),subtab(P,{'uid','dmask'}));
T.protein_position = cellstr2num(T.protein_position,'celloutput',true);
T.protein_position = cellfun(@min,T.protein_position);
T.dolost = cellfun(@(x,y)sum(x(y:end)),T.dmask,num2cell(T.protein_position));
Ttrunc = mergetabl(Ttrunc,subtab(T,{'pmi','dolost'},1));
clear T;

% loss of lcs from truncation
T = mergetabi(subtab(Ttrunc,{'uid','pmi','protein_position'},1),subtab(lcsmask,{'uid','rmask'}));
T.protein_position = cellstr2num(T.protein_position,'celloutput',true);
T.protein_position = cellfun(@min,T.protein_position);
T.lcslost = cellfun(@(x,y)sum(x(y:end)),T.rmask,num2cell(T.protein_position));
Ttrunc = mergetabl(Ttrunc,subtab(T,{'pmi','lcslost'},1));
clear T;

% add mid-related metrics to truncation mutations
T = mergetabi(subtab(Ttrunc,{'uid','pmi','protein_position'},1),subtab(D,{'uid','ipro','start','stop'},1));
T.protein_position = cellstr2num(T.protein_position,'celloutput',true);
T.protein_position = cellfun(@min,T.protein_position);
T.domlost = T.protein_position<=T.stop;
T = tvf2(T,'pmi','domlost',@sum,'colname','domlost');
Ttrunc = mergetabl(Ttrunc,T);
Ttrunc.domlost(isnan(Ttrunc.domlost)) = 0;
clear T;

% loss of mid from mutation
T = mergetabi(subtab(Ttrunc,{'uid','pmi','protein_position'},1),subtab(mid,{'uid','ipro','start','stop'},1));
T.protein_position = cellstr2num(T.protein_position,'celloutput',true);
T.protein_position = cellfun(@min,T.protein_position);
T.midlost = T.protein_position<=T.stop;
T = tvf2(T,'pmi','midlost',@sum,'colname','midlost');
Ttrunc = mergetabl(Ttrunc,T);
Ttrunc.midlost(isnan(Ttrunc.midlost)) = 0;
clear T;

% add lcs-related metrics to missense/indel mutations
stcoords = mask2coords(P,'smask');
T = mergetabi(subtab(Tpsa,{'uid','pmi','protein_position'},1),subtab(stcoords,{'uid','start','stop'},1));
T.protein_position = cellstr2num(T.protein_position,'celloutput',true);
T.ppmin = cellfun(@min,T.protein_position);
T.ppmax = cellfun(@max,T.protein_position);
T.inst = cellfun(@(w,x,y,z)(w>=y&w<=z)|(x>=y&x<=z),num2cell(T.ppmin),num2cell(T.ppmax),num2cell(T.start),num2cell(T.stop));
T = T(T.inst,:);
Tpsa = mergetabl(Tpsa,subtab(T,{'pmi','inst'},1));

% affects on idr
docoords = mask2coords(P,'dmask');
T = mergetabi(subtab(Tpsa,{'uid','pmi','protein_position'},1),subtab(docoords,{'uid','start','stop'},1));
T.protein_position = cellstr2num(T.protein_position,'celloutput',true);
T.ppmin = cellfun(@min,T.protein_position);
T.ppmax = cellfun(@max,T.protein_position);
T.indo = cellfun(@(w,x,y,z)(w>=y&w<=z)|(x>=y&x<=z),num2cell(T.ppmin),num2cell(T.ppmax),num2cell(T.start),num2cell(T.stop));

T = T(T.indo,:);
Tpsa = mergetabl(Tpsa,subtab(T,{'pmi','indo'},1));
clear T docoords;

% affects on lcs
lcscoords = mask2coords(lcsmask,'rmask');
T = mergetabi(subtab(Tpsa,{'uid','pmi','protein_position'},1),subtab(lcscoords,{'uid','start','stop'},1));
T.protein_position = cellstr2num(T.protein_position,'celloutput',true);
T.ppmin = cellfun(@min,T.protein_position);
T.ppmax = cellfun(@max,T.protein_position);
T.inlcs = cellfun(@(w,x,y,z)(w>=y&w<=z)|(x>=y&x<=z),num2cell(T.ppmin),num2cell(T.ppmax),num2cell(T.start),num2cell(T.stop));
T = T(T.inlcs,:);
Tpsa = mergetabl(Tpsa,subtab(T,{'pmi','inlcs'},1));
clear T lcscoords;

% affects on domains
T = mergetabi(subtab(Tpsa,{'uid','pmi','protein_position'},1),subtab(D,{'uid','ipro','start','stop'},1));
T.protein_position = cellstr2num(T.protein_position,'celloutput',true);

T.ppmin = cellfun(@min,T.protein_position);
T.ppmax = cellfun(@max,T.protein_position);
T.indom = cellfun(@(w,x,y,z)(w>=y&w<=z)|(x>=y&x<=z),num2cell(T.ppmin),num2cell(T.ppmax),num2cell(T.start),num2cell(T.stop));
T = T(T.indom,:);
Tpsa = mergetabl(Tpsa,subtab(T,{'pmi','indom'},1));

T = mergetabi(subtab(Tpsa,{'uid','pmi','protein_position'},1),subtab(mid,{'uid','ipro','start','stop'},1));
T.protein_position = cellstr2num(T.protein_position,'celloutput',true);
T.ppmin = cellfun(@min,T.protein_position);
T.ppmax = cellfun(@max,T.protein_position);
T.inmid = cellfun(@(w,x,y,z)(w>=y&w<=z)|(x>=y&x<=z),num2cell(T.ppmin),num2cell(T.ppmax),num2cell(T.start),num2cell(T.stop));

T = T(T.inmid,:);
Tpsa = mergetabl(Tpsa,subtab(T,{'pmi','inmid'},1));

% add metrics to all pathogenic mutation working table
Tpall = mergetabl(Tpall,subtab(Tpsa,{'pmi','inlcs','inmid','indom','inst','indo'},1));
Tpall = mergetabl(Tpall,subtab(Ttrunc,{'pmi','lcslost','midlost','domlost','stlost','dolost'},1));

% compute additional metrics for valency
Ttrunc = mergetabl(Ttrunc,subtab(P,{'uid','lcsval','midval','llps'}));
Tpsa = mergetabl(Tpsa,subtab(P,{'uid','lcsval','midval','llps'}));
Tpall = mergetabl(Tpall,subtab(P,{'uid','lcsval','midval','llps'}));

% fractional valency loss for truncation mutations
Ttrunc.flcslost = Ttrunc.lcslost./Ttrunc.lcsval;
Ttrunc.fmidlost = Ttrunc.midlost./Ttrunc.midval;
Tpall.flcslost = Tpall.lcslost./Tpall.lcsval;
Tpall.fmidlost = Tpall.midlost./Tpall.midval;

% cutoff for fractional valency
fc = 0.25;

% define mutations that 'affect' condensate promoting features
Tpall.cpf = (Tpall.psa&(Tpall.inmid|Tpall.inlcs))|(Tpall.trunc&~Tpall.nmd&(Tpall.fmidlost>=fc|Tpall.flcslost>=fc));
Tpall.av = Tpall.psa|Tpall.trunc; % analyze missense/indel/truncation mutations

% define whether mutation affects mids, lcss, or both classes of
% condensate-promoting features
Tpall.cpfclass = repmat({''},[height(Tpall),1]);

f = (Tpall.psa&(Tpall.inmid&~Tpall.inlcs))|(Tpall.trunc&~Tpall.nmd&(Tpall.fmidlost>=fc&~(Tpall.flcslost>=fc)));
Tpall(f,:).cpfclass = repmat({'mid'},[sum(f),1]);

f = (Tpall.psa&(~Tpall.inmid&Tpall.inlcs))|(Tpall.trunc&~Tpall.nmd&(~(Tpall.fmidlost>=fc)&Tpall.flcslost>=fc));
Tpall(f,:).cpfclass = repmat({'lcs'},[sum(f),1]);

f = (Tpall.psa&(Tpall.inmid&Tpall.inlcs))|(Tpall.trunc&~Tpall.nmd&(Tpall.fmidlost>=fc&Tpall.flcslost>=fc));
Tpall(f,:).cpfclass = repmat({'both'},[sum(f),1]);