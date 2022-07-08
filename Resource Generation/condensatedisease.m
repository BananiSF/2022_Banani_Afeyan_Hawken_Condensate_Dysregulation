% requires GO object and data variables with mendelian and cancer disease
% lists and with condensate associated GO terms in workspace

goids = unique(condgo.go);
decs = num2goid(getdescendants(go,cellfun(@goid2num,goids)));
annos = unique(filtertable(goaexp,'go',decs).upkb);

menddx = menddx(ismember(menddx.upkb,annos),:);
cancdx = cancdx(ismember(cancdx.upkb,annos),:);

clear goids decs annos;

tic;
dx = unique(menddx.dx);
cores = unique(condgo.core);

for i=1:length(cores)
	
	goids = filtertable(condgo,'core',cores{i}).go;
	decs = num2goid(getdescendants(go,cellfun(@goid2num,goids)));
	annos = unique(filtertable(goaexp,'go',decs).upkb);
	
	x1 = ismember(F.upkb,annos);
	totcond = sum(x1);
	
	for j=1:length(dx)
		
		x2 = ismember(F.upkb,filtertable(menddx,'dx',dx{j}).upkb);
		n = sum(x1&x2);
		totdx = sum(x2);
		
		x = crosstab(x1,x2);

		if numel(x)==2
			p = nan;
		else
			[~,p] = fishertest(x,'Tail','right');
		end
		
		mendstats = [mendstats;...
			table({cores{i}},{dx{j}},n,totcond,totdx,p,'VariableNames',{'cond','dx','n','totcond','totdx','p'})];
	end
end
toc;
clear goids decs annos x1 x2 x n totcond totdx;

tic;
dx = unique(cancdx.dx);
cores = unique(condgo.core);

for i=1:length(cores)
	
	goids = filtertable(condgo,'core',cores{i}).go;
	decs = num2goid(getdescendants(go,cellfun(@goid2num,goids)));
	annos = unique(filtertable(goaexp,'go',decs).upkb);
	
	x1 = ismember(F.upkb,annos);
	totcond = sum(x1);
	
	for j=1:length(dx)
		
		x2 = ismember(F.upkb,filtertable(cancdx,'dx',dx{j}).upkb);
		n = sum(x1&x2);
		totdx = sum(x2);
		
		x = crosstab(x1,x2);

		if numel(x)==2
			p = nan;
		else
			[~,p] = fishertest(x,'Tail','right');
		end
		
		cancstats = [cancstats;...
			table({cores{i}},{dx{j}},n,totcond,totdx,p,'VariableNames',{'cond','dx','n','totcond','totdx','p'})];
	end

end
toc;
mendstats.pct = mendstats.n./mendstats.totdx;
mendstats.fdr = bhfdr(mendstats.p);
mendstats.sig = mendstats.fdr<0.05;
cancstats.pct = cancstats.n./cancstats.totdx;
cancstats.fdr = bhfdr(cancstats.p);
cancstats.sig = cancstats.fdr<0.05;

clear goids decs annos x1 x2 x n totcond totdx i j dx cores p;