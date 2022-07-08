% requires GO object and data variables with Mendelian and cancer disease
% lists and with GO annotations in workspace

allgoterms = unique(goaexp.go);

for i=1:length(allgoterms)
	warning off;
	goid = allgoterms{i};
	allgostats.term{i} = allgoterms{i};
	decs = num2goid(getdescendants(go,goid2num(goid)));
	annos = unique(filtertable(goaexp,'go',decs).upkb);
	x1 = F.ismv&ismember(F.upkb,dzg.upkb);
	x2 = ismember(F.upkb,annos);
	x = crosstab(x1,x2);
	
	if numel(x)==2
		x = [x [0;0]];
		allgostats.pval(i) = nan;
	else
		allgostats.n(i) = x(4);
		allgostats.tot(i) = sum(x(3:4));
		[~,allgostats.pval(i)] = fishertest(x,'Tail','right');
	end

	progressf(i,length(allgoterms),10);
	
end
warning on;
toc;
allgostats.fdr = bhfdr(allgostats.pval);
allgostats.sig = allgostats.fdr<0.05;