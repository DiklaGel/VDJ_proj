% Read table
tableFile = 'Z:/diklag/output/result.csv';
f = fopen(tableFile);
header = fgetl(f);
header = regexp(header,',','split');
fields = containers.Map();
for i=1:length(header)
    fields(header{i}) = i;
end

Data = {};
i=1;
while ~feof(f)
    l = fgetl(f);
    l =  regexp(l,',','split');
    if (length(l{fields('CDR3 first')})>0)
        Data.CDR3{i} = l{fields('CDR3 first')};
        Data.well{i} = l{fields('well_id')};
        Data.Patient{i} = l{fields('Patient')};
        i=i+1;
    end
end
fclose(f);

% Calculate the pairwise distance between all CDR3s
CDR3_dist = seqpdist(Data.CDR3);
Z = linkage(CDR3_dist,'average');
[H,T,row_perm] = dendrogram(Z,0,'ColorThreshold','default');