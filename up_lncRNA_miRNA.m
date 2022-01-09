

clear;
[allempty,all_up_miRNA]= xlsread('all_up_mRNA.xlsx');
[empty,up_miRNA_unique]= xlsread('up_miRNA_unique.xlsx');

num=length(up_miRNA_unique)
for n=1:num
    xls_name=['ENCORI_hg19_degradome-seq',up_miRNA_unique(n),'all.xlsx'];
    file=join(xls_name,'_')
    [empty1,miRNA_name]= xlsread(string(file));
    select_mRNA=miRNA_name(:,3);
    select_mRNA=unique(select_mRNA);
    count=0;
    for i=1:length(select_mRNA)
        [idd,liang]= find(ismember(all_up_miRNA,select_mRNA{i}));
        if isempty(liang)
             continue;
             
        else
            count=count+1;
            select_up_mRNA(n,count)=select_mRNA(i);
        end
    end
    
    
end

xlswrite('select_up_lnc_up_mRNA.xlsx',select_up_mRNA);





































