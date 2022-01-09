

clear;
[allempty,all_up_down_mRNA]= xlsread('all_down_mRNA.xlsx');
[empty,down_miRNA_unique]= xlsread('down_miRNA_unique.xlsx');

num=length(down_miRNA_unique)
for n=1:num
    xls_name=['ENCORI_hg19_degradome-seq',down_miRNA_unique(n),'all.xlsx'];
    file=join(xls_name,'_')
    [empty1,miRNA_name]= xlsread(string(file));
    select_mRNA=miRNA_name(:,3);
    select_mRNA=unique(select_mRNA);
    count=0;
    for i=1:length(select_mRNA)
        [idd,liang]= find(ismember(all_up_down_mRNA,select_mRNA{i}));
        if isempty(liang)
             continue;
             
        else
            count=count+1;
            select_up_down_mRNA(n,count)=select_mRNA(i);
        end
    end
    
    
end
xlswrite('select_down_lnc_down_mRNA.xlsx',select_up_down_mRNA);
















