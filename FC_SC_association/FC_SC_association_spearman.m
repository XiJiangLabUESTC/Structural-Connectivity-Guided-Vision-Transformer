clear
clc
load('/media/D/fred/dhcp/birth_scan_below_real.mat');
load('/media/D/fred/dhcp/birth_below_scan_above_real.mat');
load('/media/D/fred/dhcp/birth_scan_above_real.mat');
for j=1:86
    load(['/media/D/fred/dhcp/fc/birth_scan_below/' birth_scan_below_real(j).name '_ses-' num2str(birth_scan_below_real(j).session) '_fc.mat']);
    load(['/media/D/fred/dhcp/sc/birth_scan_below/' birth_scan_below_real(j).name '_ses-' num2str(birth_scan_below_real(j).session) '_sc.mat']);
    sum_sc=sum(sc);   
    i=10;
    topk_sum=maxk(sum_sc,i);
    topk_roi=zeros(1,i);
    for k=1:i
        topk_roi(k)=find(sum_sc==topk_sum(k),1);
    end
    fc_masked=fc(:,topk_roi);
    sc_masked=sc(:,topk_roi);
    coa(j)=corr(fc_masked(:),sc_masked(:),'type','Spearman');
    coa_2(j)=corr(sum(fc_masked)',sum(sc_masked)','type','Spearman');
end

for j=1:80
    load(['/media/D/fred/dhcp/fc/birth_below_scan_above/' birth_below_scan_above_real(j).name '_ses-' num2str(birth_below_scan_above_real(j).session) '_fc.mat']);
    load(['/media/D/fred/dhcp/sc/birth_below_scan_above/' birth_below_scan_above_real(j).name '_ses-' num2str(birth_below_scan_above_real(j).session) '_sc.mat']);
    sum_sc=sum(sc);   
    i=10;
    topk_sum=maxk(sum_sc,i);
    topk_roi=zeros(1,i);
    for k=1:i
        topk_roi(k)=find(sum_sc==topk_sum(k),1);
    end
    fc_masked=fc(:,topk_roi);
    sc_masked=sc(:,topk_roi);
    cob(j)=corr(fc_masked(:),sc_masked(:),'type','Spearman');
    cob_2(j)=corr(sum(fc_masked)',sum(sc_masked)','type','Spearman');
end

for j=1:271
    load(['/media/D/fred/dhcp/fc/birth_scan_above/' birth_scan_above_real(j).name '_ses-' num2str(birth_scan_above_real(j).session) '_fc.mat']);
    load(['/media/D/fred/dhcp/sc/birth_scan_above/' birth_scan_above_real(j).name '_ses-' num2str(birth_scan_above_real(j).session) '_sc.mat']);
    sum_sc=sum(sc);   
    i=10;
    topk_sum=maxk(sum_sc,i);
    topk_roi=zeros(1,i);
    for k=1:i
        topk_roi(k)=find(sum_sc==topk_sum(k),1);
    end
    fc_masked=fc(:,topk_roi);
    sc_masked=sc(:,topk_roi);
    coc(j)=corr(fc_masked(:),sc_masked(:),'type','Spearman');
    coc_2(j)=corr(sum(fc_masked)',sum(sc_masked)','type','Spearman');
end

mean_p=[mean(coa),std(coa),mean(cob),std(cob),mean(coc),std(coc)];
mean_p_2=[mean(coa_2),std(coa_2),mean(cob_2),std(cob_2),mean(coc_2),std(coc_2)];

[~,p_connections(1,1)]=ttest2(coa,cob,'Alpha',0.05,'Vartype','unequal');
[~,p_connections(1,2)]=ttest2(cob,coc,'Alpha',0.05,'Vartype','unequal');
[~,p_connections(1,3)]=ttest2(coa,coc,'Alpha',0.05,'Vartype','unequal');
[~,p_rois(1,1)]=ttest2(coa_2,cob_2,'Alpha',0.05,'Vartype','unequal');
[~,p_rois(1,2)]=ttest2(cob_2,coc_2,'Alpha',0.05,'Vartype','unequal');
[~,p_rois(1,3)]=ttest2(coa_2,coc_2,'Alpha',0.05,'Vartype','unequal');

corr_connection=zeros(271,3);
corr_connection(1:86,1)=coa';
corr_connection(1:80,2)=cob';
corr_connection(1:271,3)=coc';
corr_roi=zeros(271,3);
corr_roi(1:86,1)=coa_2';
corr_roi(1:80,2)=cob_2';
corr_roi(1:271,3)=coc_2';
