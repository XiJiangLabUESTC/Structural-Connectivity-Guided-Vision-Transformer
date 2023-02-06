
clear
clc
load('PP_list.mat')
load('PT_list.mat')
load('TT_list.mat')
load('preterm_gene_atlas.mat');
load('term_gene_atlas.mat');

gyri_roi=[7,8,30,29,32,31,9,10,30,29,30,29,30,29,32,31,7,8,14,15,30,29];
sulci_roi=gyri_roi+32;
cortical_roi=[gyri_roi,sulci_roi];

for gene_id=1:120
    for i=1:length(pp)
        load(pp(i).name)  %load i-th subject's FC
        cortex=sum(fc(cortical_roi,:),2);
        cortex=(cortex-mean(cortex))/std(cortex);
        tmp_atlas=[preterm_gene_atlas(gene_id,:),preterm_gene_atlas(gene_id,:)]';
        tmp_atlas=(tmp_atlas-mean(tmp_atlas))./sqrt(tmp_atlas); %normalized
        pp_gene_roi_corr(gene_id,i)=corr(tmp_atlas,cortex,'type','Spearman');
    end
end

for gene_id=1:120
    for i=1:length(pt)
        load(pt(i).name)
        cortex=sum(fc(cortical_roi,:),2);
        cortex=(cortex-mean(cortex))/std(cortex);
        tmp_atlas=[term_gene_atlas(gene_id,:),term_gene_atlas(gene_id,:)]';
        tmp_atlas=(tmp_atlas-mean(tmp_atlas))./sqrt(tmp_atlas);
        pt_gene_roi_corr(gene_id,i)=corr(tmp_atlas,cortex,'type','Spearman');
    end
end

for gene_id=1:120
    for i=1:length(tt)
        load(tt(i).name)
        cortex=sum(fc(cortical_roi,:),2);
        cortex=(cortex-mean(cortex))/std(cortex);
        tmp_atlas=[term_gene_atlas(gene_id,:),term_gene_atlas(gene_id,:)]';
        tmp_atlas=(tmp_atlas-mean(tmp_atlas))./sqrt(tmp_atlas);
        tt_gene_roi_corr(gene_id,i)=corr(tmp_atlas,cortex,'type','Spearman');
    end
end

% gene_corr_sum_sum=sum([pp_gene_roi_corr,pt_gene_roi_corr,tt_gene_roi_corr],2);
% tmp=gene_corr_sum_sum;
% for i=1:10
%     max_gene(i)=find(abs(tmp)==max(abs(tmp)));
%     tmp(max_gene(i))=0;
% end
% sort(max_gene)

%%
%anova
label=[zeros(1,86),ones(1,80),ones(1,271)+1];
for i=1:120
    anova_p_all(i,1)=anova1([pp_gene_roi_corr(i,:),pt_gene_roi_corr(i,:),tt_gene_roi_corr(i,:)],label,'off');
end
close all



%%
gene_corr_mean=[mean(pp_gene_roi_corr,2),mean(pt_gene_roi_corr,2),mean(tt_gene_roi_corr,2)];
significant_p=find(anova_p_all*120<=0.05);
sign_gene_corr_mean=gene_corr_mean(significant_p,:);
for i=1:89
    [~,sorted_index(i,:)]=sort(sign_gene_corr_mean(i,:));
end
situation1=sorted_index(find(sorted_index(:,3)==1),:);
situation2=sorted_index(find(sorted_index(:,3)==2),:);
situation3=sorted_index(find(sorted_index(:,3)==3),:);
% ttest
for i=1:120
    [~,p(i,1)]=ttest2(pp_gene_roi_corr(i,:),pt_gene_roi_corr(i,:),'Alpha',0.05,'Vartype','unequal');
    [~,p(i,2)]=ttest2(pp_gene_roi_corr(i,:),tt_gene_roi_corr(i,:),'Alpha',0.05,'Vartype','unequal');
    [~,p(i,3)]=ttest2(pt_gene_roi_corr(i,:),tt_gene_roi_corr(i,:),'Alpha',0.05,'Vartype','unequal');
end

