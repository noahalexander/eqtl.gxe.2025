#crude data organization for the time being 

most segregant single-cell data is in /u/home/n/nalexand/project-kruglyak/sc.data.2303
byxrm salt: /u/home/n/nalexand/project-kruglyak/sc.data.2303/salt_20230424/salt_tc_20230424
byxrm q: /u/project/kruglyak/nalexand/sc.data.2303/segregant/FASTQ_Generation_2023-08-30_09_21_20Z-690101418/
cbsxyjm salt: /u/home/n/nalexand/project-kruglyak/cellranger.testing/3004.salt.seg/

#####processing from fastq

#in t0_seg dir in cellranger.testing ...
qrsh -l highp,h_rt=23:00:00,h_data=8G -pe shared 8

#login4:
cellranger count --id=sample-I --fastqs /u/project/kruglyak/nalexand/sc.data.2303/segregant/FASTQ_Generation_2023-08-30_09_21_20Z-690101418/sample_I --sample=sample-I --transcriptome=/u/project/kruglyak/smilefre/sceqtl/2021/ref/yeast

#login3:
#in t10_seg dir in cellranger.testing ...
cellranger count --id=sample-II --fastqs /u/project/kruglyak/nalexand/sc.data.2303/segregant/FASTQ_Generation_2023-08-30_09_21_20Z-690101418/sample_II --sample=sample-II --transcriptome=/u/project/kruglyak/smilefre/sceqtl/2021/ref/yeast

# data in /u/home/n/nalexand/project-kruglyak/cellranger.testing/q.seg/t0.r3/sample-I/outs and corresponding dir for t10




#tmux, qrsh, tmux again:

#module load bcftools
#python /u/home/n/nalexand/project-kruglyak/vcf_parents.t30/extract_parents_hoffman.navartrixpath.20230822.py.cpp --vcf /u/home/n/nalexand/parents.nostar.vcf --crosses 3004 --threads 8 --cellranger-outdir /u/home/n/nalexand/project-kruglyak/cellranger.testing/q.seg/t0.r3/sample-I/ --parent-list /u/project/kruglyak/smilefre/sceqtl/2021/ref/parent_list.txt 
#login4:
#tmux, qrsh, tmux again:
#module load bcftools
#python /u/home/n/nalexand/project-kruglyak/vcf_parents.t30/extract_parents_hoffman.navartrixpath.20230822.py.cpp --vcf /u/home/n/nalexand/parents.nostar.vcf --crosses 3004 --threads 8 --cellranger-outdir /u/home/n/nalexand/project-kruglyak/cellranger.testing/q.seg/t10/sample-II/ --parent-list /u/project/kruglyak/smilefre/sceqtl/2021/ref/parent_list.txt 


####----------20230905

3051:

#login2:
#tmux, qrsh, tmux again:

module load bcftools
python /u/home/n/nalexand/project-kruglyak/vcf_parents.t30/extract_parents_hoffman.navartrixpath.20230822.py.cpp --vcf /u/home/n/nalexand/parents.nostar.vcf --crosses A --threads 8 --cellranger-outdir /u/home/n/nalexand/project-kruglyak/cellranger.testing/q.seg/t0.r3/sample-I --parent-list /u/project/kruglyak/smilefre/sceqtl/2021/ref/parent_list.txt 

#login3:
#tmux, qrsh, tmux again:
module load bcftools
python /u/home/n/nalexand/project-kruglyak/vcf_parents.t30/extract_parents_hoffman.navartrixpath.20230822.py.cpp --vcf /u/home/n/nalexand/parents.nostar.vcf --crosses A --threads 8 --cellranger-outdir /u/home/n/nalexand/project-kruglyak/cellranger.testing/q.seg/t10/sample-II --parent-list /u/project/kruglyak/smilefre/sceqtl/2021/ref/parent_list.txt 


#################################################################################### CBSxYJM salt time course processing 
qrsh -l highp,h_rt=23:00:00,h_data=8G -pe shared 8

login4:
cellranger count --id=sample-I --fastqs /u/home/n/nalexand/project-kruglyak/sc.data.2303/salt.3004/sample_I --sample=sample-I --transcriptome=/u/project/kruglyak/smilefre/sceqtl/2021/ref/yeast

cellranger count --id=sample-II --fastqs /u/home/n/nalexand/project-kruglyak/sc.data.2303/salt.3004/sample_II --sample=sample-II --transcriptome=/u/project/kruglyak/smilefre/sceqtl/2021/ref/yeast


/u/home/n/nalexand/project-kruglyak/cellranger.testing/3004.salt.seg/t0/sample-I/outs/filtered_feature_bc_matrix



login1
module load bcftools
python /u/home/n/nalexand/project-kruglyak/vcf_parents.t30/extract_parents_hoffman.navartrixpath.20230822.py.cpp --vcf /u/home/n/nalexand/parents.nostar.vcf --crosses 3004 --threads 8 --cellranger-outdir /u/home/n/nalexand/project-kruglyak/cellranger.testing/3004.salt.seg/t0/sample-I/ --parent-list /u/project/kruglyak/smilefre/sceqtl/2021/ref/parent_list.txt 


login2
module load bcftools
python /u/home/n/nalexand/project-kruglyak/vcf_parents.t30/extract_parents_hoffman.navartrixpath.20230822.py.cpp --vcf /u/home/n/nalexand/parents.nostar.vcf --crosses 3004 --threads 8 --cellranger-outdir /u/home/n/nalexand/project-kruglyak/cellranger.testing/3004.salt.seg/t30/sample-II/ --parent-list /u/project/kruglyak/smilefre/sceqtl/2021/ref/parent_list.txt
