# 16S_from_fasta
Process 16S data from a fasta file.


### trying this https://docs.qiime2.org/2023.2/tutorials/otu-clustering/
### made new fasta in R with same formatting as example in above link 

(activate qiime)

qiime tools import \
  --input-path new_fasta.fa \
  --output-path seqs.qza \
  --type 'SampleData[Sequences]'


qiime vsearch dereplicate-sequences \
  --i-sequences seqs.qza \
  --o-dereplicated-table table.qza \
  --o-dereplicated-sequences rep-seqs.qza

#running in slum script denovo_cluster.sh

qiime vsearch cluster-features-de-novo \
  --i-table table.qza \
  --i-sequences rep-seqs.qza \
  --p-perc-identity 0.99 \
  --o-clustered-table table-dn-99.qza \
  --o-clustered-sequences rep-seqs-dn-99.qza


### visualize feature table

qiime feature-table summarize \
 --i-table table-dn-99.qza  \
 --o-visualization denovo_99_table.qzv \
 --m-sample-metadata-file metadata_v1.txt


### visualize representative sequences

qiime feature-table tabulate-seqs \
--i-data rep-seqs-dn-99.qza \
--o-visualization denovo-99-rep-seqs.qzv

### assign taxonomy

qiime feature-classifier classify-sklearn \
--i-classifier /home/Database/qiime2/classifiers/qiime2-2023.9/GTDBclassifier214.1_EMP.qza \
--i-reads rep-seqs-dn-99.qza \
--o-classification taxonomy_gtdb_214.qza

qiime metadata tabulate \
--m-input-file taxonomy_gtdb_214.qza \
--o-visualization taxonomy_gtdb_214.qzv
