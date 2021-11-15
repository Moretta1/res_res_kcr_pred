# res_res_kcr_pred 

## Sample code for research "Residue-residue contact can be a potential feature for the prediction of lysine crotonylation sites"

### subscrpt.py: 

#### contains some comment functions used in this research. For instance: reading files of different format(FASTA, tsv and SVM), find the whole peptide sequence with given ID, saving prediction result.

### rrc_generator.py: 

#### feature encoding scheme of Residue-Residue Composition (RRC). The output contains one .fasta file with peptides with all contacted residues expanded at the tail of original peptides (original peptides refer to the cutted sequence with length 2 * window-size + 1 ), and one .tsv format RRC feature file.

### rrpc_generator.py: 

#### feature encoding scheme of Residue-Residue Pair Composition (RRC). The output contains one .txt format RRPC feature file.

### classifiers.py: 

#### classifiers involved in this research. Three machine learning methods are shown here: SVM, random forest and logistic regression.

### classification.py: 

#### the main function in this demo. It takes feature file as input(for example, here RRC_10.txt means RRC features of window size ten), by taking different classifier, it shows performance matrix file and ROC curve figures as outputs.

### data: 

#### in this folder, some demo data were provided.

##### MP007031_results folder: contains the residue-residue contact prediction result obtained by MapPred tool. 

##### mapped_segments_length_10.txt: shows the original peptide sequence(here window size equals to 10 were used for instance, so the sequence length is 21 here), sample type(1 for positive and 0 for negative), residue_id, two indices from the whole-length peptide.

##### residues.fasta: whole-length peptide file, which consists of lable '>residue_id', and whole peptide sequence. 

##### id_mapping.csv: mapping between job_id and residue_id, which ensures that we can find the residue contact prediction job and its corresponding peptide.





