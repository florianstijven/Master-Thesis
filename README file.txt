README file – Analysing the GASTRIC datasets

Reference:

Two datasets are made available alongside the following paper:
Buyse M, Molenberghs G, Paoletti X, Oba K, Alonso A, Van der Elst W, Burzykowski T (2015). 
Statistical evaluation of surrogate endpoints with examples from cancer trials. Biometrical Journal. 

Conditions of use:

The GASTRIC Group data can be used for research purposes, under the conditions that 
(1) the research be scientifically appropriate, 
(2) the confidentiality of individual patient data be protected, 
(3) the results of the analyses be shared with the GASTRIC Group prior to public communication, 
(4) the source of data be fully acknowledged as above, 
(5) resulting data and results be further shared with the research community.

Permission to publish:

Prior to publication of results based on the GASTRIC Group data, permission to publish must be sought by contacting the GASTRIC Secretariat: 
- xavier.paoletti@curie.fr
- oba@epistat.m.u-tokyo.ac.jp

Acknowledgments:

If GASTRIC Group data are used for a publication, the origin of the data must be acknowledged as follows: 
“The authors thank the GASTRIC (Global Advanced/Adjuvant Stomach Tumor Research International Collaboration) Group for permission to use their data. 
The investigators who contributed to GASTRIC are listed in references (1,2).” 
(1) The GASTRIC (Global Advanced/Adjuvant Stomach Tumor Research International Collaboration) Group (2010). 
Benefit of adjuvant chemotherapy for resectable gastric cancer: a meta-analysis. Journal of the American Medical Association 303:1729-37.
(2) The GASTRIC (Global Advanced/Adjuvant Stomach Tumor Research International Collaboration) Group (2013). 
Role of chemotherapy for advanced / recurrent gastric cancer: an individual-patient-data meta-analysis. European Journal of Cancer 49:1565-77.

Datasets:

- GASTRIC adjuvant data.csv: data on 3288 patients randomized in 14 clinical trials for patients with resectable gastric cancer 
- GASTRIC advanced data.csv: data on 4069 patients randomized in 20 clinical trials for patients with advanced gastric cancer 

Data structure of GASTRIC adjuvant data.csv:
- OS, overall survival in days
- last_status_0Cens_1Death, survival status (0=alive,1=dead)
- dfs, disease-free survival in days
- dfs_status_0cens_1relapsedeath, disease-free survival status (0=alive without relapse, 1=relapsed or dead)
- arm_txt, treatment arm (Surgery alone=control, Adj CT=experimental)
- study, study identifier
- id-pat, patient identifier

Data structure of GASTRIC advanced data.csv:
- OS, overall survival in days
- last_status_0cens_1death, survival status (0=alive,1=dead)
- PFS_status_0cens_1ProgDeath, progression-free survival status (0=alive without progression, 1=progressed or dead)
- pfs, progression-free survival in days
- id-study, study identifier
- arm, treatment arm (1=control, 2=experimental)
- id-pat, patient identifier

Software:

SAS Software to perform surrogacy analyses is freely available from the following website: www.ibiostat.be/software. 
In particular, the macro for the Plackett copula, "plactad8withtauandrho.sas", have been used to perform the analyses of surrogacy. 
The macro is contained in the .zip file available at http://ibiostat.be/software/uploads/surr1.zip. 
The macro contains a preamble providing instructions for its use.

