Network-based Antiviral ranking

AntiViralPred.py runs the ranking algorithm as descirbed in
A systems-based approach to identify potential antivirals with a COVID-19 showcase

by Mengran Wang, Johanna Withers, Piero Ricchiuto, Michael McAnally, 
Helia Sanchez, Alif Saleh, Slava Akmaev and Dina Ghiassian 


Instruction to use the source code:

Download the code.
Make sure to make the code executable by chmod +x AntiViralPred.py
Make sure you replace all the input data and path to them according to the desired data
Data directory (data_dir) where all input data are stored must be set prior to running the code
Run the code (./AntiViralPred.py in command line or “run AntiViralPred.py" in ipython).
This will generate the final ranking matrix, given the compatible input data (See below)
-------------------

Directory "Data"

contains five input files:

Symbol2Synonyms.pcl (A dictionary generated from NCBI to convert gene symbol to all its synonyms)
Entrez2Symbol.pcl (A dictionary generated form NCBI to convert Entrez ID to gene symbol)
DrugBank.csv (Drug target information extracted from DrugBank)
bait_prey_high_confidence.xlsx (High confidence interactions between SARS-CoV-2 viral proteins and human proteins published by Gordon et. al)
HIV-1_physical_interactions.tsv (Physical interations among HIV proteins and human proteins extracted from NCBI)

The sixth input file is needed to run the code and it should be provided by the used depending on what network they want to use
Human_Interactome.txt (Protein-protein interaction network. note that gene IDs should be consistent in the two input files and in entrez ID)
