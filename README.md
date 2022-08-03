README

AirPollution_DistributedLag_ADRD

These files provide the code for the Alzheimer’s Disease and Related Dementia distributed lag data. By following Medicare beneficiaries for at least 10 years, starting in the year 2000, the yearly average air pollution exposure of the past decade is measured based on their place of residence and followed through 2016. Using about 8 million people that remained in Medicare in 2010, Medicare part A data is used to ascertain the first ADRD hospitalization for each beneficiary. 

Note that these files are very large, and it is necessary to run this code with a large memory allocation. The needs of each file vary and are indicated below, the largest of the requirements being 96 GB. Running the scripts with insufficient memory will result in errors. 

The files needed to run each script are specified below. Please ensure that the file path you use matches your file management as some paths are user specific and RCE or FASSE specific. Be aware that scripts may take multiple hours to run due to their computational intensity.


Summary of files in code folder:

1.	QID-year gridded exposure data at place of residence.R (96 GB to run)
Code: match exposures to FFS beneficiary zipcode of residence               
Inputs: denominator files, QD exposure data, meteorological data              
Outputs: exposures in grid by qid (rows) and year (cols)
      
2.	First hospitalization.R (32 GB to run)
Code: extract first hospitalization for AD/ADRD based on CCW algorithm      
Inputs: ADRD hospitalization data (2000-2016)                                 
Outputs: first AD/ADRD primary and secondary hospitalization     

3.	Denominator.R (32 GB to run)
Code: create individual denominator data for years 2009-2016     
Inputs: denominator data (2000-2016)                                 
Outputs: cleaned yearly qid enrollment and individual data   

4.	Year-zipcode confounders.R (32 GB to run)
Code: gather census and BRFSS data for every year/zipcode     
Inputs: denominator data (2000-2016), BRFSS data (exclude b/c no zip after 2012)                             
Outputs: year/zip confounders   

5.	Single state analysis data.R (32 GB to run)
Code: single state analysis data creation  
Inputs: denominator data (2000-2016)                                 
Outputs: cleaned yearly qid enrollment and individual data

6.	Region analysis data.R (64 GB)
Code: region analysis data creation  
Inputs: denominator data (2000-2016)                                 
Outputs: cleaned yearly qid enrollment and individual data       

7.	Cont US analysis data.R (64 GB to run)
Code: continental US analysis data creation  
Inputs: denominator data (2000-2016)                                 
Outputs: cleaned yearly qid enrollment and individual data
