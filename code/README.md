This code is split into 3 sections: UD generation, overlap calculation, and mark-recapture model fitting. 

The UD generation goes through the steps of fitting the crawl, making the kde using just one individual (due to privacy constraints) and then the kde_data_all_old file has all the kdes for all individuals used in the paper such that users can average them and get the social group UDs. 

original code used in paper ran the crawl in parallel using custom functions available in "helper", scripts using these functions can be made available upon request. 


After unzipping "kde_data_all_old" you have to coerce the file extension to ".RData"
