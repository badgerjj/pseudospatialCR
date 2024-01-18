This code is split into 3 sections: UD generation, overlap calculation, and mark-recapture model fitting. 

The UD generation goes through the steps of fitting the movement model in crawl and forming the kernel density estimate using just one tagged individual (due to data privacy constraints). The kde_data_all_old file has the kernel density estimates for all individuals used in the manuscript, such that users can average them and obtain the social group UDs. Importantly, after unzipping "kde_data_all_old", users must manually coerce the file extension to ".RData". Original code (used in the paper) ran crawl in parallel using custom functions available in helper.R, and scripts using these functions can be made available on request.

