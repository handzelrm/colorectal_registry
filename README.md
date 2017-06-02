# colorectal_registry
Organize a colorectal registry database

Steps to running or updating data

1. Sign into redcap and export the most recent set of data. This may require a custom export in multiple chuncks
2. Export an updated database dictionary (Project Home -> Download the current Data Dictionary -> Codebook)
3. run create_data_dict function as it will update the pickle file that is used.
