# This script is used to download the alphafold2 data for model
# organisms and global health proteomes from the html page

# 0. Importing packages and functions----------------------------------------
import os
from beautifulsoup4 import BeautifulSoup
import subprocess
import multiprocessing as mp


# Defining a function to download a tar file
def download_tar_file(x):
    print(x)
    cmd = ["wget", x]
    subprocess.run(cmd)


# 1. Defining variables------------------------------------------------------

# Defining the number of workers for the multiprocessing
NUM_WORKERS = 20

# Defining the output path
output_path = "/scratch/alphafold_model_organisms/"

# 2. Downloading tar files---------------------------------------------------

# Extracting the urls from the html page using beutiful soup
with open("Index of _pub_databases_alphafold_latest.html") as fp:
    soup = BeautifulSoup(fp, "html.parser")

soup = soup.find_all("a")
soup_filtered = [x["href"] for x in soup if "UP" in x["href"]]

# Change directory to download the tar files in
os.chdir(output_path)

# Initialising the mprocess pool and number of workers
with mp.Pool(processes=NUM_WORKERS) as process_pool:
    process_pool.map(download_tar_file, soup_filtered)
