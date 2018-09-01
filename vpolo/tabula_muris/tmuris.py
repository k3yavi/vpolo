import os
import click
import sys
from urllib.request import urlretrieve
import subprocess

cellranger_base_path = "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3040nnn/GSM"
gsm_start = 3040890
gsm_end =  3040917

organs = ["Tongue","Tongue","Liver","Bladder","Bladder","Kidney", \
"Kidney","Spleen","Liver","Liver","Marrow","Marrow","Heart","Kidney", \
"Spleen","Bladder","Lung","Lung","Tongue","Thymus","Mammary","Mammary", \
"Muscle","Muscle","Lung","Lung","Trachea","Trachea"]

geo_ids = ["P4_0","P4_1","P4_2","P4_3","P4_4","P4_5","P4_6","P4_7", \
"P7_0","P7_1","P7_2","P7_3","P7_4","P7_5","P7_6","P7_7","P7_8","P7_9", \
"P7_10","P7_11","P7_12","P7_13","P7_14","P7_15","P8_12","P8_13","P8_14", \
"P8_15"]

bam_base_path = "https://sra-download.ncbi.nlm.nih.gov/traces/sra61/SRZ/006835/SRR"
bam_start = 6835844
bam_end = 6835871

def reporthook(blocknum, blocksize, totalsize):
# https://stackoverflow.com/questions/13881092/download-progressbar-for-python-3
    readsofar = blocknum * blocksize
    if totalsize > 0:
        percent = readsofar * 1e2 / totalsize
        s = "\r%5.1f%% %*d / %d" % (
            percent, len(str(totalsize)), readsofar, totalsize)
        sys.stderr.write(s)
        if readsofar >= totalsize: # near the end
            sys.stderr.write("\n")
    else: # total size is unknown
        sys.stderr.write("read %d\n" % (readsofar,))

@click.command()
@click.option('-t', '--tool', default="tenx", help='Currently only tenx.')
@click.option('--out', help='Path to the folder where to save the data; make sure you have write permission')
@click.option('-o', '--organ', help='Choose a type of organ',
    type=click.Choice(['Tongue', 'Liver', 'Bladder', 'Kidney', 'Spleen', 'Marrow', 'Heart', 'Lung', 'Trachea', 'Thymus', 'Mammary', 'Muscle']))
@click.option('-d', '--datatype', help='type of data to download, fastq assumes 10x tool bamtofastq is available',
    type=click.Choice(['bam', 'mtx', 'fastq']))
def download(tool, organ, datatype, out):
    if not os.path.exists(out):
        os.makedirs(out)

    if organ == None:
        print("Please provide a type of organ to download")
    elif datatype == None:
        print("Please provide a type of data to download")
    elif tool != "tenx":
        print("Unfortunately this tools works with tenx data only")
    else:
        organ_indices = [i for i, x in enumerate(organs) if x == organ]
        print("Found " + str(len(organ_indices)) + " data for " + organ)
        print("Downloading: ", [geo_ids[x] for x in organ_indices])
        if datatype in ["bam", "fastq"]:
            for index in organ_indices:
                file_name = os.path.join(out, geo_ids[index]+".bam")
                if os.path.isfile(file_name):
                    print("{} file already exist".format(file_name))
                else:
                    print("Downloading: "+geo_ids[index])
                    url = bam_base_path + str(bam_start+index) + "/10X_" + geo_ids[index] + ".bam"

                    urlretrieve(url, file_name, reporthook)
                    print("Done Downloading: "+geo_ids[index])
                if datatype == "fastq":
                    print("Making fastq: "+geo_ids[index])
                    fastq_path = os.path.join(out, "fastq")
                    stdoutdata = subprocess.getoutput("bamtofastq {} {}".format(file_name, fastq_path))
                    print("stdoutdata: {}".format(stdoutdata.split()))

        elif datatype == "mtx":
            for index in organ_indices:
                print("Downloading: "+geo_ids[index])
                url = cellranger_base_path + str(gsm_start+index) + "/suppl/GSM" + str(gsm_start+index) + "_" \
                    + organ + "-10X_" + geo_ids[index] + ".tar.gz"
                print(url)
                file_name = os.path.join(out, geo_ids[index]+".tar.gz")

                urlretrieve(url, file_name, reporthook)
                print("Done Downloading: "+geo_ids[index])

if __name__ == '__main__':
    '''
    Tool to download Tabule Muris data from SRA
    '''
    download()
