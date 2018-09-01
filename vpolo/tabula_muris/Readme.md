# Requirements:
Python 3 
click

⇒  python tmuris.py --help
Usage: tmuris.py [OPTIONS]

Options:
  -t, --tool TEXT                 Currently only tenx.
  --out TEXT                      Path to the folder where to save the data;
                                  make sure you have write permission
  -o, --organ [Tongue|Liver|Bladder|Kidney|Spleen|Marrow|Heart|Lung|Trachea|Thymus|Mammary|Muscle]
                                  Choose a type of organ
  -d, --datatype [bam|mtx|fastq]  type of data to download, if fastq make sure
                                  you have bamtofastq in path
  --help                          Show this message and exit.
