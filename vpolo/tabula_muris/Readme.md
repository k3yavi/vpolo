# Motication
quick and dirty way to download tabula muris data from [czBiohub/Geo](https://github.com/czbiohub/tabula-muris)

# Install
`pip install vpolo`

# Requirements:
Python 3  
click  

# Use
```
git clone https://github.com/k3yavi/vpolo.git
```

```
⇒  python tmuris.py --help  
Usage: tmuris.py [OPTIONS]  

Options:
  -t, --tool TEXT                 Currently only tenx.
  --out TEXT                      Path to the folder where to save the data;
                                  make sure you have write permission
  -o, --organ [Tongue|Liver|Bladder|Kidney|Spleen|Marrow|Heart|Lung|Trachea|Thymus|Mammary|Muscle]
                                  Choose a type of organ
  -d, --datatype [bam|mtx|fastq]  type of data to download, fastq is WIP
  --help                          Show this message and exit.
```
