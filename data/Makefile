
convert:
	xlsx2csv naive-percentage-britanova.xlsx > naive-percentage-britanova.csv
	xlsx2csv metadata-emerson.xlsx -s 1 > metadata-emerson.csv
	xlsx2csv metadata-emerson.xlsx -s 2 > metadata-emerson2.csv

decompress: 
	unzip vdjdb.zip
	unzip emerson-2017-natgen.zip -d raw/emerson/
	gzip raw/emerson/*.tsv
