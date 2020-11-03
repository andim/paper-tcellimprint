import os

libdir = os.path.dirname(os.path.realpath(__file__))
root_directory = os.path.abspath(libdir + '/' + os.path.pardir + '/' + os.path.pardir)
data_directory = root_directory + '/data/'
figure_directory = root_directory + '/figures/raw/'

chu_rawdata_directory = data_directory + 'raw/chu/'
emerson_rawdata_directory = data_directory + 'raw/emerson/'
lindau_rawdata_directory = data_directory + 'raw/lindau/'
britanova_rawdata_directory = data_directory + 'raw/britanova/'
emerson_processeddata_directory = data_directory + 'processed/emerson/'
britanova_processeddata_directory = data_directory + 'processed/britanova/'
