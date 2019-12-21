# HumanGWAS_Query_within_MSSQL_Database

# Analysis_human_gwas_lipid_From_SQL_server
This python script is ran from the command line taking an input txt file, a base pair margin, and a pvalue cutoff. The margin and pvalue are optional arguments, so if not entered they will defer to default values. The script will take these arguments and locate all of the genes in the hg19 table. Once we get our information from hg19, it is cross referenced with all of the lipid gwas tables. The result is a pandas data frame, sorted by ascending pvalue, that is saved as a text file in an output folder. 

# INPUT ARGUMENTS 
In order to run this script naviagate to the directory containing script. From here, enter the following commands:
"py LIPID_GWAS_SCRIPT.py [-h] [-m] [-p] -i" Brackets indicate optional arguments and should not be included when running the script from the command line. Every arguments will come in a short and long version. Use one dash infront of short version and use a double dash for the long version. Arguments do not have to be in any specific order. Use "py" and "pip3" do use the updated version of the technologies.

 - ARGUMENT LEGEND:
   * -h, --help     use this flag, by itself, to bring up the help menu for the arguments from
                   the command line
   * -p, --pval     enter this flag, then your desired value of the pvalue cutoff   
   * -m, --margin   enter this flag, then your desired Base Pair margin value interval 
   * -i, --input    enter this flag, then the name of your .txt file to run on the script
                    
   * Additional Info
        * If [-p] not entered default value of .00005 will be used
        * If [-m] not entered default value 200000 will be used
        * input.txt should be in the form "geneName\n"
 
