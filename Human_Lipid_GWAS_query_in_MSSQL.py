from sqlalchemy import create_engine
import pandas as pd
import urllib
import pyodbc
import os
import pyodbc
import csv
import os
from functools import reduce
import argparse
from functools import reduce
import argparse


def error_messeage():
    print("""
          ***********************************
          *****Please enter valid input!*****
          ***********************************
          """)


def get_table():
# This function is to get a 'list' containing indices of tables to be included from user

    while True:
        # It prompts to get input asking which table to include
        choose_table = input("""
Which tables do you want to include?\n
1. Lipid_EHR_Hoffmann_NG
2. Lipid_GLGC_Teslovich_N
3. Lipid_Korean_Exome_Moon_SR
4. Lipid_Engage_Surakka_NG
5. Lipid_GLGC_Exome_East_Asian_Lu_NG
6. Lipid_GLGC_Exome_European_and_East_Asian_Lu_NG
7. Lipid_GLGC_Willer_NG
8. Lipid_Japanese_lipid_trait_Kanai_NG
9. Lipid_MVP_Klarin_NG
10.Lipid_EastAsian_Spracklen_Hum_Mol_Genetics
11.Lipid_UKBB_high_cholesterol_ukbb_Connor_alkesgroup
12.Lipid_UKBB_lipid_trait_Neale
13.Lipid_UKBB_statin_usage_Neale
14.Lipid_GLGC_Exome_Liu_NG
15.Lipid_UK10K_prins_SR
16.Lipid_Bivariate_GLGC_CAD_Siewert_CGPM
17.Lipid_Finnish_Exome_Locke
18.Lipid_Chinese_exome_Lu_HMG

Press 'ENTER' to include all tables.
Otherwise, enter the table numbers that you want to include seperataed by comma. e.g.) 1,3,5
""")

        #Obtain index of chosen tables as a list
        try:
            #select all
            if choose_table == '':
                chosen_table = [i for i in range(1,19)]
                break

            #select multiple tables
            elif ',' in choose_table:
                # The sets module provides classes for constructing and
                # manipulating unordered collections of unique elements.
                # Then cast to list again.
                chosen_table = list(set([int(i) for i in choose_table.split(",")]))

                if False not in [i in range(1,19) for i in chosen_table]:
                    break
                else:
                    error_messeage()
                    continue

            #select a single table
            elif int(choose_table) in [i for i in range(1,19)]:
                chosen_table = [int(choose_table)]
                break

            else:
                error_messeage()
                continue

        except ValueError:
            error_messeage()
            continue

    return(chosen_table)


def get_chr():
# This function is to get a 'list' of chromosome to be included from user

    while True:
        choose_chr = input("""
Press "ENTER' to include all chromosomes.
Otherwise, enter chromosomes you want to include separate by comma.:\n""")

        # Valid input list for chromosome
        chr_list   = list(range(1,24))

        #Obtain index of chosen tables as a list
        try:
            #select all
            if choose_chr == '':
                chosen_chr = chr_list
                break

            #select multiple tables
            elif ',' in choose_chr:
                # The sets module provides classes for constructing and
                # manipulating unordered collections of unique elements.
                # Then cast to list again.
                chosen_chr = list(set([int(i) for i in choose_chr.split(",")]))

                if False not in [i in chr_list for i in chosen_chr]:
                    break
                else:
                    error_messeage()
                    continue

            #select a single table
            elif int(choose_chr) in [i for i in chr_list]:
                chosen_chr = [int(choose_chr)]
                break

            else:
                error_messeage()
                continue

        except ValueError:
            error_messeage()
            continue


    return([str(i) for i in chosen_chr])


def get_margin():

    while True:
        try:
            margin_input = input("""
Press just 'ENTER' to use default value of 200,000.
Otherwise, specify your margin:""")

            if margin_input == '':
                return(200000)

            else:
                margin = int(margin_input)
                break

        except ValueError:
            print('Enter a valid input.\n')

    return(margin)


def get_genename():
    gene_name    = input("Type gene name: , or press ENTER to import a txt file : ").lower()

    if gene_name == "":
        filepath = input("input your file path: ")
        file = open(filepath, 'r')
        gene_name_list = file.read()
        genes = gene_name_list.split(",")
        print ("input genes: ")
        print (genes)
        file.close()
        return genes

    else:
        if ',' in gene_name:
            print('yes')
            return(list(set([i.replace(' ','') for i in gene_name.split(",")])))

        else:
            return([gene_name])


def get_pvalue():

    cutoff = input("""
Press just 'ENTER' to use the default value of 0.05.
Otherwise, specify the cutoff for p-value:\n""")

    if cutoff == '':
        return(0.05)


    while True:

        try:
            cutoff = float(cutoff)
            while (cutoff < 0) or (cutoff > 1):
                print('Enter a valid input.\n')
                cutoff = float(input("""
Press just 'ENTER' to use the default value of 0.05.\n
Otherwise, specify the cutoff for p-value:\n"""))

        except ValueError:
            print('Enter a valid input.\n')
            cutoff = input("""
Press just 'ENTER' to use the default value of 0.05.\n
Otherwise, specify the cutoff for p-value:\n""")

        return(cutoff)


def where(hg19,cutoff):
    chr_list = list(set(hg19['chr'].tolist()))
    where = "("
    for chr in chr_list:
        where += "chr = \'%s\'" %(chr)

        if chr_list.index(chr) != len(chr_list)-1:
            where += " or "

    where += ") and ("
    for i in range(len(hg19)):
        where += "(bp between %d and %d)" %(hg19['adj_chr_start'][i], hg19['adj_chr_end'][i])

        if i != len(hg19)-1:
            where += ' or '

    where += ") and (p_value < %f)" %(cutoff)

    return(where)


def where_varchar(hg19, cutoff):
# This function is exactly same as 'def where'
# But for the table whose p_value column is stored as varchar
    chr_list = list(set(hg19['chr'].tolist()))
    where = "("
    for chr in chr_list:
        where += "chr = \'%s\' " %(chr)

        if chr_list.index(chr) != len(chr_list)-1:
            where += "or "

    where += ") and ("
    for i in range(len(hg19)):
        where += "(bp between %d and %d)" %(hg19['adj_chr_start'][i], hg19['adj_chr_end'][i])

        if i != len(hg19)-1:
            where += ' or '

    where += ") and (p_value < \'%f\')" %(cutoff)

    return(where)


def get_hg19(sql_conn, gene_name, margin):
#      This will spits out 'hg19' table with chosen chromosome and gene_name

#     for chr in chosen_chr:
#         where += "chr = \'%s\' " %(chr)

#         if chosen_chr.index(chr) != len(chosen_chr)-1:
#             where += "or "

    # where += ") and ("

    # concatenating a condition regarding 'gene_name'
    where = "("
#     for gene_i in gene_name:
    where      += "gene_name = \'%s\'" %(gene_name)

#         if gene_name.index(gene_i) != len(gene_name)-1:
#             where += "or "

    where += ")"

    # concatenating 'where' statement and modify the interval of 'chr_start' and 'chr_end'
    query_hg19 = """
                 SELECT chr, gene_name, chr_start - %d as 'adj_chr_start',
                        chr_end + %d as 'adj_chr_end', 'hg19' as table_name
                 FROM hg19
                 WHERE """ %(margin, margin) + where


    hg19 = pd.read_sql(query_hg19, sql_conn)
    hg19 = hg19.astype({"chr": "category","adj_chr_start": "int64", "adj_chr_end": "int64", "table_name": "category"})

    return(hg19)

def load_to_sql():
    choice = input("what name you want to save this file in sql server ?\n")
    tablename = choice
    server = 'PARKSLAB'
    create = False
    path = "E:\\cross_ref\\load"
    database = "Human_GWAS"

    print (database)
    # connect to database

    cnxn = pyodbc.connect('DRIVER={SQL Server}' + \
                        ';SERVER=' + server + \
                        ';DATABASE=' + database + \
                        ';Trusted_Connection= Yes')

    print('connected to the database: %s successfully!' % database)
    cursor = cnxn.cursor()

    # if create is specified in the command line
    # if create:
    query = "create table {!s} ".format(tablename) + "(" + \
            "beta varchar (100)," \
            "bp varchar (100)," \
            "chr varchar (100)," \
            "gene_name varchar (100)," \
            "p_value varchar (100),"\
            "rsid varchar (100)," \
            "table_name varchar (100)," \
            "trait varchar (100));"


    print(query)
    cursor.execute(query)
    cursor.commit()
    print("table %s successfully created in database %s" % (tablename, database))



    # check if path exists
    if not os.path.isdir(path):
        print('path does not exist')
        exit(1)

        # change directory
    os.chdir(path)

    fileNames = []

    # for each of the files in the dir ending with txt, add to the list
    for file in os.listdir(path):
        if file.endswith(".txt"):
            fileNames.append(file)

    for fileName in fileNames:

        print("Reading from: " + str(fileName))

        with open(fileName, 'r') as txtFile:

            txtReader = csv.reader(txtFile, delimiter='\t')

            fileFormat = next(txtReader)

            errorCounter = 0

            index = 0
            # Resort each row to resemble the database format
            for rows in txtReader:

                index = index + 1

                list = []
                wanted = [0,1, 2, 3, 4, 5, 6, 7]

                listCounter = 0

                for cols in rows[0:]:
                    if listCounter in wanted:
                        list.append(cols)
                    listCounter += 1

                try:
                    query = "insert into dbo.{!s}".format(tablename) + \
                            " values ( "

                    for j in range(0, (len(list)-1)):
                        query += ("{!r}".format(list[j]) + ", ")

                    query += "{!r}".format(list[len(list) - 1])
                    query += ");"

                    cursor.execute(query)
                    cursor.commit()

                # write errmsg if file I/O exception
                except Exception as eex:
                    errorCounter += 1

                    if errorCounter == 1:
                        f = open("GC_{!r}_err.txt".format(tablename), "w")
                        f.write(str(list) + "\n")
                        print("ERROR in index " + str(index) + "!")
                    else:
                        print("ERROR in index " + str(index) + "!")
                        f.write(str(list) + "\n")

                else:
                    print("Insert " + str(index) + " was successful!")

    cursor.commit()
    print("File Read Done!" + str(fileName))
    cnxn.close()

def get_df(sql_conn, hg19, chosen_table, margin, cutoff):
    df = pd.DataFrame(columns=["chr", "bp", "beta", "p_value", "trait", "rsid","table_name"])
    for i in chosen_table:
        ct = 1
        if ct == i:
            query = "SELECT chr,rsid, bp, beta, p_value, trait, 'Lipid_EHR_Hoffmann_NG' as table_name FROM Lipid_EHR_Hoffmann_NG WHERE "
            query += where(hg19,cutoff)
            record_base = pd.read_sql(query, sql_conn)
            record_base = record_base.astype({"chr": "category","rsid": "category", "bp": "int64", "beta": "float64",  "p_value": "float64", "table_name": "category"})
            df = df.append(record_base)
            continue

        ct += 1
        if ct == i:
            query = "SELECT chr,rsid, bp, beta, p_value, trait, 'Lipid_Teslovich_Nature' as table_name FROM Lipid_Teslovich_Nature WHERE "
            query +=where(hg19,cutoff)
            tesl = pd.read_sql(query, sql_conn)
            tesl = tesl.astype({"chr": "category", "rsid": "category", "bp": "int64", "beta": "float64",  "p_value": "float64", "table_name": "category"})
            df = df.append(tesl)
            continue
        ct += 1

        if ct == i:
            query = "SELECT chr,rsid, bp, beta, p_value, trait, 'Lipid_Korean_Exome_Moon_SR' as table_name FROM Lipid_Korean_Exome_Moon_SR WHERE "
            query +=where(hg19,cutoff)
            Korean_exome = pd.read_sql(query, sql_conn)
            Korean_exome = Korean_exome.astype({"chr": "category", "rsid": "category", "bp": "int64", "beta": "float64",  "p_value": "float64", "table_name": "category"})
            df = df.append(Korean_exome)
            continue
        ct += 1

        if ct == i:
            query = "SELECT chr, rsid, bp, beta, p_value, trait, 'Lipid_Engage_Surakka_NG' as table_name FROM Lipid_Engage_Surakka_NG WHERE "
            query += where_varchar(hg19, cutoff)
            surakka = pd.read_sql(query, sql_conn)
            surakka = surakka.astype({"chr": "category", "rsid": "category", "bp": "int64", "beta": "float64",  "p_value": "float64", "table_name": "category"})
            df = df.append(surakka)
            continue
        ct += 1

        if ct == i:
            query = "SELECT chr, rsid, bp, beta, p_value, trait, 'Lipid_Exome_Lu_East_Asian_NG' as table_name FROM Lipid_Exome_Lu_East_Asian_NG WHERE "
            query +=where(hg19,cutoff)
            east_asian = pd.read_sql(query, sql_conn)
            east_asian = east_asian.astype({"chr": "category", "rsid": "category", "bp": "int64", "beta": "float64",  "p_value": "float64", "table_name": "category"})
            df = df.append(east_asian)
            continue
        ct += 1

        if ct == i:
            query = "SELECT chr, rsid, bp, beta, p_value, trait, 'Lipid_Exome_Lu_European_and_East_Asian_NG' as table_name FROM Lipid_Exome_Lu_European_and_East_Asian_NG WHERE "
            query +=where(hg19,cutoff)
            european = pd.read_sql(query, sql_conn)
            european = european.astype({"chr": "category", "rsid": "category", "bp": "int64", "beta": "float64",  "p_value": "float64", "table_name": "category"})
            df = df.append(european)
            continue
        ct += 1

        if ct == i:
            query = "SELECT chr, rsid, bp, beta, p_value, trait, 'Lipid_GLGC_Willer_NG' as table_name FROM Lipid_GLGC_Willer_NG WHERE "
            query +=where(hg19,cutoff)
            glgc = pd.read_sql(query, sql_conn)
            glgc = glgc.astype({"chr": "category", "rsid": "category", "bp": "int64", "beta": "float64",  "p_value": "float64", "table_name": "category"})
            df = df.append(glgc)
            continue
        ct += 1

        if ct == i:
            query = "SELECT chr, rsid, bp, beta, p_value, trait, 'Lipid_Japanese_lipid_trait_Kanai_NG' as table_name FROM Lipid_Japanese_lipid_trait_Kanai_NG WHERE "
            query +=where_varchar(hg19,cutoff)
            lipid_japanese = pd.read_sql(query, sql_conn)
            lipid_japanese = lipid_japanese.astype({"chr": "category", "rsid": "category", "bp": "int64", "beta": "float64",  "p_value": "float64", "table_name": "category"})
            df = df.append(lipid_japanese)
            continue
        ct += 1

        if ct == i:
            query = "SELECT chr, rsid, bp, beta, p_value, trait, 'Lipid_MVP_Klarin_NG' as table_name FROM Lipid_MVP_Klarin_NG WHERE "
            query +=where(hg19,cutoff)
            lipid_mvp = pd.read_sql(query, sql_conn)
            lipid_mvp = lipid_mvp.astype({"chr": "category", "rsid": "category", "bp": "int64", "beta": "float64",  "p_value": "float64", "table_name": "category"})
            df = df.append(lipid_mvp)
            continue
        ct += 1

        if ct == i:
            query = "SELECT chr, rsid, bp, beta, p_value, trait, 'Lipid_Spracklen_Hum_Mol_Genetics' as table_name FROM Lipid_Spracklen_Hum_Mol_Genetics WHERE "
            query +=where_varchar(hg19,cutoff)
            lipid_spracklen = pd.read_sql(query, sql_conn)
            lipid_spracklen = lipid_spracklen.astype({"chr": "category", "rsid": "category", "bp": "int64", "beta": "float64",  "p_value": "float64", "table_name": "category"})
            df = df.append(lipid_spracklen)
            continue
        ct += 1

        if ct == i:
            query = "SELECT chr, rsid, bp, beta, p_value, trait, 'Lipid_UKBB_high_cholesterol_ukbb_Connor_alkesgroup' as table_name FROM Lipid_UKBB_high_cholesterol_ukbb_Connor_alkesgroup WHERE "
            query +=where_varchar(hg19,cutoff)
            high_chol = pd.read_sql(query, sql_conn)
            high_chol = high_chol.astype({"chr": "category", "rsid": "category", "bp": "int64", "beta": "float64",  "p_value": "float64", "table_name": "category"})
            df = df.append(high_chol)
            continue
        ct += 1

        if ct == i:
            query = "SELECT chr, rsid, bp, beta, p_value, trait, 'Lipid_UKBB_lipid_trait_Neale' as table_name FROM Lipid_UKBB_lipid_trait_Neale WHERE "
            query +=where(hg19,cutoff)
            ukbb_lipid_trait = pd.read_sql(query, sql_conn)
            ukbb_lipid_trait = ukbb_lipid_trait.astype({"chr": "category", "rsid": "category", "bp": "int64", "beta": "float64",  "p_value": "float64", "table_name": "category"})
            df = df.append(ukbb_lipid_trait)
            continue
        ct += 1

        if ct == i:
            query = "SELECT chr, rsid, bp, beta, p_value, trait, 'Lipid_UKBB_statin_usage_Neale' as table_name FROM Lipid_UKBB_statin_usage_Neale WHERE "
            query +=where(hg19,cutoff)
            ukbb_statin = pd.read_sql(query, sql_conn)
            ukbb_statin = ukbb_statin.astype({"chr": "category", "rsid": "category", "bp": "int64", "beta": "float64",  "p_value": "float64", "table_name": "category"})
            df = df.append(ukbb_statin)
            continue

        ct += 1

        if ct == i:
            query = "SELECT chr, rsid, bp, beta, p_value, trait, 'Lipid_GLGC_Exome_Liu_NG' as table_name FROM Lipid_GLGC_Exome_Liu_NG WHERE "
            query +=where(hg19,cutoff)
            GLGC_exome = pd.read_sql(query, sql_conn)
            GLGC_exome = GLGC_exome.astype({"chr": "category", "rsid": "category", "bp": "int64", "beta": "float64",  "p_value": "float64", "table_name": "category"})
            df = df.append(GLGC_exome)
            continue
        ct += 1

        if ct == i:
            query = "SELECT chr, rsid, bp, beta, p_value, trait, 'Lipid_UK10K_prins_SR' as table_name FROM Lipid_UK10K_prins_SR WHERE "
            query += where(hg19, cutoff)
            uk10k = pd.read_sql(query, sql_conn)
            uk10k = uk10k.astype(
                    {"chr": "category", "rsid": "category", "bp": "int64", "beta": "float64", "p_value": "float64",
                     "table_name": "category"})
            df = df.append(uk10k)
            continue

        ct += 1

        if ct == i:
            query = "SELECT chr, rsid, bp, beta, p_value, trait, 'Lipid_Bivariate_Siewert_CGPM' as table_name FROM Lipid_Bivariate_Siewert_CGPM WHERE "
            query += where(hg19, cutoff)
            biavi = pd.read_sql(query, sql_conn)
            biavi = biavi.astype(
                    {"chr": "category", "rsid": "category", "bp": "int64", "beta": "float64", "p_value": "float64",
                     "table_name": "category"})
            df = df.append(biavi)
            continue
        ct += 1

        if ct == i:
            query = "SELECT chr, rsid, bp, beta, p_value, trait, 'Lipid_Finnish_Locke' as table_name FROM Lipid_Finnish_Locke WHERE "
            query += where(hg19, cutoff)
            finnish = pd.read_sql(query, sql_conn)
            finnish = finnish.astype(
                    {"chr": "category", "rsid": "category", "bp": "int64", "beta": "float64", "p_value": "float64",
                     "table_name": "category"})
            df = df.append(finnish)
            continue


        ct += 1

        if ct == i:
            query = "SELECT chr, rsid, bp, beta, p_value, trait, 'Lipid_Chinese_exome_Lu_HMG' as table_name FROM Lipid_Chinese_exome_Lu_HMG WHERE "
            query += where(hg19, cutoff)
            chinese = pd.read_sql(query, sql_conn)
            chinese = chinese.astype(
                    {"chr": "category", "rsid": "category", "bp": "int64", "beta": "float64", "p_value": "float64",
                     "table_name": "category"})
            df = df.append(chinese)
            continue
    df = df.sort_values(by=['p_value'], ascending=True)
    return(df)

def save_option(df, sql_conn):
    option = input("""How do you want to save this result?
1. Save as txt file Only
2. Save and Push to SQL Server
3. Plot data
4. Nothing:\n""")

    while option not in ['1','2','3','4']:
        print("Enter a valid input\n")
        option = input("""How do you want to save this result?
1. Save as txt file Only
2. Save and Push to SQL Server
3. Plot data
4. Nothing:\n""")


    try:
        if option == '1':
            path = 'E:/cross_ref/output/save_file.txt'
            df.to_csv(path, encoding = 'utf-8', sep='\t',index = False)
            # os.startfile(path)
            choice = input ("Name your txt file to save:\n")
            os.rename ('E:/cross_ref/output/save_file.txt','E:/cross_ref/output/%s.txt' % choice)

        elif option == '2':
            path = 'E:/cross_ref/output/save_file.txt'
            df.to_csv(path, encoding = 'utf-8', sep='\t',index = False)
            # os.startfile(path)
            choice = input ("Name your txt file to save:\n")
            os.rename ('E:/cross_ref/output/save_file.txt','E:/cross_ref/output/%s.txt' % choice)

            newpath = 'E:/cross_ref/load/cross_ref.txt'
            df.to_csv(newpath, encoding = 'utf-8', sep= '\t', index = None)
            # os.startfile(path)
            print("load this table to server")
            load_to_sql()

        elif option == '3':
            # setgroup = set();
            # nameoftable = df [['table_name']]
            # setgroup.add(nameoftable)
            grouped = df.groupby('table_name')
            print (grouped)
            print(grouped ['trait'])
            # print (df[['p_value','trait','table_name']])

    except PermissionError:
        print('Close the csv file currently opened and try again')
        pass

def main():
    # Connect to SQL Server
    sql_conn =  pyodbc.connect("""
                               DRIVER={ODBC Driver 11 for SQL Server};
                               SERVER=PARKSLAB;
                               DATABASE=Human_GWAS;
                               Trusted_Connection=Yes;
                               """)

    while True:
        # Input for which tables&chromosome to include and gene name and margin(+/-)
        chosen_table = get_table()
        print('')

        # chosen_chr   = get_chr()
        # print('')

        gene_name_list = get_genename()
        print('')

        margin         = get_margin()
        print('')

        cutoff         = get_pvalue()

        df = pd.DataFrame(columns=["chr", "bp", "beta", "p_value", "trait", "table_name"])
        for gene_name in gene_name_list:
            hg19 = get_hg19(sql_conn, gene_name, margin)
            df_temp   = get_df(sql_conn, hg19, chosen_table, margin, cutoff)
            df_temp['gene_name'] = gene_name
            df = df.append(df_temp, sort = True)

        print(df.to_string())
        save_option(df,sql_conn)


if __name__ == "__main__":
    main()
