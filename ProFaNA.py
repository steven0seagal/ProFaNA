import pandas
from collections import Counter
import pandas as pd
from statsmodels.stats.multitest import multipletests as correction  #
import datetime
import json
import math
from scipy.stats import poisson

""" Functions common for all classes """


def pfam_to_go(domains):
    """
    Function that takes list of domains, parse GO database and match GO terms for each domain.

    Args:
        domains (list): list of pfam domains

    Returns:
        biological_process (dict): dictionary with GO terms for biological process
    """

    # open GO-PFAM map
    with open('/usr/src/app/important_files/go_pfam', 'r') as handler:
        go_pfams = [x.strip().split('\t') for x in handler]
    # open GO database
    with open('/usr/src/app/important_files/go_terms', 'r') as handler:
        go_terms = [x.strip().split('\t') for x in handler]



    # found GO terms for each domain
    founded_gos = []
    for domain in domains:
        for go_pfam in go_pfams:
            if domain == go_pfam[0]:
                founded_gos.append((domain, go_pfam[3]))

    # create dicts for each ontology
    biological_process = {}
    molecular_function = {}
    cellular_component = {}
    for go in founded_gos:
        for go_term in go_terms:
            if go[1] == go_term[0]:
                if go_term[2] == 'biological_process':
                    biological_process[go[1]] = {'term': go_term[1], 'domains': []}
                elif go_term[2] == 'molecular_function':
                    molecular_function[go[1]] = {'term': go_term[1], 'domains': []}
                elif go_term[2] == 'cellular_component':
                    cellular_component[go[1]] = {'term': go_term[1], 'domains': []}

    # add description for each go
    for go in founded_gos:
        if go[1] in biological_process.keys():
            biological_process[go[1]]['domains'].append(go[0])
        elif go[1] in molecular_function.keys():
            molecular_function[go[1]]['domains'].append(go[0])
        elif go[1] in cellular_component.keys():
            cellular_component[go[1]]['domains'].append(go[0])

    return biological_process, molecular_function, cellular_component


def parse_go_classification(dict_data):
    """
    Function that change dict GO classification to dataframe with GO terms and domains

    Args:
        dict_data (dict): dictionary with GO terms and domains

    Returns:
        df (dataframe): dataframe with GO terms and domains

    """
    df = pd.DataFrame(columns=['GO_id', 'term', 'domains'])
    for k, v in dict_data.items():
        data_to_append = [k, v['term'], ','.join(v['domains'])]
        df.loc[len(df)] = data_to_append
    return df


def match_go_to_pfam_and_save(dataframe, file_name):
    """
    Function that collect results and GO terms to save it to file.

    Args:
        dataframe (pd.Dataframe): dataframe with results
        file_name (str): name of file with results
    """
    # filter dataframe to show only significant results
    try:
        filter_poisson = dataframe.loc[dataframe['Poisson_correct'].astype(float) < 0.00001]
    except KeyError:
        filter_poisson = dataframe.loc[dataframe['Poisson'].astype(float) < 0.00001]

    # filter dataframe to show only results with percentage of family > 5%
    filter_poisson['In what percentage'] = filter_poisson['In what percentage'].str.rstrip('%').str.replace(',',
                                                                                                            '').astype(
        'float')
    filter_percentage = filter_poisson.loc[filter_poisson['In what percentage'] > 0.05]

    # extract domains
    domains = list(filter_percentage.index.values)

    # match GO terms to domains
    biological_process, molecular_function, cellular_component = pfam_to_go(domains)
    df_bp = parse_go_classification(biological_process)
    df_mf = parse_go_classification(molecular_function)
    df_cp = parse_go_classification(cellular_component)

    # save
    new_file_name = file_name.split('.')[0]
    with open(f'/usr/src/app{new_file_name}_ontology.txt', 'w') as f:
        f.write('BIOLOGICAL PROCESS')
        f.write("\n")
        f.write('Index')
        df_bp.to_csv(f, sep='\t')
        f.write("\n")
        f.write('MOLECULAR FUNCTION')
        f.write("\n")
        f.write('Index')
        df_mf.to_csv(f, sep='\t')
        f.write("\n")
        f.write('CELLULAR COMPONENT')
        f.write('Index')
        f.write("\n")
        df_cp.to_csv(f, sep='\t')


def match_go_to_pfam_and_save_average(dataframe, file_name):
    """
    Function that collect results and GO terms to save it to file.

    Args:
        dataframe (pd.Dataframe): dataframe with results
        file_name (str): name of file with results

    """
    # filter dataframe to show only significant results
    try:
        filter_poisson = dataframe.loc[dataframe['Poisson_correct'].astype(float) < 0.00001]
    except KeyError:
        filter_poisson = dataframe.loc[dataframe['Poisson'].astype(float) < 0.00001]

    # filter dataframe to show only results with percentage of family > 5%
    filter_poisson['Family percentage'] = filter_poisson['Family percentage'].astype('float')
    filter_percentage = filter_poisson.loc[filter_poisson['Family percentage'] > 0.05]

    # extract domains
    domains = list(filter_percentage.index.values)

    # match GO terms to domains
    biological_process, molecular_function, cellular_component = pfam_to_go(domains)
    df_bp = parse_go_classification(biological_process)
    df_mf = parse_go_classification(molecular_function)
    df_cp = parse_go_classification(cellular_component)
    # save
    new_file_name = file_name.split('.')[0]
    with open(f'/usr/src/app/media/results/{new_file_name}_ontology.txt', 'w') as f:
        f.write('BIOLOGICAL PROCESS')
        f.write("\n")
        f.write('Index')
        df_bp.to_csv(f, sep='\t')
        f.write("\n")
        f.write('MOLECULAR FUNCTION')
        f.write("\n")
        f.write('Index')
        df_mf.to_csv(f, sep='\t')
        f.write("\n")
        f.write('CELLULAR COMPONENT')
        f.write('Index')
        f.write("\n")
        df_cp.to_csv(f, sep='\t')


def open_singleline(path_to_file):
    """
    Opens file that contain 1 column and strip it by space.

    Args:
        path_to_file (str): path to file to open

    Returns:
        file_names (list): list with file names
    """

    with open(path_to_file, 'r') as f:
        file_names = [line.strip() for line in f]
    return file_names


def open_singleline_num(path_to_file):
    """
    Opens file that contain 1 column and strip it by space.

    Args:
        path_to_file (str): path to file to open

    Returns:
        file_names (list): list with file names
    """

    with open(path_to_file, 'r') as f:
        file_names = [float(line.strip()) for line in f]
    return file_names


def open_multiple_line(path_to_file):
    """
    Opens file that contains data in multiple columns and strip it by space.

    Args:
        path_to_file (str): path to file to open

    Returns:
        file_names (list): list with file names

    """

    data = []
    with open(path_to_file) as input_file:
        for line in input_file:
            data.append(line.strip().split())
    return data


def save_data(file_name, complete_data):
    """
    Saves all stuff together as one file

    Args:
        file_name (str): name of file to save
        complete_data (pd.Dataframe): list with data to save
    """

    complete_data.index.name = 'Domain'
    complete_data.to_csv('/usr/src/app/' + file_name, sep='\t', mode='w')


def create_six_list(file):
    """
    Open single genome file and chew it and return 6 lists -> GENE, START_COORD,
    END_COORDS ,ORIENTATION, DOMAINS

    Args:
        file (str): path to file to open

    Returns:
        genes (list): list with genes
        start_coords (list): list with start coordinates
        end_coords (list): list with end coordinates
        orientation (list): list with orientation
        domains (list): list with domains
        contig (list): list with contigs
    """

    data = open_multiple_line("/usr/src/app/important_files/nowa_baza_danych/" + file)
    start_coord = []
    end_coord = []
    orientation = []
    domains = []
    genes = []
    contig = []
    for bit in data:
        genes.append(int(bit[0]))
        start_coord.append(int(bit[1]))
        end_coord.append(int(bit[2]))
        orientation.append(bit[3])
        domains.append(bit[4])
        contig.append(bit[5])
    return genes, start_coord, end_coord, orientation, domains, contig


def open_database(tax):
    """
    Function that open database and return it as list of genomes

    Args:
        tax (str): taxonomic group

    Returns:
        genomes (list): list of genomes
    """
    tax = tax.replace("_", " ")
    # all genomes
    if tax == "all genomes":
        with open("/home/bartek/PycharmProjects/ProFaNA-small/modules/data/all_genomes", "r") as handler:
            genomes = [x.strip() for x in handler]
        return genomes

    # just entero
    elif tax == "entero":
        with open("/usr/src/app/important_files/entero", "r") as handler:
            genomes = [x.strip() for x in handler]
        return genomes

    # everything else
    else:
        with open("/usr/src/app/important_files/database.json", "r") as handler:
            data = json.load(handler)
        if tax in data.keys():
            genomes = []
            for genus, genome in data[tax].items():
                genomes += genome
            return genomes
        else:
            for family, genera in data.items():
                if tax in genera:
                    genomes = data[family][tax]
                    return genomes


def genome_size_in_gene(genome_id, list_of_genome, list_of_size):
    """
    Takes genome's id, list of genome and lists with data about how much genes in genome

    Args:
        genome_id (str): genome's id
        list_of_genome (list): list with genome's id
        list_of_size (list): list with data about how much genes in genome

    Returns:
        size (int): how much genes in genome
    """

    index_genome = list_of_genome.index(genome_id)
    genome_size_in_gene = list_of_size[index_genome]
    return genome_size_in_gene


def size_in_genes(genome_ids):
    """
    Takes list of genes and returns list of uniqe genes

    Args:
        genome_ids (list): list with genome's id

    Returns:
        uniqe_genome_ids (list): list with uniqe genome's id
    """

    list_of_genes = []
    for gene in genome_ids:
        list_of_genes.append(gene)
    return len(set(list_of_genes))


def list_of_genes_in_neigh(list_of_genes, index_list):
    """
    Takes overall gene's list in genome and creates complete list of genes in neighbourhood

    Args:
        list_of_genes (list): list with genes in genome
        index_list (list): list with indexes of genes in neighbourhood

    Returns:
        list_of_genes_in_neigh (list): list with genes in neighbourhood
    """

    list_of_genes_in_neigh = []
    for pfam_domain in index_list:
        list_of_genes_in_neigh.append(list_of_genes[pfam_domain])
    return list_of_genes_in_neigh


def searching_for_domain_in_genome(pfam, start_coord, end_coord, orient, domains,contig, genes):
    """
    Takes user's pfam domain searches through lists contains coordinates,
    orientation, domains, contigs and returns complete data about all
    users' domains in genome

    Args:
        pfam (str): pfam domain
        start_coord (list): list with start coordinates
        end_coord (list): list with end coordinates
        orient (list): list with orientation
        domains (list): list with domains
        contig (list): list with contigs
        genes (list): list with genes

    Returns:
        coords (list): list with coordinates
    """
    coords = []
    coords_counter = 0
    for domain in domains:
        one_coords_part = []
        if domain == pfam:
            orientation_pfam_domain = orient[coords_counter]
            one_coords_part.append(start_coord[coords_counter])
            one_coords_part.append(end_coord[coords_counter])
            one_coords_part.append(orientation_pfam_domain)
            one_coords_part.append(contig[coords_counter])
            one_coords_part.append(genes[coords_counter])
            coords.append(one_coords_part)
            coords_counter += 1
        else:
            coords_counter += 1
            continue
    return coords


def presence_confirmation(coords, file, pfam):
    """
    Just prints out information how many pfam domains is in each genome
    during analysis probably I will have to remove it before putting it to website ...

    Args:
        coords (list): list with coordinates
        file (str): path to file
        pfam (str): pfam domain

    Returns:
        None
    """

    if len(coords) == 0:
        print("In genome  " + file + " pfam domain you have been searching do not exist")
    elif len(coords) == 1:
        print("In genome  " + file + " there is " + str(len(coords)) + " " + pfam + " domain")
    else:
        print("In genome  " + file + " there are " + str(len(coords)) + " " + pfam + " domains")


def get_range_coordinates(target_pfam, distance):
    """
    Based on how large neighbourhood user wants to analyze creates
    points FROM and TO, additionally shows where main user's pfam begins and ENDS
    with orientation on strands

    Args:
        target_pfam (str): pfam domain
        distance (str): distance

    Returns:
        range_coordinates (list): list with coordinates
    """

    last_coordinate = target_pfam[1] + int(distance)
    first_coordinate = target_pfam[0] - int(distance)
    pfam_beg = target_pfam[0]
    pfam_end = target_pfam[1]
    searched_pfam_orientation = target_pfam[2]
    return last_coordinate, first_coordinate, pfam_beg, pfam_end, searched_pfam_orientation


def get_domain_both_zeros(ref_gene, q_gene):
    """
    Takes list of genes in neighbourhood and gene from query genome and returns index where genes are placed where
    distance is zero.

    Args:
        ref_gene (list): list with genes in genome
        q_gene (str): query gene

    Returns:
        index_list_gene (list): list with indexes of query gene
    """
    index_list_gene = []
    counter = 0
    for gene in ref_gene:
        if gene == int(q_gene):
            index_list_gene.append(counter)
        counter += 1
    return index_list_gene


def get_domains_both_strands(start_coord, end_coord, last_coordinate, first_coordinate,
                             pfam_beg, pfam_end, contig, point):
    """
    Takes list of coordinates, list of contigs, list of genes and coordinates of neighbourhood and returns list of
    indexes of domains in neighbourhood.

    Args:
        start_coord (list): list with start coordinates
        end_coord (list): list with end coordinates
        last_coordinate (int): last coordinate
        first_coordinate (int): first coordinate
        pfam_beg (int): pfam domain start
        pfam_end (int): pfam domain end
        contig (list): list with contigs
        point (list): list with genes

    Returns:
        pfam_index_to_neigh (list): list with indexes of pfam domains
    """
    pfam_index_to_neigh = []

    # Start coordinates
    s_coord_counter = 0
    for s_coord in start_coord:
        if last_coordinate >= s_coord >= pfam_beg and contig[s_coord_counter] == point[3] and s_coord_counter not in \
                pfam_index_to_neigh:
            pfam_index_to_neigh.append(s_coord_counter)
            s_coord_counter += 1
        else:
            s_coord_counter += 1
            continue
    # End coordinates
    e_coord_counter = 0
    for e_coord in end_coord:
        if first_coordinate <= e_coord <= pfam_end and contig[e_coord_counter] == point[3] and e_coord_counter not in \
                pfam_index_to_neigh:
            pfam_index_to_neigh.append(e_coord_counter)
            e_coord_counter += 1
        else:
            e_coord_counter += 1
            continue
    return pfam_index_to_neigh


def collect_pfam_domain_description(file):
    """
    Opens file only to gather information about pfam domains and returns dictionary with pfam domains as keys and
    description as values

    Args:
        file (str): path to file

    Returns:
        data (dict): dictionary with pfam domains as keys and description as values
    """

    data = {}
    with open(file) as input_file:
        for line in input_file:
            try:
                one_line = line.strip().split("@")
                domain = one_line[0]
                pre_family = one_line[1][8:]
                delete_this = "(" + domain.upper() + ")"
                family = pre_family.replace(delete_this, "")
                summary = one_line[2][9:]
                data[domain] = (family, summary)

            except IndexError:
                continue
    return data


def get_domains_plus_minus_strand(start_coord, end_coord, last_coordinate, first_coordinate,
                                  pfam_beg, pfam_end, contig, point, orientation):
    """
    Takes list of coordinates, list of contigs, list of genes and coordinates of neighbourhood and returns list of
    indexes of domains in neighbourhood for opposite strand and same strand.

    Args:
        start_coord (list): list with start coordinates
        end_coord (list): list with end coordinates
        last_coordinate (int): last coordinate
        first_coordinate (int): first coordinate
        pfam_beg (int): pfam domain start
        pfam_end (int): pfam domain end
        contig (list): list with contigs
        point (list): list with genes
        orientation (list): list with orientation of pfam domains

    Returns:
        pfam_index_to_neigh_same_strand (list): list with indexes of pfam domains on the same strand
        pfam_index_to_neigh_opposite_strand (list): list with indexes of pfam domains on the opposite strand
    """

    pfam_index_to_neigh_same_strand = []
    pfam_index_to_neigh_oposite_strand = []
    s_coord_counter = 0
    for s_coord in start_coord:
        if last_coordinate >= s_coord >= pfam_beg and orientation[s_coord_counter] == point[2] and \
                contig[s_coord_counter] == point[3] and s_coord_counter not in pfam_index_to_neigh_same_strand:
            pfam_index_to_neigh_same_strand.append(s_coord_counter)
            s_coord_counter += 1
        elif last_coordinate >= s_coord >= pfam_beg and orientation[s_coord_counter] != point[2] and \
                contig[s_coord_counter] == point[3] and s_coord_counter not in pfam_index_to_neigh_oposite_strand:
            pfam_index_to_neigh_oposite_strand.append(s_coord_counter)
            s_coord_counter += 1
        else:
            s_coord_counter += 1
            continue
    e_coord_counter = 0
    for e_coord in end_coord:
        if first_coordinate <= e_coord <= pfam_end and orientation[e_coord_counter] == point[2] and \
                contig[e_coord_counter] == point[3] and e_coord_counter not in pfam_index_to_neigh_same_strand:
            pfam_index_to_neigh_same_strand.append(e_coord_counter)
            e_coord_counter += 1
        elif first_coordinate <= e_coord <= pfam_end and orientation[e_coord_counter] != point[2] and \
                contig[e_coord_counter] == point[3] and e_coord_counter not in pfam_index_to_neigh_oposite_strand:
            pfam_index_to_neigh_oposite_strand.append(e_coord_counter)
            e_coord_counter += 1
        else:
            e_coord_counter += 1
            continue
    return pfam_index_to_neigh_same_strand, pfam_index_to_neigh_oposite_strand


def get_list_of_domain_in_neigh(pfam_index, file, domains):
    """
    Takes list of indexes of pfam domains and returns list of pfam domains in neighbourhood.
    Args:
        pfam_index (list): list of indexes of pfam domains
        file (str): path to file
        domains (list): list of pfam domains
    Returns:
        pfam_domains_in_neigh (list): list pf pfam domains in neighbourhood
    """
    to_counter = []
    party = []
    party.append(file)
    for part in pfam_index:
        party.append(domains[part])
        to_counter.append(domains[part])

    return to_counter


def multiple_test_correction(some_dataframe, correction_met):
    """
    Takes dataframe and correction method and returns dataframe with corrected p-values

    Args:
        some_dataframe (pd.Dataframe): dataframe with p-values
        correction_met (str): correction method

    Returns:
          some_dataframe (pd.Dataframe): dataframe with corrected p-values
    """
    if correction_met == 'none':
        return some_dataframe
    else:
        value_to_correct = [float(x) for x in some_dataframe.PVALUE.tolist()]
        reject, pvals_corrected, alphaSidak, alphaBonf = correction(pvals=value_to_correct,
                                                                    method=correction_met,
                                                                    is_sorted=False, returnsorted=False)
        pvals = pvals_corrected.tolist()
        some_dataframe.PVALUE = pvals
        return some_dataframe


def multiple_test_correction_list(some_tuple_list, correction_met):
    """
    Takes list of tuples and correction method and returns list of tuples with corrected p-values

    Args:
        some_tuple_list (list): list of tuples of values
        correction_met (str): correction method

    Returns:
        outout (list): list of tuples of values with corrected p-values
    """
    pvalues = [float(x) for x in some_tuple_list]

    if correction_met == 'none':
        output = pvalues

    else:
        try:

            reject, pvals_corrected, alphaSidak, alphaBonf = correction(pvals=pvalues,
                                                                        method=correction_met,
                                                                        is_sorted=False, returnsorted=False)
            pvals_after_corr = pvals_corrected.tolist()
            output = pvals_after_corr
        except ZeroDivisionError:
            output = pvalues
    return output


def cutoff_value(some_dataframe, cutoff):
    """
    Takes dataframe and cutoff value and returns dataframe with p-values lower than cutoff value

    Args:
        some_dataframe (pd.Dataframe): dataframe with p-values
        cutoff (str): cutoff value

    Returns:
        some_dataframe (pd.Dataframe): dataframe with p-values lower than cutoff value
    """
    domain_list = list(some_dataframe.index)
    if cutoff == 'none':
        return some_dataframe
    elif cutoff == '0':
        for i in domain_list:
            diff = float(some_dataframe.loc[i, 'Poisson_correct'])
            if diff < 0 < diff:
                some_dataframe = some_dataframe.drop([i])
        return some_dataframe
    else:
        cutoff = float(cutoff)
        for i in domain_list:
            diff = float(some_dataframe.loc[i, 'Poisson_correct'])
            if diff > cutoff:
                some_dataframe = some_dataframe.drop([i])
        return some_dataframe


def sort_table(some_dataframe,):
    """
    Takes dataframe and sorts it by p-value

    Args:
        some_dataframe (pd.Dataframe): dataframe with p-values

    Returns:
        sorted_dataframe (pd.Dataframe): sorted dataframe with p-values
    """
    sorted_data = some_dataframe.sort_values('Poisson', ascending=True)
    return sorted_data


def add_information(some_dataframe, dictionary):
    """
    Takes dataframe and dictionary and adds domain information from dictionary to dataframe

    Args:
        some_dataframe (pd.Dataframe): dataframe with p-values
        dictionary (dict): dictionary with domain information

    Returns:
        some_dataframe (pd.Dataframe): dataframe with p-values and domain information
    """
    indeksy = list(some_dataframe.index)
    for i in indeksy:
        try:

            pfam = i[0:2] + i[4:]
            family = dictionary[pfam][0]
            summary = dictionary[pfam][1]

            some_dataframe.at[i, 'Family'] = family
            some_dataframe.at[i, 'Summary'] = summary
        except KeyError:
            continue
    return some_dataframe



def poisson_distribution(p, n, domain, percentage_counter):
    """
    Return value of poisson distribution

    Args:
        p (int): number of domain in all database
        n (int): number of genomes
        domain (str): domain name
        percentage_counter (dict): dictionary with domain occurrence in all database

    Returns:
        float: value of poisson distribution

    """
    if int(percentage_counter[domain]) < p * n:
        return 1
    else:
        sums = 0
        k = 0
        number_of_samples = int(percentage_counter[domain])
        while k < number_of_samples:
            sums += poisson.pmf(k=k, mu=n * p)
            k += 1
        if 1 - sums >= 0.0:
            return 1 - sums
        else:
            return 0.0


def pfam_for_pf(dataframe, correction):
    """
    Final customization of dataframe
        --> PF02696 instead of pfam02696
        --> PVALUE in 10e-3 format
        --> In what percentage % format
        --> avg occurences and density  have now 3 digits after coma
        --> drop NOam
        --> drop columns with 0 values
    Args:
        dataframe (pd.Dataframe): dataframe with p-values
        correction (str): correction method
    Returns:
        dataframe (pd.Dataframe): final dataframe
    """

    indexy = dataframe.index
    for i in indexy:
        dataframe = dataframe.rename(index={i: i[:2].upper() + i[4:]})
    dataframe['Poisson'] = dataframe['Poisson'].map('{:.3e}'.format)
    dataframe['Poisson_correct'] = dataframe['Poisson_correct'].map('{:.3e}'.format)

    dataframe['In what percentage'] = dataframe['In what percentage'].map('{:.3%}'.format)
    dataframe['average occurence in neighborhood'] = dataframe['average occurence in neighborhood'].map(
        '{:.3}'.format)
    dataframe['average occurence in genome'] = dataframe['average occurence in genome'].map('{:.3}'.format)

    if "NOam" in indexy:
        dataframe = dataframe.drop(index="NOam")
    # dataframe = dataframe.style.format({'In what percentage': '{:,.3%}'.format})
    dataframe.drop('average occurence in neighborhood', axis=1, inplace=True)
    dataframe.drop('average occurence in genome', axis=1, inplace=True)
    dataframe.drop('avg neigh size', axis=1, inplace=True)
    dataframe.drop('domain in whole database', axis=1, inplace=True)
    if correction == 'none':
        dataframe.drop('Poisson_correct', axis=1, inplace=True)
    else:
        dataframe.drop('Poisson', axis=1, inplace=True)
    try:
        dataframe = dataframe.drop(index="NOam")
    except KeyError:
        return dataframe
    return dataframe



def open_single_column_file(file_name):
    """
    Opens file that contains one column of data
    Args:
        file_name (str): name of file
    Returns:
        list: list of data
            """
    plik = []
    with open(file_name) as inputfile:
        for line in inputfile:
            plik.append(line.strip())
    return plik


def open_multiple_column_file(file_name, split_mark=" "):
    """
    Opens file that contain more than one column and split it by space.
    Args:
        file_name (str): name of file
        split_mark (str): mark that split columns
    Returns:
        list: list of data
    """

    plik = []
    with open(file_name) as inputfile:
        for line in inputfile:
            plik.append(line.strip().split(split_mark))
    return plik


def open_json_file(file_name):
    """
    Opens json file
    Args:
        file_name (str): name of file

    Returns:
        dict: dictionary of data
    """
    with open(file_name, 'r') as fp:
        data = json.load(fp)
    all_data = []
    for key, values in data.items():
        all_data.append([key, values[0], values[1]])
    return all_data


def open_query(path_to_file):
    """
    Opens file with query genes
    Args:
        path_to_file (str): name of file

    Returns:
        file_names (str): list of genes

    """
    with open(path_to_file, 'r') as f:
        file_names = [line.strip() for line in f]
    return file_names


def find_correct_genome(gene, reference_data):
    """
    Finds genome in which gene is located
    Args:
        gene (str): gene name
        reference_data (list): list of data

    Returns:
        genome (str): genome name
    """
    for line in reference_data:
        genome = line[0]
        starts = [x for x in line[1::2]]
        ends = [x for x in line[2::2]]
        counter = 0
        for start in starts:
            end = ends[counter]
            counter += 1
            if int(gene) in range(int(start), int(end) + 1):
                return genome


def search_for_gene_in_genome(q_gene, genes, start_coords, end_coords, orients,
                              contig):
    """
    Searches for gene in data
    Args:
        q_gene (str): query gene
        genes (list): genes in genome
        start_coords (list): start coordinates of genes
        end_coords (list): end coordinates of genes
        orients (list): orientation of genes
        contig (list): contig of genes

    Returns:
        result (list): gene name
    """
    result = []
    main_index = genes.index(int(q_gene))
    result.append(start_coords[main_index])
    result.append(end_coords[main_index])
    result.append(orients[main_index])
    result.append(contig[main_index])
    result.append(q_gene)
    return result


def collect_gene_without_genome(gene):
    """
    Collects gene without genome
    Args:
        gene (str): gene name

    """
    with open('/usr/src/app/important_files/missing_genes/data_to_download.txt',
              'a+') as output:
        output.write(gene)
        output.write("\n")


class SuperSpeedAnalysisFromDomain:
    """
    This class is used to analyze neighborhood based on domain.
    """

    def __init__(self, user_pfam, user_distance, user_organisms, user_cutoff, user_correction,
                 user_strand, user_output, skip_negative):
        """ Initializing input data
            Args:
                user_pfam:  str (from pfam00000 to pfam99999)
                user_distance: int non negative integer 1-20000
                user_organisms: str
                user_cutoff: str (none, 0, 0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005,
                                 0.01, 0.05, 0.1, 0.2, 0.5)
                user_correction: str (none, bonferroni, fdr_bh)
                user_strand: str (both)
                user_output: str
         """
        self.user_pfam = user_pfam
        self.user_distance = user_distance
        self.user_organisms = user_organisms
        self.user_cutoff = user_cutoff
        self.user_correction = user_correction
        self.user_strand = user_strand
        self.user_output = user_output
        self.skip_negative = skip_negative

    def go(self):

        """
        Zipping all functions and execute them
        """
        # 1. Get all important data from files
        domain_information = collect_pfam_domain_description("/usr/src/app/important_files/domain_information")
        genome_id_size_in_gene = open_multiple_line("/usr/src/app/important_files/GENOME_ID_SIZE_IN_GENE.txt")
        with open('/usr/src/app/important_files/counter_domains', 'r') as handler:
            pfam_counter_domains = json.load(handler)
        file_names = open_database(self.user_organisms)

        genome_id = [x[0] for x in genome_id_size_in_gene]
        size_in_gene = [x[1] for x in genome_id_size_in_gene]

        # 2. Initialize all important variables
        genome_number_overall = 0
        genome_number_to_stat = 0
        just_neigh_data = []
        genomes_with_domain = []
        neigh_genome_size = []
        counter = 0
        counter_db = 0
        sum_of_genes_in_neighborhoods = 0
        number_of_neighborhoods = 0
        percentage_counter = Counter()

        # 3. Go through all files
        for file in file_names:
            counter_db += 1
            counter += 1
            # 3.1 Get data from single file
            try:
                genes, start_coord, end_coord, orientation, domains, contig = create_six_list(file)
            except FileNotFoundError:
                continue
            except ValueError:
                continue
            file_name_raw = file.split('/')
            tax_name = file_name_raw[-1]
            genome_number_overall += 1

            # 3.2 Get number of genes in genome
            number_of_genes_in_genome = genome_size_in_gene(tax_name, genome_id, size_in_gene)
            coords = searching_for_domain_in_genome(self.user_pfam, start_coord, end_coord, orientation,
                                                    domains, contig, genes)

            # 3.3 If there is anything to analyze add genome_number_to_stat
            if len(coords) > 0:
                genome_number_to_stat += 1

            # 3.4 Go through all founded domains in genome
            for point in coords:
                number_of_neighborhoods += 1

                # 3.4.1 Get data from neighborhood if distance = 0
                if self.user_distance != 0:

                    last_coordinate, first_coordinate, pfam_beg, pfam_end, searched_pfam_orientation = \
                        get_range_coordinates(point, self.user_distance)
                    # 3.4.2 Get data from neighborhood for both strands
                    if self.user_strand == 'both':
                        pfam_index_to_neigh = get_domains_both_strands(start_coord, end_coord, last_coordinate,
                                                                       first_coordinate, pfam_beg, pfam_end,
                                                                       contig,
                                                                       point)
                    # 3.4.3 Get data from neighborhood for positive strand and negative strand
                    else:
                        pfam_index_to_neigh_same, pfam_index_to_neigh_oposite = get_domains_plus_minus_strand(
                            start_coord, end_coord, last_coordinate,
                            first_coordinate, pfam_beg, pfam_end,
                            contig, point, orientation)
                        if self.user_strand == 'same':
                            pfam_index_to_neigh = pfam_index_to_neigh_same
                        elif self.user_strand == 'opposite':
                            pfam_index_to_neigh = pfam_index_to_neigh_oposite
                else:
                    pfam_index_to_neigh = get_domain_both_zeros(genes, point[4])
                # 3.5 Collect data for each domain in one genome
                genes_in_neigh = list_of_genes_in_neigh(genes, pfam_index_to_neigh)
                number_of_genes_in_neigh = size_in_genes(genes_in_neigh)
                whole = get_list_of_domain_in_neigh(pfam_index_to_neigh, file, domains)
                percentage_counter += Counter(set(whole))
                just_neigh_data.append(whole)
                genomes_with_domain.append(file)
                neigh_genome_size.append((number_of_genes_in_neigh, number_of_genes_in_genome))
                sum_of_genes_in_neighborhoods += int(number_of_genes_in_neigh)

        # 4. Gather all neighborhood data
        alls = []
        neigh_counter = Counter()

        for i in just_neigh_data:
            for j in i:
                if j not in alls:
                    alls.append(j)
        genome_counter = Counter()

        # 5. Create matrix of data to calculate
        counter = 0
        for domain_list, genome_in_domain, ng_size in zip(just_neigh_data, genomes_with_domain, neigh_genome_size):
            counter += 1
            genes, start_coord, end_coord, orientation, domains, contig = create_six_list(genome_in_domain)
            genome_domains_counter = Counter(domains)
            genome_counter += genome_domains_counter
            neigh_domains_counter = Counter(domain_list)
            neigh_counter += neigh_domains_counter

        # 6. Create pandas dataframe for storing results
        najlepsze_dane = pd.DataFrame(columns=['occurence in neighborhoods',
                                               'average occurence in neighborhood', 'occurence genomes',
                                               'average occurence in genome',
                                               'In what percentage', 'Poisson', 'Poisson_correct', 'P', 'N', 'K',
                                               'avg neigh size', 'domain in whole database',
                                               'Family', 'Summary', ])

        # 7. Calculate all domains
        scores_poisson = []
        counter = 0
        for i in alls:
            counter += 1

            p_KP = ((sum_of_genes_in_neighborhoods / genome_number_to_stat) * pfam_counter_domains[i]) / 326628067

            najlepsze_dane.at[i, 'P'] = p_KP
            najlepsze_dane.at[i, 'avg neigh size'] = sum_of_genes_in_neighborhoods / genome_number_to_stat
            najlepsze_dane.at[i, 'domain in whole database'] = pfam_counter_domains[i]
            poisson_test_result = poisson_distribution(p_KP, number_of_neighborhoods, i, percentage_counter)
            najlepsze_dane.at[i, 'K'] = int(percentage_counter[i])
            najlepsze_dane.at[i, 'N'] = number_of_neighborhoods
            najlepsze_dane.at[i, 'P'] = p_KP
            najlepsze_dane.at[i, 'Poisson'] = poisson_test_result
            scores_poisson.append(poisson_test_result)

        poisson_after_correction = multiple_test_correction_list(scores_poisson, self.user_correction)

        # 8. Combine all results in one dataframe
        for domena, poisson_corr in zip(alls, poisson_after_correction):
            najlepsze_dane.at[domena, 'Poisson_correct'] = poisson_corr
            najlepsze_dane.at[domena, 'In what percentage'] = percentage_counter[domena] / len(just_neigh_data)
            najlepsze_dane.at[domena, 'occurence in neighborhoods'] = neigh_counter[domena]
            najlepsze_dane.at[domena, 'average occurence in neighborhood'] = neigh_counter[domena] / len(
                just_neigh_data)
            najlepsze_dane.at[domena, 'occurence genomes'] = genome_counter[domena]
            najlepsze_dane.at[domena, 'average occurence in genome'] = genome_counter[domena] / len(
                set(genomes_with_domain))

        # 9. Customize dataframe
        after_cutoff = cutoff_value(najlepsze_dane, self.user_cutoff)
        sorted_table = sort_table(after_cutoff)
        added_information = add_information(sorted_table, domain_information)
        final_data = pfam_for_pf(added_information, self.user_correction)

        # 10. Save results
        match_go_to_pfam_and_save(final_data, self.user_output)
        save_data(self.user_output, final_data)


class SuperSpeedAnalysisFromDomainAll:
    """
    Class for analysis of all domains in all database
    """

    def __init__(self, user_pfam, user_distance, user_organism, user_cutoff, user_correction,
                 user_strand, user_output, user_level, skip_negative):
        """ Initializing input data
            Args:
                user_pfam:  str (from pfam00000 to pfam99999)
                user_distance: str non negative integer 1-20000
                user_level: str
                user_cutoff: str (none, 0, 0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005,
                                 0.01, 0.05, 0.1, 0.2, 0.5)
                user_correction: str (none, bonferroni, fdr_bh)
                user_strand: str (both)
                        user_output: str
         """
        self.user_pfam = user_pfam  # pfam02696
        self.user_distance = user_distance  # 5000
        self.user_organism = user_organism  #
        self.user_cutoff = user_cutoff  # none
        self.user_correction = user_correction  # none
        self.user_strand = user_strand  # both
        self.user_output = user_output  # place
        self.user_level = user_level  # order
        self.skip_negative = skip_negative

    def open_database(self, level, organism):
        """
        Function for opening database and creating list of genomes
        Args:
            level (str): level of taxonomy
            organism (str): name of organism

        Returns:
            taxonomy (dict): dictionary with taxonomy
        """
        if level == "genus":
            with open('/usr/src/app/important_files/genus_level.json', "r") as handler:
                taxonomy = json.load(handler)
        elif level == "family":
            with open('/usr/src/app/important_files/family_level.json', "r") as handler:
                taxonomy = json.load(handler)
        elif level == "order":
            with open('/usr/src/app/important_files/order_level.json', "r") as handler:
                taxonomy = json.load(handler)
        elif level == "class":
            with open('/usr/src/app/important_files/class_level.json', "r") as handler:
                taxonomy = json.load(handler)
        elif level == "phylum":
            with open('/usr/src/app/important_files/phylum_level.json', "r") as handler:
                taxonomy = json.load(handler)

        return taxonomy

    def cutoff_value(self, some_dataframe, cutoff):
        """
        Function for cutting off values
        Args:
            some_dataframe (pd.Dataframe): dataframe with results
            cutoff (str): cutoff value

        Returns:
            some_dataframe (pd.Dataframe): dataframe with results
        """
        domain_list = list(some_dataframe.index)
        if cutoff == 'none':
            return some_dataframe
        elif cutoff == '0':
            for i in domain_list:

                diff = float(some_dataframe.loc[i, 'PVALUE'])
                if diff < 0 and diff > 0:
                    some_dataframe = some_dataframe.drop([i])
            return some_dataframe
        else:
            cutoff = float(cutoff)
            for i in domain_list:
                diff = float(some_dataframe.loc[i, 'PVALUE'])
                if diff > cutoff:
                    some_dataframe = some_dataframe.drop([i])
            return some_dataframe

    def pfam_for_pf(self, dataframe):
        """Final customization of dataframe
            --> PF02696 instead of pfam02696
            --> PVALUE in 10e-3 format
            --> In what percentage % format
            --> avg occurences and density  have now 3 digits after coma
            --> drop NOam
        """

        indexy = dataframe.index
        for i in indexy:
            dataframe = dataframe.rename(index={i: i[:2].upper() + i[4:]})
        dataframe['Poisson'] = dataframe['Poisson'].map('{:.3e}'.format)
        dataframe['Poisson_correct'] = dataframe['Poisson_correct'].map('{:.3e}'.format)

        dataframe['In what percentage'] = dataframe['In what percentage'].map('{:.3%}'.format)
        dataframe['average occurence in neighborhood'] = dataframe['average occurence in neighborhood'].map(
            '{:.3}'.format)
        dataframe['average occurence in genome'] = dataframe['average occurence in genome'].map('{:.3}'.format)

        # dataframe = dataframe.round({"average occurence in neighborhood": 3,"average occurence in genome": 3,
        #                              "Density difference": 3})
        if "NOam" in indexy:
            dataframe = dataframe.drop(index="NOam")
        # dataframe = dataframe.style.format({'In what percentage': '{:,.3%}'.format})
        dataframe.drop('average occurence in neighborhood', axis=1, inplace=True)
        dataframe.drop('average occurence in genome', axis=1, inplace=True)
        dataframe.drop('avg neigh size', axis=1, inplace=True)
        dataframe.drop('domain in whole database', axis=1, inplace=True)
        if self.user_correction == 'none':
            dataframe.drop('Poisson_correct', axis=1, inplace=True)
        else:
            dataframe.drop('Poisson', axis=1, inplace=True)
        dataframe = dataframe.drop(index="NOam")
        return dataframe

    def poisson_distribution(self, p, n, domain, percentage_counter):
        """
        Function for calculating poisson distribution
        Args:
            p (int): number of domains in neighborhood
            n (int): number of domains in genome
            domain (str): name of domain
            percentage_counter (dict): dictionary with percentage of domains in genome

        Returns:
            poisson (float): poisson distribution
        """
        if int(percentage_counter[domain]) < p * n:
            return 1
        else:
            sums = 0
            k = 0
            number_of_samples = int(percentage_counter[domain])
            while k < number_of_samples:
                sums += poisson.pmf(k=k, mu=n * p)
                k += 1
            if sums >= 1:
                return float(0.0)
            else:
                return 1 - sums

    def go(self):

        """
        Zipping all functions and execute them
        """
        # 1. open import files
        domain_information = collect_pfam_domain_description("/usr/src/app/important_files/domain_information")
        genome_id_size_in_gene = open_multiple_line("/usr/src/app/important_files/GENOME_ID_SIZE_IN_GENE.txt")
        with open('/usr/src/app/important_files/counter_domains', 'r') as handler:
            pfam_counter_domains = json.load(handler)

        genome_id = [x[0] for x in genome_id_size_in_gene]
        size_in_gene = [x[1] for x in genome_id_size_in_gene]
        # 2. create variables
        genome_number_overall = 0
        genome_number_to_stat = 0

        # 3. Get all genomes from database
        level = self.open_database(self.user_level, self.user_organism)

        # 4. Go through all taxonomy levels.
        for key in level:
            just_neigh_data = []
            genomes_with_domain = []
            neigh_genome_size = []
            counter = 0
            counter_db = 0
            sum_of_genes_in_neighborhoods = 0
            number_of_neighborhoods = 0
            percentage_counter = Counter()

            # 5. Go through all genomes in taxonomy level
            for file in level[key]:
                counter_db += 1
                counter += 1
                try:
                    genes, start_coord, end_coord, orientation, domains, contig = create_six_list(file)
                except FileNotFoundError:
                    continue
                except ValueError:
                    continue
                file_name_raw = file.split('/')
                tax_name = file_name_raw[-1]
                genome_number_overall += 1
                number_of_genes_in_genome = genome_size_in_gene(tax_name, genome_id, size_in_gene)
                coords = searching_for_domain_in_genome(self.user_pfam, start_coord, end_coord, orientation,
                                                        domains, contig, genes)
                if len(coords) > 0:
                    genome_number_to_stat += 1

                for point in coords:
                    number_of_neighborhoods += 1

                    if self.user_distance != 0:
                        last_coordinate, first_coordinate, pfam_beg, pfam_end, searched_pfam_orientation = \
                            get_range_coordinates(point, self.user_distance)
                        pfam_index_to_neigh = get_domains_both_strands(start_coord, end_coord, last_coordinate,
                                                                       first_coordinate, pfam_beg, pfam_end,
                                                                       contig,
                                                                       point)
                    else:
                        pfam_index_to_neigh = get_domain_both_zeros(genes, point[4])

                    genes_in_neigh = list_of_genes_in_neigh(genes, pfam_index_to_neigh)
                    number_of_genes_in_neigh = size_in_genes(genes_in_neigh)
                    whole = get_list_of_domain_in_neigh(pfam_index_to_neigh, file, domains)
                    percentage_counter += Counter(set(whole))
                    just_neigh_data.append(whole)
                    genomes_with_domain.append(file)
                    neigh_genome_size.append((number_of_genes_in_neigh, number_of_genes_in_genome))
                    sum_of_genes_in_neighborhoods += int(number_of_genes_in_neigh)

            # 6. Gather all data from all genomes in taxonomy level
            alls = []
            neigh_counter = Counter()

            for i in just_neigh_data:
                for j in i:
                    if j not in alls:
                        alls.append(j)
            genome_counter = Counter()

            # 7. Create matrix of data to calculate statistics
            counter = 0
            for domain_list, genome_in_domain, ng_size in zip(just_neigh_data, genomes_with_domain, neigh_genome_size):
                counter += 1
                genes, start_coord, end_coord, orientation, domains, contig = create_six_list(genome_in_domain)
                genome_domains_counter = Counter(domains)
                genome_counter += genome_domains_counter
                neigh_domains_counter = Counter(domain_list)
                neigh_counter += neigh_domains_counter

            # 8. Create pandas dataframe for storring results
            najlepsze_dane = pd.DataFrame(columns=['occurence in neighborhoods',
                                                   'average occurence in neighborhood', 'occurence genomes',
                                                   'average occurence in genome',
                                                   'In what percentage', 'Poisson', 'Poisson_correct', 'P', 'N', 'K',
                                                   'avg neigh size', 'domain in whole database',
                                                   'Family', 'Summary', ])

            # 9. Calculate statistics
            scores_poisson = []
            counter = 0
            for i in alls:
                counter += 1

                p_KP = ((sum_of_genes_in_neighborhoods / genome_number_to_stat) * pfam_counter_domains[i]) / 326628067
                najlepsze_dane.at[i, 'P'] = p_KP
                najlepsze_dane.at[i, 'avg neigh size'] = sum_of_genes_in_neighborhoods / genome_number_to_stat
                najlepsze_dane.at[i, 'domain in whole database'] = pfam_counter_domains[i]
                # najlepsze_dane.at[i,'sum of database len'] = 326628067
                poisson_test_result = self.poisson_distribution(p_KP, number_of_neighborhoods, i, percentage_counter)
                najlepsze_dane.at[i, 'K'] = int(percentage_counter[i])
                najlepsze_dane.at[i, 'N'] = number_of_neighborhoods
                najlepsze_dane.at[i, 'P'] = p_KP
                najlepsze_dane.at[i, 'Poisson'] = poisson_test_result
                scores_poisson.append(poisson_test_result)

            poisson_after_correction = multiple_test_correction_list(scores_poisson, self.user_correction)

            # 10. Collect data for database
            for domena, poisson_corr, in zip(alls, poisson_after_correction):
                najlepsze_dane.at[domena, 'Poisson_correct'] = poisson_corr
                najlepsze_dane.at[domena, 'In what percentage'] = percentage_counter[domena] / len(just_neigh_data)
                najlepsze_dane.at[domena, 'occurence in neighborhoods'] = neigh_counter[domena]
                najlepsze_dane.at[domena, 'average occurence in neighborhood'] = neigh_counter[domena] / len(
                    just_neigh_data)
                najlepsze_dane.at[domena, 'occurence genomes'] = genome_counter[domena]
                najlepsze_dane.at[domena, 'average occurence in genome'] = genome_counter[domena] / len(
                    set(genomes_with_domain))
            # 11. Customize dataframe

            after_cutoff = self.cutoff_value(najlepsze_dane, self.user_cutoff)
            sorted_table = sort_table(after_cutoff)
            added_information = add_information(sorted_table, domain_information)
            final_data = self.pfam_for_pf(added_information)

            message_down = " In selected database there is {} genomes and in {} searched domain " \
                           "was found".format(str(counter_db), str(genome_number_to_stat))
            match_go_to_pfam_and_save(final_data, self.user_output)
            save_data('{}_{}'.format(self.user_level, key), final_data)


class SuperSpeedAnalysisFromDomainAllAvg:

    def __init__(self, user_pfam, user_distance, user_organism, user_cutoff, user_correction,
                 user_strand, user_output, user_level):
        """ Initializing input data
            Inputs:
                    user_pfam:  str (from pfam00000 to pfam99999)
                    user_distance: str non negative integer 1-20000
                    user_level: str
                    user_cutoff: str (none, 0, 0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005,
                                 0.01, 0.05, 0.1, 0.2, 0.5)
                    user_correction: str (none, bonferroni, fdr_bh)
                    user_strand: str (both)
                    user_output: str
         """
        self.user_pfam = user_pfam  # pfam02696
        self.user_distance = user_distance  # 5000
        self.user_organism = user_organism  # none
        self.user_cutoff = user_cutoff  # none
        self.user_correction = user_correction  # none
        self.user_strand = user_strand  # both
        self.user_output = user_output  # place
        self.user_level = user_level  # order

    def save_data(self, file_name, complete_data,):
        """
        Saves all stuff together as one file
        Args:
            file_name: str
            complete_data: str
        """

        complete_data.index.name = 'Domain'
        complete_data.to_csv('/usr/src/app/media/results/temp/' + file_name, sep='\t', mode='w')

    def open_database(self, level, organism):
        """
        Function for opening database and creating dictionary of genomes
        Args:
            level (str): level of taxonomy
            organism (str): name of organism

        Returns:
            taxonomy (dict): dictionary with taxonomy
        """
        tax = organism.replace("_", " ")
        if level == "genus":
            with open('/usr/src/app/important_files/genus_level.json', "r") as handler:
                taxonomy = json.load(handler)
        elif level == "family":
            with open('/usr/src/app/important_files/family_level.json', "r") as handler:
                taxonomy = json.load(handler)
        elif level == "order":
            with open('/usr/src/app/important_files/order_level.json', "r") as handler:
                taxonomy = json.load(handler)
        elif level == "class":
            with open('/usr/src/app/important_files/class_level.json', "r") as handler:
                taxonomy = json.load(handler)
        elif level == "phylum":
            with open('/usr/src/app/important_files/phylum_level.json', "r") as handler:
                taxonomy = json.load(handler)

        return taxonomy

    def cutoff_value(self, some_dataframe, cutoff):
        """
        Function for cutting off values
        Args:
            some_dataframe (pd.Dataframe): dataframe with results
            cutoff (str): cutoff value

        Returns:
            some_dataframe (pd.Dataframe): dataframe with results
        """
        domain_list = list(some_dataframe.index)
        if cutoff == 'none':
            return some_dataframe
        elif cutoff == '0':
            for i in domain_list:

                diff = float(some_dataframe.loc[i, 'PVALUE'])
                if diff < 0 and diff > 0:
                    some_dataframe = some_dataframe.drop([i])
            return some_dataframe
        else:
            cutoff = float(cutoff)
            for i in domain_list:
                diff = float(some_dataframe.loc[i, 'PVALUE'])
                if diff > cutoff:
                    some_dataframe = some_dataframe.drop([i])
            return some_dataframe

    def pfam_for_pf(self, dataframe):
        """Final customization of dataframe
            --> PF02696 instead of pfam02696
            --> PVALUE in 10e-3 format
            --> In what percentage % format
            --> avg occurences and density  have now 3 digits after coma
            --> drop NOam
        """

        indexy = dataframe.index
        for i in indexy:
            dataframe = dataframe.rename(index={i: i[:2].upper() + i[4:]})

        # if "NOam" in indexy:
        # dataframe = dataframe.drop(index="NOam")
        # dataframe = dataframe.style.format({'In what percentage': '{:,.3%}'.format})
        try:
            dataframe = dataframe.drop(index="NOam")
        except KeyError:
            return dataframe
        return dataframe

    def poisson_distribution(self, p, n, K):
        """
        Function for calculating poisson distribution
        Args:
            p (int): number of genes in genome
            n (int): number of genes in domain
            K (int): number of genes in domain in genome

        Returns:
            poisson (float): poisson distribution
        """
        if int(K) < p * n:
            return 1
        else:
            sums = 0
            k = 0
            number_of_samples = int(K)
            while k < number_of_samples:
                sums += poisson.pmf(k=k, mu=n * p)
                k += 1
            if sums >= 1:
                return float(0.0)
            else:
                return 1 - sums

    def magic_number_creator(self, p0, p1, wie):
        """
        Function for calculating magic number based on values distribution
        Args:
            p0 (float): first number
            p1 (float): second number
            wie (float): how many numbers

        Returns:
            magic_number (float): magic number
        """
        import numpy as np

        np.random.choice([0, 1], p=[p0, p1])

        result = []
        for i in range(wie):
            result.append(np.random.choice([0, 1], p=[p0, p1]))
        return result

    def analysis(self, groups):
        """
        Function for analysis of data for averaging
        Args:
            groups (list): list of groups to analyze

        Returns:
            results (pd.Dataframe): dataframe with results
            all_dictionary (dict): dictionary with all results
        """
        max_neigh = 0

        for org in groups:
            data = pd.read_csv('/usr/src/app/media/results/temp/' + org, sep='\t', index_col=0)

            domain_query = 'PF' + self.user_pfam[4:]
            max_neigh = max(max_neigh, data.at[domain_query, 'N'])

        all_dictionary = {}
        all_domains = []
        group = {}
        neigh_c = {}
        genome_c = {}

        for org in groups:

            data = pd.read_csv('/usr/src/app/media/results/temp/' + org, sep='\t', index_col=0)

            domains = data.index.tolist()

            neigh_counter = list(data['occurence in neighborhoods'])
            genome_counter = list(data['occurence genomes'])
            k = list(data['K'])
            n = list(data['N'])
            p = list(data['P'])

            # domains
            for domain in domains:
                if domain not in group.keys():
                    group[domain] = 0
                group[domain] += 1
            # K N P
            for i in domains:
                if i not in all_domains:
                    all_domains.append(i)
                    all_dictionary[i] = {'k': [], 'n': [], 'p': [], 'neigh': [], 'genome': []}

            for domain, _k, _n, _p, neigh, genome in zip(domains, k, n, p, neigh_counter, genome_counter):
                all_dictionary[domain]['k'].append(_k)
                all_dictionary[domain]['p'].append(_p)
                all_dictionary[domain]['n'].append(_n)
                all_dictionary[domain]['neigh'].append(neigh)
                all_dictionary[domain]['genome'].append(genome)

        new_dictionary = {}
        for domain in list(all_dictionary.keys()):
            print(domain)
            small_data = all_dictionary[domain]
            new_dictionary[domain] = {'k': [], 'n': [], 'p': [], 'neigh': [], 'genome': []}

            for _k, _n, _p in zip(small_data['k'], small_data['n'], small_data['p']):
                new_dictionary[domain]['p'].append(float(_p))
                rate_of_one = int(_k) / int(_n)
                missing_values = max_neigh - int(_n)
                magic_result = self.magic_number_creator(1 - rate_of_one, rate_of_one, missing_values)
                new_dictionary[domain]['k'].append(int(_k) + sum(magic_result))
                new_dictionary[domain]['n'].append(int(_n) + len(magic_result))
        results = {}
        for domain in list(new_dictionary.keys()):
            k = sum(new_dictionary[domain]['k'])
            n = sum(new_dictionary[domain]['n'])
            p = sum(new_dictionary[domain]['p']) / len(new_dictionary[domain]['p'])
            print(domain)
            results[domain] = self.poisson_distribution(p, n, k)

        return results, all_dictionary

    def merge_all_data_and_save(self, results, all_dict, size_of_group, domain_information,
                                save, overall_num_of_genomes, overall_num_of_neigh):
        """
        Function for merging all data and saving it to file
        Args:
            results (dict): dictionary with results
            all_dict (dict): dictionary with all data
            size_of_group (int): size of group
            domain_information (dict): dictionary with domain information
            save (str): path to save file
            overall_num_of_genomes (int): number of genomes
            overall_num_of_neigh (int): number of neighborhoods

        Returns:
            dataframe (pd.Dataframe): dataframe with results that wa saved
        """

        complete = pd.DataFrame(columns=['Occurence in genome', 'Avg occurence in genome', 'Occurence in neighborhood',
                                         'Avg occurence in neigh', 'Family percentage',
                                         'Poisson', 'Summary', 'Family'])

        domains = list(results.keys())
        values_to_correct = list(results.values())
        if self.user_correction != 'none':
            reject, pvals_corrected, alphaSidak, alphaBonf = correction(pvals=values_to_correct,
                                                                        method=self.user_correction,
                                                                        is_sorted=False, returnsorted=False)
            pvals = pvals_corrected.tolist()
            for domain, res in zip(domains, pvals):
                complete.at[domain, 'Poisson'] = res
        else:
            for domain, res in zip(domains, values_to_correct):
                complete.at[domain, 'Poisson'] = res

        for domain in domains:
            complete.at[domain, 'Occurence in genome'] = sum(all_dict[domain]['genome'])
            complete.at[domain, 'Avg occurence in genome'] = sum(all_dict[domain]['genome'])
            complete.at[domain, 'Occurence in neighborhood'] = sum(all_dict[domain]['neigh']) / overall_num_of_neigh
            complete.at[domain, 'Avg occurence in neigh'] = sum(all_dict[domain]['neigh']) / overall_num_of_genomes
            complete.at[domain, 'Family percentage'] = len(all_dict[domain]['neigh']) / size_of_group
        dataframe = add_information(complete, domain_information)
        dataframe.index.name = 'Domain'
        dataframe.to_csv('/usr/src/app/media/results/' + save, sep='\t', mode='w')
        return dataframe

    def go(self):
        """
        Zipping all functions and execute them
        """
        # 1. Open correct important files
        groups_calculated = []
        if self.user_level == "genus":
            with open('/usr/src/app/important_files/genus_level.json', "r") as handler:
                taxonomy = json.load(handler)

        elif self.user_level == "family":
            with open('/usr/src/app/important_files/family_level.json', "r") as handler:
                taxonomy = json.load(handler)

        elif self.user_level == "order":
            with open('/usr/src/app/important_files/order_level.json', "r") as handler:
                taxonomy = json.load(handler)

        elif self.user_level == "class":
            with open('/usr/src/app/important_files/class_level.json', "r") as handler:
                taxonomy = json.load(handler)

        elif self.user_level == "phylum":
            with open('/usr/src/app/important_files/phylum_level.json', "r") as handler:
                taxonomy = json.load(handler)

        # genus_database = list(taxonomy.keys())

        # 2. open important files
        domain_information = collect_pfam_domain_description(
            "/usr/src/app/important_files/domain_information")
        genome_id_size_in_gene = open_multiple_line(
            "/usr/src/app/important_files/GENOME_ID_SIZE_IN_GENE.txt")
        with open('/usr/src/app/important_files/counter_domains', 'r') as handler:
            pfam_counter_domains = json.load(handler)

        genome_id = [x[0] for x in genome_id_size_in_gene]
        size_in_gene = [x[1] for x in genome_id_size_in_gene]

        # 3. create variables
        genome_number_overall = 0
        genome_number_to_stat = 0
        overall_num_of_genomes = 0
        overall_num_of_neigh = 0
        level = self.open_database(self.user_level, self.user_organism)

        # 4. Go through all taxonomy levels.
        for key in list(level.keys()):
            just_neigh_data = []
            genomes_with_domain = []
            neigh_genome_size = []
            counter = 0
            counter_db = 0
            sum_of_genes_in_neighborhoods = 0
            number_of_neighborhoods = 0
            percentage_counter = Counter()
            if len(level[key]) > 25:
            # 5. Go through all genomes in taxonomy level
                for file in level[key]:
                    counter_db += 1
                    counter += 1
                    try:
                        genes, start_coord, end_coord, orientation, domains, contig = create_six_list(file)
                    except FileNotFoundError:
                        continue
                    except ValueError:
                        continue
                    file_name_raw = file.split('/')
                    tax_name = file_name_raw[-1]
                    genome_number_overall += 1
                    number_of_genes_in_genome = genome_size_in_gene(tax_name, genome_id, size_in_gene)
                    coords = searching_for_domain_in_genome(self.user_pfam, start_coord, end_coord, orientation,
                                                            domains, contig, genes)
                    if len(coords) > 0:
                        genome_number_to_stat += 1
                        overall_num_of_genomes += 1
                    for point in coords:
                        number_of_neighborhoods += 1
                        overall_num_of_neigh += 1

                        if self.user_distance != 0:
                            last_coordinate, first_coordinate, pfam_beg, pfam_end, searched_pfam_orientation = \
                                get_range_coordinates(point, self.user_distance)
                            pfam_index_to_neigh = get_domains_both_strands(start_coord, end_coord, last_coordinate,
                                                                           first_coordinate, pfam_beg, pfam_end,
                                                                           contig,
                                                                           point)
                        else:
                            pfam_index_to_neigh = get_domain_both_zeros(genes, point[4])

                        genes_in_neigh = list_of_genes_in_neigh(genes, pfam_index_to_neigh)
                        number_of_genes_in_neigh = size_in_genes(genes_in_neigh)
                        whole = get_list_of_domain_in_neigh(pfam_index_to_neigh, file, domains)
                        percentage_counter += Counter(set(whole))
                        just_neigh_data.append(whole)
                        genomes_with_domain.append(file)
                        neigh_genome_size.append((number_of_genes_in_neigh, number_of_genes_in_genome))
                        sum_of_genes_in_neighborhoods += int(number_of_genes_in_neigh)

                # 6. Gather all data from all genomes in taxonomy level
                alls = []
                neigh_counter = Counter()

                for i in just_neigh_data:
                    for j in i:
                        if j not in alls:
                            alls.append(j)
                genome_counter = Counter()

                # 7. Create matrix of data to calculate statistics
                counter = 0
                for domain_list, genome_in_domain, ng_size in zip(just_neigh_data, genomes_with_domain, neigh_genome_size):
                    counter += 1
                    genes, start_coord, end_coord, orientation, domains, contig = create_six_list(genome_in_domain)
                    genome_domains_counter = Counter(domains)
                    genome_counter += genome_domains_counter
                    neigh_domains_counter = Counter(domain_list)
                    neigh_counter += neigh_domains_counter

                # 8. Create pandas dataframe for storring results
                najlepsze_dane = pd.DataFrame(columns=['occurence in neighborhoods',
                                                       'occurence genomes',
                                                       'P', 'N', 'K',
                                                       'domain in whole database', 'In what percentage'])
                # 9. Calculate statistics
                counter = 0
                for i in alls:
                    counter += 1
                    p_KP = ((sum_of_genes_in_neighborhoods / genome_number_to_stat) * pfam_counter_domains[i]) / 326628067
                    najlepsze_dane.at[i, 'P'] = p_KP
                    najlepsze_dane.at[i, 'domain in whole database'] = pfam_counter_domains[i]
                    najlepsze_dane.at[i, 'K'] = int(percentage_counter[i])
                    najlepsze_dane.at[i, 'N'] = number_of_neighborhoods
                    najlepsze_dane.at[i, 'occurence in neighborhoods'] = neigh_counter[i]
                    najlepsze_dane.at[i, 'occurence genomes'] = genome_counter[i]
                    najlepsze_dane.at[i, 'In what percentage'] = neigh_counter[i] / number_of_neighborhoods

                added_information = add_information(najlepsze_dane, domain_information)
                final_data = self.pfam_for_pf(added_information)
                if len(najlepsze_dane) > 0:
                    groups_calculated.append(key)
                    self.save_data(key, final_data)

        complete_results, additional = self.analysis(groups_calculated)

        dataframe = self.merge_all_data_and_save(complete_results, additional, len(list(level.keys())),
                                                 domain_information, self.user_output, overall_num_of_genomes,
                                                 overall_num_of_neigh)
        if self.user_level == 'yes':
            match_go_to_pfam_and_save(dataframe, self.user_output)
        else:
            match_go_to_pfam_and_save_average(dataframe, self.user_output)


class SuperSpeedAnalysisFromDomainFamilyAvg:
    """
    Class for analysis of domain families using averageing method
    """

    def __init__(self, user_pfam, user_distance, user_organism, user_cutoff, user_correction,
                 user_strand, user_output):

        """ Initializing input data
            Args:
                user_pfam:  str (from pfam00000 to pfam99999)
                user_distance: str non negative integer 1-20000
                user_level: str
                user_cutoff: str (none, 0, 0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005,
                                 0.01, 0.05, 0.1, 0.2, 0.5)
                user_correction: str (none, bonferroni, fdr_bh)
                user_strand: str (both)
                user_output: str
         """

        self.user_pfam = user_pfam
        self.user_distance = user_distance
        self.user_organism = user_organism
        self.user_cutoff = user_cutoff
        self.user_correction = user_correction
        self.user_strand = user_strand
        self.user_output = user_output

    def save_data(self, file_name, complete_data,):
        """
        Saves all stuff together as one file
        Args:
            file_name (str): name of file
            complete_data (pd.DataFrame): data to save
        """

        complete_data.index.name = 'Domain'
        complete_data.to_csv('/usr/src/app/media/results/temp/' + file_name, sep='\t', mode='w')

    def open_database(self, family):
        """
        Opens database with domain families
        Args:
            family (str): name of family

        Returns:
            family_data (pd.DataFrame): data from database
        """

        # json_family = json.load(open('/home/djangoadmin/final_site-project/important_files/database_file/family.json'))
        json_family = json.load(open('/usr/src/app/important_files/families.json'))
        family_data = json_family[family]
        return family_data

    def pfam_for_pf(self, dataframe):
        """Final customization of dataframe
            --> PF02696 instead of pfam02696
            --> PVALUE in 10e-3 format
            --> In what percentage % format
            --> avg occurences and density  have now 3 digits after coma
            --> drop NOam
        """

        indexy = dataframe.index
        for i in indexy:
            dataframe = dataframe.rename(index={i: i[:2].upper() + i[4:]})

        if "NOam" in indexy:
            dataframe = dataframe.drop(index="NOam")
        # dataframe = dataframe.style.format({'In what percentage': '{:,.3%}'.format})

        return dataframe

    def poisson_distribution(self, p, n, K):
        """
        Calculates poisson distribution
        Args:
            p (int): number of genes in neighborhood
            n (int): number of genes in genome
            K (int): number of genes in domain family

        Returns:
            poisson (float): poisson distribution
        """
        if int(K) < p * n:
            return 1
        else:
            sums = 0
            k = 0
            number_of_samples = int(K)
            while k < number_of_samples:
                sums += poisson.pmf(k=k, mu=n * p)
                k += 1
            if sums >= 1:
                return float(0.0)
            else:
                return 1 - sums

    def magic_number_creator(self, p0, p1, wie):
        """
        Function for calculating magic number based on values distribution
        Args:
            p0 (float): first number
            p1 (float): second number
            wie (float): how many numbers

        Returns:
            magic_number (float): magic number
        """
        import numpy as np

        np.random.choice([0, 1], p=[p0, p1])

        result = []
        for i in range(wie):
            result.append(np.random.choice([0, 1], p=[p0, p1]))
        return result

    def analysis(self, groups):
        """
        Function for analysis of data for averaging
        Args:
            groups (list): list of groups to analyze

        Returns:
            results (pd.Dataframe): dataframe with results
            all_dictionary (dict): dictionary with all results
        """
        max_neigh = 0

        for org in groups:
            data = pd.read_csv('temp/' + org, sep='\t', index_col=0)

            domain_query = 'PF' + self.user_pfam[4:]
            max_neigh = max(max_neigh, data.at[domain_query, 'N'])

        all_dictionary = {}
        all_domains = []
        group = {}
        neigh_c = {}
        genome_c = {}

        for org in groups:

            data = pd.read_csv('temp/' + org, sep='\t', index_col=0)

            domains = data.index.tolist()

            neigh_counter = list(data['occurence in neighborhoods'])
            genome_counter = list(data['occurence genomes'])
            k = list(data['K'])
            n = list(data['N'])
            p = list(data['P'])

            # domains
            for domain in domains:
                if domain not in group.keys():
                    group[domain] = 0
                group[domain] += 1
            # K N P
            for i in domains:
                if i not in all_domains:
                    all_domains.append(i)
                    all_dictionary[i] = {'k': [], 'n': [], 'p': [], 'neigh': [], 'genome': []}

            for domain, _k, _n, _p, neigh, genome in zip(domains, k, n, p, neigh_counter, genome_counter):
                all_dictionary[domain]['k'].append(_k)
                all_dictionary[domain]['p'].append(_p)
                all_dictionary[domain]['n'].append(_n)
                all_dictionary[domain]['neigh'].append(neigh)
                all_dictionary[domain]['genome'].append(genome)

        new_dictionary = {}
        for domain in list(all_dictionary.keys()):
            print(domain)
            small_data = all_dictionary[domain]
            new_dictionary[domain] = {'k': [], 'n': [], 'p': [], 'neigh': [], 'genome': []}

            for _k, _n, _p in zip(small_data['k'], small_data['n'], small_data['p']):
                new_dictionary[domain]['p'].append(float(_p))
                rate_of_one = int(_k) / int(_n)
                missing_values = max_neigh - int(_n)
                magic_result = self.magic_number_creator(1 - rate_of_one, rate_of_one, missing_values)
                new_dictionary[domain]['k'].append(int(_k) + sum(magic_result))
                new_dictionary[domain]['n'].append(int(_n) + len(magic_result))
        results = {}
        for domain in list(new_dictionary.keys()):
            k = sum(new_dictionary[domain]['k'])
            n = sum(new_dictionary[domain]['n'])
            p = sum(new_dictionary[domain]['p']) / len(new_dictionary[domain]['p'])
            print(domain)
            results[domain] = self.poisson_distribution(p, n, k)

        return results, all_dictionary

    def merge_all_data_and_save(self, results, all_dict, size_of_group, domain_information,
                                save, overall_num_of_genomes, overall_num_of_neigh):
        """
        Function for merging all data and saving it to file
        Args:
            results (dict): dictionary with results
            all_dict (dict): dictionary with all data
            size_of_group (int): size of group
            domain_information (dict): dictionary with domain information
            save (str): path to save file
            overall_num_of_genomes (int): number of genomes
            overall_num_of_neigh (int): number of neighborhoods

        Returns:
            dataframe (pd.Dataframe): dataframe with results that wa saved
        """
        complete = pd.DataFrame(columns=['Occurence in genome', 'Avg occurence in genome', 'Occurence in neighborhood',
                                         'Avg occurence in neigh', 'Family percentage',
                                         'Poisson', 'Summary', 'Family'])

        domains = list(results.keys())
        values_to_correct = list(results.values())
        if self.user_correction != 'none':
            reject, pvals_corrected, alphaSidak, alphaBonf = correction(pvals=values_to_correct,
                                                                        method=self.user_correction,
                                                                        is_sorted=False, returnsorted=False)
            pvals = pvals_corrected.tolist()
            for domain, res in zip(domains, pvals):
                complete.at[domain, 'Poisson'] = res
        else:
            for domain, res in zip(domains, values_to_correct):
                complete.at[domain, 'Poisson'] = res

        for domain in domains:
            complete.at[domain, 'Occurence in genome'] = sum(all_dict[domain]['genome'])
            complete.at[domain, 'Avg occurence in genome'] = sum(all_dict[domain]['genome'])

            complete.at[domain, 'Occurence in neighborhood'] = sum(all_dict[domain]['neigh']) / overall_num_of_neigh
            complete.at[domain, 'Avg occurence in neigh'] = sum(all_dict[domain]['neigh']) / overall_num_of_genomes
            complete.at[domain, 'Family percentage'] = len(all_dict[domain]['neigh']) / size_of_group
        dataframe = add_information(complete, domain_information)
        dataframe.index.name = 'Domain'
        dataframe.to_csv('temp/' + save, sep='\t', mode='w')

    def go(self):
        """
        Zipping all functions and execute them
        """
        # 1. Open correct database

        groups_calculated = []
        genera = self.open_database(self.user_organism)

        # 2. open important files
        genome_id_size_in_gene = open_multiple_line('/usr/src/app/important_files/GENOME_ID_SIZE_IN_GENE_NEW.txt')
        domain_information = collect_pfam_domain_description('/usr/src/app/important_files/domain_information')
        with open('/usr/src/app/important_files/counter_domains', 'r') as handler:
            pfam_counter_domains = json.load(handler)

        genome_id = [x[0] for x in genome_id_size_in_gene]
        size_in_gene = [x[1] for x in genome_id_size_in_gene]

        # 3. create variables
        genome_number_overall = 0
        genome_number_to_stat = 0
        overall_num_of_genomes = 0
        overall_num_of_neigh = 0
        print("checkpoint 2")

        # 4. Go through all taxonomy levels.
        for genus, genomes in genera.items():
            just_neigh_data = []
            genomes_with_domain = []
            neigh_genome_size = []
            counter = 0
            counter_db = 0
            sum_of_genes_in_neighborhoods = 0
            number_of_neighborhoods = 0
            percentage_counter = Counter()
            if len(genomes) > 25:

                # 5. Go through all genomes in taxonomy level
                for genome in genomes:
                    counter_db += 1
                    counter += 1
                    try:
                        genes, start_coord, end_coord, orientation, domains, contig = create_six_list(genome)
                    except FileNotFoundError:
                        continue
                    except ValueError:
                        continue
                    number_of_genes_in_genome = self.genome_size_in_gene(genome, genome_id, size_in_gene)
                    coords = searching_for_domain_in_genome(self.user_pfam, start_coord, end_coord, orientation,
                                                            domains, contig, genes)
                    if len(coords) > 0:
                        genome_number_to_stat += 1
                        overall_num_of_genomes += 1
                    for point in coords:
                        number_of_neighborhoods += 1
                        overall_num_of_neigh += 1

                        if self.user_distance != 0:
                            last_coordinate, first_coordinate, pfam_beg, pfam_end, searched_pfam_orientation = \
                                get_range_coordinates(point, self.user_distance)
                            pfam_index_to_neigh = get_domains_both_strands(start_coord, end_coord, last_coordinate,
                                                                           first_coordinate, pfam_beg, pfam_end,
                                                                           contig,
                                                                           point)
                        else:
                            pfam_index_to_neigh = get_domain_both_zeros(genes, point[4])

                        genes_in_neigh = list_of_genes_in_neigh(genes, pfam_index_to_neigh)
                        number_of_genes_in_neigh = size_in_genes(genes_in_neigh)
                        whole = get_list_of_domain_in_neigh(pfam_index_to_neigh, genome, domains)
                        percentage_counter += Counter(set(whole))
                        just_neigh_data.append(whole)
                        genomes_with_domain.append(genome)
                        neigh_genome_size.append((number_of_genes_in_neigh, number_of_genes_in_genome))
                        sum_of_genes_in_neighborhoods += int(number_of_genes_in_neigh)

                # 6. Gather all data from all genomes in taxonomy level
                alls = []
                neigh_counter = Counter()

                for i in just_neigh_data:
                    for j in i:
                        if j not in alls:
                            alls.append(j)
                genome_counter = Counter()

                # 7. Create matrix of data to calculate statistics
                counter = 0
                for domain_list, genome_in_domain, ng_size in zip(just_neigh_data, genomes_with_domain,
                                                                  neigh_genome_size):
                    counter += 1
                    genes, start_coord, end_coord, orientation, domains, contig = create_six_list(genome_in_domain)
                    genome_domains_counter = Counter(domains)
                    genome_counter += genome_domains_counter
                    neigh_domains_counter = Counter(domain_list)
                    neigh_counter += neigh_domains_counter

                # 8. Create pandas dataframe for storring results
                najlepsze_dane = pd.DataFrame(columns=['occurence in neighborhoods',
                                                       'occurence genomes',
                                                       'P', 'N', 'K',
                                                       'domain in whole database'])

                # 9. Calculate statistics
                counter = 0
                for i in alls:
                    counter += 1

                    p_KP = ((sum_of_genes_in_neighborhoods / genome_number_to_stat) * pfam_counter_domains[
                        i]) / 326628067

                    najlepsze_dane.at[i, 'P'] = p_KP
                    najlepsze_dane.at[i, 'domain in whole database'] = pfam_counter_domains[i]
                    najlepsze_dane.at[i, 'K'] = int(percentage_counter[i])
                    najlepsze_dane.at[i, 'N'] = number_of_neighborhoods
                    najlepsze_dane.at[i, 'occurence in neighborhoods'] = neigh_counter[i]
                    najlepsze_dane.at[i, 'occurence genomes'] = genome_counter[i]

                added_information = add_information(najlepsze_dane, domain_information)
                final_data = self.pfam_for_pf(added_information)

                if len(najlepsze_dane) > 0:
                    groups_calculated.append(genus)
                    self.save_data(genus, final_data)

            complete_results, additional = self.analysis(groups_calculated)

            self.merge_all_data_and_save(complete_results, additional, len(groups_calculated), domain_information,
                                         self.user_output, overall_num_of_genomes, overall_num_of_neigh)


class NeighborhoodAnalyzerFromGene:
    """
    Class for analyzing neighborhoods from gene list.
    """

    def __init__(self, user_list_of_genes, user_distance_value, user_database, user_cutoff, user_correction,
                 user_strand_value, user_output):
        """ Initializing input data
            Args:
                user_list_of_genes:  str (from pfam00000 to pfam99999)
                user_distance_value: str non negative integer 1-20000
                user_database: str
                user_cutoff: str (none, 0, 0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005,
                                 0.01, 0.05, 0.1, 0.2, 0.5)
                user_correction: str (none, bonferroni, fdr_bh)
                user_strand_value: str (both)
                user_output: str
         """
        self.user_list_of_genes = user_list_of_genes
        self.user_distance_value = user_distance_value
        self.user_database = user_database
        self.user_cutoff = user_cutoff
        self.user_correction = user_correction
        self.user_strand_value = user_strand_value
        self.user_output = user_output

    def save_complete(self, file, complete_data,):
        """
        Saves all stuff together as one file
        Args:
            file (str): path_of_file
            complete_data (pd.DataFrame): results to save

        """
        complete_data.index.name = 'Domain'
        complete_data.to_csv('/usr/src/app' + file, sep='\t', mode='w')


    def get_size_in_genes(self, genome_or_gene_list):
        """ Takes list of genes and returns size of neighbourhood """

        list_of_genes = []
        for gene in genome_or_gene_list:
            list_of_genes.append(gene)
        return len(set(list_of_genes))



    def pfam_for_pf(self, dataframe,):
        """Final customization of dataframe
            --> PF02696 instead of pfam02696
            --> PVALUE in 10e-3 format
            --> In what percentage % format
            --> avg occurences and density  have now 3 digits after coma
            --> drop NOam
            -->
        """

        indexy = dataframe.index
        for i in indexy:
            dataframe = dataframe.rename(index={i: i[:2].upper() + i[4:]})
        dataframe['Poisson'] = dataframe['Poisson'].map('{:.3e}'.format)
        dataframe['In what percentage'] = dataframe['In what percentage'].map('{:.3%}'.format)
        dataframe = dataframe.round({"average occurence in neighbourhood": 3, "average occurence in genome": 3,
                                     "Density difference": 3})
        if "NOam" in indexy:
            dataframe = dataframe.drop(index="NOam")
        # dataframe = dataframe.style.format({'In what percentage': '{:,.3%}'.format})
        dataframe.drop('average occurence in neighborhood', axis=1, inplace=True)
        dataframe.drop('average occurence in genome', axis=1, inplace=True)
        dataframe.drop('avg neigh size', axis=1, inplace=True)
        # dataframe.drop('domain in whole database', axis=1, inplace=True)
        # if self.user_correction == 'none':
        #     # dataframe.drop('Poisson_correct', axis=1, inplace=True)
        # elif :
        #     dataframe.drop('Poisson', axis=1, inplace=True)
        try:
            dataframe = dataframe.drop(index="NOam")
        except KeyError:
            pass
        return dataframe

    def make_genera_statistic(self, genome_list):
        """
        Makes statistic of genera in neighbourhood
        Args:
            genome_list (list): list of genomes in neighbourhood
        Returns:
        """

        with open("/usr/src/app/important_files/genera_statistics", "r") as handler:
            data = [x.strip().split() for x in handler]  # data to lista list [[genom,genus],[genome,genus]]
        genomes = [x[0] for x in data]
        genera = [x[1] for x in data]

        just_genera = []

        for genome in genome_list:
            indeks = genomes.index(genome)
            just_genera.append(genera[indeks])
        counter_genera = Counter(just_genera)

        zliczanie = []
        for genus, value in counter_genera.items():
            zliczanie.append((genus, value / len(genome_list)))
        return zliczanie



    def go(self):
        # 1. Open correct important files
        query_data = open_query("/usr/src/app" + self.user_list_of_genes)
        domain_information = collect_pfam_domain_description(
            "/usr/src/app/important_files/domain_information")
        genome_gene_reference = open_multiple_column_file(
            "/usr/src/app/important_files/genomes_map")
        genome_size_in_gene = open_multiple_column_file(
            "/usr/src/app/important_files/GENOME_ID_SIZE_IN_GENE.txt")
        with open('/usr/src/app/important_files/counter_domains', 'r') as handler:
            pfam_counter_domains = json.load(handler)

        # 2. create variables
        genome_ids = [x[0] for x in genome_size_in_gene]
        genome_size = [x[1] for x in genome_size_in_gene]
        genome_size_in_gene = None
        just_neigh_data = []
        genomes_with_domains = []
        neigh_genome_size = []
        counter = 0
        counter_db = 0
        sum_of_genes_in_neighborhoods = 0
        number_of_neighborhoods = 0
        genome_number_to_stat = 0

        percentage_counter = Counter()

        # 3. Go through all genes in query
        for gene in query_data:
            counter_db += 1
            counter += 1
            try:

                correct_genome = find_correct_genome(gene, genome_gene_reference)
                correct_size = genome_size_in_gene(correct_genome, genome_ids, genome_size)
                ref_gene, ref_start_coord, ref_end_coord, ref_orient, ref_domain, ref_contig = create_six_list(
                    correct_genome)
                search_result = search_for_gene_in_genome(gene, ref_gene, ref_start_coord, ref_end_coord,
                                                               ref_orient,
                                                               ref_contig)
                if self.user_distance_value != 0:
                    last_coordiate, first_coordinate, gene_beg, gene_end, searched_gene_orientation, = \
                        get_range_coordinates(search_result, self.user_distance_value)
                    pfam_index_to_neigh = get_domains_both_strands(ref_start_coord, ref_end_coord, last_coordiate,
                                                                   first_coordinate, gene_beg, gene_end, ref_contig,
                                                                   search_result)

                else:
                    pfam_index_to_neigh = get_domain_both_zeros(ref_gene, search_result[4])

                genes_in_neigh = list_of_genes_in_neigh(ref_gene, pfam_index_to_neigh)
                number_of_genes_in_neigh = self.get_size_in_genes(genes_in_neigh)
                whole_neigh = get_list_of_domain_in_neigh(pfam_index_to_neigh, gene, ref_domain)
                just_neigh_data.append(whole_neigh)
                genomes_with_domains.append(correct_genome)
                neigh_genome_size.append((number_of_genes_in_neigh, correct_size))
                percentage_counter += Counter(whole_neigh)
                genome_number_to_stat += 1
                sum_of_genes_in_neighborhoods += int(number_of_genes_in_neigh)
                number_of_neighborhoods += 1
            except ValueError:
                collect_gene_without_genome(gene)
                continue
        alls = []
        neigh_counter = Counter()

        # 4. Gather all data from all genomes in taxonomy level
        for i in just_neigh_data:
            for j in i:
                if j not in alls:
                    alls.append(j)
        genome_counter = Counter()
        # 5. Create matrix of data to calculate statistics
        counter = 0
        for domain_list, genome_in_domain, ng_size in zip(just_neigh_data, genomes_with_domains, neigh_genome_size):
            counter += 1
            genes, start_coord, end_coord, orientation, domains, contig = create_six_list(genome_in_domain)
            genome_domains_counter = Counter(domains)
            genome_counter += genome_domains_counter
            neigh_domains_counter = Counter(domain_list)
            neigh_counter += neigh_domains_counter

        #6. Create dataframe for all data
        najlepsze_dane = pd.DataFrame(columns=['occurence in neighborhoods',
                                               'average occurence in neighborhood', 'occurence genomes',
                                               'average occurence in genome', 'In what percentage', 'Poisson',
                                               'avg neigh size', 'Family', 'Summary'])

        # 7. Calculate statistics
        counter = 0
        scores_poisson = []
        for i in alls:
            counter += 1
            p_KP = ((sum_of_genes_in_neighborhoods / genome_number_to_stat) * pfam_counter_domains[i]) / 326628067

            najlepsze_dane.at[i, 'avg neigh size'] = sum_of_genes_in_neighborhoods / genome_number_to_stat
            najlepsze_dane.at[i, 'domain in whole database'] = pfam_counter_domains[i]
            poisson_test_result = poisson_distribution(p_KP, number_of_neighborhoods, i, percentage_counter)
            najlepsze_dane.at[i, 'Poisson'] = poisson_test_result
            scores_poisson.append(poisson_test_result)
        scores_after_correction = multiple_test_correction_list(scores_poisson, self.user_correction)

        for domain, wynik in zip(najlepsze_dane.index.to_list(), scores_after_correction):
            najlepsze_dane.at[domain, 'Poisson_correct'] = wynik
            najlepsze_dane.at[domain, 'In what percentage'] = percentage_counter[domain] / len(just_neigh_data)
            najlepsze_dane.at[domain, 'occurence genomes'] = genome_counter[domain]
            najlepsze_dane.at[domain, 'average occurence in genome'] = genome_counter[domain] / len(
                set(genomes_with_domains))
            najlepsze_dane.at[domain, 'occurence in neighborhoods'] = neigh_counter[domain]
            najlepsze_dane.at[domain, 'average occurence in neighborhood'] = neigh_counter[domain] / len(
                just_neigh_data)


        after_cutoff = cutoff_value(najlepsze_dane, self.user_cutoff)
        sorted_table = sort_table(after_cutoff)
        added_information = add_information(sorted_table, domain_information)
        final_data = self.pfam_for_pf(added_information)

        match_go_to_pfam_and_save(final_data, self.user_output)

        self.save_complete(self.user_output, final_data)
