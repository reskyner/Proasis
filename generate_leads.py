#!/dls/science/groups/i04-1/software/pandda-update/ccp4-linux64-2016-10-22-0115/libexec/python2.7

import sys, getopt, os, subprocess, re
import sqlite3

from Bio.PDB import NeighborSearch, PDBParser, Atom, Residue, PDBIO
import pandas as pd
import numpy as np

def get_id_string(out):
    """
    Regex function for finding proasis strucid
    """
    try:
        strucidstr = re.search(r"strucid='.....'", out)
        strucidstr = strucidstr.group()
        strucidstr = strucidstr.replace('strucid=', '')
        strucidstr = strucidstr.replace("'", '')
    except:
        strucidstr = ''
    return strucidstr


def main(argv):
    reference_structure = ''
    pandda_analyse_centroids = ''
    name = ''
    processing_directory=''
    database_file = 'database/soakDBDataFile.sqlite'
    try:
       opts, args = getopt.getopt(argv,"hr:p:n:d:",["rfile=","pfile=","name="])
    except getopt.GetoptError:
       print 'usage: generate_leads.py -r <reference_pdb.pdb> -p <pandda_analyse_sites.csv> -n <name>'
       sys.exit(2)
    except:
       print 'usage: generate_leads.py -r <reference_pdb.pdb> -p <pandda_analyse_sites.csv> -n <name>'
       sys.exit(2)
    for opt, arg in opts:
       if opt == '-h':
          print 'usage: generate_leads.py -r <reference_pdb.pdb> -p <pandda_analyse_sites.csv> -n <name>'
          sys.exit()
       elif opt in ("-r", "--rfile"):
          reference_structure = arg
       elif opt in ("-p", "--pfile"):
          pandda_analyse_centroids = arg
       elif opt in ("-n", "--name"):
           name = arg
       elif opt in ("-d", "--dir"):
           processing_directory = arg


    print('\n The apo structure ' + str(name) + ' in ' + str(reference_structure) +
          ' will be searched for leads defined by the site centroids in ' + str(pandda_analyse_centroids) + '\n')

    # read structure from reference file with BioPy PDBParser
    structure = PDBParser(PERMISSIVE=0).get_structure(str(name), str(reference_structure))

    # read centroids from pandda analyse sites list (native)
    site_list = pd.read_csv(str(pandda_analyse_centroids))['native_centroid']
    print(' Searching for residue atoms for ' + str(len(site_list)) + ' site centroids \n')
    print(' NOTE: 3 residue atoms are required for each site centroid \n')

    print site_list

    no = 0
    for centroid in site_list:
        print('next centroid')
        structure = PDBParser(PERMISSIVE=0).get_structure(str(name), str(reference_structure))
        no += 1
        res_list = []

        # initial distance for nearest neighbor (NN) search is 2A
        neighbor_distance = 20

        centroid_coordinates = centroid.replace('(', '[')
        centroid_coordinates = centroid_coordinates.replace(')', ']')
        centroid_coordinates = eval(str(centroid_coordinates))

        # define centroid as an atom object for NN search
        centroid_atom = Atom.Atom('CEN', centroid_coordinates, 0, 0, 0, 0, 9999, 'C')
        atoms = list(structure.get_atoms())
        center = np.array(centroid_atom.get_coord())
        ns = NeighborSearch(atoms)

        # calculate NN list
        neighbors = ns.search(center, neighbor_distance)
        res_list = []

        # for each atom in the NN list
        for neighbor in neighbors:
            try:
                # get the residue that the neighbor belongs to
                parent = Atom.Atom.get_parent(neighbor)
                # if the residue is not a water etc. (amino acids have blank)
                if parent.get_id()[0] == ' ':
                    # get the chain that the residue belongs to
                    chain = Residue.Residue.get_parent(parent)
                    # if statements for fussy proasis formatting
                if len(str(parent.get_id()[1]))==3:
                    # residue string = 'RES CHAIN NUMBER :...'
                    res = (str(parent.get_resname()) + ' ' + str(chain.get_id()) + ' ' + str(parent.get_id()[1]))
                    res_list.append(res)
                if len(str(parent.get_id()[1]))==2:
                    res = (str(parent.get_resname()) + ' ' + str(chain.get_id()) + '  ' + str(parent.get_id()[1]))
                    res_list.append(res)
            except:
                break
    res_list = (list(set(res_list)))
    print res_list
    lig1 = str("'" + str(res_list[0]) + ' :' + str(res_list[1]) + ' :'
                                        + str(res_list[2]) + " ' ")
    print lig1

    res_string = "-o '"

    for i in range(3,len(res_list)-1):
        res_string += str(res_list[i] + ' ,')
        res_string += str(res_list[i+1] + ' ')
    #print str(res_string[i+1])
    submit_to_proasis = str('/usr/local/Proasis2/utils/submitStructure.py -p ' + str(name) + ' -t ' + str(name + '_lead -d admin -f '+ str(reference_structure) + ' -l ' + str(lig1)) +  str(res_string) + "' -x XRAY -n")
    process = subprocess.Popen(submit_to_proasis, stdout=subprocess.PIPE, shell=True)
    out, err = process.communicate()
    print out
    strucidstr = re.search(r"strucid='.....'", out)
    strucidstr = strucidstr.group()
    strucidstr = strucidstr.replace('strucid=','')
    strucidstr = strucidstr.replace("'",'')
    add_lead = str('/usr/local/Proasis2/utils/addnewlead.py -p ' + str(name) + ' -s ' + str(strucidstr))
    os.system(add_lead)

    #process = subprocess.Popen(add_lead, stdout=subprocess.PIPE, shell=True)
    #out, err = process.communicate()
    #strucidstr = get_id_string(out)
    if strucidstr != '':
        print('Success... ' + name + ' submitted. ProasisID: ' + str(strucidstr))
        conn = sqlite3.connect(os.path.join(processing_directory, database_file))
        c = conn.cursor()
        c.execute("CREATE TABLE IF NOT EXISTS 'proasisLead' ('protein' TEXT, "
                  "'leadDir' TEXT, 'proasisID');")
        c.execute(str("INSERT INTO proasisLead (protein, leadDir, proasisID) VALUES ('"
                            + str(name) + "','" +
                            str(reference_structure) + "','" +
                            str(strucidstr) + "');"))
        conn.commit()
        c.close()
    else:
        print('Error: ' + str(err))

if __name__ == "__main__":
   main(sys.argv[1:])
