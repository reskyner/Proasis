#!/dls/science/groups/i04-1/software/pandda-update/ccp4-linux64-2016-10-22-0115/libexec/python2.7

import sqlite3
import sys, os, getopt, subprocess, re
from datetime import datetime
from rdkit import Chem
from rdkit.Chem import AllChem

USAGE = """
populate_hits.py - Autodetect and add refined structures to proasis

SYNOPSIS: 

    usage: populate_hits.py -d <processing_directory> [-db <database file>]

OPTIONS:

    -d      processing directory (the same one used for XCE)
    -db     custom database file link (part of filepath after processing directory)

NOTES:

    This script uses the proasis autodetect script to add new ligand refinements to proasis.
    For this to work, refine.bound.pdb files must exist, and the ligand should be specified
    within the pdb file by the residue name "LIG", which is the default for XCE. Additionally,
    soakDB will be updated with a 'proasis' table, containing each of the succesful uploads 
    information: crystal name, filepath for the refine.bound.pdb file, the compound smiles 
    (to generate 2D sdf files) and the proasis strucid for each successful upload. 

"""


def create_sd_file(name, smiles, save_directory):
    """
    Create a 2D sdf file in the proasis project directory for succesfully detected ligands
    """
    # create sdf file for ligand and save to hit directory
    canon_smiles = Chem.CanonSmiles(smiles)
    mol = Chem.MolFromSmiles(canon_smiles)
    AllChem.Compute2DCoords(mol)
    print('Generating sdf file and saving to ' + name + ' directory...\n')
    sd_file = Chem.SDWriter(save_directory)
    sd_file.write(mol)


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


def submit_hits(bound_pdb, hit_directory, crystal, protein_name, smiles):
    print('Copying refine.bound.pdb...')
    os.system(str('cp ' + bound_pdb + ' ' + hit_directory))

    # create 2D sdf files for all ligands from SMILES string
    create_sd_file(crystal, smiles, str(os.path.join(hit_directory, crystal + '.sdf')))

    print('detecting ligand for ' + crystal)
    pdb_file = open(bound_pdb, 'r')
    ligands = []
    for line in pdb_file:
        if "LIG" in line:
            lig_string = re.search(r"LIG.......", line).group()
            ligands.append(str(lig_string))

    ligands = list(set(ligands))


    if len(ligands)==1:
        print('submission string:\n')
        submit_to_proasis1 = str("/usr/local/Proasis2/utils/submitStructure.py -d 'admin' -f " + hit_directory +
                                "/refine.bound.pdb -l '" + lig_string + "' -m " +
                                str(os.path.join(hit_directory, crystal + '.sdf')) +
                                " -p " + str(protein_name) + " -t " + str(crystal) + " -x XRAY -N")
	submit_to_proasis2 = str("/usr/local/Proasis2/utils/submitStructure.py -d 'admin' -f " + hit_directory +
                                "/refine.output.bound-state.pdb -l '" + lig_string + "' -m " +
                                str(os.path.join(hit_directory, crystal + '.sdf')) +
                                " -p " + str(protein_name) + " -t " + str(crystal) + " -x XRAY -N")
        #print(submit_to_proasis + '\n')
        process = subprocess.Popen(submit_to_proasis1, stdout=subprocess.PIPE, shell=True)
        out, err = process.communicate()
        strucidstr = get_id_string(out)

	process = subprocess.Popen(submit_to_proasis2, stdout=subprocess.PIPE, shell=True)
        out, err = process.communicate()
        strucidstr = get_id_string(out)
        if strucidstr != '':
            print('Success... ' + crystal + ' submitted. ProasisID: ' + str(strucidstr) + '\n *** \n')
        else:
            print('Error: ' + str(err))

    elif len(ligands)>1:
        lig1 = ligands[0]
        lign = " -o '"
        for i in range(1,len(ligands)-1):
            lign += str(ligands[i] + ',')
        lign += str(ligands[len(ligands)-1] + "'")

        submit_to_proasis = str("/usr/local/Proasis2/utils/submitStructure.py -d 'admin' -f " + hit_directory +
                                "/refine.bound.pdb -l '" + lig1 + "' " + lign + " -m " +
                                str(os.path.join(hit_directory, crystal + '.sdf')) +
                                " -p " + str(protein_name) + " -t " + str(crystal) + " -x XRAY -N")

        print(submit_to_proasis + '\n')
        process = subprocess.Popen(submit_to_proasis, stdout=subprocess.PIPE, shell=True)
        out, err = process.communicate()
        strucidstr = get_id_string(out)
        if strucidstr != '':
            print('Success... ' + crystal + ' submitted. ProasisID: ' + str(strucidstr) + '\n *** \n')
        else:
            print('Error: ' + str(err))


    return strucidstr


def main(argv):
    # set
    processing_directory = ''
    database_file = 'database/soakDBDataFile.sqlite'
    initial_model_directory = 'analysis/initial_model'
    proasis_directory = '/dls/science/groups/proasis/LabXChem'
    try:
       opts, args = getopt.getopt(argv,"hd:bi:",["dir=", "dbfile=", "initialmodel"])
    except getopt.GetoptError:
       print USAGE
       sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print USAGE
            sys.exit()
        elif opt in ("-d", "--dir"):
            processing_directory = arg
        elif opt in ("-b", "--dbfile"):
            database_file = arg
            print('Custom database file selected as: ' + str(database_file))
	elif opt in ("-i", "--initialmodel"):
	    initial_model_directory = arg

    print('**********************************************************************************\n')
    print('Checking datafile exists: ' + str(os.path.isfile(os.path.join(processing_directory, database_file))))
    print('Searching the datafile ' + str(os.path.join(processing_directory, database_file) +
                                          ' for hits\n'))
    print('**********************************************************************************\n\n')

    # connect to soakDB and open a cursor
    conn = sqlite3.connect(os.path.join(processing_directory, database_file))
    c = conn.cursor()
    index = -1
    pro_no = 0
    for row in c.execute("SELECT Protein FROM soakDB;"):
        protein_name = row[0]
        print('**********************************************************************************\n')
        print('Protein name (from soakDB): ' + protein_name + '\n')
        print('**********************************************************************************\n\n')
        pro_no += 1
        if pro_no > 1:
            print('ERROR: more than one protein name found!')
            exit()

    # set dictionary keys for storing existing and found (table_dict) data
    proasis_table_dict = {'crystal': [], 'SMILES': [], 'boundDir': [], 'modified': [], 'proasisID':[]}
    proasis_existing_dict = {'crystal': [], 'SMILES': [], 'boundDir': [], 'modified': [], 'proasisID':[]}

    # set a list up to store structures which should be skipped (no ligand found)
    skip_list=[]
    
    # select data to search for refinements from mainTable in soakDB
    for row in c.execute("SELECT CrystalName, CompoundSMILES, RefinementPDB_latest, RefinementOutcome FROM mainTable;"):
        # If a refinement is found and refinement outcome is greater than or equal to 3
        if str(row[2]) != 'None' and int(re.search(r"\d+", row[3]).group()) >= 3:
            index+=1
            # search the refinement directory for refine.bound.pdb
            pdb = [os.path.join(processing_directory, initial_model_directory, row[0], 'refine.bound.pdb'), os.path.join(processing_directory, initial_model_directory, row[0], 'refine.output.bound-state.pdb')]
            for bound_pdb in pdb:
		    print('checking for ' + str(bound_pdb))
		    if os.path.isfile(bound_pdb):
			# create time string
			last_modified = str(datetime.fromtimestamp(os.path.getmtime(bound_pdb)))
			last_modified = last_modified.replace(":",'')
			last_modified = last_modified.replace('.','')
			last_modified = last_modified.replace('-','')
			last_modified = last_modified.replace(' ','')
			last_modified = last_modified[:12]

			# add data to proasis table dictionary
			proasis_table_dict['crystal'].append(str(row[0]))
			proasis_table_dict['SMILES'].append(str(row[1]))
			proasis_table_dict['boundDir'].append(bound_pdb)
			proasis_table_dict['modified'].append(str(last_modified))
			proasis_table_dict['proasisID'].append('')

			# if a directory for the proasis project does not exist
			if not os.path.isdir(os.path.join(proasis_directory, protein_name, str(row[0]))):
			    # make directory for hit and copy refine.bound.pdb
			    print('Creating directory for ' + str(row[0]) + '...')
			    hit_directory = os.path.join(proasis_directory, protein_name, str(row[0]) + '/')
			    os.system('mkdir ' + hit_directory)

			    crystal = str(row[0])
			    smiles = str(row[1])

			    strucidstr = submit_hits(bound_pdb, hit_directory, crystal, protein_name, smiles)

			    proasis_table_dict['proasisID'][index]=strucidstr
			else:
			    print('Files already exist in proasis directories! Nothing to do here.\n *** \n')
		    else:
			print('Warning: No refine.bound.pdb file found for ' + str(row[0]) + ', refinement may have failed!\n *** \n')
			index-=1
    # If a proasis table doesn't exist in soakDB, create it and add info for uploaded structures
    c.execute("CREATE TABLE IF NOT EXISTS 'proasis' ('crystal' TEXT, "
              "'smiles' TEXT, 'boundDir' TEXT, 'modifiedStamp' TEXT, 'proasisID' TEXT);")
    conn.commit()
    # a list to contain queries to be run
    query_list = []

    print('**********************************************************************************\n')

    # for all of the found crystal id's with refinements
    for i in range(0, len(proasis_table_dict['crystal'])):
        # if no strucid (not currently in proasis)
        if str(proasis_table_dict['proasisID'][i]) != '':
            try:
                # check whether there has already an entry in the proasis table
                counter = 0
                check = str("SELECT crystal FROM proasis WHERE crystal GLOB '%s'"
                            %str(proasis_table_dict['crystal'][i]))

                for row in c.execute(check):
                    counter+=1
                # yes - no autoupdate
                update = 0
            except:
                # otherwise - update required
                update = 1

            if counter<1:
                # add the relevant data into the proasis table
                query = str("INSERT INTO proasis (crystal, smiles, boundDir, modifiedStamp, proasisID) VALUES ('"
                            + str(proasis_table_dict['crystal'][i]) + "','" +
                            str(proasis_table_dict['SMILES'][i]) + "','" +
                            str(proasis_table_dict['boundDir'][i]) + "','" +
                            str(proasis_table_dict['modified'][i]) + "','" +
                            str(proasis_table_dict['proasisID'][i]) + "');")
                # add query to list to be executed later
                query_list.append(query)

    # perform all queries and commit changes to table
    if len(query_list) > 1:
        print('Updating soakDB with new proasis data...\n')
    else:
        print('No new entries...')
    for i in range(0, len(query_list)):
        print(str(query_list[i]))
        c.execute(query_list[i])
    conn.commit()
    try:

        print('**********************************************************************************\n')
        print('Will now check for updates to existing refinements...')
        print('**********************************************************************************\n\n')
        # retrieve last proasis update info from proasis table
        for row in c.execute("SELECT crystal, smiles, boundDir, modifiedStamp, proasisID FROM proasis"):
            proasis_existing_dict['crystal'].append(str(row[0]))
            proasis_existing_dict['SMILES'].append(str(row[1]))
            proasis_existing_dict['boundDir'].append(str(row[2]))
            proasis_existing_dict['modified'].append(str(row[3]))
            proasis_existing_dict['proasisID'].append(str(row[4]))

        # list for new entries
        new_entries = []
        remove_unrefined = []

        for i in range(0,len(proasis_table_dict['crystal'])):
            # check for duplicates
            if proasis_table_dict['crystal'][i] in proasis_existing_dict['crystal']:
                print('Ligand ' + str(proasis_table_dict['crystal'][i]) + ' already in db...')
                print('Checking if files have been modified...')
                # check whether the timestamps for the refine.bound.pdb files have changed
                for j, k in enumerate(proasis_existing_dict['crystal']):
                    if k == proasis_table_dict['crystal'][i]:
                        if proasis_existing_dict['modified'][j] != proasis_table_dict['modified'][i]:
                            print('The refine.bound.pdb file for ' + str(proasis_table_dict['crystal'][i])
                                  + 'has changed')
                            print('Updating file in proasis directories...\n')
                            bound_pdb = proasis_table_dict['boundDir'][i]
                            hit_directory = os.path.join(proasis_directory, protein_name,
                                                         proasis_table_dict['crystal'][i] + '/')
                            # if changed, copy over updated file
                            os.system(str('cp ' + bound_pdb + ' ' + hit_directory))
                        else:
                            print('No modification... do nothing\n')
                print('checking for changes to refinement outcome...')
                for row in c.execute("SELECT CrystalName, RefinementOutcome FROM mainTable WHERE CrystalName GLOB '%s'"
                                         %proasis_table_dict['crystal'][i]):
                    print(str("SELECT CrystalName, RefinementOutcome FROM mainTable WHERE CrystalName GLOB '%s'"
                                         %proasis_table_dict['crystal'][i]))
                    print row
                    if int(re.search(r"\d+", row[1]).group()) < 3:
                        remove_unrefined = str(row[0])
                        for struc in c.execute("SELECT proasisID FROM proasis WHERE crystal GLOB '%s'" %remove_unrefined):
                            strucstr = str(struc[0])
                            c.execute("DELETE FROM proasis WHERE crystal GLOB '%s'" %remove_unrefined)
                            conn.commit()
                        print('Deleting old entry: ' + strucstr)
                        delete = str('/usr/local/bin/python /usr/local/Proasis2/utils/removestruc.py -s ' + strucstr)
                        os.system(delete)

                print(remove_unrefined)
            else:
                new_entries.append(str(proasis_table_dict['crystal'][i]))


    except:
        print('Problem at update refinement state')

    remove_skips = list(set(skip_list) - set(new_entries))

    try:
        # if there is nothing to be done, carry on
        if len(remove_skips) == 0:
            print('**********************************************************************************\n')
            print('No modifications...')
            print('**********************************************************************************\n')

        for i in range(0,len(remove_skips)):

            try:
                # remove old structures from proasis
                for row in c.execute("SELECT crystal, smiles, boundDir, proasisID FROM proasis WHERE crystal GLOB '%s'" %remove_skips[i]):
                    print('Deleting old entry: ' + str(row[3]))
                    delete = str('/usr/local/bin/python /usr/local/Proasis2/utils/removestruc.py -s ' + str(row[3]))
                    os.system(delete)
                    crystal = str(row[0])
                    smiles = str(row[1])
                    hit_directory = os.path.join(proasis_directory, protein_name, str(remove_skips[i]) + '/')
                    bound_pdb = str(row[2])
                    print('Adding new entry for ' + crystal)
                    strucidstr = submit_hits(bound_pdb, hit_directory, crystal, protein_name, smiles)
                    print('New proasisID for ' + crystal + ': ' + strucidstr)
                    print('Updating soakDB...\n')
                    update_id = str("UPDATE proasis SET proasisID GLOB '%s' WHERE crystal GLOB '%s'" %strucidstr %remove_skips[i])
                    c.execute(update_id)
                    conn.commit()


            except:
                print('There was an error whilst attempting to update ' + crystal + ', but the old entry has probably been removed from proasis. Please check proasis, and the refine.bound.pdb file!\n')

        print('**********************************************************************************\n\n')

    except:
        print('Something else went wrong!')

    print('DONE! :)\n\n')
    conn.close()

if __name__ == "__main__":
    main(sys.argv[1:])
