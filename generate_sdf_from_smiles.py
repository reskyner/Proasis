#!/dls/science/groups/i04-1/software/pandda-update/ccp4-linux64-2016-10-22-0115/libexec/python2.7

from rdkit import Chem
from rdkit.Chem import AllChem
import getopt, sys

def create_sd_file(smiles, save_directory):
    """
    Create a 2D sdf file in the proasis project directory for succesfully detected ligands
    """
    # create sdf file for ligand and save to hit directory
    canon_smiles = Chem.CanonSmiles(smiles)
    mol = Chem.MolFromSmiles(canon_smiles)
    AllChem.Compute2DCoords(mol)
    #print('Generating sdf file and saving to ' + name + ' directory...\n')
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

def main(argv):
    #struc_name=''
    smiles=''
    save_directory=''
    #pdb_file=''
    try:
	opts, args = getopt.getopt(argv, "s:d:", ["smiles=","directory="])
    except:
	sys.exit(2)
    #if opt in ("-n", "--name"):
#	struc_name = arg
    for opt, arg in opts:
	if opt in ("-s", "--smiles"):
		smiles = arg
    	elif opt in ("-d", "--directory"):
		save_directory = arg
    #elif opt in ("-f", "--file"):
#	pdb_file = arg

    create_sd_file(smiles, save_directory)
    
if __name__ == "__main__":
    main(sys.argv[1:])


