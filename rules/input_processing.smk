import os, yaml
from snakemake.utils import validate


# Limitting Step for the run of Snakemake, creating wildcards
# If the input is a text file containing the folder structure to the fastqs
if os.path.isfile(config['select_fastqs']) and not config['select_fastqs'].endswith('.yaml') and not config['select_fastqs'].endswith('.yml'):

    #Lists that will contain wildcards
    pool=[] # wildcard 'pool'
    # If fastq files are called Sample_xxxxx-cDNA/xxxxx-cDNA_R1.fastq.gz

    with open(config['select_fastqs']) as fq:
        for line in fq:
            if not line.startswith('#'):
                line_sp = line.split('/')
                pool.append(line_sp[0].strip().replace('-cDNA', ''))

    # Create a dict of wildcards
    wildcards_list={'pool':pool}

    # EXAMPLE: 2 WILDCARDS for FASTQs
    # round_num=[] # wildcard 'num'
    # sample_name=[] # wildcard 'id1'
    # for eg. round5/Sample_xxxxx-cDNA/xxxxx-cDNA
    # round_num = round5
    # sample_name = xxxxx
    # 
    # with open(config['select_fastqs']) as fq:
    #     for line in fq:
    #         if not line.startswith('#'):
    #             line_sp = line.split('/')
    #             round_num.append(line_sp[0])
    #             sample_name.append(line_sp[1].strip().replace('-cDNA', ''))

    # Create a dict of wildcards
    # wildcards_list={'num':round_num, 'id1':sample_name}


# If the input is yaml file then validate it and process it
# The example in here is for more than one wildcards
elif os.path.isfile(config['select_fastqs']) and (config['select_fastqs'].endswith('.yaml') or config['select_fastqs'].endswith('.yml')):
    # Validating config['select_fastqs']
    # with open(config['select_fastqs']) as fout:
    #     samples_df = pd.json_normalize(yaml.load(fout, Loader=yaml.SafeLoader))
    # validate(samples_df, "samples.schema.json")

    # Value of last_step should be regulated by the json schema
    if not isinstance(config['last_step'], str) and not (config['last_step'].endswith('.yaml') and not config['last_step'].endswith('.yml')):
        raise ValueError("Value in last_step should be one of the modules names!")


    # File exists or not
    # if not os.path.isfile(config['select_fastqs']):
    #     raise OSError("The provided file doesn't exists! Check the path.")

    # Parse the samples containing yaml file
    with open(config['select_fastqs']) as fout:
        sample_dict = yaml.load(fout, Loader=yaml.SafeLoader)

    # contains the keys that categorizes multiple
    # modules (same names should be present in the modules_yaml file)
    # sample_set_names = [] 
    # files_list = [] # List of list
    # modules_list = [] # List of list


    # Yaml file can't be a mix of dirs and folder structures
    # Need refinement (with this code it will be checking each set of 
    # samples but we need a code which checks ALL samples present in 
    # ALL separate sets to be similar)
    # Checking with respect to fastq folders
    if isinstance(sample_dict, dict):
        for k, v in sample_dict.items():
            if all([ isfolder_struct(config['cDNA_fastqs_dir'], x) for x in v ]):
                continue
            elif all(list(map(os.path.isdir, v))):
                continue
            else:
                raise ValueError("Expected either all folder stuctures for the samples or a dir containing the fastqs but neither was provided")

    elif isinstance(sample_dict, list):
        # If it contains folder_structure for all samples
        if all([ isfolder_struct(config['cDNA_fastqs_dir'], x) for x in sample_dict ]):
            FOLD_STRUCT = True
            FASTQ_DIR = False

        # cDNA fastqs dir and HTO fastqs dir need to be none when providing fastq dir and inferring wildcards from there
        # During this process cDNA and HTO files are expected to be segregated by a dir called 'cDNA' or 'HTO' and not containing
        # these as suffixes in the filename.
        elif os.path.isdir(sample_dict[0]):
            FOLD_STRUCT = False
            FASTQ_DIR = True

        else:
            raise ValueError("Expected either all folder stuctures for the samples or a dir containing the fastqs but neither was provided")

    # If a string is provided, it is expected to be a folder structure with wildcards to extracted
    elif isinstance(sample_dict, str):
        samp_type, proj_type, sample_name = glob_wildcards(sample_dict) # Make sure there are same number of wildcards on LHS and RHS

    else:
        raise ValueError("Unexpected Type of input in samples containing yaml file!")

    # Check if the yaml file has list of file names (or folder structure
    # for each sample) or list of dirs each containing samples that 
    # need separate modules.
    # Similar to previous 'for loop' needs refinement
    # for k, vals in sample_dict.items():
    #     sample_set_names.append(k)
    #     files_list.append(vals)


    wildcards_list = dict.fromkeys(sample_dict.keys()) # dict of list
    if FOLD_STRUCT and not FASTQ_DIR:
        # How to parse the values of wildcards (also separate wildcards 
        # for each module)
        for k, f_l in sample_dict.items():
            samp_type = [] # wildcard 'samp_type'
            proj_type = [] # wildcard 'proj'
            sample_name = [] # wildcard 'samp_name'
            for f in f_l:
                line_sp = f.split('/')
                samp_type.append(line_sp[0].strip())
                proj_type.append(line_sp[1].strip())
                sample_name.append(line_sp[2].strip())

            wildcards_list[k] = {'samp_type':samp_type, 'proj':proj_type, 'samp_name':sample_name }
            # if k == 'multip_samples':
            #     wildcards_list[k] = {'samp_type':samp_type, 'proj':proj_type, 'samp_name':sample_name }
            # else:
            #     wildcards_list[k] = {'samp_type':samp_type, 'proj':proj_type}

    # Create a dict of wildcards
    # For multi_modules create list of dict of wildcards (or dict of 
    # wildcards dicts with the same values use to differentiate each 
    # multi_module set as key)


# If the input is a dir containing all the fastqs or files
# Extract the required wildcards here
elif os.path.isdir(config['select_fastqs']):
    assert config['wildcards_select'] is not None, "if 'select_fastqs' option in the yaml file is a dir then 'wildcards_select' can't be 'None'!"
    # Example for one wildcard extraction
    # ID1, = glob_wildcards(config['select_fastqs'])
    ID1, DONOR = glob_wildcards(os.path.join(config['select_fastqs'], config['wildcards_select']))
    
    # Example of how to filter out a list of files
    # In here some of the values in ID1 (wildcard)
    # are to be removed. So to remove it's corresponding 
    # values from DONOR

    # wants=["Ch-Plexus-11-C1", "Ch-Plexus-11-C2", 
    #        "Ch-Plexus-12-C1", "Ch-Plexus-12-C2", 
    #        "Ch-Plexus-16-C1", "Ch-Plexus-16-C2", 
    #        "Ch-Plexus-21-C1", "Ch-Plexus-21-C2", 
    #        "Ch-Plexus-24-A1", "Ch-Plexus-24-A2", 
    #        "Ch-Plexus-28-A1", "Ch-Plexus-28-A2", 
    #        "Ch-Plexus-4-A1", "Ch-Plexus-4-A2", 
    #        "Ch-Plexus-7-C1", "Ch-Plexus-7-C2",
    #        ]
    
    # to_rem=[]
    # ID1 = [ j if j in wants else to_rem.append(i) for i, j in enumerate(ID1) ]
    # ID1 = [ j for j in ID1 if j is not None ]
    # DONOR = [ j for i, j in enumerate(DONOR) if i not in to_rem]


    # Create a dict of wildcards
    wildcards_list={'pool': ID1, 'donor': DONOR}

else:
    raise ValueError("Unrecognized input file! Can't go ahead with the pipeline!")


# For multiruns of cellSNP and vireoSNP
# VCF_TYPE=config['gt_demux_pipeline']['vcf_info_columns'][2:]
