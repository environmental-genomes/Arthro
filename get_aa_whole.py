#python3
import os
import glob
from Bio import SeqIO


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '../../../Desktop/get_aa_test'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


# define input file
faa_file_folder = '/mnt/d/work/arthrobacter/AEM/revisionForAEM/Arthro_faaTable1.3Houskeeping/grouped32faa/'


amino_acid_str = 'ARNDBCEQZGHILKMFPSTWYV'

faa_file_re = '%s/*.faa' % faa_file_folder
faa_file_list = [os.path.basename(file_name) for file_name in glob.glob(faa_file_re)]


print('Genome\t%s' % '\t'.join([i for i in amino_acid_str]))
for faa_file in faa_file_list:

    pwd_faa = '%s/%s' % (faa_file_folder, faa_file)
    faa_file_path, faa_file_basename, faa_file_extension = sep_path_basename_ext(pwd_faa)

    aa_count_dict = {}
    aa_count_total = 0
    for seq_record in SeqIO.parse(pwd_faa, 'fasta'):
        seq_record_seq = str(seq_record.seq)
        aa_count_total += len(seq_record_seq)
        for aa in seq_record_seq:
            if aa not in aa_count_dict:
                aa_count_dict[aa] = 1
            else:
                aa_count_dict[aa] += 1

    aa_count_dict_pct = {}
    for aa in aa_count_dict:
        aa_count_dict_pct[aa] = float("{0:.4f}".format(aa_count_dict[aa]*100/aa_count_total))

    amino_acid_pct_list = []
    for amino_acid in amino_acid_str:

        amino_acid_pct = 0
        if amino_acid in aa_count_dict_pct:
            amino_acid_pct = aa_count_dict_pct[amino_acid]

        amino_acid_pct_list.append(str(amino_acid_pct))

    print('%s\t%s' % (faa_file_basename, '\t'.join(amino_acid_pct_list)))
