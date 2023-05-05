import os
import sys
import argparse
import subprocess
from Bio import SeqIO
import shutil
zol_main_directory = '/'.join(os.path.realpath(__file__).split('/')[:-2]) + '/'

def create_parser():
    """ Parse arguments """
    parser = argparse.ArgumentParser(description="""
	Program: setup_annotation_dbs.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology
        
    Downloads annotation databases for KO, PGAP, 
	""", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-p', '--download_path', help='Path to where the databases should be downloaded. Default is /path/to/zol_Github_clone/db/', required=False, default=zol_main_directory + 'db/')
    parser.add_argument('-c', '--cpus', type=int, help="Number of cpus/threads to use.", required=False, default=4)
    parser.add_argument('-m', '--minimal', action='store_true', help="Minimal mode - will only download PGAP.", required=False, default=False)

    args = parser.parse_args()
    return args


def setup_annot_dbs():
    myargs = create_parser()

    download_path = os.path.abspath(myargs.download_path) + '/'
    cpus = myargs.cpus
    minimal_mode = myargs.minimal

    try:
        assert(os.path.isdir(download_path))
    except:
        sys.stderr.write('Error: Provided directory for downloading annotation files does not exist!\n')

    if minimal_mode:
        sys.stdout.write('Minimal mode requested, will only be downloading the PGAP database.\n')

    try:
        shutil.rmtree(download_path)
        os.mkdir(download_path)
        readme_outf = open(download_path + 'README.txt', 'w')
        readme_outf.write('Default space for downloading KOFam annotation databases.\n')
        readme_outf.close()
    except Exception as e:
        sys.stderr.write('Issues clearing contents of db/ to re-try downloads.\n')
        sys.stderr.write(str(e) + '\n')
        sys.exit(1)

    listing_file = zol_main_directory + 'db/database_location_paths.txt'
    listing_handle = open(listing_file, 'w')

    if minimal_mode:
        # Final annotation files
        pgap_info_file = download_path + 'hmm_PGAP.tsv'
        pgap_phmm_file = download_path + 'PGAP.hmm'

        download_links = ['https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM.tgz',
                          'https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.tsv']

        # Download
        print('Starting download of files!')
        os.chdir(download_path)
        try:
            for dl in download_links:
                axel_download_dbs_cmd = ['axel', '-a', '-n', str(cpus), dl]
                os.system(' '.join(axel_download_dbs_cmd))
        except Exception as e:
            sys.stderr.write('Error occurred during downloading!\n')
            sys.stderr.write(str(e) + '\n')
            sys.exit(1)

        try:
            print('Setting up PGAP database!')
            os.system(' '.join(['tar', '-zxvf', 'hmm_PGAP.HMM.tgz']))
            assert (os.path.isfile(pgap_info_file))
            assert (os.path.isdir(download_path + 'hmm_PGAP.HMM/'))
            for f in os.listdir(download_path + 'hmm_PGAP.HMM/'):
                os.system(' '.join(['cat', download_path + 'hmm_PGAP.HMM/' + f, '>>', pgap_phmm_file]))
            pgap_descriptions_file = download_path + 'pgap_descriptions.txt'
            pdf_handle = open(pgap_descriptions_file, 'w')
            with open(pgap_info_file) as opil:
                for i, line in enumerate(opil):
                    if i == 0: continue
                    line = line.strip()
                    ls = line.split('\t')
                    label = ls[2]
                    description = ls[10]
                    pdf_handle.write(label + '\t' + description + '\n')
            pdf_handle.close()
            assert (os.path.isfile(pgap_phmm_file))
            os.system(' '.join(['hmmpress', pgap_phmm_file]))
            listing_handle.write('pgap\t' + pgap_descriptions_file + '\t' + pgap_phmm_file + '\n')
            os.system(' '.join(
                ['rm', '-rf', download_path + 'hmm_PGAP.HMM/', download_path + 'hmm_PGAP.HMM.tgz', pgap_info_file]))
        except Exception as e:
            sys.stderr.write('Issues setting up PGAP database.\n')
            sys.stderr.write(str(e) + '\n')
            sys.exit(1)
    else:
        # Final annotation files
        pfam_phmm_file = download_path + 'Pfam-A.hmm'
        ko_annot_info_file = download_path + 'ko_list'
        ko_phmm_file = download_path + 'profile.hmm'
        pgap_info_file = download_path + 'hmm_PGAP.tsv'
        pgap_phmm_file = download_path + 'PGAP.hmm'
        pb_faa_file = download_path + 'paperblast.dmnd'
        #pb_sql_file = download_path + 'litsearch.db'
        vog_phmm_file = download_path + 'vog.hmm'
        vog_info_file = download_path + 'vog.annotations.tsv'
        is_faa_file = download_path + 'isfinder.dmnd'
        mb_faa_file = download_path + 'mibig.dmnd'
        card_faa_file = download_path + 'card.dmnd'
        vfdb_faa_file = download_path + 'vfdb.dmnd'

        download_links = ['https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz',
                          'https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM.tgz',
                          'ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz',
                          'http://papers.genomics.lbl.gov/data/uniq.faa',
                          'http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz',
                          'https://dl.secondarymetabolites.org/mibig/mibig_prot_seqs_3.1.fasta',
                          'http://fileshare.csb.univie.ac.at/vog/latest/vog.hmm.tar.gz',
                          #'http://papers.genomics.lbl.gov/data/litsearch.db'
                          'https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.tsv',
                          'http://fileshare.csb.univie.ac.at/vog/latest/vog.annotations.tsv.gz',
                          'ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz',
                          'https://raw.githubusercontent.com/thanhleviet/ISfinder-sequences/master/IS.faa',
                          'https://card.mcmaster.ca/download/0/broadstreet-v3.2.5.tar.bz2']

        # Download
        print('Starting download of files!')
        os.chdir(download_path)
        try:
            for dl in download_links:
                axel_download_dbs_cmd = ['axel', '-a', '-n', str(cpus), dl]
                os.system(' '.join(axel_download_dbs_cmd))
        except Exception as e:
            sys.stderr.write('Error occurred during downloading!\n')
            sys.stderr.write(str(e) + '\n')
            sys.exit(1)

        try:
            print('Setting up Pfam database!')
            os.system(' '.join(['gunzip', 'Pfam-A.hmm.gz']))
            assert(os.path.isfile(pfam_phmm_file))
            name = None
            desc = None
            pfam_descriptions_file = download_path + 'pfam_descriptions.txt'
            pdf_handle = open(pfam_descriptions_file, 'w')
            with open(pfam_phmm_file) as oppf:
                for line in oppf:
                    line = line.strip()
                    ls = line.split()
                    if ls[0].strip() == 'NAME':
                        name = ' '.join(ls[1:]).strip()
                    elif ls[0].strip() == 'DESC':
                        desc = ' '.join(ls[1:]).strip()
                        pdf_handle.write(name + '\t' + desc + '\n')
            pdf_handle.close()
            os.system(' '.join(['hmmpress', pfam_phmm_file]))
            listing_handle.write('pfam\t' + pfam_descriptions_file + '\t' + pfam_phmm_file + '\n')
        except Exception as e:
            sys.stderr.write('Issues setting up Pfam database.\n')
            sys.stderr.write(str(e) + '\n')
            sys.exit(1)

        try:
            print('Setting up KO database!')
            os.system(' '.join(['gunzip', 'ko_list.gz']))
            os.system(' '.join(['tar', '-zxvf', download_path + 'profiles.tar.gz']))

            assert(os.path.isfile(ko_annot_info_file))
            assert(os.path.isdir(download_path + 'profiles/'))

            if os.path.isfile(ko_phmm_file):
                os.system(' '.join(['rm', '-f', ko_phmm_file]))
            for f in os.listdir(download_path + 'profiles/'):
                if not f.endswith('.hmm'): continue
                os.system(' '.join(['cat', download_path + 'profiles/' + f, '>>', ko_phmm_file]))
            assert(os.path.isfile(ko_phmm_file))
            os.system(' '.join(['hmmpress', ko_phmm_file]))

            ko_descriptions_file = download_path + 'ko_descriptions.txt'
            kdf_handle = open(ko_descriptions_file, 'w')
            with open(ko_annot_info_file) as oaif:
                for i, line in enumerate(oaif):
                    if i == 0: continue
                    line = line.strip()
                    ls = line.split('\t')
                    ko = ls[0]
                    if ls[1] == '-': continue
                    description = ls[-1]
                    kdf_handle.write(ko + '\t' + description + '\n')
            kdf_handle.close()
            listing_handle.write('ko\t' + ko_descriptions_file + '\t' + ko_phmm_file + '\n')
            os.system(' '.join(['rm', '-rf', download_path + 'profiles/', download_path + 'profiles.tar.gz', ko_annot_info_file]))
        except Exception as e:
            sys.stderr.write('Issues setting up KO database.\n')
            sys.stderr.write(str(e) + '\n')
            sys.exit(1)

        try:
            print('Setting up PGAP database!')
            os.system(' '.join(['tar', '-zxvf', 'hmm_PGAP.HMM.tgz']))
            assert(os.path.isfile(pgap_info_file))
            assert(os.path.isdir(download_path + 'hmm_PGAP/'))
            for f in os.listdir(download_path + 'hmm_PGAP/'):
                os.system(' '.join(['cat', download_path + 'hmm_PGAP/' + f, '>>', pgap_phmm_file]))
            pgap_descriptions_file = download_path + 'pgap_descriptions.txt'
            pdf_handle = open(pgap_descriptions_file, 'w')
            with open(pgap_info_file) as opil:
                for i, line in enumerate(opil):
                    if i == 0: continue
                    line = line.strip()
                    ls = line.split('\t')
                    label = ls[2]
                    description = ls[10]
                    pdf_handle.write(label + '\t' + description + '\n')
            pdf_handle.close()
            assert(os.path.isfile(pgap_phmm_file))
            os.system(' '.join(['hmmpress', pgap_phmm_file]))
            listing_handle.write('pgap\t' + pgap_descriptions_file + '\t' + pgap_phmm_file + '\n')
            os.system(' '.join(['rm', '-rf', download_path + 'hmm_PGAP/', download_path + 'hmm_PGAP.HMM.tgz', pgap_info_file]))
        except Exception as e:
            sys.stderr.write('Issues setting up PGAP database.\n')
            sys.stderr.write(str(e) + '\n')
            sys.exit(1)

        try:
            print('Setting up MI-BiG database!')
            os.system(' '.join(['diamond', 'makedb', '--in', 'mibig_prot_seqs_3.1.fasta', '-d', mb_faa_file]))
            assert(os.path.isfile(mb_faa_file))
            mibig_descriptions_file = download_path + 'mibig_descriptions.txt'
            mdf_handle = open(mibig_descriptions_file, 'w')
            with open(download_path + 'mibig_prot_seqs_3.1.fasta') as omf:
                for rec in SeqIO.parse(omf, 'fasta'):
                    mdf_handle.write(rec.id + '\t' + rec.description + '\n')
            mdf_handle.close()
            listing_handle.write('mibig\t' + mibig_descriptions_file + '\t' + mb_faa_file + '\n')
            os.system(' '.join(['rm', '-rf', download_path + 'mibig_prot_seqs_3.1.fasta']))
        except Exception as e:
            sys.stderr.write('Issues setting up MI-BiG database.\n')
            sys.stderr.write(str(e) + '\n')
            sys.exit(1)

        try:
            print('Setting up VFDB database!')
            os.system(' '.join(['gunzip', 'VFDB_setB_pro.fas.gz']))
            os.system(' '.join(['diamond', 'makedb', '--in', 'VFDB_setB_pro.fas', '-d', vfdb_faa_file]))
            assert(os.path.isfile(mb_faa_file))
            vfdb_descriptions_file = download_path + 'vfdb_descriptions.txt'
            vdf_handle = open(vfdb_descriptions_file, 'w')
            with open(download_path + 'VFDB_setB_pro.fas') as ovf:
                for rec in SeqIO.parse(ovf, 'fasta'):
                    vdf_handle.write(rec.id + '\t' + rec.description + '\n')
            vdf_handle.close()
            listing_handle.write('vfdb\t' + vfdb_descriptions_file + '\t' + vfdb_faa_file + '\n')
            os.system(' '.join(['rm', '-rf', download_path + 'VFDB_setB_pro.fas']))
        except Exception as e:
            sys.stderr.write('Issues setting up VFDB database.\n')
            sys.stderr.write(str(e) + '\n')
            sys.exit(1)

        # have note in program to recommend folks check out ARTS webserver for more detailed analysis
        # in finding antibiotic BGCs
        try:
            print('Setting up CARD database!')
            os.mkdir(download_path + 'CARD_DB_Files/')
            os.system(' '.join(['tar', '-xf', download_path + 'card-data.tar.bz2', '-C', download_path + 'CARD_DB_Files/']))
            os.system(' '.join(['mv', download_path + 'CARD_DB_Files/protein_fasta_protein_homolog_model.fasta', download_path]))
            os.system(' '.join(['diamond', 'makedb', '--in', download_path + 'protein_fasta_protein_homolog_model.fasta', '-d', card_faa_file]))
            card_descriptions_file = download_path + 'card_descriptions.txt'
            cdf_handle = open(card_descriptions_file, 'w')
            with open(download_path + 'protein_fasta_protein_homolog_model.fasta') as ocf:
                for rec in SeqIO.parse(ocf, 'fasta'):
                    cdf_handle.write(rec.id + '\t' + rec.description + '\n')
            cdf_handle.close()
            os.system(' '.join(['rm', '-rf', download_path + 'CARD_DB_Files/', download_path + 'protein_fasta_protein_homolog_model.fasta', download_path + 'card-data.tar.bz2']))
            assert(os.path.isfile(card_faa_file))
            listing_handle.write('card\t' + card_descriptions_file + '\t' + card_faa_file + '\n')
        except Exception as e:
            sys.stderr.write('Issues setting up CARD database.\n')
            sys.stderr.write(str(e) + '\n')
            sys.exit(1)

        try:
            print('Setting up VOG database!')
            os.mkdir(download_path + 'VOG_HMM_Files/')
            os.system(' '.join(['tar', '-xzvf', 'vog.hmm.tar.gz', '-C', download_path + 'VOG_HMM_Files/']))
            os.system(' '.join(['gunzip', download_path + 'vog.annotations.tsv.gz']))
            vog_db_dir = download_path + 'VOG_HMM_Files/'
            for f in os.listdir(vog_db_dir):
                os.system('cat %s >> %s' % (vog_db_dir + f, vog_phmm_file))
            vog_descriptions_file = download_path + 'vog_descriptions.txt'
            vdf_handle = open(vog_descriptions_file, 'w')
            with open(download_path + 'vog.annotations.tsv') as ovf:
                for i, line in enumerate(ovf):
                    line = line.rstrip('\n')
                    ls = line.split('\t')
                    if i == 0: continue
                    vdf_handle.write(ls[0] + '\t' + ls[4] + '[' + ls[3] + ']\n')
            vdf_handle.close()
            os.system(' '.join(['rm', '-rf', download_path + 'VOG_HMM_Files/', vog_info_file]))
            assert(os.path.isfile(vog_phmm_file))
            assert(os.path.isfile(vog_descriptions_file))
            os.system(' '.join(['hmmpress', vog_phmm_file]))
            listing_handle.write('vog\t' + vog_descriptions_file + '\t' + vog_phmm_file + '\n')
        except Exception as e:
            sys.stderr.write('Issues setting up VOG database.\n')
            sys.stderr.write(str(e) + '\n')
            sys.exit(1)

        try:
            print('Setting up PaperBlast database!')
            pb_faa_path = download_path + 'uniq.faa'
            os.system(' '.join(['diamond', 'makedb', '--in', pb_faa_path, '-d', pb_faa_file]))
            paperblast_descriptions_file = download_path + 'paperblast_descriptions.txt'
            pbdf_handle = open(paperblast_descriptions_file, 'w')
            with open(download_path + 'uniq.faa') as opdf:
                for rec in SeqIO.parse(opdf, 'fasta'):
                    pbdf_handle.write(rec.id + '\t' + rec.description + '\n')
            pbdf_handle.close()
            os.system(' '.join(['rm', '-rf', pb_faa_path]))
            assert(os.path.isfile(paperblast_descriptions_file))
            assert(os.path.isfile(pb_faa_file))
            listing_handle.write('paperblast\t' + paperblast_descriptions_file + '\t' + pb_faa_file + '\n')

        except Exception as e:
            sys.stderr.write('Issues setting up PaperBlast database.\n')
            sys.stderr.write(str(e) + '\n')
            sys.exit(1)

        try:
            print('Setting up ISFinder database!')
            is_faa_path = download_path + 'IS.faa'
            os.system(' '.join(['diamond', 'makedb', '--in', is_faa_path, '-d', is_faa_file]))
            is_descriptions_file = download_path + 'is_descriptions.txt'
            idf_handle = open(is_descriptions_file, 'w')
            with open(is_faa_path) as oif:
                for rec in SeqIO.parse(oif, 'fasta'):
                    idf_handle.write(rec.id + '\t' + rec.description + '\n')
            idf_handle.close()
            assert(os.path.isfile(is_faa_path))
            assert(os.path.isfile(is_descriptions_file))
            listing_handle.write('isfinder\t' + is_descriptions_file + '\t' + is_faa_file + '\n')
            os.system(' '.join(['rm', '-rf', is_faa_path]))
        except Exception as e:
            sys.stderr.write('Issues setting up ISFinder database.\n')
            sys.stderr.write(str(e) + '\n')
            sys.exit(1)

    listing_handle.close()

    print("Done setting up annotation databases!")
    print("Information on final files used by zol can be found at:\n%s" % listing_file)
    sys.exit(0)
setup_annot_dbs()
