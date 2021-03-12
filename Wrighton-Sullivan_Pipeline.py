#!/usr/bin/env python

"""
#PBS -N my_job
#PBS -j oe
#PBS -m abe
set -x
trap "cd $PBS_O_WORKDIR;mkdir $PBS_JOBID;cp -R $TMPDIR/* $PBS_JOBID" TERM
cd $PBS_O_WORKDIR
cp $ff $TMPDIR
cd $TMPDIR
/users/PAS1018/osu9681/bin/pullseq.py -i $ff -m $len -o contigs_$len.fa
/users/PAS1018/osu9681/bin/prodigal -i contigs_$len.fa -o contigs_$len.genes -a contigs_$len.genes.faa -d contigs_$len.genes.fna -p meta -m
/users/PAS1018/osu9681/bin/ANNOTATION_PIPELINE_IPER_OPTION_OSC.sh contigs_$len.genes.faa contigs_$len.genes.fna NOT_USED ORIG $iper 12
cd $PBS_O_WORKDIR
mkdir $PBS_JOBID
cp -R $TMPDIR/* $PBS_JOBID

----- SCRIPTS STILL REQUIRED -----
"/users/PAS1018/osu9681/bin/convert_pfam_to_iperscan.py"
"/users/PAS1018/osu9681/bin/perl4_NEW.pl"
-----

prodigal
hmmer3
pfam_scan = bioperl
usearch

"""
import sys
import os
import subprocess
import multiprocessing as mp
import psutil
import shutil
import argparse
from Bio import SeqIO
import pandas as pd
from pprint import pprint

parser = argparse.ArgumentParser(description="Wrighton pipeline with Sullivan improvements.")

inputs = parser.add_argument_group('Input')

inputs.add_argument('-i', '--input-fasta', dest='fasta_fp', metavar='FILEPATH',
                    help="'Full' path of file to annotate")

inputs.add_argument('-t', '--na-type', type=str, choices=['nucl', 'prot'], default='nucl', dest='na_type',
                    help='Which to use.')

inputs.add_argument('-m', '--min-size', dest='min_size', metavar='INT', type=int, default=5000,
                    help='Minimum sequence length size. Only applies to nucleotide input sequences.')

inputs.add_argument('-o', '--output-dir', dest='output_fp', metavar='FILEPATH',
                    help="Output directory")

options = parser.add_argument_group('Annotations')

options.add_argument('--tool', type=str, choices=['pfam', 'iperscan'], default='pfam', dest='tool',
                     help='Which to use.')

results = parser.parse_args()


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def distributed_chunks(l, n):
    """
    http://stackoverflow.com/a/6856593
    Splits list l into n chunks with approximately equals sum of values
    see http://stackoverflow.com/questions/6855394/splitting-list-in-chunks-of-balanced-weight
    """
    result = [[] for i in range(n)]
    sums = {i: 0 for i in range(n)}
    c = 0
    for e in l:
        for i in sums:
            if c == sums[i]:
                result[i].append(e)
                break
        sums[i] += e
        c = min(sums.values())
    return result


def run(cmd):
    sys.stdout.write('Executing: {}\n'.format(cmd))
    subprocess.check_call(cmd, shell=True)


def build_db_dicts(db_header_fp, separator=' '):
    print('Processing {}'.format(os.path.basename(db_header_fp)))
    db_header_dict = {}
    with open(db_header_fp, 'rU') as db_header_fh:
        for lines in db_header_fh:
            key, value = lines.strip().split(separator, 1)
            db_header_dict[key] = value

    return db_header_dict


def process_db_dfs(search_fp, bitscore, map_dict):
    """
    Import raw output from usearch, format, filter 'excess', filter on bitscore, add DB-annotation
    @param search_fp:
    @param bitscore:
    @param map_dict:
    @return:
    """
    search_df = pd.read_csv(search_fp, delimiter='\t', header=None, index_col=False, dtype='object')

    # Remove extra ORF annotation information (start, stop, strand, etc...)
    search_df[0] = search_df[0].str.split(' ', 1).str[0]

    search_df[11] = pd.to_numeric(search_df[11])  # Pandas tries to be smart and import, but those int "roundings"...
    search_df = search_df[search_df[11] >= bitscore]

    search_df[12] = search_df[1].map(map_dict)

    return search_df


def process_rbh(forward_fp, reverse_fp, map_dict):
    forward_df = pd.read_csv(forward_fp, delimiter='\t', header=None, index_col=0, dtype='object')
    forward_df[11] = pd.to_numeric(forward_df[11])

    reverse_df = pd.read_csv(reverse_fp, delimiter='\t', header=None, index_col=0, dtype='object')
    reverse_df[11] = pd.to_numeric(reverse_df[11])

    rbh_dict = {}

    for index, series in forward_df.iterrows():  # Go through each match
        target, bitscore = series[1], series[11]  # Grab

        try:
            reverse_target = reverse_df.loc[target, 1]  # Get "reciprocal" target's match
            reverse_bitscore = reverse_df.loc[target, 11]  # Get "reciprocal" target's match

            if index == reverse_target:
                rbh_dict[index] = {
                    'target': target,  # Not sure I like using target here, it'd 'feel' better to have reverse_df
                    'bitscore': reverse_bitscore
                }

        except KeyError:
            continue

    rbh_df = pd.DataFrame.from_dict(rbh_dict, orient='index')
    rbh_df['desc'] = rbh_df['target'].map(map_dict)

    # NOT a fan of dropping data, but later perl/python scripts (not mine) can't handle additional columns
    rbh_df.drop('bitscore', axis=1, inplace=True)

    return rbh_df


def annotate_sequence(input_sequence_fp, annotation_df):
    records_out = []
    unknown_headers = []
    with open(input_sequence_fp, 'rU') as input_sequence_fh:
        for record in SeqIO.parse(input_sequence_fh, 'fasta'):

            record.description = ''  # Junk from prodigal
            try:
                annotations = annotation_df.loc[record.id]

                if annotations['rank'] == 'A':
                    record.id = record.id + '_{}'.format(annotations['RBH:KEGG'])
                if annotations['rank'] == 'B':
                    record.id = record.id + '_{}'.format(annotations['RBH:UNIREF'])
                if annotations['rank'] == 'C':
                    record.id = record.id + '_{}'.format(annotations['BLAST:KEGG'])
                if annotations['rank'] == 'D':
                    record.id = record.id + '_{}'.format(annotations['BLAST:UNIREF'])
                if annotations['rank'] == 'E':
                    record.id = record.id + '_{}'.format(annotations['IPRSCAN'])
            except KeyError:
                unknown_headers.append(record.id)  # Don't want "unknown_function" added to keys of dataframe
                record.id = record.id + '_Unknown_Function'

            records_out.append(record)

    # Should probably direct it to a user-defined output file
    prodigal_annotations_fp = input_sequence_fp.rsplit('.', 1)[0] + '.annotated.' + input_sequence_fp.rsplit('.', 1)[1]
    with open(prodigal_annotations_fp, 'w') as prodigal_annotations_fh:
        SeqIO.write(records_out, prodigal_annotations_fh, 'fasta')

    return prodigal_annotations_fp, unknown_headers


if __name__ == '__main__':

    # Need to ensure everything's working
    which_programs = ['prodigal', 'hmmscan', 'hmmalign', 'convert_pfam_to_iperscan.py']  # 'usearch' available linked
    for program in which_programs:
        if not program:
            print('Unable to find {}, are you sure it\'s available in your PATH?'.format(program))
            sys.exit(1)
    cpu_count = psutil.cpu_count(logical=False)

    # Gather only real "input"
    fasta_fp = results.fasta_fp
    fasta_fn = os.path.basename(fasta_fp)
    fasta_bn = fasta_fn.rsplit('.', 1)[0]
    min_len = results.min_size

    output_dir = results.output_fp
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    os.chdir(output_dir)

    if results.na_type == 'nucl':

        # Filter sequence length(s)
        sys.stdout.write('Checking nucleotide fasta file for minimum sequence length ({})\n'.format(min_len))
        len_filtered_contigs_fn = '{}_{}.fa'.format(fasta_bn, min_len)
        len_filtered_contigs_fp = os.path.join(output_dir, len_filtered_contigs_fn)
        if not os.path.exists(len_filtered_contigs_fp):

            # Get minimum sequence lengths
            len_filtered_contigs = []
            with open(fasta_fp, 'rU') as fasta_fh:
                for records in SeqIO.parse(fasta_fh, 'fasta'):
                    if len(records.seq) >= results.min_size:
                        len_filtered_contigs.append(records)

            # Write length-filtered sequences to file
            with open(len_filtered_contigs_fp, 'w') as len_filtered_contigs_fh:
                wrt = SeqIO.write(len_filtered_contigs, len_filtered_contigs_fh, 'fasta')
                sys.stdout.write('Wrote {} records to {}\n'.format(wrt, len_filtered_contigs_fp))

            # Why here? Need it to apply back to main basename
            fasta_bn = len_filtered_contigs_fn.rsplit('.', 1)[0]
            fasta_fp = len_filtered_contigs_fp

        # Run prodigal
        prodigal_bn = fasta_bn  # As above, now can keep track of what's called what
        prodigal_genes_fn = '{}.genes'.format(prodigal_bn)
        prodigal_genes_fp = os.path.join(output_dir, prodigal_genes_fn)
        prodigal_faa_fn = '{}.faa'.format(prodigal_bn)
        prodigal_faa_fp = os.path.join(output_dir, prodigal_faa_fn)
        prodigal_fna_fn = '{}.fna'.format(prodigal_bn)
        prodigal_fna_fp = os.path.join(output_dir, prodigal_fna_fn)

        if not os.path.exists(prodigal_faa_fp):

            prodigal_cmd = 'prodigal -i {} -o {} -a {} -d {} -p meta -m'.format(
                len_filtered_contigs_fp, prodigal_genes_fp, prodigal_faa_fp, prodigal_fna_fp)
            run(prodigal_cmd)

            fasta_fp = prodigal_faa_fp  # Now a user can select nucl or prot and both will be applied to fasta_fp

    elif results.na_type == 'prot':
        print('Identified existing prodigal/amino-acid file: {}'.format(fasta_fp))
    else:
        print('Unable to identify nucleotide or amino-acid')
        sys.exit(1)


    """
    Interproscan block in original script 1) determined number of lines & sequences, created 2-line versions of each FASTA
    sequence, wrote those to X number of files to distribute them evenly on a multi-threaded system, then ran each as a job.
    Let's not do that.
    """

    # This WHOLE SECTION (150+ lines) deals with running either interproscan or pfam_scan. If neither is run, then the
    # original protein sequence must be split up and processed individually for max speed
    iperscan_or_pfam_fp = fasta_fp.replace('.faa', '.scan.txt')
    if not os.path.exists(iperscan_or_pfam_fp):

        # Interproscan is not natively multi-threaded, so ad-hoc it
        with open(fasta_fp, 'rU') as prodigal_faa_fh:
            records = [record for record in SeqIO.parse(prodigal_faa_fh, 'fasta')]

        if len(records) % cpu_count == 0:
            chunk_size = len(records) // cpu_count
        else:
            chunk_size = (len(records) + 1) // cpu_count

        chunk_group_fps = []
        sys.stdout.write('Splitting {} into {} chunks\n'.format(len(records), chunk_size))
        for i, group in enumerate(chunks(records, chunk_size)): # SeqIO.parse returns iterator
            # sys.stdout.write('\rProcessing batch group: {}'.format(i))
            # sys.stdout.flush()
            output_fn = '{}.{}_{}.faa'.format(fasta_bn, i, chunk_size)
            output_fp = os.path.join(output_dir, output_fn)
            with open(output_fp, 'w') as output_fh:
                SeqIO.write(group, output_fh, 'fasta')
            chunk_group_fps.append(output_fp)

        sys.stdout.write('Splitting complete. Running {} on each chunk.\n'.format(results.tool))

        if results.tool == 'iperscan':  # Run interproscan

            """
            Not quite sure why there's 2 tools here, interproscan and pfam scan. Original code had one OR the other,
            and while  it's true that they can somewhat cover the same function, why the need for both?
            """

            if not os.path.exists(iperscan_or_pfam_fp):

                pool = mp.Pool(processes=cpu_count)

                interpro = 'interproscan.sh'
                interpro_apps = '--applications TIGRFAM,Pfam,ProSiteProfiles,ProSitePatterns'
                interpro_opts = '--formats TSV --disable-precalc -goterms'  # -goterms implies iprlookup
                for chunk_group_fp in chunk_group_fps:
                    interpro_cmd = interpro
                    interpro_out = chunk_group_fp.replace('.faa', '.iprscan')
                    interpro_log = chunk_group_fp.replace('.faa', '.iprscan.log')
                    interpro_cmd += ' -i {} -o {} {} {}'.format(
                        chunk_group_fp, interpro_out, interpro_apps, interpro_opts)
                    interpro_cmd += ' > {}'.format(interpro_log)

                    pool.apply_async(run, args=(interpro_cmd,))

                pool.close()
                pool.join()

                """
                A lot of the original code was devoted (~50 bash lines) to ensuring interproscan ran on 1) all the
                sequences 2) didn't return too many sequences 3) didn't return too few sequences... and counted input
                sequences + lines and output sequences and + lines.

                Does interproscan a problem running? Let's *assume* not and get there if it becomes necessary
                """

                # Merge results of interproscan
                interpro_group_fps = [fp.replace('.faa', '.iprscan') for fp in chunk_group_fps]

                # Apparently you can't cat that many files together...
                with open(iperscan_or_pfam_fp, 'wb') as iperscan_or_pfam_fh:

                    for fn in interpro_group_fps:
                        shutil.copyfileobj(open(fn, 'rb'), iperscan_or_pfam_fh)
                        os.remove(fn)

        if results.tool == 'pfam':

            """
            For the moment will need to leverage the libraries already on OSC until suitable replacements can be created

            DOUBLE-CHECK. Does Pfam scan require EACH FASTA FILE to be a singular sequence?
            """

            # The 2 outputs from the next few blocks of code
            pfamscan_merged_fn = '{}.Pfam.txt'.format(fasta_bn)
            pfamscan_merged_fp = os.path.join(output_dir, pfamscan_merged_fn)

            if not os.path.exists(pfamscan_merged_fp):

                pfam_scan = 'perl -I/fs/project/PAS1117/modules/perl/lib/perl5 /fs/project/PAS1117/modules/PfamScan/pfam_scan.pl'
                pfam_dir = '/fs/project/PAS1117/modules/hmmer3/databases/v31'

                pool = mp.Pool(processes=cpu_count)

                # So if this is cpu 1, why was this ever split into individuals?
                pfamscan_opts = '-dir {} -cpu 1'.format(pfam_dir)
                for chunk_group_fp in chunk_group_fps:
                    pfamscan_out = chunk_group_fp.replace('.faa', '.Pfam.txt')

                    pfamscan_cmd = '{} -fasta {}'.format(pfam_scan, chunk_group_fp)  # Tool and Input
                    pfamscan_cmd += ' {} -outfile {}'.format(pfamscan_opts, pfamscan_out)  # Options and Outfile

                    pool.apply_async(run, args=(pfamscan_cmd,))

                pool.close()
                pool.join()

                # Collect and merge results of pfam scan
                pfamscan_group_fps = [fp.replace('.faa', '.Pfam.txt') for fp in chunk_group_fps]

                # Apparently you can't cat that many files together...
                with open(pfamscan_merged_fp, 'wb') as interpro_merged_fh:

                    for fn in pfamscan_group_fps:
                        shutil.copyfileobj(open(fn, 'rb'), interpro_merged_fh)
                        os.remove(fn)

            """
            Also leveraging existing scripts until replacements can be generated.
            """

            if not os.path.exists(iperscan_or_pfam_fp):

                # Convert pfam output to useable format
                convert_pfam = "/fs/project/PAS1117/modules/bin/convert_pfam_to_iperscan.py"
                pfam_hmm = '/fs/project/PAS1117/modules/hmmer3/databases/v31/Pfam-A.hmm.dat'

                # Basically dropping columns and "zeroing" out others to force format
                make_iper_cmd = 'python {} -i {} -o {} -p {}'.format(
                    convert_pfam, pfamscan_merged_fp, iperscan_or_pfam_fp, pfam_hmm)
                run(make_iper_cmd)

        sys.stdout.write('Removing temporary chunks.\n')
        for chunk_group_fp in chunk_group_fps:
            os.remove(chunk_group_fp)

    # Reverse best BLAST

    usearch = '/users/PAS1018/osu9681/bin/usearch'
    uniref90_db = '/users/PAS1018/osu9681/bin/uniref90.udb'
    uniref90_fna = '/users/PAS1018/osu9681/bin/uniref90.fasta'
    uniref90_headers = '/users/PAS1018/osu9681/bin/uniref90_HEADERS.fasta'
    kegg_db = '/users/PAS1018/osu9681/bin/kegg-all-orgs_12102013.pep.udb'
    kegg_fna = '/users/PAS1018/osu9681/bin/kegg-all-orgs_12102013.pep'
    kegg_headers = '/users/PAS1018/osu9681/bin/kegg-all-orgs_12102013_HEADERS_1.pep'
    universal_opt = '-maxhits 1 -evalue 0.001 -threads {}'.format(cpu_count)

    # Build DB first so everything can be automated
    # Sequences need their database!
    seq_uclust_db = '{}.udb'.format(fasta_bn)
    if not os.path.exists(seq_uclust_db):  # Fasta_fp is either 1) new aa from fna, or 2) aa file directly
        seq_uclust_db_cmd = '{} -makeudb_ublast {} -output {}'.format(usearch, fasta_fp, seq_uclust_db)
        run(seq_uclust_db_cmd)

    usearch_template = '{} -ublast {} -db {} {} -blast6out {}'

    queries = [fasta_fp, fasta_fp, uniref90_fna, kegg_fna]
    targets = [uniref90_db, kegg_db, seq_uclust_db, seq_uclust_db]
    # Might as well set this up now, because (at this point) can't see another way to constantly re-reference them
    query_vs_uniref = os.path.join(output_dir, '{}_ublast_uniref90.b6'.format(fasta_bn))
    query_vs_kegg = os.path.join(output_dir, '{}_ublast_kegg.b6'.format(fasta_bn))
    uniref_vs_query = os.path.join(output_dir, 'uniref_ublast_{}.b6'.format(fasta_bn))
    kegg_vs_query = os.path.join(output_dir, 'kegg_ublast_{}.b6'.format(fasta_bn))

    outputs = [query_vs_uniref, query_vs_kegg, uniref_vs_query, kegg_vs_query]

    for query, target, output in zip(queries, targets, outputs):

        output_fp = os.path.join(output_dir, output)

        if not os.path.exists(output):
            usearch_cmd = usearch_template.format(usearch, query, target, universal_opt, output_fp)
            run(usearch_cmd)
        else:
            if os.stat(output).st_size == 0:
                usearch_cmd = usearch_template.format(usearch, query, target, universal_opt, output_fp)
                run(usearch_cmd)

    # Sequences vs UniRef
    # Sequences vs KEGG
    # Uniref vs Sequences
    # KEGG vs Sequences

    uniref_dict = build_db_dicts(uniref90_headers, separator=' ')
    kegg_dict = build_db_dicts(kegg_headers, separator='  ')

    for query in [query_vs_uniref, query_vs_kegg, uniref_vs_query, kegg_vs_query]:

        if query in [query_vs_uniref, query_vs_kegg]:
            query_out = query.replace('.b6', '.BIT_SCORE60.b6')

            if not os.path.exists(query_out):
                if query == query_vs_uniref:
                    query_df = process_db_dfs(query, bitscore=60, map_dict=uniref_dict)
                else:   # if query == query_vs_kegg:
                    query_df = process_db_dfs(query, bitscore=60, map_dict=kegg_dict)

                query_df.to_csv(query_out, sep='\t', header=None, index_label=False, index=False)

        if query in [uniref_vs_query, kegg_vs_query]:
            query_out = query.replace('.b6', '.BIT_SCORE300.b6')

            if not os.path.exists(query_out):
                if query == query_vs_uniref:
                    query_df = process_db_dfs(query, bitscore=300, map_dict=uniref_dict)
                else:  # if query == query_vs_kegg:
                    query_df = process_db_dfs(query, bitscore=300, map_dict=kegg_dict)

                query_df.to_csv(query_out, sep='\t', header=None, index_label=False, index=False)

    # Get best reciprocal
    uniref_rbh_fn = '{}-UniRef.rbh.txt'.format(fasta_fn)
    uniref_rbh_fp = os.path.join(output_dir, uniref_rbh_fn)
    kegg_rbh_fn = '{}-KEGG.rbh.txt'.format(fasta_fn)
    kegg_rbh_fp = os.path.join(output_dir, kegg_rbh_fn)

    # Double check to ensure that files exist/don't !!!!!!!!!!!!!!!!!!!!

    seq_vs_uniref_rbh_df = process_rbh(query_vs_uniref.replace('.b6', '.BIT_SCORE60.b6'), uniref_vs_query.replace('.b6', '.BIT_SCORE300.b6'), uniref_dict)
    seq_vs_uniref_rbh_df.to_csv(uniref_rbh_fp, sep='\t', header=None, index_label=False)

    seq_vs_kegg_rbh_df = process_rbh(query_vs_kegg.replace('.b6', '.BIT_SCORE60.b6'), kegg_vs_query.replace('.b6', '.BIT_SCORE300.b6'), kegg_dict)
    seq_vs_kegg_rbh_df.to_csv(kegg_rbh_fp, sep='\t', header=None, index_label=False)

    # Okay, at this point you need to generate a unified naming scheme that doesn't involve copy-and-paste

    # Combine all the annotations.... Why aren't we just copying all of them to the same file?
    perl_merge = '/users/PAS1018/osu9681/bin/perl4_NEW.pl'

    merged_annotations_fn = '{}.merged-annotations.txt'.format(fasta_fn)
    merged_annotations_fp = os.path.join(output_dir, merged_annotations_fn)
    merge_annot_cmd = 'perl {0} {1} {2} {3} {4} {5} > {6}'.format(
        perl_merge,
        query_vs_uniref.replace('.b6', '.BIT_SCORE60.b6'),
        query_vs_kegg.replace('.b6', '.BIT_SCORE60.b6'),
        uniref_rbh_fp,
        kegg_rbh_fp,
        iperscan_or_pfam_fp,
        merged_annotations_fp)
    run(merge_annot_cmd)

    """
    Arrange annotations from above output according to database preferences:
        RBH-KEGG, RBH-UNIREF, BLAST-KEGG, BLAST-UNIREF, IPRSCAN
    """

    final_annotations_df = pd.read_csv(merged_annotations_fp, delimiter='\t', header=None, index_col=False)

    class_a_df = final_annotations_df[
        (final_annotations_df[1].str.contains('RBH:')) & (final_annotations_df[1].str.contains('db=KEGG'))]
    class_b_df = final_annotations_df[
        (final_annotations_df[1].str.contains('RBH:')) & (final_annotations_df[1].str.contains('db=UNIREF'))]
    class_c_df = final_annotations_df[
        (final_annotations_df[1].str.contains('BLAST:')) & (final_annotations_df[1].str.contains('db=KEGG'))]
    class_d_df = final_annotations_df[
        (final_annotations_df[1].str.contains('BLAST:')) & (final_annotations_df[1].str.contains('db=UNIREF'))]
    class_e_df = final_annotations_df[final_annotations_df[1].str.contains('IPRSCAN:')]

    ranked_annotations_fn = '{}.merged-annotations-ordered.txt'.format(fasta_fn)
    ranked_annotations_fp = os.path.join(output_dir, merged_annotations_fn)

    all_classes_df = pd.concat([class_a_df, class_b_df, class_c_df, class_d_df, class_e_df])
    all_classes_df.to_csv(ranked_annotations_fp, sep='\t', header=None, index_label=False, index=False)

    # Can't think of a way to do this in pandas, but it's probably easy
    organized_annotations_fn = '{}.merged_annotations_singleLine.txt'.format(fasta_fn)
    organized_annotations_fp = os.path.join(output_dir, organized_annotations_fn)
    singleLineDict = {}
    with open(ranked_annotations_fp, 'rU') as ranked_annotations_fh:
        for lines in ranked_annotations_fh:
            line = lines.strip().split('\t')
            seq_id, hit_data = line[0], line[1]

            class_id = hit_data.split(':', 1)[0]

            # Two classes of RBH
            if 'db=UNIREF' in hit_data:
                class_id = class_id + ':' + 'UNIREF'
            if 'db=KEGG' in hit_data:
                class_id = class_id + ':' + 'KEGG'

            if seq_id not in singleLineDict:
                singleLineDict[seq_id] = {
                    class_id: hit_data
                }
            if class_id not in singleLineDict[seq_id]:
                singleLineDict[seq_id][class_id] = hit_data

    singleLineDF = pd.DataFrame.from_dict(singleLineDict, orient='index')

    # Replaces perl6.pl, which effectively ranked each ORF by how many databases they matched against
    # The perl script assumed the order of the merged annotations was in order of priority, the pandas version (below)
    # does not.
    # Ranking system
    # All databases - A
    # 4 databases(One RBH atleast) - B
    # 3 databases(One BLAST atleast) - C
    # 2 databases(One BLAST atleast) - D
    # 1 database(IPRSCAN only) - E
    for index, series in singleLineDF.iterrows():

        if not pd.isnull(series['RBH:KEGG']):
            singleLineDF.loc[index, 'rank'] = 'A'
        elif not pd.isnull(series['RBH:UNIREF']):
            singleLineDF.loc[index, 'rank'] = 'B'
        elif not pd.isnull(series['BLAST:KEGG']):
            singleLineDF.loc[index, 'rank'] = 'C'
        elif not pd.isnull(series['BLAST:UNIREF']):
            singleLineDF.loc[index, 'rank'] = 'D'
        elif not pd.isnull(series['IPRSCAN']):
            singleLineDF.loc[index, 'rank'] = 'E'
        else:
            sys.stderr.write('Appropriate header not found when parsing singleLine dataframe.')

    # Add annotations to protein
    prodigal_annotation_faa, unknown_proteins = annotate_sequence(fasta_fp, singleLineDF)

    if results.na_type == 'nucl':
        # Add annotations to nucleotide
        prodigal_annotation_fna, unknown_nucleic = annotate_sequence(results.fasta_fp, singleLineDF)

    # Add unknowns to dataframe (and thus extend annotations from sequences to summary file
    append_dict = {key: {'rank': 'F', 'MISC': 'Unknown_function'} for key in unknown_proteins}
    singleLineDF = singleLineDF.append(pd.DataFrame.from_dict(append_dict, orient='index'))

    # Write everything to file
    try:
        singleLineDF = singleLineDF[['rank', 'RBH:KEGG', 'RBH:UNIREF', 'BLAST:KEGG', 'BLAST:UNIREF', 'IPRSCAN', 'MISC']]
    except KeyError:
        singleLineDF = singleLineDF[['rank', 'RBH:KEGG', 'RBH:UNIREF', 'BLAST:KEGG', 'BLAST:UNIREF', 'IPRSCAN']]

    summary_fn = '{}.summary.csv'.format(fasta_fn)
    summary_fp = os.path.join(output_dir, summary_fn)

    singleLineDF.to_csv(summary_fp, sep=',', quotechar='"')

    # Clean up?
