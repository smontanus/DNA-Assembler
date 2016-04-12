#!/usr/bin/env python

"""
Flask microframework methods. Provides views that serve as web service endpoints and
the methods associated with them.
"""

# Utilizes the Flask library for web services and the Boto library for AWS integration.
import os
from flask import Flask, flash, render_template, redirect, url_for, request, send_from_directory
from flask_jsglue import JSGlue
from forms import FastQSubmit
from werkzeug.utils import secure_filename
from boto.s3.connection import S3Connection
import dna_assembler

BASE_DIR = os.path.dirname(__file__)
UPLOAD_FOLDER = os.path.normpath(os.path.join(BASE_DIR, 'data/upload'))  # Local destination for storage of FASTQ files.
DOWNLOAD_FOLDER = os.path.normpath(os.path.join(BASE_DIR, 'data/download'))  # Local destination for storage of FASTA files.
AWS_KEY = ""  # Must include AWS Key ID.
AWS_SECRET = ""  # Must include AWS Secret Key.
KEY_BASE_UPLOAD = "upload/"  # AWS S3 destination for storage of FASTQ files.
KEY_BASE_DOWNLOAD = "download/"  # AWS S3 destination for storage of FASTA files.

application = Flask(__name__)
jsglue = JSGlue(application)
application.config.from_object('config') # Access data in config.py
application.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
application.config['DOWNLOAD_FOLDER'] = DOWNLOAD_FOLDER


#  ####### VIEWS #########
#  #######################

@application.route('/')
@application.route('/index')
def index():
    """
    Flask microframework root endpoint. Redirects navigation
    to the fileinput endpoint.
    """
    form = FastQSubmit()
    return redirect(url_for('fileinput',
        title='DNA Assembly - Input FASTQ File',
        form=form))


@application.route('/fileinput', methods = ['GET', 'POST'])
def fileinput():
    """
    Flask microframework endpoint that provides the page for
    selecting and uploading a FASTQ file. Stores the file in
    AWS S3 and redirects navigation to the assembly endpoint
    with correct file upload.
    """
    form = FastQSubmit()
    if request.method == 'POST':
        try:
            fqfile = request.files['fqfile']  # Extract file contents from request.
        except Exception as e:
            render_template("500_error.html", error = str(e))
        else:
            if fqfile and dna_assembler.allowed_file(fqfile.filename):
                filename = secure_filename(fqfile.filename)  # Create secure file name.
                ul_fullpath = os.path.join(application.config['UPLOAD_FOLDER'], filename)
                fqfile.save(ul_fullpath)
                file_valid, seqs_valid, reads_valid = dna_assembler.validfastq(ul_fullpath)
                if not file_valid:
                    os.remove(ul_fullpath)
                    msgtxt = "The file you have submitted is not in the FASTQ format or is corrupt. Please check the file and try again."
                    flash(msgtxt, 'error')
                    return redirect(url_for('index'))
                if not seqs_valid:
                    os.remove(ul_fullpath)
                    msgtxt = "The submitted FASTQ file contains sequence reads that exceed the maximum length allowed by the program. Please select a FASTQ file with reads shorter than 300 nucleotides."
                    flash(msgtxt, 'error')
                    return redirect(url_for('index'))
                if not reads_valid:
                    os.remove(ul_fullpath)
                    msgtxt = "The submitted FASTQ file contains more reads than allowed by the program. Please select a FASTQ file with a maximum of 500 reads."
                    flash(msgtxt, 'error')
                    return redirect(url_for('index'))

                # ##### AWS S3 CONNECTION - Write FASTQ File ######
                conn = S3Connection(AWS_KEY, AWS_SECRET)  # Establish AWS S3 connection.
                bucket = conn.get_bucket('dna-assembly-data')  # Select S3 storage bucket.
                keyname = KEY_BASE_UPLOAD + filename  # Create key name for S3 file storage.
                key = bucket.new_key(keyname)  # Create empty file in S3 for storage.
                with open(ul_fullpath, 'r') as filehandle:
                    fData=filehandle.read()
                key.set_contents_from_string(fData)  # Write FASTQ file data to S3 storage.
                conn.close()  # Close S3 connection.
                # ##### END AWS S3 CONNECTION ######################

                return redirect(url_for('analyze',
                    title='DNA Assembly - Analysis Results',
                    fqFile=filename))
    return render_template('fileinput.html',
        title='DNA Assembly - Input FASTQ File',
        form=form)


@application.route('/analyze&<fqFile>')
def analyze(fqFile):
    """
    Flask microframework endpoint that provides the page for
    displaying FASTQ file contents and analysis. Also provides
    selection of assembly option.
    """
    filename = os.path.join(application.config['UPLOAD_FOLDER'], fqFile)
    with open(filename, 'r') as filehandle:  # Dump FASTQ file contents to string for display.
        fqData=filehandle.read()
    reads = dna_assembler.readfastq(filename)  # Convert FASTQ data to reads for processing.
    num, alen = dna_assembler.analyzereads(reads)  # Simple file contents analysis.
    return render_template('results.html',
        title='DNA Assembly - Analysis Results',
        num=num,
        alen=alen,
        file=fqFile,
        data=fqData)


@application.route('/assemble&<fqFile>&<blast>')
def assemble(fqFile, blast):
    """
    Flask microframework endpoint that provides the page for
    displaying assembly and Blast alignment results. Saves the
    assembled FASTA file AWS S3 storage and provides a link for
    downloading the assembled FASTA file.
    """
    filename = os.path.join(application.config['UPLOAD_FOLDER'], fqFile)
    reads = dna_assembler.readfastq(filename)  # Convert FASTQ data to reads for processing.
    k = 30  # Minimum length of read overlap for assembly.
    sequence = dna_assembler.greedy_scs(reads, k)  # Perform the assembly of the reads.
    dl_file = fqFile.split('.', 1 )[0] + ".FASTA"  # Build the FASTA file name from the FASTQ file name.
    dl_fullpath = os.path.join(application.config['DOWNLOAD_FOLDER'], dl_file)
    dna_assembler.writefasta(dl_fullpath, sequence)

    # ##### AWS S3 CONNECTION - Write FASTA File ######
    conn = S3Connection(AWS_KEY, AWS_SECRET)  # Establish AWS S3 connection.
    bucket = conn.get_bucket('dna-assembly-data')  # Select S3 storage bucket.
    keyname = KEY_BASE_DOWNLOAD + dl_file  # Create key name for S3 file storage.
    key = bucket.new_key(keyname)  # Create empty file in S3 for storage.
    with open(dl_fullpath, 'r') as filehandle:
        fData=filehandle.read()
    key.set_contents_from_string(fData)  # Write FASTA file data to S3 storage.
    conn.close()  # Close S3 connection.
    # ##### END AWS S3 CONNECTION ######################

    if blast == "True":
        e = 0.01  # Blast alignment filter value.
        blast_results = dna_assembler.searchBlast(dl_fullpath, e)  # Perform blast alignment.
        search = blast_results[0][0]  # Extract species data from results.
        species = search.split("|")[-1] # Extract species data from results.
    else:
        blast_results = "BLAST Alignment Not Performed"
        species = "BLAST Alignment Not Performed"
    return render_template('assembly.html',
        title='DNA Assembly - Complete',
        file=fqFile,
        seq=sequence,
        blast_r=blast_results,
        species=species,
        dfile=dl_file,
        dlink=dl_fullpath)


@application.route('/data/downloads/<path:filename>')
def download_file(filename):
    """
    Create the download link for the assembled FASTA file.
    """
    return send_from_directory(application.config['DOWNLOAD_FOLDER'], filename)


@application.errorhandler(500)
def server_error(error):
    """
    Flask endpoint for 500 error responses.
    """
    return render_template('500_error.html'), 500


@application.errorhandler(400)
def page_not_found(error):
    """
    Flask endpoint for 400 error responses.
    """
    return render_template('400_error.html'), 400

@application.errorhandler(Exception)
def unhandled_exception(e):
    """
    Flask endpoint for unhandled application errors.
    """
    return render_template('500_error.html'), 500

#  ###### END VIEWS #######
#  ########################
                
if __name__ == '__main__':
    application.run(debug=True)