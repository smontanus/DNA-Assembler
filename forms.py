from flask_wtf import Form
from wtforms import SubmitField, FileField, HiddenField, BooleanField


class FastQSubmit(Form):
    fqfile = FileField('FASTQ File')
    submitButton = SubmitField('Upload & Analyze')
