from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField, SelectField
from wtforms.validators import DataRequired

class FindSimilarsForm(FlaskForm):
    smiles=StringField('SMILES', validators=[DataRequired()])
    select=SelectField('Choose a dataset to query', choices=[], coerce=int)
    select_n_to_keep=SelectField('Numer of similars to keep', choices=[], coerce=int)
    submit=SubmitField('Find 3D similars')

class PredictTargetsForm(FlaskForm):
    smiles=StringField('SMILES', validators=[DataRequired()])
    submit=SubmitField('Predict targets')

class PredictLogP(FlaskForm):
    smiles=StringField('SMILES', validators=[DataRequired()])
    submit=SubmitField('Predict logP')