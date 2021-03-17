from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField, SelectField
from wtforms.validators import DataRequired

class FindSimilarsForm(FlaskForm):
    dataset_choices=[
        (1,"eMolecules 2021-03 complete"),
        (2,"eMolecules 2021-03 druglike score > 0.9"),
    ]


    smiles=StringField('SMILES', validators=[DataRequired()])
    select=SelectField('Choose a dataset to query', choices=dataset_choices)
    submit=SubmitField('Find 3D similars')

