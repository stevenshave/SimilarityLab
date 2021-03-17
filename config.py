import os

class Config(object):
    SECRET_KEY = os.environ.get('FLASK_SECRET_KEY') or '3Sy63TSRyVBaq6G37Ch7lXYMLQPuYgkUfRAN1Av4ltc'
