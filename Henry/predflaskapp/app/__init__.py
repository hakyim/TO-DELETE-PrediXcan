import os
import sys
sys.path.append("/usr/local/lib/python2.7/dist-packages")
sys.path.append(".")
from flask.ext.login import LoginManager
#from flask.ext.openid import OpenID
from flask import Flask
from config import basedir, ADMINS, MAIL_SERVER, MAIL_PORT, MAIL_USERNAME, MAIL_PASSWORD 
from flask.ext.sqlalchemy import SQLAlchemy
from werkzeug import secure_filename
from celery import Celery


app = Flask(__name__)
app.config['CELERY_BROKER_URL'] = 'redis://localhost:6379/0'
app.config['CELERY_RESULT_BACKEND'] = 'redis://localhost:6379/0'
celery = Celery(app.name, broker=app.config['CELERY_BROKER_URL'])
print app.name
#app.celery = celery
celery.conf.update(app.config)


app.static_folder = "static/"
app.config.from_object('config')

#celery stuff

db = SQLAlchemy(app)
ALLOWED_EXTENSIONS = set(['txt', 'pdf', 'png', 'jpg', 'jpeg', 'gif','gz'])

lm = LoginManager()
lm.init_app(app)
lm.login_view = 'login'

from app import views, models

#email log
"""
if not app.debug:
	import logging 
	from logging.handlers import SMTPHandler
	credentials = None
	if MAIL_USERNAME or MAIL_PASSWORD:
		credentials = (MAIL_USERNAME, MAIL_PASSWORD)
	mail_handler = SMTPHandler((MAIL_SERVER, MAIL_PORT), 'no-reply@' + MAIL_SERVER, ADMINS, 'microblog failure', credentials)
	mail_handler.setLevel(logging.ERROR)
	app.logger.addHandler(mail_handler)	
"""
"""
#file log
if not app.debug:
    import logging
    from logging.handlers import RotatingFileHandler
    file_handler = RotatingFileHandler('tmp/microblog.log', 'a', 1 * 1024 * 1024, 10)
    file_handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'))
    app.logger.setLevel(logging.INFO)
    file_handler.setLevel(logging.INFO)
    app.logger.addHandler(file_handler)
    app.logger.info('microblog startup')
"""
