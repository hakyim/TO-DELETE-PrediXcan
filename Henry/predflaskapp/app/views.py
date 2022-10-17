import os
import sys
from forms import LoginForm, EditForm, PostForm, CommandGenForm, predictForm
from flask import render_template, flash, redirect, session, url_for, request, g, send_from_directory, make_response
from flask.ext.login import login_user, logout_user, current_user, login_required
from app import app, db, lm
from models import User, Post
from datetime import datetime
import MySQLdb as db 	
from werkzeug import secure_filename
import prediction
import tarfile 
import zipfile 
import prediction as px
import helpfuncs 
import config
import celery

"""temporary globals"""
ALLOWED_EXTENSIONS = set(['txt', 'pdf', 'png', 'jpg', 'jpeg', 'gif','gz','tar'])
UPLOAD_FOLDER = 'puploads/'


"""
@app.route('/login',methods=['GET','POST'])
@oid.loginhandler
def login():
	if g.user is not None and g.user.is_authenticated():
		return redirect(url_for('index'))
	form = LoginForm()
	if form.validate_on_submit():
		session['remember_me'] = form.remember_me.data
		return oid.try_login(form.openid.data, ask_for=['nickname', 'email'])
		
		flash('Login requested for jankyID="%s", remember_me=%s' %	
			(form.openid.data, str(form.remember_me.data)))
		return redirect('/index') 
		
	return render_template('login.html',title='Sign In',form=form,providers=app.config['OPENID_PROVIDERS'])
"""

@app.route('/edit', methods=['GET','POST'])
@login_required
def edit():
	form = EditForm(g.user.nickname)
	if form.validate_on_submit():
		g.user.nickname = form.nickname.data
		g.user.about_me = form.about_me.data
		db.session.add(g.user)
		db.session.commit()
		flash("Ch-ch-ch changes saved.")
		#return redirect(url_for('edit'))
		return redirect('user/' + g.user.nickname)
	else:
		form.nickname.data = g.user.nickname
		form.about_me.data = g.user.about_me
	return render_template('edit.html',form=form)

@app.route('/newpost', methods=['GET','POST'])
@login_required 
def post():
	form = PostForm()
	if form.validate_on_submit(): 
		body = form.post_text.data 
		p = Post(body=body, timestamp = datetime.utcnow(), author = g.user)
		db.session.add(p)
		db.session.commit()
		flash("Successfully posted!")
		return redirect('user/' + g.user.nickname)
	
	return render_template('newpost.html',form=form) # jump to the actual url


"""
@oid.after_login
def after_login(resp):
	if resp.email is None or resp.email == "":
		flash('invalid login. please try again.')
		return redirect(url_for('login'))
	user = User.query.filter_by(email=resp.email).first()
	if user is None:
		nickname = resp.nickname
		if nickname is None or nickname == "":
			nickname = resp.email.split('@')[0]
		nickname = User.make_unique_nickname(nickname)
		user = User(nickname=nickname,email=resp.email)
		db.session.add(user)
		db.session.commit()
	remember_me = False

	if 'remember_me' in session:
		remember_me = session['remember_me']
		session.pop('remember_me',None)
	login_user(user, remember = remember_me)
	return redirect(request.args.get('next') or url_for('index'))

"""
@app.route('/')
@app.route('/index')
def index():
	return render_template('index.html',
							title='home')


"""here lie user profiles, leaving this code in case profiles are implemented"""
@app.route('/user/<nickname>')
@login_required
def user(nickname):
	user = User.query.filter_by(nickname=nickname).first()
	if user == None: 
		flash('user %s not found.' % nickname)
		return redirect(url_for('index'))
	posts = user.posts.all()
	return render_template('user.html', user=user, posts=posts)


@app.before_request
def before_request():
	g.user = current_user
	if g.user.is_authenticated():
		g.user.last_seen = datetime.utcnow()
		db.session.add(g.user)
		db.session.commit()

@lm.user_loader
def load_user(id):
	return User.query.get(int(id))

@app.route('/logout')
def logout():
	logout_user()
	return redirect(url_for('index'))


"""Predixcan functionality"""
@app.route('/cmdgen', methods=['GET','POST'])
def gen_command():
	form = CommandGenForm()
	database = db.connect(host="192.170.232.66", # your host 
					 user='public', # your username
					  passwd='foobar', # your password
					  db="mysql",port=3306) # name of the data base
	
	form.tissuetype.choices = helpfuncs._getTissueTypes(database) #fetch tissue types from DB
	form.study.choices = helpfuncs._getStudyNames(database) #fetch study names from DB

	if form.validate_on_submit():
		pfp = form.phenofilepath.data
		gdfp = form.genedatafilepath.data
		gth = form.genotypeheader.data 
		gtt = form.genotypetail.data
		tissue = form.tissuetyp
		e.data
		study = form.study.data		
		output = helpfuncs._generateCommand(pfp,gdfp,gth,gtt,tissue,study)
		flash (output)
		return redirect(url_for('gen_command'))
			
	return render_template('cmdgen.html',form=form) #TBA


def allowed_file(filename):
	return '.' in filename and \
		   filename.rsplit('.', 1)[1] in app.config['ALLOWED_EXTENSIONS']

"""file upload functions/tests"""
@app.route('/fileupload',methods=["POST","GET"])
def upload_file():
	if request.method == 'POST':
		file = request.files['file']
		print file.filename
		if file and allowed_file(file.filename):
			filename = secure_filename(file.filename)
			file.save(os.path.join(UPLOAD_FOLDER,filename))
			for thing in os.listdir(UPLOAD_FOLDER):
				print thing
			return render_template('uploadfile.html',filename=filename,filetext=open("/tmp/"+filename,"r").read())
		else:
			flash("Invalid file type, please try again. Allowed extensions are:")
			flash(ALLOWED_EXTENSIONS)

	return render_template('uploadfile.html')

@app.route('/tarupload',methods=["POST","GET"])
def tar_upload():
	#REQUIRES tar contain folder with same FILENAME as tarfile
	if request.method == 'POST':
		uploaded_tar = request.files["file"]
		if uploaded_tar and allowed_file(uploaded_tar.filename):
			tarname = secure_filename(uploaded_tar.filename)
			path = UPLOAD_FOLDER + str(tarname.rsplit('.',1)[0]) + '/'
			#TODO: Check if file is aleady there
			uploaded_tar.save(UPLOAD_FOLDER+tarname)
			tar = tarfile.open(UPLOAD_FOLDER + tarname,'r')
			tar.extractall(UPLOAD_FOLDER)
			files=[f for f in os.listdir(path) if f.rsplit('.',1)[1] != "tar"]
			#TODO: remove files 
			return render_template("tarupload.html",files=files)

	return render_template("tarupload.html")

@app.route('/uploads/<filename>')
def uploaded_file(filename):
	return send_from_directory(app.config['UPLOAD_FOLDER'],filename)

"""helper function - saves files in tar and returns list of files in it"""
def _save_tar(tarfile_):
	filename = tarfile_.data.filename
	print "file name in _save_tar is %s" % filename
	if filename and allowed_file(filename):	
		tarname = secure_filename(filename)
		path = UPLOAD_FOLDER + str(tarname.rsplit('.',1)[0]) + '/'
		#TODO: Check if file is aleady there
		tarfile_.data.save(UPLOAD_FOLDER+tarname)
		tar = tarfile.open(UPLOAD_FOLDER + tarname,'r')
		tar.extractall(UPLOAD_FOLDER)
		files=[f for f in os.listdir(path) if f.rsplit('.',1)[1] != "tar"]
		return (files,tarname)
	else:
		print "error, could not verify tarfile"
		return (None,None)

@app.route('/predict',methods=["POST","GET"])
def run_predict():
	print "predict one"
	form = predictForm()
	if request.method == 'POST':
		if form.validate_on_submit():
			print "validated form"
			uploaded_tar = form.tarfile
			prefix = form.prefix.data
			files,tarname = _save_tar(uploaded_tar) #returns tuple of files and tarname	
			if files:
				path = "./puploads/" + str(tarname.rsplit('.',1)[0]) + '/'
				if form.genelist.data:
						glistname = form.genelist.data.filename
						form.genelist.data.save(UPLOAD_FOLDER+glistname)
				else:
						glistname = None
				predictor = px.prediction_maker(gene_list=glistname,dosage_dir=path,dosage_prefix=prefix)				
				#queue up a job?
				fname = doPredictionJob.delay(predictor)
				#predictor.do_predictions()
				return render_template("predict.html",filename=fname,form=form)
				#TODO: delete PredXResult after delivery
				#TODO: encrypt,email
		else:
			flash("bad file upload, for some godforsaken reason")

	return render_template("predict.html",form=form)

#takes in a precitor object, runs it
@celery.task
def doPredictionJob(predictor):
	print "starting a predict job"
	fname = predictor.do_predictions()
	return fname

@app.route('/<path:filename>',methods=['GET','POST'])
def download(filename):
	f = open(os.path.join(app.static_folder,filename))
	contents = f.read()
	response = make_response(contents)
	#set it to be a download response
	response.headers["Content-Disposition"] = "attachment; filename=%s" % filename
	return response


"""Here lie error handlers"""
@app.errorhandler(404)
def not_found_error(error):
	return render_template('404.html'), 404

@app.errorhandler(500)
def internal_error(error):
	return render_template('500.html'), 500
