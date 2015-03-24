import os
import helpfuncs 
from flask import render_template, flash, redirect, session, url_for, request, g, send_from_directory
from flask.ext.login import login_user, logout_user, current_user, login_required
from app import app, db, lm, config
from forms import LoginForm, EditForm, PostForm, CommandGenForm
from models import User, Post
from datetime import datetime
import MySQLdb as db 
from werkzeug import secure_filename
import prediction
import tarfile 
import zipfile 

"""temporary globals"""
ALLOWED_EXTENSIONS = set(['txt', 'pdf', 'png', 'jpg', 'jpeg', 'gif','gz','tar'])
UPLOAD_FOLDER = 'puploads'


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
		p = Post(body=body, timestamp = datetime.utcnow(), author = g.user) #correct author ?
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
		tissue = form.tissuetype.data
		study = form.study.data
		output = helpfuncs._generateCommand(pfp,gdfp,gth,gtt,tissue,study)
		
		flash ("data i got: %s %s %s %s %s %s" % (pfp,gdfp,gth,gtt,tissue,study))
		flash (output)
		return redirect(url_for('gen_command'))
		#printCommandline(stuff from forms)
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
			flash(ALLOWED_EXTENSIONS	)
	return render_template('uploadfile.html')

@app.route('/tarupload',methods=["POST","GET"])
def tar_upload():
	if request.method == 'POST':
		uploaded_tar = request.files["file"]
		if uploaded_tar and allowed_file(uploaded_tar.filename):
			tarname = secure_filename(uploaded_tar.filename)

			path = "./puploads/" + str(tarname.rsplit('.',1)[0]) + '/'
			if not os.path.isdir(path):
				os.mkdir(path)
			uploaded_tar.save(os.path.join(path,tarname))
			print tarname
			print path
			tar = tarfile.open(path + tarname,'r')
			tar.extractall(path)
			os.listdir(UPLOAD_FOLDER)
			os.listdir(path)
			for f in os.listdir(path):
				print f
				print f.rsplit(".",1)
			files=[f for f in os.listdir(path) if f.rsplit('.',1)[1] != "tar"]
			print files
			return render_template("tarupload.html",files=files)
	return render_template("tarupload.html")

@app.route('/uploads/<filename>')
def uploaded_file(filename):
	return send_from_directory(app.config['UPLOAD_FOLDER'],filename)

"""helper function - saves files in tar and returns list of files in it"""
def _save_tar(tarfile):
	if tarfile and allowed_file(tarfile.filename):	
		tarname = secure_filename(tarfile.filename)
		path = "./puploads/" + str(tarname.rsplit('.',1)[0]) + '/'
		if not os.path.isdir(path):
			os.mkdir(path)
		uploaded_tar.save(os.path.join(path,tarname))
		print tarname
		print path
		tar = tarfile.open(path + tarname,'r')
		tar.extractall(path)
		os.listdir(UPLOAD_FOLDER)
		os.listdir(path)
		files=[f for f in os.listdir(path) if f.rsplit('.',1)[1] != "tar"]
		print files
		return (files,tarname)
	else:
		print "error, could not verify tarfile"
		return None



"""dont use this yet! still need to hook in main model code"""
@app.route('/predict',methods=["POST","GET"])
def predict_test():
	form = predictForm()
	if request.method == 'POST':
		uploaded_tar = form.tarfile#["file"]
		files,tarname = _save_tar(uploaded_tar) #what should return be?	
		if files:
			path = "./puploads/" + str(tarname.rsplit('.',1)[0]) + '/'
			prefix = form.dosageprefix
			#more stuff happens
			#fill in...
			predictor = px.prediction_maker(gene_list=None,dosage_dir=path,dosage_prefix=prefix)
			#then check the file
	return render_template("predict.html",form=form)



"""Here lie error handlers"""
@app.errorhandler(404)
def not_found_error(error):
	return render_template('404.html'), 404

@app.errorhandler(500)
def internal_error(error):
	db.session.rollback()
	return render_template('500.html'), 500
