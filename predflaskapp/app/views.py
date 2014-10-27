from flask import render_template, flash, redirect
from app import app
from .forms import LoginForm

@app.route('/login',methods=['GET','POST'])
def login():
	form = LoginForm()
	if form.validate_on_submit():
		flash('Login requested for jankyID="%s", remember_me=%s' %	
			(form.openid.data, str(form.remember_me.data)))
		return redirect('/index') 
	return render_template('login.html',title='Sign In',form=form,providers=app.config['OPENID_PROVIDERS'])

@app.route('/')
@app.route('/index')
def index():
    user = {'nickname': 'Henry'}  
    posts = [ 
	    {
	    	'author':{'nickname':'John'},
	    	'body': 'Beautiful day in Portland!'

	    },
	    {
	    	'author': {'nickname':'Susan'},
	    	'body': 'Whats going on?'
	    }
    ]
    return render_template('index.html',
    						title='home',
    						user=user,
    						posts=posts)		


