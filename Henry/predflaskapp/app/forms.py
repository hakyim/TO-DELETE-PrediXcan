from flask.ext.wtf import Form
from wtforms.fields import StringField, BooleanField, TextAreaField, SelectField, FileField, SubmitField
from flask.ext.wtf.file import FileAllowed, FileRequired 
from wtforms.validators import DataRequired, Length
from app.models import User


class LoginForm(Form):
	openid = StringField('openid', validators=[DataRequired()])
	remember_me = BooleanField('remember_me', default=False)

class EditForm(Form):
	nickname = StringField('nickname', validators=[DataRequired()])
	about_me = TextAreaField('about_me', validators=[Length(min=0,max=140)])

	def __init__(self,original_nickname, *args, **kwargs):
		Form.__init__(self,*args,**kwargs)
		self.original_nickname = original_nickname

	def validate(self):
		if not Form.validate(self):
			return False
		if self.nickname.data == self.original_nickname:
			return True 
		user = User.query.filter_by(nickname=self.nickname.data).first()
		if user != None:
			self.nickname.errors.append("This nickname is already in use. Pick another.")
			return False 
		return True

class PostForm(Form):
	post_text = TextAreaField('post', validators=[DataRequired(),Length(min=1,max=256)]) 	

class CommandGenForm(Form):
	phenofilepath = TextAreaField('phenofilepath',validators=[DataRequired(),Length(min=1,max=1024)])
	genedatafilepath = TextAreaField('genedatafilepath',validators=[DataRequired(),Length(min=1,max=1024)])
	genotypeheader = TextAreaField('genotypeheader',validators=[DataRequired(),Length(min=1,max=1024)])
	genotypetail = TextAreaField('genotypetail', validators=[DataRequired(),Length(min=1,max=1024)])
	#genotypefilepath = TextAreaField('genotypefilepath',validators=[DataRequired(),Length(min=1,max=1024)])
	tissuetype = SelectField(u'Tissue')
	study = SelectField(u'Study')
	 
class predictForm(Form):
	tarfile = FileField('tarfile')
	genelist = FileField('genelist',default=None)
	prefix = TextAreaField('dosageprefix')
	submit = SubmitField("submit")
