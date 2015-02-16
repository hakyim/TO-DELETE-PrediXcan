#!flask/bin/python
from app import app
from OpenSSL import SSL

"""
context = SSL.Context(SSL.TLSv1_METHOD)
context.use_privatekey_file('server.key')
context.use_certificate_file('server.crt')
"""
app.run(debug=True,port=80,host='0.0.0.0.')
