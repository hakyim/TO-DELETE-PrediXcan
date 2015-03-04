#!flask/bin/python
from app import app
from OpenSSL import SSL
import ssl
import traceback

"""
context = SSL.Context(SSL.TLSv1_METHOD)
context.use_privatekey_file('server.key')
context.use_certificate_file('server.crt')
"""

context = ssl.SSLContext(ssl.PROTOCOL_SSLv23)

context.load_cert_chain('server.crt','server.key')
app.run(debug=True,ssl_context=context,port=443,host='0.0.0.0')
