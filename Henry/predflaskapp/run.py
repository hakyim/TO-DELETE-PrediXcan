#!flask/bin/python
from app import app
from OpenSSL import SSL
import ssl
import traceback
from tornado.wsgi import WSGIContainer
from tornado.httpserver import HTTPServer
from tornado.ioloop import IOLoop
import os

"""
context = SSL.Context(SSL.TLSv1_METHOD)
context.use_privatekey_file('server.key')
context.use_certificate_file('server.crt')
"""

app.debug = True 

context = ssl.SSLContext(ssl.PROTOCOL_SSLv23)
context.load_cert_chain('server.crt','server.key')
#app.run(debug=True,ssl_context=context,port=443,host='0.0.0.0')


settings = { "certfile": os.path.join("server.crt"), "keyfile": os.path.join("server.key") }

http_server = HTTPServer(WSGIContainer(app),ssl_options=settings)
http_server.listen(443)
IOLoop.instance().start()
