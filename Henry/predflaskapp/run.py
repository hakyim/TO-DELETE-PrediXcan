#!flask/bin/python
from app import app
from OpenSSL import SSL
import ssl
import traceback
from tornado.web import FallbackHandler, RequestHandler, Application
from tornado.wsgi import WSGIContainer
from tornado.httpserver import HTTPServer
from tornado.ioloop import IOLoop
from tornado.log import enable_pretty_logging
import os


class MainHandler(RequestHandler):
    def get(self):
        self.write("This message comes from Tornado!")

context = ssl.SSLContext(ssl.PROTOCOL_SSLv23)
context.load_cert_chain('server.crt','server.key')
#app.run(debug=True,ssl_context=context,port=443,host='0.0.0.0')

app_wrap = WSGIContainer(app)

application = Application([
	(r"/tornado",MainHandler),
	(r".*",FallbackHandler,dict(fallback=app_wrap))
],debug=True)

settings = { "certfile": os.path.join("server.crt"), "keyfile": os.path.join("server.key") }

enable_pretty_logging()
http_server = HTTPServer(application,ssl_options=settings,max_buffer_size=2**30)
http_server.listen(443, address="0.0.0.0")
IOLoop.instance().start()
